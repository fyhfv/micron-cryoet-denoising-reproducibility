"""Frame-disjoint, held-out-region diagnostic; never treats a processed map as truth."""
from pathlib import Path
import os,sys,json,time,hashlib
os.environ.setdefault('OMP_NUM_THREADS','4')
import numpy as np,pandas as pd,mrcfile
from scipy.ndimage import gaussian_filter
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'code'))
from metrics import fsc,fsc_crossing
def dump(p,x):p.write_text(json.dumps(x,indent=2),encoding='utf-8')
def fullfft_fsc(a,b):
    a=a.astype(float);b=b.astype(float);fa=np.fft.fftn(a);fb=np.fft.fftn(b)
    grids=np.meshgrid(*[np.fft.fftfreq(n) for n in a.shape],indexing='ij')
    rad=np.sqrt(sum(g*g for g in grids));n=min(a.shape)//2
    valid=(rad>0)&(rad<=.5);idx=np.minimum((rad[valid]/(.5/n)).astype(int),n-1)
    num=np.bincount(idx,weights=(fa[valid]*fb[valid].conj()).real,minlength=n)
    aa=np.bincount(idx,weights=abs(fa[valid])**2,minlength=n);bb=np.bincount(idx,weights=abs(fb[valid])**2,minlength=n)
    out=np.full(n,np.nan);ok=(aa>aa.max()*1e-12)&(bb>bb.max()*1e-12)
    out[ok]=num[ok]/np.sqrt(aa[ok]*bb[ok]);return out
def baselines():
    import torch,topaz,topaz.denoising.models as tdm
    torch.set_num_threads(4)
    def load(pkg,path,map_location='cpu'):
        hits=list(Path(topaz.__file__).parent.rglob(Path(path).name));assert len(hits)==1
        return torch.load(hits[0],map_location=map_location,weights_only=False)
    tdm.load_state_dict_from_pkg=load
    from topaz.denoise import Denoise3D
    model=Denoise3D('unet-3d-10a',use_cuda=torch.cuda.is_available(),dims=3)
    spec=json.loads((ROOT/'cryocare_protocol.json').read_text(encoding='utf-8'))
    cfg=dict(centers_zyx=spec['centers_zyx'],context=96,scoring_center=64,gaussian_sigma=1.5,topaz='unet-3d-10a',patch_size=96,padding=24,
             summary_band=[.05,.20],normalization='cross-half MSE in original units divided by raw even-odd MSE; no test-target intensity fitting',
             limitations='Shared alignment/reconstruction and same fixed fitted network can induce common signal/bias; FSC is not gold-standard STA or verified physical resolution. Topaz pretraining overlap with Tomo110 cannot be excluded.')
    dump(ROOT/'half_metrics_protocol.json',cfg)
    out=ROOT/'half_evaluation';out.mkdir(exist_ok=True)
    a=mrcfile.mmap(ROOT/'data/tutorial_data/tomo_even_frames.rec',permissive=True)
    b=mrcfile.mmap(ROOT/'data/tutorial_data/tomo_odd_frames.rec',permissive=True)
    allmap=mrcfile.mmap(ROOT/'data/tutorial_data/tomo_all_frames.rec',permissive=True)
    sourcechecks=[]
    for i,center in enumerate(cfg['centers_zyx']):
        s=tuple(slice(c-48,c+48) for c in center);x=np.array(a.data[s]);y=np.array(b.data[s]);ts=time.time()
        p=out/f'baseline_{i}.npz'
        if not p.exists():
            ga=gaussian_filter(x,1.5);gb=gaussian_filter(y,1.5)
            ta=model.denoise(x,patch_size=96,padding=24,verbose=False);tb=model.denoise(y,patch_size=96,padding=24,verbose=False)
            np.savez_compressed(p,even=x,odd=y,gaussian_even=ga,gaussian_odd=gb,topaz_even=ta,topaz_odd=tb)
        sourcechecks.append(dict(crop=i,even_odd_identical=bool(np.array_equal(x,y)),even_odd_mse=float(np.mean((x.astype(float)-y)**2)),average_to_all_mse=float(np.mean(((x.astype(float)+y)/2-allmap.data[s])**2)),seconds=time.time()-ts))
        print('[half baseline]',i,flush=True)
    a.close();b.close();allmap.close();dump(ROOT/'half_source_checks.json',sourcechecks)
def evaluate():
    assert (ROOT/'cryocare_complete.json').exists()
    cfg=json.loads((ROOT/'half_metrics_protocol.json').read_text());rows=[];curves=[];checks=[]
    for i in range(len(cfg['centers_zyx'])):
        base=np.load(ROOT/f'half_evaluation/baseline_{i}.npz');care=np.load(ROOT/f'half_evaluation/crop_{i}.npz')
        assert np.array_equal(base['even'],care['even']) and np.array_equal(base['odd'],care['odd'])
        sl=(slice(16,80),)*3;a=base['even'][sl];b=base['odd'][sl]
        rawmse=float(np.mean((a.astype(float)-b)**2))
        for method in ['raw','gaussian','topaz','cryocare']:
            if method=='raw':pa,pb=a,b
            else:
                source=care if method=='cryocare' else base
                pa,pb=source[method+'_even'][sl],source[method+'_odd'][sl]
            assert np.isfinite(pa).all() and np.isfinite(pb).all()
            mse=.5*(np.mean((pa.astype(float)-b)**2)+np.mean((pb.astype(float)-a)**2))
            f,c=fsc(pa-pa.mean(),pb-pb.mean());direct=fullfft_fsc(pa-pa.mean(),pb-pb.mean())
            assert np.array_equal(np.isfinite(c),np.isfinite(direct))
            err=float(np.nanmax(abs(c-direct)));assert err<1e-5
            checks.append(dict(crop=i,method=method,max_abs_fsc_error=err))
            band=(f>=.05)&(f<=.20);cross=fsc_crossing(f,c,.143)
            row=dict(crop=i,method=method,cross_half_mse=float(mse),raw_pair_mse=rawmse,cross_mse_ratio=float(mse/rawmse),fsc_band_mean=float(np.nanmean(c[band])),fsc_cross_status=cross['status'],crossing_cycles_voxel=cross['frequency'],frequency_bound=cross['frequency_bound'])
            rows.append(row)
            curves.extend(dict(crop=i,method=method,window='none',frequency=float(ff),fsc=float(cc)) for ff,cc in zip(f,c))
            # Boundary-sensitivity amendment after inspecting the initial hard-crop
            # curves: fixed separable Hann taper, with no threshold/band selection.
            h=np.hanning(64);window=h[:,None,None]*h[None,:,None]*h[None,None,:]
            # Float64 avoids roundoff dominating strongly attenuated high-frequency
            # shells after Gaussian filtering plus tapering.
            wa=(pa.astype(float)-float(pa.mean()))*window;wb=(pb.astype(float)-float(pb.mean()))*window
            hf,hc=fsc(wa,wb);hd=fullfft_fsc(wa,wb)
            assert np.array_equal(np.isfinite(hc),np.isfinite(hd))
            he=float(np.nanmax(abs(hc-hd)));assert he<1e-5
            checks.append(dict(crop=i,method=method,window='hann',max_abs_fsc_error=he))
            hx=fsc_crossing(hf,hc,.143)
            row.update(hann_fsc_band_mean=float(np.nanmean(hc[band])),hann_fsc_cross_status=hx['status'],hann_crossing_cycles_voxel=hx['frequency'],hann_frequency_bound=hx['frequency_bound'])
            curves.extend(dict(crop=i,method=method,window='hann',frequency=float(ff),fsc=float(cc)) for ff,cc in zip(hf,hc))
    df=pd.DataFrame(rows);df.to_csv(ROOT/'half_metrics.csv',index=False)
    pd.DataFrame(curves).to_csv(ROOT/'half_fsc_curves.csv',index=False)
    summary=df.groupby('method')[['cross_mse_ratio','fsc_band_mean','hann_fsc_band_mean']].agg(['mean','std'])
    summary.to_csv(ROOT/'half_summary.csv')
    dump(ROOT/'half_verification.json',dict(curves_checked=len(checks),max_abs_fullfft_error=max(x['max_abs_fsc_error'] for x in checks),checks=checks,raw_ratio_all_one=bool((df.loc[df.method=='raw','cross_mse_ratio']==1).all()),note='Six spatial crops of ONE tomogram; descriptive SD only, not biological inference.'))
    print(summary.to_string(),flush=True)
if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='evaluate':evaluate()
    else:baselines()
