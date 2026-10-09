"""Frozen-protocol extension: five unused DS10440 runs. No original output writes."""
from pathlib import Path
import os, sys, json, time, hashlib, glob, traceback
os.environ.setdefault('OMP_NUM_THREADS','4')
import numpy as np
import pandas as pd
import mrcfile
from scipy.ndimage import gaussian_filter, binary_dilation
from scipy.spatial.distance import cdist
from skimage.feature import match_template, peak_local_max
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
sys.path.insert(0,str(BASE/'code'))
from metrics import pick_prf, edge_sharpness
RUNS=json.loads((ROOT/'protocol.json').read_text())['new_CZII_runs']
META=json.loads((BASE/'data/czii_TS_5_4_meta.json').read_text())
RADII=META['radii_vox']
def dump(p,obj): p.write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def masks(shape,points,r):
    sig=np.zeros(shape,bool)
    for xyz in points:
        center=np.array(xyz)[::-1]
        lo=np.maximum(0,np.floor(center-r).astype(int))
        hi=np.minimum(shape,np.ceil(center+r).astype(int)+1)
        z,y,x=np.ogrid[lo[0]:hi[0],lo[1]:hi[1],lo[2]:hi[2]]
        sig[tuple(slice(a,b) for a,b in zip(lo,hi))] |= (z-center[0])**2+(y-center[1])**2+(x-center[2])**2<=r*r
    return sig,binary_dilation(sig,iterations=int(r))&~binary_dilation(sig,iterations=2)
def selftest():
    shape=(20,21,22);pts=[[.4,3,5],[15.2,19,18]];r=3.2
    z,y,x=np.indices(shape);ref=np.zeros(shape,bool)
    for xx,yy,zz in pts:ref |= (z-zz)**2+(y-yy)**2+(x-xx)**2<=r*r
    sig,bg=masks(shape,pts,r)
    assert np.array_equal(ref,sig)
    assert np.array_equal(bg,binary_dilation(ref,iterations=int(r))&~binary_dilation(ref,iterations=2))
def freeze():
    source=BASE/'revision_20261003_peng_followup/picking_threshold_sweep.csv'
    curves=pd.read_csv(source);thresholds={}
    for cls in RADII:
        thresholds[cls]={}
        for method in ['raw','topaz']:
            rows=curves[(curves.run=='TS_5_4')&(curves.method==method)&(curves.target==cls)].sort_values('threshold')
            thresholds[cls][method]=float(rows.loc[rows.f1.idxmax(),'threshold'])
    spec=dict(calibration_run='TS_5_4',source=str(source),source_sha256=sha(source),thresholds=thresholds,radii_vox=RADII,
              selection='max calibration F1, lowest threshold on ties; no new-run labels used',protocols=['raw_fixed','per_input'])
    p=ROOT/'frozen_thresholds.json'
    if p.exists():assert json.loads(p.read_text())==spec
    else:dump(p,spec)
    return thresholds
def preprocess():
    import torch, topaz, topaz.denoising.models as tdm
    from importlib.metadata import version
    torch.set_num_threads(4);weight_info=[]
    def load(pkg,path,map_location='cpu'):
        hits=list(Path(topaz.__file__).parent.rglob(Path(path).name));assert len(hits)==1
        weight_info.append(dict(path=str(hits[0]),sha256=sha(hits[0])))
        return torch.load(hits[0],map_location=map_location,weights_only=False)
    tdm.load_state_dict_from_pkg=load
    from topaz.denoise import Denoise3D,denoise_tomogram
    model=Denoise3D('unet-3d-10a',use_cuda=torch.cuda.is_available(),dims=3)
    dump(ROOT/'preprocess_config.json',dict(model='unet-3d-10a',weights=weight_info,torch=torch.__version__,topaz=version('topaz-em'),cuda=torch.cuda.is_available(),patch_size=96,padding=24,gaussian_sigma=1.5,input='global z-score; Topaz also performs its standard normalization',code_sha256=sha(__file__)))
    allrows=[]
    for run in RUNS:
        d=ROOT/'data'/run;out=d/'processed';out.mkdir(exist_ok=True)
        if (out/'metrics.json').exists():allrows.extend(json.loads((out/'metrics.json').read_text()));continue
        with mrcfile.open(d/'raw.mrc',permissive=True) as m:raw=m.data.copy().astype('float32')
        raw=(raw-raw.mean())/(raw.std()+1e-8)
        pts=json.loads((d/'points.json').read_text());sig,bg=masks(raw.shape,pts['cytosolic ribosome'],RADII['cytosolic ribosome'])
        pin=out/'normalized.mrc'
        with mrcfile.new(pin,overwrite=True) as m:m.set_data(raw);m.voxel_size=10.012
        ts=time.time();g=gaussian_filter(raw,1.5);tg=time.time()-ts
        tout=out/'topaz';tout.mkdir(exist_ok=True)
        ts=time.time();denoise_tomogram(str(pin),model,outdir=str(tout),patch_size=96,padding=24,verbose=False);tt=time.time()-ts
        files=list(tout.glob('*.mrc'));assert len(files)==1
        with mrcfile.open(files[0],permissive=True) as m:top=m.data.copy()
        assert top.shape==raw.shape and np.isfinite(top).all()
        np.save(out/'topaz.npy',top)
        rows=[]
        for name,v,t in [('raw',raw,0),('gaussian',g,tg),('topaz',top,tt)]:
            s,b=v[sig],v[bg];delta=abs(float(s.mean()-b.mean()))
            row=dict(run=run,method=name,n_ribosomes=len(pts['cytosolic ribosome']),SNR_region=delta/(float(b.std())+1e-8),CNR=delta/(float(np.sqrt(s.var()+b.var()))+1e-8),sharpness=edge_sharpness(v),runtime_s=t)
            rows.append(row);print('[sample]',row,flush=True)
        dump(out/'metrics.json',rows);allrows.extend(rows)
        pd.DataFrame(allrows).to_csv(ROOT/'extra_sample_metrics.csv',index=False)
        del raw,g,top,sig,bg
    pd.DataFrame(allrows).to_csv(ROOT/'extra_sample_metrics.csv',index=False)
    print('PREPROCESS COMPLETE',flush=True)
def independent_tp(pred,gt,tol):
    used=set();tp=0
    for row in cdist(pred,gt):
        for j in np.argsort(row,kind='stable'):
            if row[j]>tol:break
            if int(j) not in used:used.add(int(j));tp+=1;break
    return tp
def picking():
    thresholds=freeze();rows=[]
    for run in RUNS:
        d=ROOT/'data'/run;out=d/'processed';cache=out/'picking';cache.mkdir(exist_ok=True)
        pts=json.loads((d/'points.json').read_text())
        for method in ['raw','topaz']:
            if method=='raw':
                with mrcfile.open(out/'normalized.mrc') as m:v=m.data.copy()
            else:v=np.load(out/'topaz.npy')
            for cls,r in RADII.items():
                rr=round(r);dist=max(1,int(r));p=cache/f'{method}_r{rr}_d{dist}.npz'
                if p.exists():
                    c=np.load(p);xyz,scores,mx=c['xyz'],c['scores'],float(c['maximum'])
                else:
                    start=time.time();z,y,x=np.mgrid[-rr:rr+1,-rr:rr+1,-rr:rr+1]
                    template=-((z*z+y*y+x*x)<=rr*rr).astype('float32');template-=template.mean()
                    response=match_template(v,template,pad_input=True);mx=float(response.max())
                    peaks=peak_local_max(response,min_distance=dist,threshold_rel=.30,exclude_border=False)
                    scores=response[tuple(peaks.T)];xyz=peaks[:,::-1].astype(float)
                    np.savez_compressed(p,xyz=xyz,scores=scores,maximum=mx);del response
                    print('[candidates]',run,method,rr,dist,len(xyz),round(time.time()-start,1),'s',flush=True)
                truth=np.asarray(pts[cls],float).reshape(-1,3)
                for protocol in ['raw_fixed','per_input']:
                    t=thresholds[cls]['raw' if protocol=='raw_fixed' else method]
                    pred=xyz[scores>t*mx];score=pick_prf(pred,truth,r)
                    check=independent_tp(pred,truth,r);assert score['tp']==check,(run,method,cls,score,check)
                    rows.append(dict(run=run,method=method,target=cls,protocol=protocol,threshold=t,n_gt=len(truth),n_pred=len(pred),independent_tp=check,**score))
            del v
            pd.DataFrame(rows).to_csv(ROOT/'extra_sample_picking.csv',index=False)
        print('[picking done]',run,flush=True)
    df=pd.DataFrame(rows)
    df.groupby(['run','protocol','method']).f1.mean().reset_index(name='macro_F1').to_csv(ROOT/'extra_sample_picking_macro.csv',index=False)
    dump(ROOT/'samples_complete.json',dict(runs=RUNS,rows=len(df),independent_matcher_all_passed=True,completed=time.strftime('%Y-%m-%dT%H:%M:%S'),script_sha256=sha(__file__)))
    print('PICKING COMPLETE',flush=True)
if __name__=='__main__':
    selftest();freeze()
    if len(sys.argv)<2 or sys.argv[1]=='preprocess':preprocess()
    if len(sys.argv)>1 and sys.argv[1]=='picking':picking()
