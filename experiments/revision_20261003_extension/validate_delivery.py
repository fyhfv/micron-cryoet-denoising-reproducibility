"""Audit new numerical claims and preserved manuscript assets; render changed pages."""
from pathlib import Path
import json,re,hashlib,difflib,shutil,zipfile
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parent;PROJECT=ROOT.parents[1]
OLD=PROJECT/'正式投稿/Micron_submission_20261003_Peng_literature'
OUT=PROJECT/'正式投稿/Micron_submission_20261003_Peng_extension'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False),encoding='utf-8')
def audit():
    summary=json.loads((ROOT/'extension_summary.json').read_text(encoding='utf-8'))
    assert sha(OLD/'manuscript.tex')==summary['source_manuscript_sha256']
    assert sha(OLD/'Supplementary_methods.tex')==summary['source_supplement_sha256']
    preserved=[p.name for p in OLD.glob('fig_*.pdf')]+['Title_page.pdf','Title_page.tex','Figure_captions.pdf','Figure_captions.tex','table_quant.tex']
    assert all(sha(OLD/f)==sha(OUT/f) for f in preserved)
    old=(OLD/'manuscript.tex').read_text(encoding='utf-8');new=(OUT/'manuscript.tex').read_text(encoding='utf-8')
    assert old.split(r'\section*{Abstract}')[0]==new.split(r'\section*{Abstract}')[0]
    assert '62175037' in new and 'Sihua.Peng@uga.edu' in new
    assert old.split(r'\begin{thebibliography}')[1]==new.split(r'\begin{thebibliography}')[1]
    keys=re.findall(r'\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}',new)
    cites={x.strip() for g in re.findall(r'\\cite\w*(?:\[[^\]]*\])?\{([^}]+)\}',new) for x in g.split(',')}
    assert len(keys)==len(set(keys)) and not cites-set(keys)
    assert re.findall(r'\\label\{fig:[^}]+\}',old)==re.findall(r'\\label\{fig:[^}]+\}',new)
    for name in ['manuscript','Supplementary_methods']:
        log=(OUT/(name+'.log')).read_text(errors='replace')
        issues=[l for l in log.splitlines() if any(w in l for w in ['Overfull','undefined','multiply defined'])]
        if issues:print(name,issues)
        assert not issues
        delta=''.join(difflib.unified_diff((OLD/(name+'.tex')).read_text(encoding='utf-8').splitlines(True),(OUT/(name+'.tex')).read_text(encoding='utf-8').splitlines(True),fromfile='prior/'+name+'.tex',tofile='extension/'+name+'.tex'))
        (OUT/(name+'_extension.diff')).write_text(delta,encoding='utf-8')
    p=pd.read_csv(ROOT/'extra_sample_picking.csv')
    assert np.allclose(p.f1,2*p.tp/(2*p.tp+p.fp+p.fn),rtol=0,atol=1e-14)
    macro=p.groupby(['run','protocol','method']).f1.mean().sort_index()
    exported=pd.read_csv(ROOT/'extra_sample_picking_macro.csv').set_index(['run','protocol','method']).macro_F1.sort_index()
    assert np.allclose(macro,exported,rtol=0,atol=1e-14)
    h=pd.read_csv(ROOT/'half_metrics.csv');errors=[]
    for i in range(6):
        b=np.load(ROOT/f'half_evaluation/baseline_{i}.npz');c=np.load(ROOT/f'half_evaluation/crop_{i}.npz');sl=(slice(16,80),)*3
        a=b['even'][sl].astype(float);z=b['odd'][sl].astype(float);den=np.mean((a-z)**2)
        for method in ['raw','gaussian','topaz','cryocare']:
            src=c if method=='cryocare' else b
            x,y=(a,z) if method=='raw' else (src[method+'_even'][sl].astype(float),src[method+'_odd'][sl].astype(float))
            check=(np.sum((x-z)**2)+np.sum((y-a)**2))/(2*a.size*den)
            val=h[(h.crop==i)&(h.method==method)].cross_mse_ratio.iloc[0]
            errors.append(abs(check-val));assert abs(check-val)<1e-12
    for n in ['train','val']:
        a=np.load(Path('C:/Users/wangz/cryoet_extension_runtime_20261003/training_data')/(n+'_data.npz'),allow_pickle=True)
        coords=a['coords'];shape=a['sample_shape']
        assert np.all(coords[:,:,1]>=0) and np.all(coords[:,:,1]+shape[1]<=768)
    for n in ['extra_sample_metrics.csv','extra_sample_picking.csv','extra_sample_picking_macro.csv','half_metrics.csv','half_fsc_curves.csv']:
        assert sha(ROOT/n)==sha(OUT/n)
    abstract=new.split(r'\section*{Abstract}')[1].split(r'\noindent')[0].strip()
    (OUT/'Abstract.txt').write_text(abstract+'\n',encoding='utf-8')
    result=dict(prior_source_unchanged=True,preserved_assets=preserved,authors_and_funding_preserved=True,bibliography_entries=len(keys),undefined_citations=[],main_figure_labels_unchanged=True,picking_rows=120,picking_formula_checks_passed=True,cross_mse_rows_recomputed=24,max_cross_mse_error=max(errors),training_validation_patches_within_y_0_768=True,half_curves_verified=json.loads((ROOT/'half_verification.json').read_text()),latex_overflow_and_reference_checks='passed')
    dump(OUT/'extension_verification.json',result)
    print('DELIVERY AUDIT PASSED',len(keys),'references',len(preserved),'preserved assets')
def render():
    import fitz
    from PIL import Image,ImageDraw
    qa=ROOT/'qa';qa.mkdir(exist_ok=True);selection={}
    for name in ['manuscript','Supplementary_methods']:
        doc=fitz.open(OUT/(name+'.pdf'));thumbs=[];selected=[]
        for i,page in enumerate(doc):
            pix=page.get_pixmap(matrix=fitz.Matrix(.3,.3));im=Image.frombytes('RGB',[pix.width,pix.height],pix.samples)
            thumb=Image.new('RGB',(195,280),'#dddddd');thumb.paste(im,(8,18));ImageDraw.Draw(thumb).text((8,2),str(i+1),fill='black');thumbs.append(thumb)
            txt=page.get_text()
            if (name=='Supplementary_methods' and i>=9) or any(t in txt for t in ['Additional acquisitions','Frozen evaluation','The original ten restoration','Data and code availability','Principal tomography methods']):
                page.get_pixmap(matrix=fitz.Matrix(1.4,1.4)).save(qa/f'{name}_{i+1}.png');selected.append(i+1)
        for j in range(0,len(thumbs),20):
            chunk=thumbs[j:j+20];sheet=Image.new('RGB',(195*5,280*((len(chunk)+4)//5)),'white')
            for k,im in enumerate(chunk):sheet.paste(im,((k%5)*195,(k//5)*280))
            sheet.save(qa/f'{name}_contact_{j//20+1}.png')
        selection[name]=dict(pages=len(doc),rendered_full_pages=selected)
    dump(qa/'rendered_pages.json',selection);print(json.dumps(selection))
if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='render':render()
    else:audit()
