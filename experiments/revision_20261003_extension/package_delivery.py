from pathlib import Path
import hashlib,json,zipfile,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;PROJECT=BASE.parent
OUT=PROJECT/'正式投稿/Micron_submission_20261003_Peng_extension'
WORK=Path('C:/Users/wangz/cryoet_extension_runtime_20261003')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False),encoding='utf-8')
def main():
    initial=json.loads((ROOT/'preprocess_config.json').read_text(encoding='utf-8'))['code_sha256']
    assert sha(ROOT/'run_samples_initial_executed.py')==initial,'Initial executed source hash not recovered'
    for key,py in [('main',BASE/'.venv/Scripts/python.exe'),('cryocare',ROOT/'.venv_cryocare/Scripts/python.exe')]:
        cmd='import importlib.metadata as m,json,sys;print(json.dumps({"python":sys.version,"packages":{d.metadata["Name"]:d.version for d in m.distributions() if d.metadata["Name"]}},indent=2))'
        r=subprocess.run([str(py),'-X','utf8','-c',cmd],capture_output=True,text=True,encoding='utf-8',check=True)
        (ROOT/f'environment_{key}.json').write_text(r.stdout,encoding='utf-8')
    for source,target in [(WORK/'official_cryocare',ROOT/'cryocare_model'),(WORK/'training_data',ROOT/'cryocare_training_data')]:
        target.mkdir(exist_ok=True)
        for p in source.iterdir():
            if p.is_file():shutil.copy2(p,target/p.name)
    claims=[dict(id='M1_half_data',status='partial',completed='Real frame-disjoint Tomo110 reconstructions and separately processed held-out regions',remaining='Shared reconstruction effects; not gold-standard experimental resolution'),dict(id='M2_official_cryocare',status='addressed_for_bounded_demonstration',completed='Official package trained and scored, saved weights/configuration/history',remaining='Single seed, 400 updates, no specialised supervised method or matched-resource ranking'),dict(id='M4_more_runs',status='partial',completed='All five previously unused DS10440 runs, frozen six-class evaluation, 120 checked count rows',remaining='Same dataset, no demonstrated biological or inter-laboratory replication or new learned-picker calibration'),dict(id='M3_restoration',status='partial',completed='Experimental half-data acquired',remaining='No new missing-wedge method, repeated DDW training, IsoNet or ICECREAM'),dict(id='M5_STA',status='pending',completed='Half-FSC boundary sensitivity retained as separate diagnostic',remaining='No gold-standard STA or before-splitting training-leakage control')]
    dump(OUT/'reviewer_response_matrix.json',claims)
    provenance=dict(executed_initial_sample_script_sha256=initial,current_sample_script_sha256=sha(ROOT/'run_samples.py'),executed_cryocare_script_sha256=sha(ROOT/'run_cryocare.py'),cryocare_saved_script_hash=json.loads((ROOT/'cryocare_training_complete.json').read_text())['script_sha256'],metrics_module_sha256=sha(BASE/'code/metrics.py'),source_hashes_preserved=True,all_new_computations_complete=True,pdf_visual_review='All page contact sheets and changed/new pages inspected; removed a nearly empty supplement page; final 49-page main text and 13-page supplement.',final_pdf_pages=dict(manuscript=49,Supplementary_methods=13))
    assert provenance['executed_cryocare_script_sha256']==provenance['cryocare_saved_script_hash']
    dump(OUT/'extension_provenance.json',provenance)
    shutil.copy2(ROOT/'README_reproducibility.txt',OUT/'README_reproducibility.txt')
    files=[p for p in ROOT.iterdir() if p.is_file() and p.suffix in ['.py','.json','.csv','.txt','.log']]
    for dirname in ['half_evaluation','cryocare_model','cryocare_training_data']:
        files.extend(p for p in (ROOT/dirname).rglob('*') if p.is_file())
    for d in (ROOT/'data').iterdir():
        if not d.is_dir():
            if d.suffix=='.json':files.append(d)
            continue
        if d.name=='tutorial_data':files.append(d/'README.txt')
        else:files.extend(p for p in d.iterdir() if p.is_file() and p.suffix in ['.json','.ndjson'])
    files += [BASE/'code/metrics.py',BASE/'data/czii_TS_5_4_meta.json',BASE/'revision_20261003_peng_followup/picking_threshold_sweep.csv']
    manifest=[dict(path=str(p.relative_to(PROJECT)).replace('\\','/'),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(set(files))]
    dump(ROOT/'reproducibility_manifest.json',manifest)
    zipname=OUT/'Peng_extension_reproducibility.zip'
    with zipfile.ZipFile(zipname,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for p in sorted(set(files)):z.write(p,p.relative_to(PROJECT))
        z.write(ROOT/'reproducibility_manifest.json',(ROOT/'reproducibility_manifest.json').relative_to(PROJECT))
    with zipfile.ZipFile(zipname) as z:assert z.testzip() is None
    with zipfile.ZipFile(OUT/'Micron_LaTeX_source.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in OUT.iterdir():
            if p.suffix=='.tex' or (p.name.startswith('fig_') and p.suffix=='.pdf'):z.write(p,p.name)
    manifest=[dict(name=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in OUT.iterdir() if p.is_file() and p.suffix in ['.pdf','.tex','.csv','.zip','.json','.html','.txt']]
    dump(OUT/'delivery_manifest.json',manifest)
    print('PACKAGED',zipname.name,round(zipname.stat().st_size/1e6,1),'MB; CRC verified')
if __name__=='__main__':main()
