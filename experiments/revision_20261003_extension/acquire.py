from pathlib import Path
import requests,json,hashlib,time,concurrent.futures,datetime,zipfile
from cryoet_data_portal import Client,Annotation
ROOT=Path(__file__).parent; DATA=ROOT/'data';DATA.mkdir(exist_ok=True)
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def get(url,dest,size=None):
    if dest.exists():
        assert size is None or dest.stat().st_size==size
        return dict(url=url,file=str(dest),bytes=dest.stat().st_size,sha256=sha(dest),cached=True)
    start=time.time();part=dest.with_suffix(dest.suffix+'.part')
    for attempt in range(3):
        try:
            with requests.get(url,stream=True,timeout=(20,90)) as r:
                r.raise_for_status();expected=int(r.headers.get('Content-Length',0));n=0;last=0
                with part.open('wb') as f:
                    for b in r.iter_content(2**20):
                        f.write(b);n+=len(b)
                        if n-last>=128*2**20:print(dest.name,round(n/1e6),'MB',flush=True);last=n
                assert (not expected or n==expected) and (size is None or n==size),(n,expected,size)
                part.replace(dest)
                result=dict(url=url,file=str(dest),bytes=n,sha256=sha(dest),completed_utc=stamp(),seconds=time.time()-start,etag=r.headers.get('ETag'),last_modified=r.headers.get('Last-Modified'))
                dest.with_suffix(dest.suffix+'.provenance.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
                print('downloaded',dest.name,round(time.time()-start,1),'s',flush=True)
                return result
        except Exception as e:
            print('retry',dest.name,attempt+1,str(e)[:180],flush=True)
            if attempt==2:raise
def tomo(row):
    t=next(t for t in row['tomograms'] if t['processing']=='filtered' and t['voxel_spacing']==10.012)
    d=DATA/row['name'];d.mkdir(exist_ok=True)
    download=get(t['https_mrc_file'],d/'raw.mrc',int(t['file_size_mrc']))
    points={};evidence=[]
    for a in Annotation.find(Client(),[Annotation.run_id==row['run_id']]):
        if not a.ground_truth_status or 'hybrid curation' not in (a.annotation_method or '').lower():continue
        for sh in a.annotation_shapes:
            if sh.shape_type!='Point':continue
            for af in sh.annotation_files:
                if af.format!='ndjson' or '/VoxelSpacing10.012/' not in af.https_path:continue
                r=requests.get(af.https_path,timeout=60);r.raise_for_status()
                raw=r.content; (d/f'annotation_{af.id}.ndjson').write_bytes(raw)
                entries=[json.loads(line) for line in r.text.splitlines() if line.strip()]
                xyz=[[e.get('location',e)[k] for k in ['x','y','z']] for e in entries]
                assert len(xyz)==a.object_count
                assert all(0<=p[0]<t['size_x'] and 0<=p[1]<t['size_y'] and 0<=p[2]<t['size_z'] for p in xyz)
                assert a.object_name not in points
                points[a.object_name]=xyz
                evidence.append(dict(annotation_id=a.id,file_id=af.id,url=af.https_path,sha256=hashlib.sha256(raw).hexdigest(),object=a.object_name,count=len(xyz),method=a.annotation_method))
    assert len(points)==6,(row['name'],points.keys())
    (d/'points.json').write_text(json.dumps(points),encoding='utf-8')
    (d/'metadata.json').write_text(json.dumps(dict(run=row['name'],run_id=row['run_id'],tomogram=t,download=download,annotations=evidence),indent=2),encoding='utf-8')
    print('ready',row['name'],sum(map(len,points.values())),'particles',flush=True)
    return row['name']
def halves():
    p=DATA/'Tomo110_tutorial.zip'
    get('https://api.figshare.com/v2/file/download/45582309',p,1947310911)
    with zipfile.ZipFile(p) as z:
        members=[a for a in z.infolist() if a.filename.startswith('tutorial_data/') and not a.is_dir()]
        print('half archive',[(a.filename,a.file_size) for a in members],flush=True)
        for a in members:
            dest=(DATA/a.filename).resolve();assert dest.is_relative_to(DATA.resolve())
            if not dest.exists():z.extract(a,DATA)
        manifest=[dict(name=a.filename,bytes=a.file_size,sha256=sha(DATA/a.filename)) for a in members]
    (ROOT/'half_archive_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    return 'halves'
if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='halves': halves()
    else:
        rows=json.loads((ROOT/'portal_inventory.json').read_text())
        todo=[r for r in rows if r['name'] not in ['TS_5_4','TS_86_3']]
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
            for r in pool.map(tomo,todo):print('DONE',r,flush=True)
        (ROOT/'extra_samples_downloaded.json').write_text(json.dumps(dict(runs=[r['name'] for r in todo],completed=stamp())))
