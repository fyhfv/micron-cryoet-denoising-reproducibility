from pathlib import Path
import json,requests,datetime
from cryoet_data_portal import Client,Run,Tomogram,Annotation
ROOT=Path(__file__).parent;ROOT.mkdir(exist_ok=True)
c=Client(); records=[]
for r in Run.find(c,[Run.dataset_id==10440]):
    ts=[]
    for t in Tomogram.find(c,[Tomogram.run_id==r.id]):
        ts.append({k:getattr(t,k,None) for k in ['id','name','processing','voxel_spacing','size_x','size_y','size_z','https_mrc_file','file_size_mrc','tomogram_version','ctf_corrected','reconstruction_software']})
    records.append(dict(run_id=r.id,name=r.name,tomograms=ts))
(ROOT/'portal_inventory.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
print(json.dumps(records,indent=2),flush=True)
u='https://api.figshare.com/v2/file/download/45582309'
r=requests.head(u,allow_redirects=True,timeout=30)
print('HALF_PAIR_ARCHIVE',r.status_code,dict(r.headers),flush=True)
(ROOT/'discovery_status.json').write_text(json.dumps(dict(access_time=datetime.datetime.now(datetime.timezone.utc).isoformat(),half_url=u,half_status=r.status_code,half_headers=dict(r.headers))),encoding='utf-8')
