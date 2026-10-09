"""Official cryoCARE 0.2.2 CPU experiment, held-out spatial regions of frame-split Tomo110."""
from pathlib import Path
import os,sys,json,time,hashlib,random,traceback
os.environ['CUDA_VISIBLE_DEVICES']='-1'
os.environ['TF_NUM_INTRAOP_THREADS']='4'
os.environ['TF_NUM_INTEROP_THREADS']='2'
os.environ['OMP_NUM_THREADS']='4'
os.environ['TF_CPP_MIN_LOG_LEVEL']='2'
ROOT=Path(__file__).resolve().parent
WORK=Path('C:/Users/wangz/cryoet_extension_runtime_20261003')
WORK.mkdir(exist_ok=True)
class Tee:
    def __init__(self,stream,file):self.stream,self.file=stream,file
    def write(self,s):self.stream.write(s);self.file.write(s);self.file.flush()
    def flush(self):self.stream.flush();self.file.flush()
log=(ROOT/'cryocare_run.log').open('a',encoding='utf-8')
sys.stdout=Tee(sys.stdout,log);sys.stderr=Tee(sys.stderr,log)
import numpy as np,mrcfile,tensorflow as tf
from importlib.metadata import version
from csbdeep.models import Config
from cryocare.internals.CryoCAREDataModule import CryoCARE_DataModule
from cryocare.internals.CryoCARE import CryoCARE,predict_direct
def dump(p,x):p.write_text(json.dumps(x,indent=2),encoding='utf-8')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for x in iter(lambda:f.read(2**20),b''):h.update(x)
    return h.hexdigest()
def main():
    seed=20261003;random.seed(seed);np.random.seed(seed);tf.random.set_seed(seed)
    spec=dict(method='official cryoCARE',package_version=version('cryocare'),tensorflow=tf.__version__,seed=seed,device='CPU',
              training_crop_zyx=[[0,208],[0,768],[0,544]],validation_fraction=.15,patch_shape=[64,64,64],n_samples=200,n_normalization_samples=32,
              epochs=20,steps_per_epoch=20,batch_size=2,unet_depth=3,unet_first=16,learning_rate=.0004,
              training_budget='400 gradient updates; bounded local demonstration, not a converged or tuned method benchmark',
              evaluation='six geometric 96-cube contexts, central 64-cube scoring; same model applied to each held-out half separately, then pair average as secondary output',
              centers_zyx=[[104,y,x] for y in [960,1088,1216] for x in [144,400]],test_y=[896,1280],guard_y=[768,896],
              independence='frame-disjoint inputs; shared acquisition, alignment and reconstruction, spatial crops are technical subsamples not biological replicates',
              voxel_spacing='MRC header is zero; report frequency in cycles/voxel, no Angstrom resolution inferred',
              wrapper='Use unchanged official DataModule, CryoCARE.train and predict_direct; custom held-out crop orchestration. Packaged pretrained tutorial checkpoint is NOT loaded.',
              script_sha256=sha(Path(__file__)))
    dump(ROOT/'cryocare_protocol.json',spec)
    source=ROOT/'data/tutorial_data';paths={}
    for half in ['even','odd']:
        p=WORK/f'train_{half}.mrc';paths[half]=str(p)
        with mrcfile.mmap(source/f'tomo_{half}_frames.rec',permissive=True) as m:
            assert m.data.shape==(209,1280,550)
            if not p.exists():
                with mrcfile.new(p) as n:n.set_data(np.ascontiguousarray(m.data[:208,:768,:544],dtype='float32'))
    modeldir=WORK/'official_cryocare'
    if (ROOT/'cryocare_training_complete.json').exists():
        info=json.loads((ROOT/'cryocare_training_complete.json').read_text());mean,std=info['mean'],info['std']
        model=CryoCARE(None,'official_cryocare',basedir=str(WORK))
    else:
        if modeldir.exists():raise RuntimeError('Partial training exists; preserve it and inspect before restarting.')
        started=time.time();dm=CryoCARE_DataModule()
        dm.setup([paths['odd']],[paths['even']],n_samples_per_tomo=200,validation_fraction=.15,sample_shape=(64,64,64),tilt_axis='Y',n_normalization_samples=32)
        dp=WORK/'training_data';dp.mkdir(exist_ok=True);dm.save(str(dp))
        mean,std=float(dm.train_dataset.mean),float(dm.train_dataset.std)
        conf=Config(axes='ZYXC',train_loss='mse',train_epochs=20,train_steps_per_epoch=20,train_batch_size=2,unet_kern_size=3,unet_n_depth=3,unet_n_first=16,train_tensorboard=False,train_learning_rate=.0004)
        model=CryoCARE(conf,'official_cryocare',basedir=str(WORK))
        hist=model.train(dm.get_train_dataset(),dm.get_val_dataset())
        dump(ROOT/'cryocare_history.json',{k:[float(x) for x in v] for k,v in hist.history.items()})
        dump(modeldir/'norm.json',dict(mean=mean,std=std))
        dm.close()
        info=dict(training_seconds=time.time()-started,mean=mean,std=std,epochs_completed=len(hist.history['loss']),weights_sha256=sha(modeldir/'weights_best.h5'),model_dir=str(modeldir),script_sha256=sha(Path(__file__)))
        assert info['epochs_completed']==20
        dump(ROOT/'cryocare_training_complete.json',info)
    evaldir=ROOT/'half_evaluation';evaldir.mkdir(exist_ok=True)
    manifest=[]
    even=mrcfile.mmap(source/'tomo_even_frames.rec',permissive=True)
    odd=mrcfile.mmap(source/'tomo_odd_frames.rec',permissive=True)
    for i,(z,y,x) in enumerate(spec['centers_zyx']):
        slices=tuple(slice(c-48,c+48) for c in [z,y,x]);start=time.time()
        a=np.array(even.data[slices],dtype='float32');b=np.array(odd.data[slices],dtype='float32')
        assert a.shape==b.shape==(96,96,96)
        # Unchanged official function; each pass sees only one noisy half.
        pa=predict_direct(model.keras_model,a[...,None],mean,std,'ZYXC',verbose=0)[...,0]
        pb=predict_direct(model.keras_model,b[...,None],mean,std,'ZYXC',verbose=0)[...,0]
        assert np.isfinite(pa).all() and np.isfinite(pb).all()
        p=evaldir/f'crop_{i}.npz'
        np.savez_compressed(p,even=a,odd=b,cryocare_even=pa,cryocare_odd=pb)
        manifest.append(dict(crop=i,center_zyx=[z,y,x],file=str(p),sha256=sha(p),inference_seconds=time.time()-start))
        print('[cryocare] held-out crop',i,'complete',flush=True)
    even.close();odd.close()
    dump(ROOT/'cryocare_complete.json',dict(training=info,crops=manifest,completed=time.strftime('%Y-%m-%dT%H:%M:%S')))
    print('OFFICIAL CRYOCARE TRAINING AND HELD-OUT INFERENCE COMPLETE',flush=True)
if __name__=='__main__':
    try:main()
    except Exception:
        traceback.print_exc();raise
