from run_samples import *
thresholds=freeze();d=ROOT/'data/TS_69_2/processed'
with mrcfile.open(d/'normalized.mrc') as m:v=m.data.copy()
rr=3;dist=2;t=thresholds['Beta-amylase']['raw']
z,y,x=np.mgrid[-rr:rr+1,-rr:rr+1,-rr:rr+1]
template=-((z*z+y*y+x*x)<=rr*rr).astype('float32');template-=template.mean()
resp=match_template(v,template,pad_input=True);mx=float(resp.max())
peaks=peak_local_max(resp,min_distance=dist,threshold_rel=t,exclude_border=False)
prior=np.load(d/'picking/raw_r3_d2.npz');old=prior['xyz'][prior['scores']>t*float(prior['maximum'])]
new=peaks[:,::-1].astype(float)
assert np.array_equal(old,new),(old.shape,new.shape)
dump(ROOT/'peak_pruning_verification.json',dict(cached_initial_threshold=.3,frozen_threshold=t,selected_candidates=len(new),exact_coordinates_and_order_equal=True,
     reason='Fixed-threshold evaluation needs no lower-score peaks; the initial slow all-threshold extraction was stopped, its cached candidates retained. No thresholds or scored outputs changed.'))
print('PRUNING EXACTLY MATCHES ORIGINAL COORDINATES AND ORDER',len(new))
