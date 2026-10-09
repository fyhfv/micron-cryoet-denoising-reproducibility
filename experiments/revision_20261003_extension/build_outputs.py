"""Export evidence-backed supplemental tables/plots, no original-file changes."""
from pathlib import Path
import json,hashlib,shutil,re,difflib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parents[1]
OUT=PROJECT/'正式投稿/Micron_submission_20261003_Peng_extension'
OLD=PROJECT/'正式投稿/Micron_submission_20261003_Peng_literature'
RUNS=json.loads((ROOT/'protocol.json').read_text())['new_CZII_runs']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,obj):p.write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')
def table(file,title,cols,head,rows,note=''):
    s='\\begin{table}[htbp]\\centering\\small\n\\caption*{'+title+'}\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{'+cols+'}\\toprule\n'+head+'\\\\\\midrule\n'
    s+='\n'.join(' & '.join(row)+r'\\' for row in rows)
    s+='\n\\bottomrule\\end{tabular}\n'
    if note:s+='\\par\\smallskip\\begin{minipage}{.98\\textwidth}\\footnotesize '+note+'\\end{minipage}\n'
    (OUT/file).write_text(s+'\\end{table}\n',encoding='utf-8')
def main():
    assert (ROOT/'samples_complete.json').exists() and (ROOT/'cryocare_complete.json').exists()
    m=pd.read_csv(ROOT/'extra_sample_metrics.csv');p=pd.read_csv(ROOT/'extra_sample_picking.csv');macro=pd.read_csv(ROOT/'extra_sample_picking_macro.csv');h=pd.read_csv(ROOT/'half_metrics.csv');curves=pd.read_csv(ROOT/'half_fsc_curves.csv')
    assert len(m)==15 and len(p)==120 and len(h)==24 and len(curves)==1536
    assert (p.tp==p.independent_tp).all() and (p.tp+p.fp==p.n_pred).all() and (p.tp+p.fn==p.n_gt).all()
    names={'raw':'Raw','gaussian':'Gaussian','topaz':'Topaz','cryocare':'cryoCARE'}
    rows=[]
    for run in RUNS:
        for method in ['raw','gaussian','topaz']:
            r=m[(m.run==run)&(m.method==method)].iloc[0]
            rows.append([run.replace('_',r'\_'),names[method],str(int(r.n_ribosomes)),f'{r.SNR_region:.3f}',f'{r.CNR:.3f}',f'{r.sharpness:.2f}'])
    table('extension_table_metrics.tex','Table S6. Additional full-volume region-metric measurements.','llrrrr','Run & Input & Ribosomes & SNR & CNR & Sobel',rows,'SNR and CNR are linear ratios, not dB. Sobel magnitude is computed after global z-scoring and is sensitive to noise as well as blur. No molecular-resolution claim is derived from these metrics.')
    rows=[];counts=[]
    for run in RUNS:
        n=int(p[(p.run==run)&(p.method=='raw')&(p.protocol=='per_input')].n_gt.sum());counts.append(n)
        values=[float(macro[(macro.run==run)&(macro.method==method)&(macro.protocol==protocol)].macro_F1.iloc[0]) for method,protocol in [('raw','per_input'),('topaz','raw_fixed'),('topaz','per_input')]]
        rows.append([run.replace('_',r'\_'),str(n)]+[f'{v:.4f}' for v in values])
    means=[float(macro[(macro.method==method)&(macro.protocol==protocol)].macro_F1.mean()) for method,protocol in [('raw','per_input'),('topaz','raw_fixed'),('topaz','per_input')]]
    rows.append(['Mean','---']+[f'{v:.4f}' for v in means])
    table('extension_table_picking.tex','Table S7. Six-class macro-F1 on five additional acquisition runs.','lrrrr','Run & Particles & Raw & Topaz/fixed & Topaz/own',rows,'All thresholds are calibrated on TS\\_5\\_4. Macro-F1 weights the six classes equally within each run; the last row weights the five runs equally. Full class-level TP/FP/FN and thresholds are retained in the CSV, including zero detections. Particle totals are not independent replicate counts.')
    rows=[]
    for method in ['raw','gaussian','topaz','cryocare']:
        d=h[h.method==method]
        rows.append([names[method]]+[f'${d[k].mean():.4f}\\pm{d[k].std():.4f}$' for k in ['cross_mse_ratio','fsc_band_mean','hann_fsc_band_mean']])
    table('extension_table_halves.tex','Table S8. Tomo110: descriptive mean $\\pm$ SD across six spatial crops.','lccc','Input & Cross-MSE ratio & FSC, no taper & FSC, Hann',rows,'FSC columns average shells in 0.05--0.20 cycles/voxel. The six crops belong to one tomogram; SD describes spatial variation, not uncertainty across biological samples. Smaller cross-MSE is better prediction of the opposite noisy half, not an accuracy estimate against truth.')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    colors={'raw':'#4b5563','gaussian':'#b65f17','topaz':'#196da1','cryocare':'#19816b'}
    fig,axs=plt.subplots(1,2,figsize=(10,3.8),layout='constrained');x=np.arange(5)
    for method in ['raw','gaussian','topaz']:
        vals=m[m.method==method].set_index('run').loc[RUNS].SNR_region
        axs[0].plot(x,vals,'o-',label=names[method],color=colors[method],ms=4)
    for method,rule,label,color in [('raw','per_input','Raw',colors['raw']),('topaz','raw_fixed','Topaz/fixed','#7bb8d6'),('topaz','per_input','Topaz/own',colors['topaz'])]:
        vals=macro[(macro.method==method)&(macro.protocol==rule)].set_index('run').loc[RUNS].macro_F1
        axs[1].plot(x,vals,'o-',label=label,color=color,ms=4)
    for ax,title,ylabel in zip(axs,['(a) Ribosome-region contrast','(b) Frozen six-class detection'],['Region SNR (linear ratio)','Macro-F1']):
        ax.set_xticks(x,RUNS,rotation=25,ha='right');ax.set_ylabel(ylabel);ax.set_title(title,loc='left',fontsize=10);ax.set_ylim(bottom=0);ax.legend(frameon=False,fontsize=8);ax.grid(axis='y',alpha=.2)
    fig.savefig(OUT/'fig_extension_samples.pdf');fig.savefig(OUT/'fig_extension_samples.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(11,3.7),layout='constrained')
    methods=['raw','gaussian','topaz','cryocare']
    for j,method in enumerate(methods):
        d=h[h.method==method];vals=d.cross_mse_ratio.to_numpy()
        axs[0].scatter(j+np.linspace(-.12,.12,6),vals,c=colors[method],s=22,alpha=.8)
        axs[0].plot([j-.25,j+.25],[vals.mean()]*2,c=colors[method],lw=2)
        for ax,window in zip(axs[1:],['none','hann']):
            dd=curves[(curves.method==method)&(curves.window==window)]
            for _,g in dd.groupby('crop'):ax.plot(g.frequency,g.fsc,c=colors[method],alpha=.15,lw=.6)
            avg=dd.groupby('frequency').fsc.mean();ax.plot(avg.index,avg.values,c=colors[method],label=names[method],lw=1.6)
    axs[0].set_xticks(range(4),[names[m] for m in methods],rotation=28,ha='right');axs[0].set_ylabel('Symmetric cross-half MSE / raw MSE');axs[0].set_ylim(.55,1.04);axs[0].set_title('(a) Noisy-target prediction',loc='left',fontsize=10)
    for ax,title in zip(axs[1:],['(b) Half-FSC: no taper','(c) Half-FSC: Hann taper']):
        ax.axhline(.143,c='#888888',ls=':',lw=1);ax.set_xlim(0,.5);ax.set_ylim(-.08,1.02);ax.set_xlabel('Spatial frequency (cycles/voxel)');ax.set_ylabel('FSC');ax.set_title(title,loc='left',fontsize=10)
    axs[2].legend(frameon=False,fontsize=8);fig.savefig(OUT/'fig_extension_halves.pdf');fig.savefig(OUT/'fig_extension_halves.png',dpi=160);plt.close(fig)
    ledger=dict(new_runs=RUNS,particles=sum(counts),ribosomes=int(m[m.method=='raw'].n_ribosomes.sum()),picking_rows=len(p),metric_rows=len(m),half_rows=len(h),fsc_curves=48,macro_F1=dict(zip(['raw','topaz_fixed','topaz_own'],means)),half_summary={k:dict(mean_ratio=float(g.cross_mse_ratio.mean()),sd_ratio=float(g.cross_mse_ratio.std())) for k,g in h.groupby('method')},source_manuscript_sha256=sha(OLD/'manuscript.tex'),source_supplement_sha256=sha(OLD/'Supplementary_methods.tex'))
    dump(ROOT/'extension_summary.json',ledger);dump(OUT/'extension_summary.json',ledger)
    for filename in ['extra_sample_metrics.csv','extra_sample_picking.csv','extra_sample_picking_macro.csv','half_metrics.csv','half_fsc_curves.csv','half_verification.json','frozen_thresholds.json','cryocare_protocol.json','half_boundary_amendment.json']:
        shutil.copy2(ROOT/filename,OUT/filename)
    print(json.dumps(ledger,indent=2))
if __name__=='__main__':main()
