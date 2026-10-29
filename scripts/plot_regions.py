from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
p=Path(__file__).resolve().parents[1]
with (p/'primer_pairs.csv').open() as f:rows=[r for r in csv.DictReader(f) if r['name'] in ['MiScle-S','MiScle-MR']]
with (p/'barcode_resolution.csv').open() as f:lengths={r['name']:r for r in csv.DictReader(f)}
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'path','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.spines.left':False})
fig,ax=plt.subplots(figsize=(12,4.5));fig.subplots_adjust(left=.16,right=.95,top=.72,bottom=.28)
fig.text(.05,.92,'Scleractinia 16S · reference design windows',fontsize=21,weight='bold',color='#174c3b');fig.text(.05,.83,'Ordinary and D primers share endpoints within each window',fontsize=11,color='#50655e')
for i,(r,col) in enumerate(zip(rows,['#087e8b','#c28b29'])):
 fs,fe,rs,re=(int(r[k]) for k in ('F_start','F_end','R_start','R_end'));ax.plot([fe,rs],[i,i],lw=3,alpha=.35,color=col)
 for a,b in ((fs,fe),(re,rs)):ax.annotate('',xy=(b,i),xytext=(a,i),arrowprops=dict(arrowstyle='-|>',lw=4,color=col,shrinkA=0,shrinkB=0))
 ax.text((fe+rs)/2,i-.19,f"Reference span: {r['reference_span']} bp",ha='center',fontsize=11);ax.text(fs,i+.25,f'F {fs}–{fe}',color=col,fontsize=10);ax.text(re,i+.25,f'R {rs}–{re}',color=col,ha='right',fontsize=10)
ax.set_yticks(range(2),['Window '+r['region'] for r in rows],fontsize=13);ax.set_ylim(1.6,-.6);ax.set_xlim(1520,1885);ax.grid(axis='x',alpha=.15);ax.tick_params(axis='y',length=0);ax.set_xlabel('NC_003522.1 · ungapped mitochondrial 16S positions (1-based)',labelpad=10)
fig.text(.05,.055,'Predicted products: S, 147–209 bp; MR/MRD, 186–298 bp; MV/MVD, 99–186 bp.\nM-R and M-V have independent F/R sequences at identical reference endpoints.',fontsize=10,color='#50655e',linespacing=1.6)
for ext in ('svg','png'):fig.savefig(p/f'figures/reference_windows.{ext}',dpi=180,facecolor='white')
