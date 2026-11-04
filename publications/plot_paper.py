from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path(__file__).resolve().parents[1];out=r/'publications/figures'
rows=list(csv.DictReader((r/'by_suborder.csv').open()));allrows=list(csv.DictReader((r/'comparison.csv').open()))
names=['MiScle-S','MiScle-SD','MiScle-MR','MiScle-MRD','MiScle-MV','MiScle-MVD'];chosen=[]
for n in names:
 chosen.append(next(x for x in (allrows if n in names[:2] else rows) if x['name']==n and (n in names[:2] or x['taxon']==('Refertina' if n in names[2:4] else 'Vacatina'))))
fig,ax=plt.subplots(figsize=(8,4))
for offset,key,col in [(-.18,'covered','#268797'),(.18,'exact','#e9aa38')]:
 vals=[100*int(x[key])/int(x['evaluable']) for x in chosen];bars=ax.bar([i+offset for i in range(6)],vals,.36,label=key.capitalize(),color=col)
 ax.bar_label(bars,labels=[f"{x[key]}/{x['evaluable']}" for x in chosen],fontsize=8,padding=3)
ax.set_xticks(range(6),[n.replace('MiScle-','') for n in names]);ax.set_ylabel('Evaluable species (%)');ax.set_ylim(0,115);ax.legend(frameon=False);ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(out/'performance.svg');plt.close(fig)
fig,ax=plt.subplots(figsize=(8,3.8));labels=['Ordinary WC','Ordinary GT','D WC','D GT']
mix=list(csv.DictReader((r/'mixture_summary.csv').open()))
for off,vals,label,col in [(-.18,[int(x['separate_pair_union']) for x in mix],'Separate intended pairs','#268797'),(.18,[int(x['mixed_all_pair_union']) for x in mix],'All potential same-tube pairs','#e9aa38')]:
 b=ax.bar([i+off for i in range(4)],vals,.36,label=label,color=col);ax.bar_label(b,padding=3,fontsize=9)
ax.set_xticks(range(4),labels);ax.set_ylabel('Background records');ax.set_ylim(0,1300);ax.legend(frameon=False,fontsize=9);ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(out/'background.svg')
