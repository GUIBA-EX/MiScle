from pathlib import Path
import csv,json,sys,itertools,collections,functools
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.IUPACData import ambiguous_dna_values as I
p=Path(__file__).parent;root=p.parents[1];sys.path.insert(0,str(root/'work/milily'));from common import deg,expand,thermal,oligo_thermal
panel=list(csv.DictReader((p/'curated/panel_complement.csv').open()));al={r.id:str(r.seq).upper() for r in SeqIO.parse(p/'aligned.fasta','fasta')};ref='NC_003522.1';cols=[i for i,c in enumerate(al[ref]) if c!='-'];names=[r['accession'] for r in panel];arr=np.array([[dict(A=1,C=2,G=4,T=8).get(c,0) for c in al[acc]] for acc in names],dtype=np.uint8);N=len(panel);rev={frozenset(v):k for k,v in I.items() if k!='X'};fambits={f:sum(1<<i for i,r in enumerate(panel) if r['family']==f) for f in set(r['family'] for r in panel)}
def terminal(q):return any(set(I[c])<=set('GC') for c in q[-3:]) and not any(set.intersection(*(set(I[c]) for c in q[j:j+5])) for j in range(len(q)-4))
def merge(q,t):return ''.join(rev[frozenset(set(I[a])|set(I[b]))] for a,b in zip(q,t))
def bits(values):return sum(1<<int(i) for i in np.flatnonzero(values))
@functools.lru_cache(None)
def diversity(cv,valid):
 fs=[(cv&b).bit_count()/(valid&b).bit_count() for b in fambits.values() if valid&b];return sum(fs)/len(fs) if fs else 0
pool=[];diagnostic=[]
for side in ('F','R'):
 for a in range(0,len(cols)-24,5):
  for L in (20,22,24):
   ix=cols[a:a+L];mat=arr[:,ix];mat=mat[:,::-1] if side=='R' else mat
   if side=='R':mat=np.array([0,8,4,0,2,0,0,0,1],dtype=np.uint8)[mat]
   complete=(mat!=0).all(axis=1)
   if complete.sum()<70:continue
   strings=[]
   for i in np.flatnonzero(complete):
    q=''.join('ACGT'[[1,2,4,8].index(int(v))] for v in mat[i]);strings.append(q)
   counts=collections.Counter(strings);observed=[q for q,count in counts.most_common(18)];seed=[q for q in observed[:5] if terminal(q)]
   consensus=''.join('ACGT'[np.argmax([(mat[:,j]==v).sum() for v in (1,2,4,8)])] for j in range(L))
   if terminal(consensus):seed.append(consensus)
   if not seed:continue
   valid=bits(complete);cache={}
   def evaluate(q):
    if q not in cache:
     masks=np.array([sum(dict(A=1,C=2,G=4,T=8)[c] for c in I[ch]) for ch in q],dtype=np.uint8);bad=(mat&masks)==0;mm=bad.sum(axis=1);cv=bits(complete&(mm<=2)&(~bad[:,-3:].any(axis=1)));ex=bits(complete&(mm==0));cache[q]=(cv,ex)
    return cache[q]
   cand=set(seed)
   for q in seed:
    for _ in range(3):
     options=[merge(q,t) for t in observed if deg(merge(q,t))<=8 and terminal(merge(q,t))];options=sorted(set(options))
     if not options:break
     q=max(options,key=lambda z:(evaluate(z)[0].bit_count(),evaluate(z)[1].bit_count(),-deg(z)));cand.add(q)
   rated=[]
   for q in sorted(cand):
    cv,ex=evaluate(q);rated.append(dict(side=side,start=a+1,end=a+L,sequence=q,degeneracy=deg(q),cv=cv,ex=ex,valid=valid,score=cv.bit_count()+N*.25*diversity(cv,valid)))
   keep=sorted(rated,key=lambda x:(x['score'],x['ex'].bit_count()),reverse=True)[:2]+sorted(rated,key=lambda x:(x['ex'].bit_count(),x['score']),reverse=True)[:1]+sorted([x for x in rated if x['degeneracy']==1],key=lambda x:x['score'],reverse=True)[:1];pool+=list({x['sequence']:x for x in keep}.values())
 print(side,'windows evaluated; pool',len(pool),flush=True)
# Keep positional diversity so highly conserved sites in one area cannot crowd out pairs.
keep=[]
for side in ('F','R'):
 bins=collections.defaultdict(list)
 for x in pool:
  if x['side']==side:bins[(x['start']-1)//20].append(x)
 for rr in bins.values():keep+=sorted(rr,key=lambda x:(x['score'],x['ex'].bit_count()),reverse=True)[:3]+sorted(rr,key=lambda x:(x['ex'].bit_count(),x['score']),reverse=True)[:2]
keep=list({(x['side'],x['start'],x['end'],x['sequence']):x for x in keep}.values());sidepool={side:[x for x in keep if x['side']==side] for side in ('F','R')};pairs=[]
for f,r in itertools.product(sidepool['F'],sidepool['R']):
 span=r['end']-f['start']+1
 if not (150<=span<=500 and r['start']-f['end']>=50):continue
 valid=f['valid']&r['valid'];n=valid.bit_count()
 if n<70:continue
 cv=f['cv']&r['cv'];ex=f['ex']&r['ex'];reg='S' if span<=250 else 'M' if span<=350 else 'L';score=cv.bit_count()+N*.25*diversity(cv,valid)
 pairs.append(dict(region=reg,F=f['sequence'],R=r['sequence'],F_start=f['start'],F_end=f['end'],R_start=r['start'],R_end=r['end'],reference_span=span,degeneracy_F=f['degeneracy'],degeneracy_R=r['degeneracy'],pred_covered=cv.bit_count(),pred_exact=ex.bit_count(),pred_evaluable=n,score=score))
(p/'complement_all_pairs.json').write_text(json.dumps(pairs));(p/'complement_all_oligos.json').write_text(json.dumps(keep))
selected=[];therm={}
for reg in ('S','M','L'):
 rr=[x for x in pairs if x['region']==reg];top=sorted(rr,key=lambda x:(x['score'],x['pred_exact']),reverse=True)[:140]+sorted(rr,key=lambda x:(x['pred_exact'],x['score']),reverse=True)[:70];top=list({(x['F'],x['R'],x['F_start'],x['R_start']):x for x in top}.values());passed=[]
 for x in top:
  f,r=x['F'],x['R']
  for q in (f,r):
   if q not in therm:therm[q]=oligo_thermal(q)
  if not all(max(a[2],b[2])<55 and abs((a[0]+a[1]-b[0]-b[1])/2)<=5 and min(a[0],b[0])>=45 for a,b in zip(therm[f],therm[r])):continue
  ok,t=thermal(f,r)
  if ok:passed.append(dict(x,thermal=t))
 selected+=sorted(passed,key=lambda x:(x['score'],x['pred_exact']),reverse=True)[:30]+sorted(passed,key=lambda x:(x['pred_exact'],x['score']),reverse=True)[:15];print(reg,'pairs',len(rr),'checked',len(top),'physical_pass',len(passed),flush=True)
selected=list({(x['F'],x['R'],x['F_start'],x['R_start']):x for x in selected}.values());(p/'complement_candidates.json').write_text(json.dumps(selected,indent=2));(p/'complement_design_manifest.json').write_text(json.dumps(dict(reference=ref,representatives=N,pool=len(pool),retained_oligos=len(keep),pairs=len(pairs),chemical_pass_retained=len(selected),site_min_evaluable=130,max_degeneracy_per_oligo=8,lengths=[20,22,24],reference_start_step=5),indent=2));print('ready',len(selected),flush=True)
