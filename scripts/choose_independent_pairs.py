from pathlib import Path
import csv,json,collections,functools,itertools,math,sys
import primer3
p=Path(sys.argv[1]).resolve();project=Path(sys.argv[2]).resolve();I=dict(A='A',C='C',G='G',T='T',Y='CT',R='AG',W='AT',K='GT',M='AC',S='CG',B='CGT',D='AGT',H='ACT',V='ACG',N='ACGT');pools=json.loads((p/'pools.json').read_text());pool={(r['side'],r['sequence']):r['id'] for r in json.loads((p/'oligos.json').read_text())};meta={r['accession']:r for r in csv.DictReader((project/'data/panel.csv').open())};groups={g:{a for a,r in meta.items() if g in r['lineage'].split(';')} for g in ['Refertina','Vacatina']};valid={r['accession'] for r in csv.DictReader((project/'results/baseline/by_species.csv').open()) if r['name']=='MiScle-MD' and r['evaluable']=='1'};data={}
for role in ['target','background']:
 h=collections.defaultdict(lambda:collections.defaultdict(list))
 for name,acc,ori,gt,pos,mm,last in csv.reader((p/role/'sites/background_sites.tsv').open(),delimiter='\t'):h[name,int(gt)][acc,ori].append((int(pos),int(mm),int(last)))
 data[role]=h
@functools.lru_cache(None)
def stats(F,R):
 answer={}
 for role in ['target','background']:
  for gt in ([0] if role=='target' else [0,1]):
   h=data[role];f,r=h[pool['F',F],gt],h[pool['R',R],gt];cv=set();ex=set();loose=set();l4=set()
   for acc,ori in f.keys()&r.keys():
    for fp,fm,fl in f[acc,ori]:
     for rp,rm,rl in r[acc,ori]:
      inner=rp-fp-len(F);size=rp+len(R)-fp
      if not ((inner>=20 and 70<=size<=500) if role=='target' else 0<=inner<=1000):continue
      if fm<=2 and rm<=2 and not fl and not rl:cv.add(acc)
      if fm==0 and rm==0:ex.add(acc)
      if fm<=(2 if gt else 3) and rm<=(2 if gt else 3):loose.add(acc)
      l4.add(acc)
   if role=='target':answer.update(cv=cv&valid,ex=ex&valid)
   else:answer['GT' if gt else 'WC']=cv;answer['GT2' if gt else 'WC3']=loose
   if role=='background' and not gt:answer['WC4']=l4
 return answer
kw=dict(mv_conc=50,dv_conc=1.5,dntp_conc=.6)
@functools.lru_cache(None)
def expand(q):return list(map(''.join,itertools.product(*(I[c] for c in q))))
@functools.lru_cache(None)
def single(q,c):return primer3.calc_tm(q,dna_conc=c,**kw),max(primer3.calc_hairpin(q,dna_conc=c,**kw).tm,primer3.calc_homodimer(q,dna_conc=c,**kw).tm)
@functools.lru_cache(None)
def hetero(a,b,c):return primer3.calc_heterodimer(a,b,dna_conc=c,**kw).tm
@functools.lru_cache(None)
def chemistry(F,R):
 reports=[]
 for scenario in ['50_per_variant','200_per_direction']:
  qs=[(s,50 if scenario=='50_per_variant' else 200/len(expand(q))) for q in [F,R] for s in expand(q)];tm=[[single(s,50 if scenario=='50_per_variant' else 200/len(expand(q)))[0] for s in expand(q)] for q in [F,R]];st=max(single(s,c)[1] for s,c in qs)
  for (a,ac),(b,bc) in itertools.combinations(qs,2):st=max(st,hetero(a,b,max(ac,bc)))
  delta=abs(sum(tm[0])/len(tm[0])-sum(tm[1])/len(tm[1]));mid=abs((min(tm[0])+max(tm[0])-min(tm[1])-max(tm[1]))/2)
  reports.append(dict(scenario=scenario,F_tm_min=min(tm[0]),F_tm_max=max(tm[0]),R_tm_min=min(tm[1]),R_tm_max=max(tm[1]),structure_max=st,F_R_midpoint_delta=mid))
 return all(r['structure_max']<55 and r['F_R_midpoint_delta']<=5 and min(r['F_tm_min'],r['R_tm_min'])>=45 for r in reports),reports
selected=[];evidence=[]
for mode in ['ordinary','D']:
 for g in groups:
  candidates=[]
  for F in pools['F_'+g]:
   for R in pools['R_'+g]:
    if mode=='ordinary' and (set(F+R)-set('ACGT')):continue
    if mode=='D' and not (set(F+R)-set('ACGT')):continue
    z=stats(F,R)
    # Lower strict backgrounds than the shared-F release (103 WC /111 GT).
    if len(z['WC'])>(50 if g=='Refertina' else 10) or len(z['GT'])>(120 if g=='Refertina' else 20):continue
    if len(z['WC3'])>(800 if g=='Refertina' else 100):continue
    ok,t=chemistry(F,R)
    delta_limit=6 if mode=='ordinary' and g=='Vacatina' else 5
    ok=all(z['structure_max']<55 and z['F_R_midpoint_delta']<=delta_limit and min(z['F_tm_min'],z['R_tm_min'])>=45 for z in t)
    row=dict(mode=mode,suborder=g,F=F,R=R,covered=len(z['cv']&groups[g]),exact=len(z['ex']&groups[g]),all_covered=len(z['cv']),all_exact=len(z['ex']),strict_WC=len(z['WC']),strict_GT=len(z['GT']),WC3=len(z['WC3']),WC4=len(z['WC4']),thermal_pass=ok,tm_delta_limit=delta_limit,thermal=t,degeneracy_F=len(expand(F)),degeneracy_R=len(expand(R)))
    evidence.append(row)
    if ok:candidates.append(row)
  if not candidates:raise RuntimeError('No passing candidates '+mode+' '+g)
  best=max(x['covered'] for x in candidates);short=[x for x in candidates if x['covered']>=best-2]
  x=max(short,key=lambda r:(r['exact'],r['covered'],-r['strict_WC'],-r['strict_GT'],-r['degeneracy_F']*r['degeneracy_R'],r['F'],r['R']))
  x=dict(x,name='MiScle-M'+('R' if g=='Refertina' else 'V')+('D' if mode=='D' else ''),region='M',intended_suborder=g,F_start=1561,F_end=1580,R_start=1741,R_end=1760,reference_span=200)
  selected.append(x);print('SELECT',json.dumps(x),flush=True)
# Retain M-only optimum as a comparator, then choose the R-D option with
# lower backgrounds if adding unchanged SD preserves both union metrics.
(p/'selected_M_only.json').write_text(json.dumps(selected,indent=2))
v=next(x for x in selected if x['name']=='MiScle-MVD');vz=stats(v['F'],v['R'])
srows=list(csv.DictReader((project/'by_species.csv').open()));scv={r['accession'] for r in srows if r['name']=='MiScle-SD' and r['covered']=='1'};sex={r['accession'] for r in srows if r['name']=='MiScle-SD' and r['exact']=='1'}
rc=[x for x in evidence if x['mode']=='D' and x['suborder']=='Refertina' and x['thermal_pass']];topcov=max(x['covered'] for x in rc);rc=[x for x in rc if x['covered']>=topcov-2]
def integrated(x):
 z=stats(x['F'],x['R']);return (len(z['cv']|vz['cv']|scv),len(z['ex']|vz['ex']|sex),-max(len(z['WC']|vz['WC']),len(z['GT']|vz['GT'])),x['exact'],-len(z['WC']|vz['WC'])-len(z['GT']|vz['GT']))
r=max(rc,key=integrated);previous=next(x for x in selected if x['name']=='MiScle-MRD');new=dict(r,**{k:previous[k] for k in ['name','region','intended_suborder','F_start','F_end','R_start','R_end','reference_span']},selection_context='Lower background at unchanged SD+M coverage and exact matching')
selected=[new if x['name']=='MiScle-MRD' else x for x in selected];print('INTEGRATED R-D',json.dumps(new),flush=True)
(p/'selected.json').write_text(json.dumps(selected,indent=2))
(p/'screened_pairs.json').write_text(json.dumps(evidence,indent=2))
