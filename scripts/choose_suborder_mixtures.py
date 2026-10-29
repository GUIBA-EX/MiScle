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
def chemistry(F,R,V):
 reports=[]
 for scenario in ['50_per_variant','200_F_100_each_R']:
  qs=[(s,50 if scenario=='50_per_variant' else (200 if i==0 else 100)/len(expand(q))) for i,q in enumerate([F,R,V]) for s in expand(q)];tm=[[single(s,50 if scenario=='50_per_variant' else (200 if i==0 else 100)/len(expand(q)))[0] for s in expand(q)] for i,q in enumerate([F,R,V])];st=max(single(s,c)[1] for s,c in qs)
  for (a,ac),(b,bc) in itertools.combinations(qs,2):st=max(st,hetero(a,b,max(ac,bc)))
  mids=[(min(t)+max(t))/2 for t in tm];reports.append(dict(scenario=scenario,F_tm=[min(tm[0]),max(tm[0])],RR_tm=[min(tm[1]),max(tm[1])],RV_tm=[min(tm[2]),max(tm[2])],max_structure_tm=st,max_F_R_midpoint_delta=max(abs(mids[0]-m) for m in mids[1:])))
 return all(r['max_structure_tm']<55 and r['max_F_R_midpoint_delta']<=5 and min(r['F_tm']+r['RR_tm']+r['RV_tm'])>=45 for r in reports),reports
allresults=[];choices=[]
for mode in ['ordinary','D']:
 fs=[q for q in pools['F_all'] if mode=='D' or set(q)<=set('ACGT')];rs={g:[q for q in pools['R_'+g] if mode=='D' or set(q)<=set('ACGT')] for g in groups};ranked=[]
 for F in fs:
  options={}
  for g in groups:
   rr=[]
   for R in rs[g]:
    z=stats(F,R)
    if len(z['WC'])>150 or len(z['GT'])>150 or len(z['WC3'])>1200:continue
    if len(z['cv']&groups[g])/len(valid&groups[g])<0.65:continue
    rr.append((len(z['cv']&groups[g]),len(z['ex']&groups[g]),len(z['cv']),len(z['ex']),R,z))
   options[g]=sorted(rr,key=lambda r:r[:4],reverse=True)[:8]
  for a,b in itertools.product(options['Refertina'],options['Vacatina']):
   z,w=a[-1],b[-1];cv=z['cv']|w['cv'];ex=z['ex']|w['ex'];wc=z['WC']|w['WC'];gt=z['GT']|w['GT'];balanced=min(len(cv&groups[g])/len(valid&groups[g]) for g in groups);score=(len(cv),len(ex),balanced,-max(len(wc),len(gt)),-sum(len(I[c])>1 for q in [F,a[4],b[4]] for c in q))
   ranked.append((score,F,a[4],b[4],z,w))
 best_coverage=max(r[0][0] for r in ranked)
 ranked=[r for r in ranked if r[0][0]>=best_coverage-2]
 ranked.sort(key=lambda r:(r[0][1],r[0][0],r[0][2:]),reverse=True)
 print(mode,'ranked',len(ranked),flush=True)
 for score,F,R,V,z,w in ranked:
  ok,t=chemistry(F,R,V)
  record=dict(mode=mode,F=F,RR=R,RV=V,covered=score[0],exact=score[1],balanced_coverage=score[2],strict_WC=len(z['WC']|w['WC']),strict_GT=len(z['GT']|w['GT']),thermal_pass=ok,thermal=t)
  allresults.append(record)
  if ok:
   choices.append(record);print('SELECT',json.dumps(record),flush=True);break
 else:raise RuntimeError('no chemical pass '+mode)
(p/'screened_combinations.json').write_text(json.dumps(allresults,indent=2));(p/'selected_mixtures.json').write_text(json.dumps(choices,indent=2));selected=[]
for x in choices:
 for clade,side in [('Refertina','RR'),('Vacatina','RV')]:selected.append(dict(name='MiScle-M'+('R' if clade=='Refertina' else 'V')+('D' if x['mode']=='D' else ''),region='M',intended_suborder=clade,F=x['F'],R=x[side],F_start=1561,F_end=1580,R_start=1741,R_end=1760,reference_span=200,degeneracy_F=len(expand(x['F'])),degeneracy_R=len(expand(x[side])),mixture_thermal=x['thermal']))
(p/'selected.json').write_text(json.dumps(selected,indent=2));verify=p/'verify';verify.mkdir(exist_ok=True);(verify/'independent_primers.tsv').write_text(''.join(f"{x['name']}\t{x['F']}\t{x['R']}\n" for x in selected))
for role,file in [('target','target_panel.fasta'),('background','background.fasta')]:
 f=verify/(role+'.fasta')
 if not f.exists():f.symlink_to(project/'data'/file)
