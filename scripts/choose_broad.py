from pathlib import Path
import csv,json,collections,sys,functools
p=Path(__file__).parent;root=p.parents[1];sys.path.insert(0,str(root/'work/milily'));from common import thermal,oligo_thermal
rows=json.loads((p/'all_pairs.json').read_text());oligos={(r['side'],r['sequence']):r['id'] for r in json.loads((p/'all_sites/oligos.json').read_text())};h=collections.defaultdict(lambda:collections.defaultdict(list))
for name,acc,ori,gt,pos,mm,last in csv.reader((p/'all_sites/background_sites.tsv').open(),delimiter='\t'):h[name,int(gt)][acc,ori].append((int(pos),int(mm),int(last)))
@functools.lru_cache(None)
def background(F,R):
 result={}
 for gt in (0,1):
  f,r=h[oligos['F',F],gt],h[oligos['R',R],gt];strict=set();loose3=set();loose4=set()
  for acc,ori in f.keys()&r.keys():
   for fp,fm,fl in f[acc,ori]:
    for rp,rm,rl in r[acc,ori]:
     if not 0<=rp-fp-len(F)<=1000:continue
     if fm<=2 and rm<=2 and not fl and not rl:strict.add(acc)
     if fm<=(2 if gt else 3) and rm<=(2 if gt else 3):loose3.add(acc)
     loose4.add(acc)
  result['strict_GT' if gt else 'strict_WC']=len(strict);result['GT2' if gt else 'WC3']=len(loose3)
  if not gt:result['WC4']=len(loose4)
 return result
for i,x in enumerate(rows):
 x.update(background(x['F'],x['R']))
 if i%2000==0:print('screened',i,len(rows),flush=True)
selected=[];diagnostic=[];cache={}
def physical(x):
 key=(x['F'],x['R'])
 if key not in cache:cache[key]=thermal(*key)
 return cache[key]
for reg in ('S','M','L'):
 rr=[x for x in rows if x['region']==reg];good=[x for x in rr if max(x['strict_WC'],x['strict_GT'])<=20 and max(x['WC3'],x['GT2'])<=100];print(reg,'background_pass',len(good),'minimumstrict',min(max(x['strict_WC'],x['strict_GT']) for x in rr),flush=True)
 # Keep both coverage-first and exact-first candidates within modest coverage loss.
 for criterion in ('coverage','exact'):
  candidates=sorted(good,key=(lambda x:(x['score'],x['pred_exact'])) if criterion=='coverage' else (lambda x:(x['pred_exact'],x['score'])),reverse=True)
  for x in candidates[:300]:
   ok,t=physical(x)
   if ok and all(min(z['F_tm_min'],z['R_tm_min'])>=45 for z in t):
    if not any((z['F'],z['R'])==(x['F'],x['R']) for z in selected):selected.append(dict(x,name=f'{reg}-{criterion}',thermal=t,background_gate=True));print('SELECT',reg,criterion,x['F'],x['R'],x['pred_covered'],x['pred_exact'],x['pred_evaluable'],x['strict_WC'],x['strict_GT'],x['WC3'],flush=True)
    break
  else:print('NO chemical pass',reg,criterion,flush=True)
 # Preserve a low-background diagnostic when conservative gates admit nothing.
 if not any(x['region']==reg for x in selected):
  for x in sorted(rr,key=lambda x:(max(x['strict_WC'],x['strict_GT']),x['WC3'],-x['score']))[:50]:
   ok,t=physical(x)
   if ok:selected.append(dict(x,name=f'{reg}-diagnostic',thermal=t,background_gate=False));print('DIAGNOSTIC',reg,x['pred_covered'],x['pred_exact'],x['strict_WC'],x['strict_GT'],flush=True);break
(p/'selected.json').write_text(json.dumps(selected,indent=2));(p/'all_screened.json').write_text(json.dumps(rows));(p/'selection_manifest.json').write_text(json.dumps(dict(pairs=len(rows),background_unique_pairs=background.cache_info().currsize,thermal_checked=len(cache),strict_WC_max=20,strict_GT_max=20,WC3_max=100,GT2_max=100),indent=2))
