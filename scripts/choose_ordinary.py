from pathlib import Path
import csv,json,collections
p=Path(__file__).parent;root=p.parents[1];rows=json.loads((p/'ordinary_pairs.json').read_text());pool={(r['side'],r['sequence']):r['id'] for r in json.loads((p/'ordinary_pool.json').read_text())};data={}
for role in ('target','background'):
 h=collections.defaultdict(lambda:collections.defaultdict(list))
 for name,acc,ori,gt,pos,mm,last in csv.reader((p/f'ordinary_{role}/sites/background_sites.tsv').open(),delimiter='\t'):h[name,int(gt)][acc,ori].append((int(pos),int(mm),int(last)))
 data[role]=h
validsets=json.loads((p/'evaluable_sets.json').read_text());selected=[];screened=[]
for i,x in enumerate(rows):
 for role in ('target','background'):
  for gt in ((0,) if role=='target' else (0,1)):
   h=data[role];f,r=h[pool['F',x['F']],gt],h[pool['R',x['R']],gt];cv=set();ex=set();loose=set()
   for acc,ori in f.keys()&r.keys():
    for fp,fm,fl in f[acc,ori]:
     for rp,rm,rl in r[acc,ori]:
      inner=rp-fp-len(x['F']);size=rp+len(x['R'])-fp
      if not ((inner>=20 and 70<=size<=500) if role=='target' else 0<=inner<=1000):continue
      if fm<=2 and rm<=2 and not fl and not rl:cv.add(acc)
      if fm==0 and rm==0:ex.add(acc)
      if fm<=(2 if gt else 3) and rm<=(2 if gt else 3):loose.add(acc)
   if role=='target':
    valid=set(validsets['S-exact' if x['region']=='S' else 'S-complement-coverage']);x['covered']=len(cv&valid);x['exact']=len(ex&valid)
   else:x['strict_GT' if gt else 'strict_WC']=len(cv);x['GT2' if gt else 'WC3']=len(loose)
 screened.append(x)
for reg in ('S','M'):
 good=[x for x in rows if x['region']==reg and max(x['strict_WC'],x['strict_GT'])<=20 and max(x['WC3'],x['GT2'])<=100];x=max(good,key=lambda x:(x['covered'],x['exact'],-x['strict_GT'],-x['WC3']));selected.append(dict(x,name=reg,degeneracy_F=1,degeneracy_R=1,background_gate=True));print('ORDINARY',reg,x['F'],x['R'],x['covered'],x['exact'],x['strict_WC'],x['strict_GT'],x['WC3'],flush=True)
deg={r['name']:r for r in json.loads((p/'verify_selected.json').read_text())}
for source,name in [('S-exact','SD'),('S-complement-coverage','MD')]:selected.append(dict(deg[source],name=name,region=name[0],source_candidate=source))
(p/'release_selected.json').write_text(json.dumps(selected,indent=2));(p/'ordinary_screened.json').write_text(json.dumps(screened));v=p/'release_verify';v.mkdir(exist_ok=True)
(v/'independent_primers.tsv').write_text(''.join(f"XS-{x['name']}\t{x['F']}\t{x['R']}\n" for x in selected));(v/'target_barcodes.tsv').write_text('name\taccession\tsequence\n')
for name,src in [('target.fasta',p/'curated/target.fasta'),('background.fasta',p/'background.fasta')]:
 f=v/name
 if not f.exists():f.symlink_to(src.resolve())
