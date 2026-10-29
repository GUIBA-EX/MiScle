from pathlib import Path
import csv,json,sys,collections,itertools
from Bio.Data.IUPACData import ambiguous_dna_values as I
p=Path(__file__).parent;root=p.parents[1];sys.path.insert(0,str(root/'work/milily'));from common import expand,thermal
sel={r['name']:r for r in json.loads((p/'verify_selected.json').read_text())};templates=list(csv.DictReader((root/'outputs/MiScle/verified_templates.csv').open()));allpool={};pairs=[]
def terminal(q):return any(c in 'GC' for c in q[-3:]) and not any(len(set(q[j:j+5]))==1 for j in range(len(q)-4))
for candidate,region in [('S-exact','S'),('S-complement-coverage','M')]:
 d=sel[candidate];tt=[r for r in templates if r['candidate']==candidate];pools={}
 for side in ('F','R'):
  seqs=set(expand(d[side]));seqs.update(r[side+'_template'] for r in tt);counts=collections.Counter(r[side+'_template'] for r in tt);L=len(d[side]);consensus=''.join(collections.Counter(q[j] for q in counts.elements()).most_common(1)[0][0] for j in range(L));seqs.add(consensus);seqs={q for q in seqs if len(q)==L and set(q)<=set('ACGT') and terminal(q)};pools[side]=sorted(seqs)
  for q in pools[side]:allpool[side,q]=None
 for F,R in itertools.product(pools['F'],pools['R']):
  ok,t=thermal(F,R)
  if ok and all(min(x['F_tm_min'],x['R_tm_min'])>=45 for x in t):pairs.append(dict(region=region,F=F,R=R,F_start=d['F_start'],F_end=d['F_end'],R_start=d['R_start'],R_end=d['R_end'],reference_span=d['reference_span'],thermal=t))
 print(region,'ordinary pools',len(pools['F']),len(pools['R']),'chemical pairs',sum(x['region']==region for x in pairs),flush=True)
for i,k in enumerate(allpool):allpool[k]=f'N{i:04d}'
(p/'ordinary_pool.json').write_text(json.dumps([dict(id=v,side=s,sequence=q) for (s,q),v in allpool.items()]));(p/'ordinary_pairs.json').write_text(json.dumps(pairs,indent=2));
for label,source in [('ordinary_target',p/'curated/panel.fasta'),('ordinary_background',p/'background.fasta')]:
 d=p/label;(d/'sites').mkdir(parents=True,exist_ok=True);f=d/'background.fasta'
 if not f.exists():f.symlink_to(source.resolve())
 (d/'sites/oligos.tsv').write_text(''.join(f'{v} {s} {q}\n' for (s,q),v in allpool.items()))
