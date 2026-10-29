"""Recheck four primer pairs against frozen raw references (Python stdlib only)."""
from pathlib import Path
import csv,collections,json
from xs_raw_helpers import fa,rc,scan
root=Path(__file__).resolve().parents[1]
selected=json.loads((root/'selected.json').read_text());ids=set(fa(root/'data/target_panel.fasta'));hits=collections.defaultdict(set)
with (root/'raw_hits.tsv').open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  for flag in ('exact','le2last3'):
   if r[flag]=='1':hits[r['name'],r['role'],r['model'],flag].add(r['accession'])
for x in selected:
 F,R=x['F'],x['R'];key='XS-'+x['trial_id']
 for role,file in [('target','target_panel.fasta'),('background','background.fasta')]:
  db=fa(root/'data'/file)
  for model in (['WC'] if role=='target' else ['WC','GT']):
   cv=set();exact=set()
   for acc,s in db.items():
    for t in (s,rc(s)):
     for fp,fm,fl in scan(t,F,False,2,model=='GT'):
      for rp,rm,rl in scan(t,R,True,2,model=='GT'):
       inner=rp-fp-len(F);size=rp+len(R)-fp
       if not ((inner>=20 and 70<=size<=500) if role=='target' else 0<=inner<=1000):continue
       if not fl and not rl:cv.add(acc)
       if fm==0 and rm==0:exact.add(acc)
   for flag,got in [('le2last3',cv),('exact',exact)]:
    expected=hits[key,role,model,flag]&(ids if role=='target' else set(db));assert got==expected,(x['name'],role,model,flag)
   print('PASS',x['name'],role,model,len(cv),len(exact),flush=True)
print('All target and strict background sets match the frozen C++ scan.')
