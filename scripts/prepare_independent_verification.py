"""Prepare six released pairs and sixteen optional four-oligo cross-pairs."""
from pathlib import Path
import sys,json,itertools
project=Path(sys.argv[1]).resolve();work=Path(sys.argv[2]).resolve();v=work/'verify';v.mkdir(exist_ok=True,parents=True)
selected=[x for x in json.loads((project/'selected.json').read_text()) if x['region']=='S']+json.loads((work/'selected.json').read_text());extras=[]
for mode in ['ordinary','D']:
 d=[x for x in selected if x.get('mode')==mode];r=next(x for x in d if x['intended_suborder']=='Refertina');t=next(x for x in d if x['intended_suborder']=='Vacatina');ol=[r['F'],r['R'],t['F'],t['R']]
 for a,b in itertools.combinations_with_replacement(range(4),2):
  if (a,b) not in [(0,1),(2,3)]:extras.append(dict(name=f'mixture_{mode}_{a}{b}',F=ol[a],R=ol[b],region='M',reference_span=200))
(work/'verification_selected.json').write_text(json.dumps(selected+extras,indent=2));(v/'independent_primers.tsv').write_text(''.join(f"{x['name']}\t{x['F']}\t{x['R']}\n" for x in selected+extras))
for role,f in [('target','target_panel.fasta'),('background','background.fasta')]:
 path=v/(role+'.fasta')
 if not path.exists():path.symlink_to(project/'data'/f)
