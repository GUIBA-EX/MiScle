from pathlib import Path
import csv,json,math,xml.etree.ElementTree as ET
from xs_raw_helpers import I,fa
p=Path(__file__).resolve().parents[1]
with (p/'primer_pairs.csv').open() as f:rows=list(csv.DictReader(f))
with (p/'comparison.csv').open() as f:metrics={r['name']:r for r in csv.DictReader(f)}
assert set(metrics)=={r['name'] for r in rows}=={'MiScle-1','MiScle-1D','MiScle-2','MiScle-2D'}
fasta=fa(p/'primers.fasta')
for r in rows:
 for side in ('F','R'):
  assert fasta[r['name']+'-'+side]==r[side]
  assert len(r[side])==int(r[side+'_end'])-int(r[side+'_start'])+1
  assert math.prod(len(I[c]) for c in r[side])==int(r['degeneracy_'+side])<=8
 assert int(r['R_end'])-int(r['F_start'])+1==int(r['reference_span'])
 assert metrics[r['name']]['strict_WC']==metrics[r['name']]['strict_GT']=='0'
 ET.parse(p/f"figures/base_matching_{r['name']}.svg")
for window in ('1','2'):
 normal=next(r for r in rows if r['name']=='MiScle-'+window);deg=next(r for r in rows if r['name']=='MiScle-'+window+'D')
 assert all(normal[k]==deg[k] for k in ('F_start','F_end','R_start','R_end','reference_span'))
assert json.loads((p/'verification.json').read_text())['disagreements']==0
print('PASS: sequences, degeneracy, shared coordinates, figures and verification records')
