"""Check released oligos, shared endpoints, metrics, figures and verification evidence."""
from pathlib import Path
import csv,json,math,xml.etree.ElementTree as ET,collections,re
from xs_raw_helpers import I,fa
p=Path(__file__).resolve().parents[1]
def read(f):return list(csv.DictReader((p/f).open()))
rows=read('primer_pairs.csv');names={r['name'] for r in rows};expected={'MiScle-S','MiScle-SD','MiScle-MR','MiScle-MV','MiScle-MRD','MiScle-MVD'};assert names==expected
metrics={r['name']:r for r in read('comparison.csv')};assert set(metrics)==names
fasta=fa(p/'primers.fasta');assert len(fasta)==12
for r in rows:
 for side in 'FR':
  assert fasta[r['name']+'-'+side]==r[side]
  assert len(r[side])==int(r[side+'_end'])-int(r[side+'_start'])+1
  assert math.prod(len(I[c]) for c in r[side])==int(r['degeneracy_'+side])<=8
 assert int(r['R_end'])-int(r['F_start'])+1==int(r['reference_span'])
 ET.parse(p/f"figures/base_matching_{r['name']}.svg")
for nn in [('MiScle-S','MiScle-SD'),('MiScle-MR','MiScle-MV','MiScle-MRD','MiScle-MVD')]:
 rr=[r for r in rows if r['name'] in nn]
 assert len({tuple(r[k] for k in ['F_start','F_end','R_start','R_end']) for r in rr})==1
for nn in [('MiScle-MR','MiScle-MV'),('MiScle-MRD','MiScle-MVD')]:assert len({r['F'] for r in rows if r['name'] in nn})==1
by=read('by_species.csv');bases=read('base_matching.csv')
for n in names:
 rr=[r for r in by if r['name']==n];assert len(rr)==307
 for field in ['evaluable','covered','exact']:assert sum(int(r[field]) for r in rr)==int(metrics[n][field])
 for side in 'FR':
  bb=[r for r in bases if r['name']==n and r['side']==side];assert len(bb)==307;assert {r['primer'] for r in bb}=={next(r for r in rows if r['name']==n)[side]}
checks=read('independent_checks.csv');v=json.loads((p/'verification.json').read_text());assert v['disagreements']==0 and len(checks)==v['independent_set_checks']==84 and all(r['disagreements']=='0' for r in checks)
for f in p.rglob('*.md'):
 if '.git' in f.parts:continue
 for link in re.findall(r'\]\(([^)]+)\)',f.read_text()):
  if not re.match(r'https?://|#|mailto:',link):assert (f.parent/link.split('#')[0]).exists(),(f,link)
a=set(re.findall(r'https://doi.org/\S+',(p/'README.md').read_text()));b=set(re.findall(r'https://doi.org/\S+',(p/'README.zh-CN.md').read_text()));assert a==b
assert {r['name'] for r in read('verified_templates.csv')}==names
print('PASS: six pairs, shared F/endpoints, per-oligo degeneracy, metrics, 3,684 base rows, 84 independent checks, links and bilingual DOI sets')
