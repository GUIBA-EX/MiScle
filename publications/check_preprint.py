from pathlib import Path
import csv,json,hashlib,subprocess,re
r=Path(__file__).resolve().parents[1];p=r/'publications';s=(p/'MiScle_preprint.md').read_text()
read=lambda f:list(csv.DictReader((r/f).open()))
for x in read('primer_pairs.csv'):
 for side in 'FR':assert x[side] in s,(x['name'],side)
for x in read('mixture_summary.csv'):
 assert f"| {x['mode']} | {x['model']} | {x['separate_pair_union']} | {x['mixed_all_pair_union']} |" in s
v=json.loads((r/'verification.json').read_text());assert v['disagreements']==0 and v['independent_set_checks']==132
for x in read('by_suborder.csv'):
 n=x['name'];tax=x['taxon']
 if n in ['MiScle-MR','MiScle-MRD'] and tax=='Refertina' or n in ['MiScle-MV','MiScle-MVD'] and tax=='Vacatina':
  assert f"| {n} | {tax} | {x['evaluable']} | {x['covered']} (" in s
files=['primer_pairs.csv','comparison.csv','by_suborder.csv','mixture_summary.csv','verification.json','barcode_resolution.csv','independent_checks.csv']
manifest=p/'source_hashes.json';hashes={f:hashlib.sha256((r/f).read_bytes()).hexdigest() for f in files}
assert json.loads(manifest.read_text())==hashes
(p/'preprint_validation.json').write_text(json.dumps(dict(status='PASS',checks=['frozen source hashes','released sequences','intended suborder counts','same-tube comparison','132 independent checks']),indent=2)+'\n')
print('PASS: MiScle paper sources and numerical comparisons')
