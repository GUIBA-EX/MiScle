"""Refresh publication hashes; retain independent checks of scientific source files."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1];c=json.loads((root/'publications/release.json').read_text())
hashfile=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
p=root/'SHA256SUMS.json'
if p.exists():
 expected=json.loads(p.read_text())
 for f in expected:
  digest=hashfile(root/f)
  if f.startswith('publications/') or f in [*c['readmes'],'CITATION.bib','CITATION.cff']:
   expected[f]=digest
  else:assert expected[f]==digest,f
 p.write_text(json.dumps(expected,indent=2,ensure_ascii=False)+'\n')
files=[root/f'publications/{c["project"]}_preprint.pdf',root/c['fasta']]
(root/'publications/SHA256SUMS.txt').write_text(''.join(hashfile(f)+'  '+f.name+'\n' for f in files))

legacy=root/'SHA256SUMS'
if legacy.exists():
 lines=[]
 for line in legacy.read_text().splitlines():
  digest,file=line.split('  ',1)
  if file in c['readmes']:digest=hashfile(root/file)
  else:assert digest==hashfile(root/file),file
  lines.append(digest+'  '+file)
 legacy.write_text('\n'.join(lines)+'\n')
