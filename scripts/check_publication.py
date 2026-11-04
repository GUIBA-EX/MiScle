"""Check publication metadata, local links and every released oligo across formats."""
from pathlib import Path
from datetime import date
import hashlib,json,re,subprocess
from pypdf import PdfReader
root=Path(__file__).resolve().parents[1];p=root/'publications';c=json.loads((p/'release.json').read_text())
name,title,version=c['project'],c['title'],c['version'];day=date.fromisoformat(c['date']);url=f'https://github.com/GUIBA-EX/{name}/releases/tag/preprint-v{version}'
md=p/f'{name}_preprint.md';paper=md.read_text();status=f'Preprint v{version} · {day.day} {day:%B %Y} · Not peer reviewed'
assert paper.startswith('# '+title+'\n') and '\n**Fei Xia**\n' in paper and status in paper
files=[root/f for f in c['readmes']]+[md,p/'README.md']
for f in files:
 s=f.read_text()
 for link in re.findall(r'\]\(([^)]+)\)',s):
  if not link.startswith(('https://','http://','#','mailto:')):assert (f.parent/link.split('#')[0]).is_file(),(f,link)
 versions=re.findall(r'preprint-v(1\.\d+)',s);assert all(v==version for v in versions),(f,versions)
for file,lang in c['readmes'].items():
 s=(root/file).read_text();assert '# '+title+'\n' in s and '**Fei Xia**' in s and url in s
 assert f'Computational preprint v{version}.' in s,(file,'citation')
 expected=status if lang=='en' else f'{"預印本" if lang=="zh-Hant" else "预印本"} v{version} · {day.year} 年 {day.month} 月 {day.day} 日'
 assert expected in s,(file,'date')
cff=(root/'CITATION.cff').read_text();bib=(root/'CITATION.bib').read_text()
assert cff.count(f'title: "{title}"')==2 and f'date-released: "{day}"' in cff
assert re.findall(r'version: "([^"]+)"',cff)==[version,version]
assert url in cff and url in bib and title in bib and f'version = {{{version}}}' in bib
pdf=p/f'{name}_preprint.pdf';reader=PdfReader(pdf)
assert reader.metadata.author=='Fei Xia' and reader.metadata.title==title
text=subprocess.check_output(['pdftotext','-layout',str(pdf),'-'],text=True);flat=re.sub(r'\s+','',text)
assert re.sub(r'\s+','',status) in flat
seqs={};key=None
for line in (root/c['fasta']).read_text().splitlines():
 if line.startswith('>'):
  key=line[1:];assert key not in seqs;seqs[key]=''
 elif line.strip():
  assert key is not None;seqs[key]+=line.strip()
for key,seq in seqs.items():assert seq and seq in flat and seq in paper,(key,'primer missing')
for file in c['readmes']:
 s=(root/file).read_text()
 for key,seq in seqs.items():assert seq in s,(file,key,'primer missing')
for i,page in enumerate(reader.pages):assert len((page.extract_text() or '').strip())>150 or len(page.images)>0,('near-empty page',i+1)
# Figure captions must be visible as text, rather than only image alt attributes.
for num in re.findall(r'\*\*Figure ([\d]+[a-z]?)\.',paper):assert f'Figure{num}.' in flat,num
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
report=dict(project=name,version=version,date=str(day),author='Fei Xia',pages=len(reader.pages),oligos=len(seqs),status='PASS',checks=['README titles, dates and release citations','CFF and BibTeX','local document links','PDF author and title metadata','all released oligos in manuscript and PDF','visible figure captions','nonempty pages'],artifacts={str(f.relative_to(root)):sha(f) for f in [md,pdf,root/c['fasta']]})
(p/'publication_check.json').write_text(json.dumps(report,indent=2)+'\n')
print(name,version,'PASS:',len(reader.pages),'pages,',len(seqs),'oligos')
