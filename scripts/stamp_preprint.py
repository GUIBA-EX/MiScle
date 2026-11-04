"""Synchronize release metadata, every README edition, and both citation formats."""
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
import json,re,sys
root=Path(__file__).resolve().parents[1]
p=root/'publications/release.json';c=json.loads(p.read_text())
version=sys.argv[1];assert re.fullmatch(r'1\.\d+',version)
day=datetime.now(ZoneInfo('Asia/Shanghai')).date()
c.update(version=version,date=str(day));p.write_text(json.dumps(c,indent=2,ensure_ascii=False)+'\n')
name,title=c['project'],c['title'];url=f'https://github.com/GUIBA-EX/{name}/releases/tag/preprint-v{version}'
date=f'{day.day} {day:%B %Y}'
citation=f'Xia, F. ({day.year}). *{title}*. Computational preprint v{version}. [GitHub Release]({url}).'
for file,lang in c['readmes'].items():
 p=root/file;s=p.read_text()
 if lang=='en':
  status=f'**Preprint v{version} · {date} · Not peer reviewed.**'
  links=f'[PDF](publications/{name}_preprint.pdf) · [Manuscript](publications/{name}_preprint.md) · [Versioned release]({url}) · [CITATION.cff](CITATION.cff) · [BibTeX](CITATION.bib)'
 elif lang=='zh-Hant':
  status=f'**預印本 v{version} · {day.year} 年 {day.month} 月 {day.day} 日 · 尚未經同儕審查。**'
  links=f'[PDF](publications/{name}_preprint.pdf) · [論文 Markdown](publications/{name}_preprint.md) · [固定版本]({url}) · [CITATION.cff](CITATION.cff) · [BibTeX](CITATION.bib)'
 else:
  status=f'**预印本 v{version} · {day.year} 年 {day.month} 月 {day.day} 日 · 未经同行评审。**'
  links=f'[PDF](publications/{name}_preprint.pdf) · [论文 Markdown](publications/{name}_preprint.md) · [固定版本]({url}) · [CITATION.cff](CITATION.cff) · [BibTeX](CITATION.bib)'
 links += ' · [Build and updates](publications/README.md)' if lang=='en' else ' · [构建与更新](publications/README.md)'
 s=re.sub(r'<!-- publication:start -->.*?<!-- publication:end -->',f'<!-- publication:start -->\n{status}\n\n{links}\n<!-- publication:end -->',s,flags=re.S)
 s=re.sub(r'^Xia, F\..*$',citation,s,flags=re.M);p.write_text(s)
p=root/f'publications/{name}_preprint.md';s=p.read_text();s=re.sub(r'^# .*', '# '+title,s,count=1)
s=re.sub(r'(?:Computational preprint · |Preprint v[\d.]+ · )\d+ [A-Za-z]+ \d{4}(?: · Not peer reviewed)?',f'Preprint v{version} · {date} · Not peer reviewed',s)
s=re.sub(r'preprint-v1\.\d+',f'preprint-v{version}',s);p.write_text(s)
(root/'CITATION.cff').write_text(f'''cff-version: 1.2.0
message: "If you use this project's primers, data or analysis code, please cite the versioned preprint."
type: software
title: "{title}"
authors:
  - family-names: "Xia"
    given-names: "Fei"
date-released: "{day}"
url: "https://github.com/GUIBA-EX/{name}"
repository-code: "https://github.com/GUIBA-EX/{name}"
version: "{version}"
preferred-citation:
  type: unpublished
  title: "{title}"
  authors:
    - family-names: Xia
      given-names: Fei
  year: {day.year}
  version: "{version}"
  notes: "Computational preprint, version {version}. Published on GitHub."
  url: "{url}"
''')
(root/'CITATION.bib').write_text(f'''@misc{{xia{day.year}{name.lower()},
  author = {{Xia, Fei}},
  title = {{{{{title}}}}},
  year = {{{day.year}}},
  version = {{{version}}},
  howpublished = {{Computational preprint v{version}}},
  url = {{{url}}}
}}
''')
