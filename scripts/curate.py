from pathlib import Path
from Bio import SeqIO
import csv,json,re,collections,hashlib
p=Path(__file__).parent;out=p/'curated';out.mkdir(exist_ok=True);paper={x.split('.')[0] for x in json.loads((p/'paper_accessions.json').read_text())};records={};sources=collections.defaultdict(set)
for path in sorted((p/'downloads').glob('*.gb')):
 for r in SeqIO.parse(path,'genbank'):
  records[r.id]=r;sources[r.id].add(path.name.split('_')[0])
rows=[];seqs={};reasons=collections.Counter()
for acc,r in records.items():
 tax=r.annotations.get('taxonomy',[]);org=r.annotations.get('organism','');is_target='Scleractinia' in tax;lineage=';'.join(tax);fam=next((x for x in reversed(tax) if x.endswith('idae')),'unresolved');sid=next((f for f in r.features if f.type=='source'),None);taxid=next((x.split(':')[1] for x in sid.qualifiers.get('db_xref',[]) if x.startswith('taxon:')),'') if sid else '';candidates=[]
 for f in r.features:
  if f.type not in ('rRNA','gene'):continue
  q=' '.join(f.qualifiers.get('product',[])+f.qualifiers.get('gene',[])).lower()
  if any(x in q for x in ('16s','large subunit','l-rrna','rrnl','rnl')) and not '12s' in q:
   candidates.append((f,str(f.extract(r.seq)).upper()))
 source='annotated_16S'
 if candidates:f,seq=max(candidates,key=lambda z:len(z[1]));location=str(f.location)
 elif re.search(r'\b16S\b',r.description,re.I) and 100<=len(r)<=4000:seq=str(r.seq).upper();location='whole_record_16S_title';source='title_16S'
 else:
  reasons['no_16S_annotation_or_title']+=1
  if acc.split('.')[0] in paper:rows.append(dict(accession=acc,organism=org,taxid=taxid,family=fam,genus=org.split()[0],target=int(is_target),paper=1,length=0,ambiguity='',source='none',location='',status='no_16S',lineage=lineage))
  continue
 ambiguous=sum(c not in 'ACGT' for c in seq)/max(1,len(seq));status='eligible' if len(seq)>=150 and ambiguous<=.05 else 'low_quality';species=bool(re.match(r'^[A-Z][a-z]+ [a-z][a-z-]+$',org));row=dict(accession=acc,organism=org,taxid=taxid,family=fam,genus=org.split()[0],target=int(is_target),paper=int(acc.split('.')[0] in paper),length=len(seq),ambiguity=ambiguous,source=source,location=location,status=status,lineage=lineage);rows.append(row)
 if status=='eligible':seqs[acc]=seq
target=[r for r in rows if r['target'] and r['status']=='eligible'];groups=collections.defaultdict(list)
for r in target:
 if re.match(r'^[A-Z][a-z]+ [a-z][a-z-]+$',r['organism']):groups[r['taxid'] or r['organism']].append(r)
panel=[]
for key,rr in groups.items():
 rr.sort(key=lambda x:(-(x['length']>=1200),-x['length'],x['ambiguity'],-x['paper'],x['accession']));panel.append(rr[0])
panel.sort(key=lambda x:(x['family'],x['genus'],x['organism']));all_target={r['accession']:seqs[r['accession']] for r in target};bgrows=[r for r in rows if not r['target'] and r['status']=='eligible'];background={r['accession']:seqs[r['accession']] for r in bgrows}
def write(name,rr):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),lineterminator='\n');w.writeheader();w.writerows(rr)
def fasta(name,d): (out/name).write_text(''.join(f'>{acc}\n{s}\n' for acc,s in d.items()))
write('audit.csv',rows);write('panel.csv',panel);fasta('panel.fasta',{r['accession']:seqs[r['accession']] for r in panel});fasta('target.fasta',all_target);fasta('background_cnidaria.fasta',background)
write('paper_audit.csv',[r for r in rows if r['paper']]);write('background_metadata.csv',bgrows)
summary=dict(downloaded_records=len(records),annotated_or_title_16S=len(seqs),target_records=len(all_target),species_representatives=len(panel),families=dict(collections.Counter(r['family'] for r in panel)),genera=len(set(r['genus'] for r in panel)),non_target_cnidarian_records=len(background),paper_accessions=len(paper),paper_16S_records=sum(r['paper'] and r['status']=='eligible' for r in rows),species_key='source GenBank taxonomy ID with binomial source names; synonyms not externally revised',rejected=dict(reasons));(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
