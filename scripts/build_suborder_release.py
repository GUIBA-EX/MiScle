from pathlib import Path
import csv,json,collections,shutil,html,math,itertools
import sys
out=Path(__file__).resolve().parents[1]
p=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else out/'results/independent_suborder'
def read(f):return list(csv.DictReader(Path(f).open()))
def write(f,rr):
 if not rr:return
 with Path(f).open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=rr[0].keys(),lineterminator='\n');w.writeheader();w.writerows(rr)
baseline=out/'results/baseline';baseline.mkdir(parents=True,exist_ok=True)
for name in ['primer_pairs.csv','primers.fasta','comparison.csv','by_species.csv','by_taxon.csv','base_matching.csv','selected.json','thermal.json','two_pair_comparison.csv','barcode_resolution.csv','verification.json']:
 if not (baseline/name).exists():shutil.copyfile(out/name,baseline/name)
sel=json.loads((p/'verification_selected.json').read_text());final=[x for x in sel if x['name'].startswith('MiScle-')];extras=[x for x in sel if not x['name'].startswith('MiScle-')];names=[x['name'] for x in final];meta={r['accession']:r for r in read(out/'data/panel.csv')};groups={a:next((g for g in ['Refertina','Vacatina'] if g in r['lineage'].split(';')),'unassigned') for a,r in meta.items()};old=read(baseline/'by_species.csv');valid={x['name']:{r['accession'] for r in old if r['name']==('MiScle-SD' if x['region']=='S' else 'MiScle-MD') and r['evaluable']=='1'} for x in final};h=collections.defaultdict(set)
for r in csv.DictReader((p/'independent_hits.tsv').open(),delimiter='\t'):
 for flag in ['exact','le2last3','le2','le3','le4']:
  if r[flag]=='1':h[r['name'],r['role'],r['model'],flag].add(r['accession'])
def hits(n,role='target',model='WC',flag='le2last3'):return h[n,role,model,flag]
cov={n:hits(n)&valid[n] for n in names};ex={n:hits(n,flag='exact')&valid[n] for n in names};pairs=[];metrics=[];by=[];tax=[];sub=[]
for x in final:
 n=x['name'];z={k:x[k] for k in ['name','region','F','R','degeneracy_F','degeneracy_R','F_start','F_end','R_start','R_end','reference_span']};z['intended_suborder']=x.get('intended_suborder','');pairs.append(z);thermal=x.get('mixture_thermal',x.get('thermal',[]));maxst=max(t.get('max_structure_tm',t.get('structure_max',0)) for t in thermal);row=dict(z,evaluable=len(valid[n]),unknown=len(meta)-len(valid[n]),covered=len(cov[n]),exact=len(ex[n]),coverage_pct=100*len(cov[n])/len(valid[n]),exact_pct=100*len(ex[n])/len(valid[n]),strict_WC=len(hits(n,'background')),strict_GT=len(hits(n,'background','GT')),WC3=len(hits(n,'background','WC','le3')),GT2=len(hits(n,'background','GT','le2')),WC4=len(hits(n,'background','WC','le4')),thermal_max=maxst);metrics.append(row)
 for a,m in meta.items():by.append(dict(name=n,accession=a,species=m['organism'],family=m['family'],genus=m['genus'],suborder=groups[a],evaluable=int(a in valid[n]),covered=int(a in cov[n]),exact=int(a in ex[n])))
 for rank in ['family','genus','suborder']:
  for t in sorted(set(groups.values() if rank=='suborder' else (r[rank] for r in meta.values()))):
   ids={a for a,m in meta.items() if (groups[a] if rank=='suborder' else m[rank])==t};v=ids&valid[n];r=dict(name=n,rank=rank,taxon=t,species_total=len(ids),evaluable=len(v),covered=len(ids&cov[n]),exact=len(ids&ex[n]));tax.append(r)
   if rank=='suborder':sub.append(r)
write(out/'primer_pairs.csv',pairs);write(out/'comparison.csv',metrics);write(out/'by_species.csv',by);write(out/'by_taxon.csv',tax);write(out/'by_suborder.csv',sub);(out/'primers.fasta').write_text(''.join(f">{r['name']}-{s}\n{r[s]}\n" for r in pairs for s in 'FR'));(out/'selected.json').write_text(json.dumps(final,indent=2)+'\n');(out/'thermal.json').write_text(json.dumps({x['name']:x.get('mixture_thermal',x.get('thermal')) for x in final},indent=2)+'\n')
combos={'M_ordinary':['MiScle-MR','MiScle-MV'],'M_D':['MiScle-MRD','MiScle-MVD'],'S_plus_M_ordinary':['MiScle-S','MiScle-MR','MiScle-MV'],'S_plus_M_D':['MiScle-SD','MiScle-MRD','MiScle-MVD']};combs=[];cby=[];csub=[]
for label,nn in combos.items():
 v=set.union(*(valid[n] for n in nn));cv=set.union(*(cov[n] for n in nn));ee=set.union(*(ex[n] for n in nn));wc=set.union(*(hits(n,'background') for n in nn));gt=set.union(*(hits(n,'background','GT') for n in nn));combs.append(dict(scheme=label,pairs='+'.join(nn),evaluable_any=len(v),unknown=len(meta)-len(v),covered=len(cv),exact=len(ee),coverage_pct=100*len(cv)/len(v),exact_pct=100*len(ee)/len(v),strict_WC=len(wc),strict_GT=len(gt)))
 for g in sorted(set(groups.values())):
  ids={a for a in meta if groups[a]==g};vv=ids&v;csub.append(dict(scheme=label,suborder=g,evaluable=len(vv),covered=len(cv&ids),exact=len(ee&ids)))
 for a,m in meta.items():cby.append(dict(scheme=label,accession=a,species=m['organism'],family=m['family'],suborder=groups[a],evaluable=int(a in v),covered=int(a in cv),exact=int(a in ee)))
write(out/'two_pair_comparison.csv',combs);write(out/'combination_by_species.csv',cby);write(out/'combination_by_suborder.csv',csub);write(out/'missed_species.csv',[r for r in cby if r['scheme']=='M_D' and r['evaluable'] and not r['covered']]);write(out/'unknown_species.csv',[r for r in cby if r['scheme']=='M_D' and not r['evaluable']])
temps=[r for r in read(p/'verified_templates.csv') if r['name'] in names and r['accession'] in valid[r['name']]];write(out/'verified_templates.csv',temps);checks=read(p/'independent_checks.csv');write(out/'independent_checks.csv',checks);shutil.copyfile(p/'verification.json',out/'verification.json');shutil.copyfile(p/'independent_hits.tsv',out/'raw_hits.tsv');shutil.copyfile(p/'verification.log',out/'verify_reproduction.log')
# Background lineage audit.
bgmeta={r['accession']:r for r in read(out/'data/background_metadata.csv')};bgr=[]
for mode,nn in [('ordinary',combos['M_ordinary']),('D',combos['M_D'])]:
 for model in ['WC','GT']:
  normal=set.union(*(hits(n,'background',model) for n in nn));extra=set.union(*(hits(x['name'],'background',model) for x in extras if x['name'].startswith('mixture_'+mode+'_')));allhits=normal|extra
  print('MIX',mode,model,'pairhits',len(normal),'extra',len(extra),'total',len(allhits))
  for a in sorted(allhits):bgr.append(dict(mode=mode,model=model,accession=a,organism=bgmeta[a]['organism'],lineage=bgmeta[a]['group'],from_intended_pairs=int(a in normal),from_other_mixture_pairs=int(a in extra)))
write(out/'mixture_background_hits.csv',bgr)
write(out/'background_hits.csv',[r for r in bgr if r['from_intended_pairs']])
mixtures=[]
for mode,nn in [('ordinary',combos['M_ordinary']),('D',combos['M_D'])]:
 for model in ['WC','GT']:
  normal=set.union(*(hits(n,'background',model) for n in nn));extra=set.union(*(hits(x['name'],'background',model) for x in extras if x['name'].startswith('mixture_'+mode+'_')));mixtures.append(dict(mode=mode,model=model,separate_pair_union=len(normal),mixed_all_pair_union=len(normal|extra),additional_background=len(extra-normal)))
write(out/'mixture_summary.csv',mixtures)
# Internal-barcode sharing, with ambiguous barcodes counted separately.
res=[];sharing=[]
for x in final:
 n=x['name'];tt=[r for r in temps if r['name']==n];b=collections.defaultdict(list)
 for r in tt:
  if set(r['internal_barcode'])<=set('ACGT'):b[r['internal_barcode']].append(r)
 affected=set();genera={meta[r['accession']]['genus'] for r in tt}
 for barcode,ls in b.items():
  gg={meta[r['accession']]['genus'] for r in ls}
  if len(gg)>1:affected|=gg
  if len(ls)>1:sharing.append(dict(name=n,species_count=len(ls),genera=';'.join(sorted(gg)),species=';'.join(meta[r['accession']]['organism'] for r in ls),barcode=barcode))
 res.append(dict(name=n,covered_species=len(tt),product_min=min(int(r['product_bp']) for r in tt),product_max=max(int(r['product_bp']) for r in tt),internal_min=min(len(r['internal_barcode']) for r in tt),internal_max=max(len(r['internal_barcode']) for r in tt),observed_unshared_species=sum(len(ls)==1 for ls in b.values()),ambiguous_internal=sum(bool(set(r['internal_barcode'])-set('ACGT')) for r in tt),covered_genera=len(genera),genera_without_observed_cross_genus_sharing=len(genera-affected)))
write(out/'barcode_resolution.csv',res);write(out/'shared_barcodes.csv',sharing)
# Base-by-base figures, preserving actual independently verified hit templates.
sysdummy=None
import sys
sys.path.insert(0,str(out/'scripts'));from xs_raw_helpers import fa,rc,I
al=fa(out/'data/alignment.fasta');cols=[i for i,c in enumerate(al['NC_003522.1']) if c!='-'];dd={(r['name'],r['accession']):r for r in by};tm={(r['name'],r['accession']):r for r in temps};basedata=[]
for oldfig in (out/'figures').glob('base_matching_*.svg'):oldfig.unlink()
for x in final:
 n=x['name'];fam=collections.defaultdict(list)
 for a,m in meta.items():fam[(groups[a],m['family'])].append(m)
 height=210+len(meta)*21+len(fam)*30+80;el=[f'<svg xmlns="http://www.w3.org/2000/svg" width="1900" height="{height}" viewBox="0 0 1900 {height}"><rect width="100%" height="100%" fill="#f7fbff"/>']
 def text(xx,y,s,size=14,mono=False):el.append(f'<text x="{xx}" y="{y}" font-family="{"monospace" if mono else "sans-serif"}" font-size="{size}" fill="#173647">{html.escape(str(s))}</text>')
 def rect(xx,y,w,h,c,stroke=None):el.append(f'<rect x="{xx}" y="{y}" width="{w}" height="{h}" fill="{c}"'+(f' stroke="{stroke}" stroke-width="2"' if stroke else '')+'/>')
 text(22,36,n+' · Primer–template base matching'+(' · '+x['intended_suborder'] if x.get('intended_suborder') else ''),24);text(22,65,'NC_003522.1 ungapped 16S | shared reference endpoints | F and R templates shown 5′→3′',17);text(22,91,'Blue: compatible | Orange: mismatch | Grey: gap / unknown | Box: final 3 nt',15);text(22,128,'Species / accession',15);text(535,128,'Hit',15)
 for side,xx in [('F',580),('R',1200)]:
  text(xx,128,f"{side} {x[side+'_start']}–{x[side+'_end']} | variants {x['degeneracy_'+side]}",15)
  for j,c in enumerate(x[side]):rect(xx+j*19,140,18,21,'#cce7f0');text(xx+j*19+4,156,c,14,True)
 y=180
 for (g,family),ms in sorted(fam.items(),key=lambda item:(item[0][0]!=x.get('intended_suborder',item[0][0]),item[0])):
  rect(16,y,1868,25,'#e2eaf3');text(22,y+18,f'{g} · {family} | {len(ms)} reference species',15);y+=30
  for i,m in enumerate(sorted(ms,key=lambda z:z['organism'])):
   if i%2==0:rect(16,y,1868,20,'#eef5f9')
   a=m['accession'];d=dd[n,a];state='?' if not d['evaluable'] else 'E' if d['exact'] else 'C' if d['covered'] else '–';text(22,y+15,m['organism']+' / '+a,12);text(540,y+15,state)
   for side,xx in [('F',580),('R',1200)]:
    if (n,a) in tm:q=tm[n,a][side+'_template'];source='verified_raw_hit'
    else:
     q=''.join(al[a][cols[j]] for j in range(x[side+'_start']-1,x[side+'_end']));q=rc(q) if side=='R' else q;source='homologous_reference_columns'
    assert len(q)==len(x[side]);basedata.append(dict(name=n,side=side,accession=a,species=m['organism'],family=m['family'],suborder=g,primer=x[side],template=q,template_source=source,state=state))
    for j,c in enumerate(q):rect(xx+j*19,y,18,19,'#d6edf4' if c in I[x[side][j]] else '#e3e8ee' if c not in 'ACGT' else '#ffd5a1');text(xx+j*19+4,y+14,c,13,True)
   y+=21
 for side,xx in [('F',580),('R',1200)]:rect(xx+(len(x[side])-3)*19-1,139,58,y-139,'none','#173647')
 text(22,y+25,'E: exact pair | C: ≤2 mismatches per primer, compatible final 3 nt | –: no predicted hit | ?: not evaluable',13);text(22,y+47,'Covered rows use independently verified raw hits; other rows show homologous reference columns.',13);el.append('</svg>');(out/f'figures/base_matching_{n}.svg').write_text(''.join(el))
write(out/'base_matching.csv',basedata)

print("Rebuilt released tables and six base-matching figures from frozen scans.")
