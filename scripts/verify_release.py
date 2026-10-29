from pathlib import Path
import csv,json,sys,collections
p=Path(__file__).parent;root=p.parents[1];out=root/'outputs/MiScle';sys.path.insert(0,str(root/'outputs/MiUrchin/scripts/species_screening'));from xs_raw_helpers import fa,rc,scan
panel=list(csv.DictReader((p/'curated/panel.csv').open()));ids={r['accession'] for r in panel};target=fa(p/'curated/panel.fasta');background=fa(p/'background.fasta');sel=json.loads((p/'release_selected.json').read_text());hits=list(csv.DictReader((p/'release_verify/independent_hits.tsv').open(),delimiter='\t'));checks=[];templates=[];h=collections.defaultdict(set);final={'S':'S','M':'M','SD':'SD','MD':'MD'}
for r in hits:
 for flag in ('exact','le2last3'):
  if r[flag]=='1':h[r['name'],r['role'],r['model'],flag].add(r['accession'])
for x in sel:
 name=x['name'];key='XS-'+name;F,R=x['F'],x['R'];roles=[('target',target)]
 if name in final:roles.append(('background',background))
 for role,db in roles:
  for model in (['WC'] if role=='target' else ['WC','GT']):
   cv=set();ex=set()
   for acc,s in db.items():
    found=[]
    for ori,t in [('+',s),('-',rc(s))]:
     for fp,fm,fl in scan(t,F,False,2,model=='GT'):
      for rp,rm,rl in scan(t,R,True,2,model=='GT'):
       inner=rp-fp-len(F);size=rp+len(R)-fp
       if not ((inner>=20 and 70<=size<=500) if role=='target' else 0<=inner<=1000):continue
       if not fl and not rl:
        cv.add(acc)
        if role=='target' and name in final:found.append((fm+rm,abs(size-x['reference_span']),fp,rp,ori,t[fp:fp+len(F)],rc(t[rp:rp+len(R)]),t[fp+len(F):rp],size))
       if fm==0 and rm==0:ex.add(acc)
    if found:
     z=min(found);templates.append(dict(name=final[name],candidate=name,accession=acc,F_template=z[5],R_template=z[6],internal_barcode=z[7],product_bp=z[8],orientation=z[4],F_start_raw=z[2]+1,R_start_raw=z[3]+1,alternative_compatible_products=len(found)))
   for flag,got in [('le2last3',cv),('exact',ex)]:
    expected=h[key,role,model,flag]&(ids if role=='target' else set(db));assert got==expected,(name,role,model,flag,len(got^expected));checks.append(dict(name=name,role=role,model=model,flag=flag,records=len(db),disagreements=0))
   print('VERIFIED',name,role,model,len(cv),len(ex),flush=True)
def write(f,rr):
 with f.open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=rr[0],lineterminator='\n');w.writeheader();w.writerows(rr)
write(out/'independent_checks.csv',checks);write(out/'verified_templates.csv',templates);(out/'verification.json').write_text(json.dumps(dict(methods=['C++ exhaustive full-sequence paired scan','Independent Python bit-parallel full-sequence paired scan'],target_pairs_checked=len(sel),strict_background_pairs_checked=4,independent_set_checks=len(checks),disagreements=0,target_product_bp=[70,500],target_min_internal_bp=20,background_internal_bp=[0,1000]),indent=2)+'\n')
