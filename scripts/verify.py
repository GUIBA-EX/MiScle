"""Independent raw-sequence checks. Python standard library only; four local workers."""
from pathlib import Path
import csv,json,sys,collections,concurrent.futures
project=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
work=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else project/'results/independent_suborder'
sys.path.insert(0,str(project/'scripts'));from xs_raw_helpers import fa,rc,scan
DB={}
def init():
 for role,f in [('target','target_panel.fasta'),('background','background.fasta')]:DB[role]=fa(project/'data'/f)
def run(job):
 x,role,model=job;F,R=x['F'],x['R'];cv=set();ex=set();templates=[]
 for acc,s in DB[role].items():
  found=[]
  for ori,t in [('+',s),('-',rc(s))]:
   for fp,fm,fl in scan(t,F,False,2,model=='GT'):
    for rp,rm,rl in scan(t,R,True,2,model=='GT'):
     inner=rp-fp-len(F);size=rp+len(R)-fp
     if not ((inner>=20 and 70<=size<=500) if role=='target' else 0<=inner<=1000):continue
     if not fl and not rl:
      cv.add(acc)
      if role=='target':found.append((fm+rm,abs(size-x['reference_span']),fp,rp,ori,t[fp:fp+len(F)],rc(t[rp:rp+len(R)]),t[fp+len(F):rp],size))
     if fm==0 and rm==0:ex.add(acc)
  if found:
   z=min(found);templates.append(dict(name=x['name'],candidate=x['name'],accession=acc,F_template=z[5],R_template=z[6],internal_barcode=z[7],product_bp=z[8],orientation=z[4],F_start_raw=z[2]+1,R_start_raw=z[3]+1,alternative_compatible_products=len(found)))
 return x['name'],role,model,sorted(cv),sorted(ex),templates
if __name__=='__main__':
 selected=json.loads((work/'verification_selected.json').read_text());h=collections.defaultdict(set)
 for r in csv.DictReader((work/'independent_hits.tsv').open(),delimiter='\t'):
  for flag in ['exact','le2last3']:
   if r[flag]=='1':h[r['name'],r['role'],r['model'],flag].add(r['accession'])
 jobs=[(x,role,model) for x in selected for role in ['target','background'] for model in (['WC'] if role=='target' else ['WC','GT'])];checks=[];templates=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=4,initializer=init) as pool:
  for name,role,model,cv,ex,tt in pool.map(run,jobs):
   for flag,got in [('le2last3',set(cv)),('exact',set(ex))]:
    expected=h[name,role,model,flag];assert got==expected,(name,role,model,flag,len(got^expected));checks.append(dict(name=name,role=role,model=model,flag=flag,disagreements=0))
   templates+=tt;print('PASS',name,role,model,len(cv),len(ex),flush=True)
 for filename,rows in [('independent_checks.csv',checks),('verified_templates.csv',templates)]:
  with (work/filename).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys(),lineterminator='\n');w.writeheader();w.writerows(rows)
 (work/'verification.json').write_text(json.dumps(dict(methods=['C++ exhaustive raw-sequence paired scan','Independent Python bit-parallel raw-sequence paired scan'],final_pairs=sum(x['name'].startswith('MiScle-') for x in selected),mixture_cross_pairs=sum(not x['name'].startswith('MiScle-') for x in selected),independent_set_checks=len(checks),disagreements=0,target_product_bp=[70,500],target_min_internal_bp=20,background_internal_bp=[0,1000]),indent=2)+'\n')
 print('All independent paired sets agree.',flush=True)
