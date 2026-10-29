from pathlib import Path
import csv,json,collections,itertools,functools,math,sys
project=Path(sys.argv[1]).resolve();work=Path(sys.argv[2]).resolve();work.mkdir(exist_ok=True,parents=True)
I=dict(A='A',C='C',G='G',T='T',Y='CT',R='AG',W='AT',K='GT',M='AC',S='CG',B='CGT',D='AGT',H='ACT',V='ACG',N='ACGT');rev={frozenset(v):k for k,v in I.items()}
def read(f):
 source=project/'results/baseline'/f if f in ['by_species.csv','base_matching.csv','primer_pairs.csv'] else project/f
 return list(csv.DictReader(source.open()))
meta={r['accession']:r for r in read('data/panel.csv')};groups={a:next((g for g in ['Refertina','Vacatina'] if g in r['lineage'].split(';')),'unassigned') for a,r in meta.items()};ev={r['accession'] for r in read('by_species.csv') if r['name']=='MiScle-MD' and r['evaluable']=='1'};base=[r for r in read('base_matching.csv') if r['name']=='MiScle-MD' and r['accession'] in ev];template={(r['side'],r['accession']):r['template'] for r in base}
def deg(q):return math.prod(len(I[c]) for c in q)
def merge(a,b):return ''.join(rev[frozenset(I[x]+I[y])] for x,y in zip(a,b))
def terminal(q):return any(set(I[c])<=set('GC') for c in q[-3:]) and not any(set.intersection(*(set(I[c]) for c in q[j:j+6])) for j in range(len(q)-5))
def score(q,side,ids):
 cv=ex=0
 for a in ids:
  t=template[side,a]
  if set(t)-set('ACGT'):continue
  mm=sum(c not in I[p] for p,c in zip(q,t));cv+=mm<=2 and all(c in I[p] for p,c in zip(q[-3:],t[-3:]));ex+=mm==0
 return cv,ex,-deg(q)
pools={};diagnostics=[]
for side,g in [('F','Refertina'),('F','Vacatina'),('R','Refertina'),('R','Vacatina')]:
 ids=sorted(a for a in ev if g=='all' or groups[a]==g);observed=collections.Counter(template[side,a] for a in ids if set(template[side,a])<=set('ACGT'));normal={q for q in observed if terminal(q)}
 consensus=''.join(collections.Counter(t[j] for t in observed.elements()).most_common(1)[0][0] for j in range(20));normal.add(consensus)
 seeds=sorted(sorted(normal),key=lambda q:score(q,side,ids),reverse=True)[:20];beam=set(seeds);allq=set(normal)
 for step in range(4):
  generated=set(beam)
  for q,t in itertools.product(beam,observed):
   z=merge(q,t)
   if deg(z)<=8 and terminal(z):generated.add(z)
  best=sorted(sorted(generated),key=lambda q:score(q,side,ids),reverse=True)[:24];exact=sorted(sorted(generated),key=lambda q:(score(q,side,ids)[1],score(q,side,ids)[0],-deg(q)),reverse=True)[:12];beam=set(best+exact);allq|=beam
 # Retain normal coverage/exact options and limited diverse degenerate options.
 pools[side,g]=sorted(set(seeds)|set(sorted(sorted(normal),key=lambda q:(score(q,side,ids)[1],score(q,side,ids)[0]),reverse=True)[:12])|set(sorted(sorted(allq),key=lambda q:score(q,side,ids),reverse=True)[:35])|set(sorted(sorted(allq),key=lambda q:(score(q,side,ids)[1],score(q,side,ids)[0],-deg(q)),reverse=True)[:20]))
 diagnostics.append(dict(side=side,clade=g,observed=len(observed),retained=len(pools[side,g]),best=[(q,score(q,side,ids)) for q in pools[side,g]][:2]))
# Always compare the released M and MD oligos.
old=read('primer_pairs.csv')
for r in old:
 if r['name'] in ['MiScle-M','MiScle-MD']:
  for side,g in [('F','Refertina'),('F','Vacatina'),('R','Refertina'),('R','Vacatina')]:pools[side,g]=sorted(set(pools[side,g])|{r[side]})
import primer3
kw=dict(mv_conc=50,dv_conc=1.5,dntp_conc=.6)
@functools.lru_cache(None)
def expand(q):return list(map(''.join,itertools.product(*(I[c] for c in q))))
@functools.lru_cache(None)
def single(q,c):return primer3.calc_tm(q,dna_conc=c,**kw),max(primer3.calc_hairpin(q,dna_conc=c,**kw).tm,primer3.calc_homodimer(q,dna_conc=c,**kw).tm)
@functools.lru_cache(None)
def oligo(q):
 vs=expand(q)
 for c in [50,200/len(vs)]:
  z=[single(v,c) for v in vs]
  if min(x[0] for x in z)<45 or max(x[1] for x in z)>=55:return False
 return True
for k in pools:pools[k]=[q for q in pools[k] if oligo(q)]
unique=sorted({(side,q) for (side,g),rr in pools.items() for q in rr});oligos=[dict(id=f'O{i:04}',side=s,sequence=q,degeneracy=deg(q)) for i,(s,q) in enumerate(unique)];(work/'oligos.json').write_text(json.dumps(oligos,indent=2));(work/'pools.json').write_text(json.dumps({s+'_'+g:v for (s,g),v in pools.items()},indent=2));(work/'generation.json').write_text(json.dumps(diagnostics,indent=2))
for role,source in [('target','target_panel.fasta'),('background','background.fasta')]:
 d=work/role;(d/'sites').mkdir(exist_ok=True,parents=True);f=d/'background.fasta'
 if not f.exists():f.symlink_to(project/'data'/source)
 (d/'sites/oligos.tsv').write_text(''.join(f"{r['id']} {r['side']} {r['sequence']}\n" for r in oligos))
print('POOLS',[(k,len(v)) for k,v in pools.items()],'UNIQUE',len(unique),flush=True)
