I=dict(A='A',C='C',G='G',T='T',R='AG',Y='CT',S='CG',W='AT',K='GT',M='AC',B='CGT',D='AGT',H='ACT',V='ACG',N='ACGT');C=dict(zip('ACGT','TGCA'))
def fa(path):
 d={};key=None
 for l in path.open():
  if l.startswith('>'):key=l[1:].split()[0];d[key]=''
  else:d[key]+=l.strip().upper()
 return d
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def sets(q,reverse=False,gt=False):
 ss=[set(I[c]) for c in q]
 if gt:
  for x in ss:
   old=set(x)
   if 'G' in old:x.add('A')
   if 'T' in old:x.add('C')
 return [{C[c] for c in s} for s in ss[::-1]] if reverse else ss
def scan(s,q,reverse=False,k=2,gt=False):
 ss=sets(q,reverse,gt);m={c:sum(1<<i for i,x in enumerate(ss) if c in x) for c in 'ACGT'};state=[0]*(k+1);end=1<<(len(q)-1);out=[]
 for j,c in enumerate(s):
  if c not in m:state=[0]*(k+1);continue
  old=state;v=m[c];state=[((old[0]<<1)|1)&v]+[(((old[e]<<1)|1)&v)|((old[e-1]<<1)|1) for e in range(1,k+1)]
  if j+1>=len(q) and state[k]&end:
   pos=j-len(q)+1;mm=next(e for e in range(k+1) if state[e]&end);jj=range(3) if reverse else range(len(q)-3,len(q));last=any(s[pos+i] not in ss[i] for i in jj);out.append((pos,mm,last))
 return out
