# Escalera por programa y clase de tamano (<=30 / >30 filas tras fusionar), una instancia por familia y clase:
# la de mas filas de la relajacion dentro de la clase. Cada peldano contra R0 (SoPlex como lo usa Ibex).
import glob,os,re,sys,collections
S='/tmp/claude-1001/-home-iaraya/2bd51c7a-9daf-4c06-b6ef-b89011baa9b5/scratchpad'
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
T=60.0; excl={'alkyl','alkylbis','bearing'}
def familia(prog,n):
    if prog=='sol': return n.split('-')[0]
    n=re.sub(r'-1$','',n)
    if n.startswith('hs'): return n                                   # Hock-Schittkowski: cada uno es un problema distinto
    n=re.sub(r'(bis|ter)$','',n)
    if re.match(r'ex\d+_\d+_\d+',n): return n                     # GLOBALLib: cada exA_B_C es un problema (ex2_*, ex6_* son solo afines)
    return re.sub(r'[_-]?\d+$','',n) or n
filas={}
for f in os.listdir(M+'/filas_m2'):
    if f.count('.')!=3: continue
    k,n,lin,_=f.split('.')
    mm=re.search(r'([0-9.]+) por linealizacion',open(M+'/filas_m2/'+f).read())
    if mm: filas[(k,n,lin)]=float(mm.group(1))
bancos=[('escf','opt','art'),('escf','sol','compo'),('escfa','sol','art')]
R=collections.defaultdict(dict)    # (prog,lin,n,seed) -> {peldano: (t,fin)}
for d,prog,lin in bancos:
    for f in glob.glob('%s/%s/%s.*.csv'%(S,d,prog)):
        n,p,s,r=open(f).read().strip().split(',',3)
        if n in excl: continue
        v=r.split(); t=float(v[1]) if r!='NA' and len(v)>1 else T
        R[(prog,lin,n,s)][p]=(min(t,T), t<T*0.999)
# eleccion: por (prog, clase, familia) la (instancia, relajacion) con mas filas, entre las que tienen datos de R0
elig={}
for (prog,lin,n,s) in R:
    if (prog,n,lin) not in filas or 'R0' not in R[(prog,lin,n,s)]: continue
    m=filas[(prog,n,lin)]; c='<=30' if m<=30 else '>30'
    key=(prog,c,familia(prog,n))
    if key not in elig or m>elig[key][0]: elig[key]=(m,lin,n)
nom=[('R0f','+ merge twin rows'),('R1','+ SoPlex tuning'),('R2','bare simplex (SPxSolver)'),('R3','SC base'),('R6','SC + limit, devex, long-step'),('R7','hybrid (SC <=30 rows)')]
for prog in ('opt','sol'):
  for c in ('<=30','>30'):
    sel=[(v[1],v[2]) for k,v in elig.items() if k[0]==prog and k[1]==c]
    if not sel: continue
    print('== %s, %s filas: %d familias'%(prog,c,len(sel)))
    for p,nm in nom:
        out=[]
        for th in (1,5):
            pr=[]
            for (prog2,lin,n,s),v in R.items():
                if prog2!=prog or (lin,n) not in sel or p not in v or 'R0' not in v: continue
                (a,fa),(b,fb)=v[p],v['R0']
                if max(a,b)>th and (fa or fb): m=max(a,b); pr.append((a/m,b/m))
            if not pr: out.append('-'); continue
            A=sum(x for x,y in pr)/len(pr); B=sum(y for x,y in pr)/len(pr)
            out.append('%.2f:%.2f %.2f (%d)'%(A,B,A/B,len(pr)))
        if any(o!='-' for o in out): print('  %-30s >1s %-20s >5s %s'%(nm,*out))
if '-v' in sys.argv:
    for k in sorted(elig): print(k, elig[k])
