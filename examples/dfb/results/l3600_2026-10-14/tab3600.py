# Tabla 3600 s: SC (l3600f, brazo x) contra caja (S/l3600, S/all3600). Excluye results/excluidas_3600.txt
import os,sys,math,collections
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
S='/tmp/claude-1001/-home-iaraya/2bd51c7a-9daf-4c06-b6ef-b89011baa9b5/scratchpad'
T=3600.0; X=sys.argv[1] if len(sys.argv)>1 else 'sc'
REFSOL=sys.argv[2] if len(sys.argv)>2 else 'caja'; CC=collections.defaultdict(dict)
EXC=set()
for l in open(os.path.expanduser('~/ibex2024/examples/dfb/results/excluidas_3600.txt')):
    if l.startswith('#') or not l.strip(): continue
    k,n=l.split()[:2]; EXC.add((k,n))
def parse(path,kind):
    if not os.path.getsize(path): return None
    n,arm,seed,r=open(path).read().strip().split(',',3)
    if r=='NA': return n,arm,seed,None
    a=r.split(); t=float(a[1]); c=max(1.0,float(a[0]))
    return n,arm,seed,(max(t,1e-3),c,t<T*0.999 and not (kind=='sol' and a[2]=='3'))
R=collections.defaultdict(dict)   # (grupo,kind,n) -> {(arm,seed):r}
for f in os.listdir(M+'/l3600f'):
    if not f.endswith('.csv'): continue
    b,kind=f.split('.')[:2]
    n,arm,seed,r=parse(M+'/l3600f/'+f,kind)
    if kind=='sol' and REFSOL!='caja' and arm==REFSOL: CC[(b,kind,n)][('caja',seed)]=r; continue
    if arm!=X: continue
    R[(b,kind,n)][('x',seed)]=r
for b,d in (('l','l3600'),('all','all3600')):
    for f in os.listdir(S+'/'+d):
        kind=f.split('.')[0]
        n,arm,seed,r=parse(S+'/'+d+'/'+f,kind)
        if arm=='caja' and (b,kind,n) in R: R[(b,kind,n)][('caja',seed)]=r
if REFSOL!='caja':
    for k,d in R.items():
        if k[1]=='sol':
            for kk in [q for q in d if q[0]=='caja']: del d[kk]
            d.update(CC.get(k,{}))
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
for g,lab in ((('l','sol'),'sistemas paper'),(('l','opt'),'opt paper'),(('all','opt'),'opt all')):
    ks=[k for k in R if k[:2]==g and (k[1],k[2]) not in EXC]
    pairs=[];fx=fr=0;nax=nar=0;sx=sr=0;both=0
    for k in ks:
        d=R[k]; seeds=sorted(s for (a,s) in d if a=='x')
        for s in seeds:
            x=d.get(('x',s)); r=d.get(('caja',s))
            nax+= x is None; nar+= r is None
            fx+= bool(x and x[2]); fr+= bool(r and r[2])
            if x and r and x[2] and r[2]: pairs.append((k[2],r[0],x[0]/r[0],x[1]/r[1]))
            elif x and x[2] and not (r and r[2]): sx+=1
            elif r and r[2] and not (x and x[2]): sr+=1
    o=[]
    for th in (0.1,1,5,60,600):
        v=[p for p in pairs if p[1]>th]; o.append('%3d c%.2f t%.2f'%(len(v),gm([p[3] for p in v]),gm([p[2] for p in v])))
    print('%-15s inst %d corr %d  terminan x %d caja %d (solo x %d, solo caja %d, NA x %d caja %d)'%(lab,len(ks),sum(1 for k in ks for a,s in R[k] if a=='x'),fx,fr,sx,sr,nax,nar))
    print('   >0.1 / >1 / >5 / >60 / >600 :', ' | '.join(o))
    fam=collections.defaultdict(list)
    for p in pairs:
        if p[1]>1: fam[p[0]].append(p[2])
    z=sorted((gm(v),n) for n,v in fam.items())
    print('   mejores:',', '.join('%s x%.2f'%(n,g) for g,n in z[:5]))
    print('   peores :',', '.join('%s x%.2f'%(n,g) for g,n in z[-5:]))
    miss=[k[2] for k in ks if not any(a=='caja' for a,s in R[k])]
    if miss: print('   sin caja:',miss)
print('--- corridas no terminadas en algun brazo')
for k in sorted(R):
    if (k[1],k[2]) in EXC: continue
    d=R[k]
    for s in sorted(s for a,s in d if a=='x'):
        x=d.get(('x',s)); r=d.get(('caja',s))
        if not (x and x[2] and r and r[2]):
            f=lambda z: 'NA' if z is None else ('%.0fs%s'%(z[0],'' if z[2] else '*'))
            print('  ',k,s,'x',f(x),'caja',f(r) if ('caja',s) in d else '-')
