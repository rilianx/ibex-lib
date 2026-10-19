# SoPlex sin presolve (d60p/cajacercap) contra con presolve (d60s/cajacercad), 60 s, mismo banco y semillas.
import os,math,collections
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
EXC=set(l.split()[1] for l in open('/home/iaraya/ibex2024/examples/dfb/results/excluidas_3600.txt') if len(l.split())>1)
def rd(p):
    n,a,s,r=open(p).read().strip().split(',',3)
    if r=='NA': return None
    x=r.split(); t=float(x[1]); return (max(t,1e-3),max(1.0,float(x[0])),t<59.94)
D=collections.defaultdict(dict)
for d,arm in (('d60p','cajacercap'),('d60s','cajacercad')):
    for f in os.listdir(M+'/'+d):
        b,n,a,s,_=f.split('.') if f.count('.')==4 else (None,)*5
        if a!=arm or n in EXC: continue
        D[(b,n,s)][arm]=rd(M+'/'+d+'/'+f)
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
for b in ('opt','sol'):
    K=[k for k in D if k[0]==b and len(D[k])==2]
    P=[];fp=fd=0;solo_p=set();solo_d=set()
    for k in K:
        p,d=D[k]['cajacercap'],D[k]['cajacercad']
        okp=bool(p and p[2]); okd=bool(d and d[2]); fp+=okp; fd+=okd
        if okp and okd: P.append((k[1],d[0],p[0]/d[0],p[1]/d[1]))
        elif okp: solo_p.add(k[1])
        elif okd: solo_d.add(k[1])
    print(f'{b}: corridas {len(K)} instancias {len(set(k[1] for k in K))} terminan sin/con presolve {fp}/{fd} pares {len(P)} t x{gm([q[2] for q in P]):.3f} celdas x{gm([q[3] for q in P]):.3f}')
    for th in (1,10):
        v=[q for q in P if q[1]>th]
        print(f'   >{th}s n={len(v)} t x{gm([q[2] for q in v]):.3f} celdas x{gm([q[3] for q in v]):.3f}  gana sin={sum(q[2]<0.95 for q in v)} gana con={sum(q[2]>1.05 for q in v)}')
    print('   solo sin presolve:',sorted(solo_p)); print('   solo con presolve:',sorted(solo_d))
    I=collections.defaultdict(list)
    for q in P:
        if q[1]>1: I[q[0]].append(q[2])
    r=sorted((gm(v),n) for n,v in I.items())
    print('   extremos:',[(n,round(x,2)) for x,n in r[:4]],'...',[(n,round(x,2)) for x,n in r[-4:]])
