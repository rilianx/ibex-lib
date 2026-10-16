# Tabla 3600 s desde l3600z (reintento en los dos motores). uso: tab3600z.py
import os,math,collections
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
T=3600.0
EXC=set(tuple(l.split()[:2]) for l in open(os.path.expanduser('~/ibex2024/examples/dfb/results/excluidas_3600.txt')) if l.strip() and not l.startswith('#'))
R=collections.defaultdict(dict)
for f in os.listdir(M+'/l3600z'):
    b,kind=f.split('.')[:2]
    n,arm,s,r=open(M+'/l3600z/'+f).read().strip().split(',',3)
    if r=='NA': x=None
    else:
        a=r.split(); t=float(a[1]); x=(max(t,1e-3),max(1.0,float(a[0])),t<T*0.999 and not (kind=='sol' and a[2]=='3'))
    R[(b,kind,n)][(arm,s)]=x
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
for g,X,lab in ((('l','sol'),'lejos','sistemas, mismo orden'),(('l','opt'),'sc','opt paper'),(('all','opt'),'sc','opt all')):
    ks=[k for k in R if k[:2]==g and (k[1],k[2]) not in EXC]
    P=[];fx=fr=0;sx=[];sr=[];na=0;nruns=0
    for k in ks:
        for (a,s),x in R[k].items():
            if a!=X: continue
            nruns+=1; r=R[k].get(('caja',s))
            if x is None or r is None: na+=1
            fx+=bool(x and x[2]); fr+=bool(r and r[2])
            if x and r and x[2] and r[2]: P.append((k[2],r[0],x[0]/r[0],x[1]/r[1]))
            elif x and x[2]: sx.append((k[2],s,round(x[0])))
            elif r and r[2]: sr.append((k[2],s,round(r[0])))
    print(f'{lab}: inst {len(ks)} corridas {nruns} terminan {fx}/{fr} NA {na} pares {len(P)} celdas x{gm([p[3] for p in P]):.3f} cpu x{gm([p[2] for p in P]):.3f}')
    for th in (1,60,600):
        v=[p for p in P if p[1]>th]; print(f'   >{th}s n={len(v)} c x{gm([p[3] for p in v]):.2f} t x{gm([p[2] for p in v]):.2f}')
    print('   solo SC:',sx); print('   solo caja:',sr)
    fam=collections.defaultdict(list)
    for p in P:
        if p[1]>1: fam[p[0]].append(p[2])
    z=sorted((round(gm(v),2),n) for n,v in fam.items()); print('   mejores',z[:4],'peores',z[-4:])
