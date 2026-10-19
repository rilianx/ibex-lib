# Caja directa: comparaciones a 60 s sobre el banco unico. uso: tabd.py dir
import os,sys,math,collections
D=sys.argv[1]; T=60.0; EXC={'alkyl','alkylbis','bearing'}
R=collections.defaultdict(dict)
for f in os.listdir(D):
    kind=f.split('.')[0]; n,arm,seed,r=open(os.path.join(D,f)).read().strip().split(',',3)
    if r=='NA': R[(kind,n)][(arm,seed)]=None; continue
    a=r.split(); t=float(a[1]); c=max(1.0,float(a[0]))
    R[(kind,n)][(arm,seed)]=(max(t,1e-3),c,t<T*0.999 and not (kind=='sol' and a[2]=='3'))
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
def line(kind,x,ref,show=False):
    P=[];fx=fr=0;solx=solr=0
    for (k,n),d in R.items():
        if k!=kind: continue
        vx=[d.get((x,s)) for (a,s) in d if a==ref]; 
        for (a,s),r in d.items():
            if a!=ref or (x,s) not in d: continue
            y=d[(x,s)]; fx+=bool(y and y[2]); fr+=bool(r and r[2])
            if kind=='opt' and n in EXC: continue
            if y and r and y[2] and r[2]: P.append((n,r[0],y[0]/r[0],y[1]/r[1]))
    o=[]
    for th in (0.1,1,5):
        v=[p for p in P if p[1]>th]; o.append(f'{gm([p[3] for p in v]):.3f}/{gm([p[2] for p in v]):.3f}(n={len(v)})')
    print(f'{kind} {x:10s}/{ref:10s} celdas/cpu >0.1 {o[0]}  >1 {o[1]}  >5 {o[2]}  terminan {fx}/{fr}')
    if show:
        fam=collections.defaultdict(list)
        for p in P:
            if p[1]>1: fam[p[0]].append(p[2])
        z=sorted((round(gm(v),2),n) for n,v in fam.items()); print('    mejores',z[:5],'\n    peores',z[-5:])
for x,r in (('sccerca','cajacercap'),('sccerca','cajacercad'),('cajacercap','cajacercad')): line('opt',x,r,x=='sccerca' and r=='cajacercap')
for x,r in (('sc','cajacercap'),('sccerca','cajacercap'),('sc','cajacercad'),('cajacercap','cajacercad'),('lejos','cajad')): line('sol',x,r,x in ('sc','sccerca') and r=='cajacercap')
