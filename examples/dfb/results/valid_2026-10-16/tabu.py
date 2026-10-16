# Validacion (102 nuevos) y banco unido (paper 139 + nuevos), SC/lejos contra caja, a 60 s y a 3600 s. uso: tabu.py
# Nuevos: h3600 + fix/h3600.* + c7200 (caja a 7200 s en los 4 que solo SC terminaba). Paper: f60 (+wall de fix) y l3600z (l.opt).
import os,math,collections
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
EXC=set(l.split()[1] for l in open('/home/iaraya/ibex2024/examples/dfb/results/excluidas_3600.txt') if len(l.split())>1)
def rd(p,T):
    n,arm,s,r=open(p).read().strip().split(',',3)
    if r=='NA': return n,arm,s,None
    a=r.split(); t=float(a[1]); return n,arm,s,(max(t,1e-3),max(1.0,float(a[0])),t<T*0.999)
NEW=collections.defaultdict(dict)
for d,pre,T in (('h3600','',3600),('fix','h3600.',3600),('c7200','',7200)):
    for f in os.listdir(M+'/'+d):
        if f.startswith(pre):
            n,a,s,x=rd(M+'/'+d+'/'+f,T); NEW[n][(a,s)]=x
PAP=collections.defaultdict(dict)
for f in os.listdir(M+'/f60'):
    if f.startswith('opt.') and f.split('.')[1] not in EXC:
        n,a,s,x=rd(M+'/f60/'+f,60.0)
        if a in('sc','caja','lejos'): PAP[n][(a,s)]=x
L=collections.defaultdict(dict)
for f in os.listdir(M+'/l3600z'):
    if f.startswith('l.opt.'):
        n,a,s,x=rd(M+'/l3600z/'+f,3600); L[n][(a,s)]=x
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
def get(B,n,a,s,lim):
    # a 3600 s el banco del paper usa l3600z donde existe, si no f60 (lo que no termino a 60 s queda sin termino)
    if B is PAP and lim>60 and (a,s) in L.get(n,{}): return L[n][(a,s)]
    return B[n].get((a,s))
def tabla(lab,bancos,X,lim):
    P=[];fx=fr=nr=0;sx=[];sr=[]
    for B in bancos:
        for n in B:
            for (a,s) in list(B[n]):
                if a!=X: continue
                x=get(B,n,X,s,lim); r=get(B,n,'caja',s,lim); nr+=1
                ok=lambda y,cap: bool(y and y[2] and y[0]<cap*0.999)
                okx=ok(x,lim); okr=ok(r,lim if lim==60 else 7200)
                fx+=okx; fr+=okr
                if okx and okr: P.append((n,r[0],x[0]/r[0],x[1]/r[1]))
                elif okx: sx.append(n)
                elif okr: sr.append(n)
    s=f'{lab:34s} {X:5s} {lim:4d}s corridas {nr} terminan {fx}/{fr} pares {len(P)} t x{gm([p[2] for p in P]):.2f} c x{gm([p[3] for p in P]):.2f}'
    for th in (1,10,60,600):
        v=[p for p in P if p[1]>th]
        if v: s+=f' | >{th}s n={len(v)} t x{gm([p[2] for p in v]):.2f} c x{gm([p[3] for p in v]):.2f}'
    print(s)
    if sx or sr: print('    solo',X,sorted(set(sx)),'solo caja',sorted(set(sr)))
    return P
for X in ('sc','lejos'):
    for lim in (60,3600):
        if X=='lejos' and lim==3600:
            tabla('validacion (102)',[NEW],X,lim); continue
        tabla('validacion (102)',[NEW],X,lim)
        tabla('paper (139-excl)',[PAP],X,lim)
        P=tabla('unido',[PAP,NEW],X,lim)
        if X=='sc' and lim==3600:
            fam=collections.defaultdict(list)
            for p in P:
                if p[1]>1: fam[p[0]].append(p[2])
            z=sorted((round(gm(v),2),n) for n,v in fam.items()); print('    peores',z[-6:])
print('instancias: nuevos',len(NEW),'paper',len(PAP))
