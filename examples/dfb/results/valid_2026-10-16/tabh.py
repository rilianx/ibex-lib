# Validacion: 102 problemas nuevos de optim/all a 3600 s (h3600), con reemplazo desde fix/ (binario con arreglo afin). uso: tabh.py
import os,math,collections
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
T=3600.0
def rd(p):
    n,arm,s,r=open(p).read().strip().split(',',3)
    if r=='NA': return n,arm,s,None
    a=r.split(); t=float(a[1]); return n,arm,s,(max(t,1e-3),max(1.0,float(a[0])),t<T*0.999)
R=collections.defaultdict(dict); src={}
for f in os.listdir(M+'/h3600'):
    n,arm,s,x=rd(M+'/h3600/'+f); R[n][(arm,s)]=x
nf=0
for f in os.listdir(M+'/fix'):
    if f.startswith('h3600.'):
        n,arm,s,x=rd(M+'/fix/'+f); R[n][(arm,s)]=x; nf+=1
print('reemplazos fix:',nf)
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
for X,lab in (('sc','SC contra caja (best)'),('lejos','lejos contra caja (same order)')):
    for lim in (60,3600):
        P=[];fx=fr=0;sx=[];sr=[];na=[];nr=0
        for n in R:
            for (a,s),x in R[n].items():
                if a!=X: continue
                nr+=1; r=R[n].get(('caja',s))
                if x is None: na.append((n,X,s))
                if r is None: na.append((n,'caja',s))
                okx=bool(x and x[2] and x[0]<lim*0.999); okr=bool(r and r[2] and r[0]<lim*0.999)
                fx+=okx; fr+=okr
                if okx and okr: P.append((n,r[0],x[0]/r[0],x[1]/r[1]))
                elif okx: sx.append((n,s,round(x[0])))
                elif okr: sr.append((n,s,round(r[0])))
        print(f'{lab} a {lim}s: inst {len(R)} corridas {nr} terminan {fx}/{fr} pares {len(P)} cpu x{gm([p[2] for p in P]):.3f} celdas x{gm([p[3] for p in P]):.3f}')
        for th in (0.1,1,10,60,600):
            v=[p for p in P if p[1]>th]
            if v: print(f'   caja>{th}s n={len(v)} t x{gm([p[2] for p in v]):.2f} c x{gm([p[3] for p in v]):.2f}')
        if lim==3600:
            print('   NA:',sorted(set(na)))
            print('   solo',X,':',sx); print('   solo caja:',sr)
            fam=collections.defaultdict(list)
            for p in P:
                if p[1]>1: fam[p[0]].append(p[2])
            z=sorted((round(gm(v),2),n) for n,v in fam.items()); print('   mejores',z[:5],'\n   peores',z[-6:])
