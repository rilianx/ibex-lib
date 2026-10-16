# Fila de optimizacion de tab:long sobre el banco unico (240): instancias donde sc o caja llega a 60 s en alguna semilla. uso: tablong.py
import os,math,collections
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
def rd(p):
    n,a,s,r=open(p).read().strip().split(',',3)
    if r=='NA': return n,a,s,None
    x=r.split(); return n,a,s,(max(float(x[1]),1e-3),max(1.0,float(x[0])))
L=collections.defaultdict(dict)   # 3600 s (c7200 aparte)
for d,pre in (('h3600',''),('fix','h3600.')):
    for f in os.listdir(M+'/'+d):
        if f.startswith(pre) and '.robot.' not in f: n,a,s,x=rd(M+'/'+d+'/'+f); L[n][(a,s)]=x
for f in os.listdir(M+'/l3600z'):
    if f.startswith('l.opt.'): n,a,s,x=rd(M+'/l3600z/'+f); L[n][(a,s)]=x
X7={}
for f in os.listdir(M+'/c7200'): n,a,s,x=rd(M+'/c7200/'+f); X7[(n,s)]=x
# instancias con alguna corrida >=60 s (sc o caja) en el banco de 60 s
U=collections.defaultdict(dict)
for f in os.listdir(M+'/u60'):
    if f.startswith('opt.'):
        n,a,s,x=rd(M+'/u60/'+f); U[n][(a,s)]=x
dif=sorted(n for n in U if any((x is None or x[0]>=59.94) for (a,s),x in U[n].items() if a in('sc','caja')))
print('dificiles',len(dif),dif)
miss=[n for n in dif if n not in L]; print('sin 3600:',miss)
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
P=[];fx=fr=fr7=nr=0
for n in dif:
    for (a,s),x in L.get(n,{}).items():
        if a!='sc': continue
        r=L[n].get(('caja',s)); nr+=1
        okx=x and x[0]<3596; okr=r and r[0]<3596
        r7=X7.get((n,s)) if not okr else r
        fx+=bool(okx); fr+=bool(okr); fr7+=bool(r7 and r7[0]<7190)
        if okx and r7 and r7[0]<7190: P.append((n,x[0]/r7[0],x[1]/r7[1]))
print(f'inst {len(dif)} runs {nr} sc fin {fx} caja fin3600 {fr} caja fin7200 {fr7} pares {len(P)} cells x{gm([p[2] for p in P]):.2f} cpu x{gm([p[1] for p in P]):.2f}')
for ex in (['bearing'],):
    v=[p for p in P if p[0] not in ex]; print(' sin',ex,len(v),f'cells x{gm([p[2] for p in v]):.2f} cpu x{gm([p[1] for p in v]):.2f}')
fam=collections.defaultdict(list)
for p in P: fam[p[0]].append(p[1])
print(sorted((round(gm(v),2),n) for n,v in fam.items()))
