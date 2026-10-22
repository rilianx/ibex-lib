# Tablas del paper con la configuracion final (SC con tope 60 = brazos c_cap_*), 2026-10-22.
import os,re,csv,math,collections
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
EXC={'alkyl','alkylbis','bearing'}
def P(k,r):
    if r=='NA': return None
    x=r.split(); t=float(x[1]); return (max(t,1e-3),max(1.,float(x[0])),t<59.94 and not(k=='sol' and x[2]=='3'))
D=collections.defaultdict(dict)
for k,n,a,s,r in list(csv.reader(open('results/final_2026-10-21.csv')))[1:]: D[(k,n)][(a,s)]=P(k,r)
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
def cmp(k,x,ref):
    pairs=[]
    for (kk,n),d in D.items():
        if kk!=k or (k=='opt' and n in EXC): continue
        for (a,s),R in d.items():
            if a!=ref or not R: continue
            X=d.get((x,s))
            if X and R[2] and X[2]: pairs.append((n,R[0],X[0]/R[0],X[1]/R[1]))
    return [(gm([p[3] for p in pairs if p[1]>th]),gm([p[2] for p in pairs if p[1]>th]),len([p for p in pairs if p[1]>th])) for th in (0.1,1,5)]
def fin(k,a):
    runs=0; solved=0
    for (kk,n),d in D.items():
        if kk!=k: continue
        v=[X for (aa,s),X in d.items() if aa==a]
        runs+=sum(1 for X in v if X and X[2])
        if v and all(X and X[2] for X in v): solved+=1
    return runs,solved
def row(k,x,ref,lab):
    o=cmp(k,x,ref); fx=fin(k,x); fr=fin(k,ref)
    cs=[c for c,_,_ in o]
    f=f'{fx[1]} / {fr[1]}' if k=='opt' else f'{fx[0]} / {fr[0]}'
    print(f'{k} {lab:28s} cells {min(cs):.2f}-{max(cs):.2f}  cpu {o[0][1]:.2f} {o[1][1]:.2f} {o[2][1]:.2f}  fin {f}  (runs {fx[0]}/{fr[0]}) n={o[0][2]}/{o[1][2]}/{o[2][2]}')
for l in ('art','compo'):
    print('==',l)
    for k in ('opt','sol'):
        row(k,f'c_cap_{l}',f'y_caja_{l}','SC / class')
        row(k,f'c_cap_{l}',f's_spx_{l}','SC / SPx')
        row(k,f'h_hib_{l}',f'y_caja_{l}','hybrid / class')
        row(k,f'h_hib_{l}',f's_spx_{l}','hybrid / SPx')
        row(k,f'h_hib_{l}',f'c_cap_{l}','hybrid / SC')
        row(k,f's_spx_{l}',f'y_caja_{l}','SPx / class')
print('== compo / afin')
for k in ('opt','sol'):
    for e in ('c_cap','y_caja','s_spx','h_hib'): row(k,f'{e}_compo',f'{e}_art',f'{e} compo/art')
    row(k,'c_cap_art','y_caja_compo','SC art / class compo'); row(k,'h_hib_compo','c_cap_art','hyb compo / SC art')
# tamano
F={}
for f in os.listdir(M+'/filas_m2'):
    k,rest=f.split('.',1); n,l=rest[:-4].rsplit('.',1)
    t=open(M+'/filas_m2/'+f).read(); m=re.search(r'filas de la relajacion: ([\d.]+)',t)
    if m: F[(k,n,l)]=float(m[1])
bins=[(0,15),(15,30),(30,60),(60,120),(120,10**9)]
print('== tamano (SC/ref cpu, instancias con ref >1 s)')
for ref in ('y_caja','s_spx'):
  for k in ('opt','sol'):
    for l in ('art','compo'):
        B=collections.defaultdict(list)
        for (kk,n),d in D.items():
            if kk!=k or (k=='opt' and n in EXC): continue
            m=F.get((k,n,l))
            if m is None: continue
            v=[d[(f'c_cap_{l}',s)][0]/d[(f'{ref}_{l}',s)][0] for s in '1234' if d.get((f'c_cap_{l}',s)) and d.get((f'{ref}_{l}',s)) and d[(f'c_cap_{l}',s)][2] and d[(f'{ref}_{l}',s)][2] and d[(f'{ref}_{l}',s)][0]>1]
            if v:
                for a,b in bins:
                    if a<=m<b: B[(a,b)].append(gm(v))
        print(f'{ref:6s} {k} {l:5s} '+' '.join(f'{gm(B[bb]):5.2f} ({len(B[bb]):3d})' for bb in bins))
# filas medianas
for k in ('opt','sol'):
    for l in ('art','compo'):
        v=sorted(m for (kk,n,ll),m in F.items() if kk==k and ll==l); print('filas',k,l,'mediana',v[len(v)//2],'cuartiles',v[len(v)//4],v[3*len(v)//4],'n',len(v))
