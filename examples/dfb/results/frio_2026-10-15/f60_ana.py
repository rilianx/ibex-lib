import os,math,sys,csv,collections
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
excl=set(l.split()[1] for l in open('/home/iaraya/ibex2024/examples/dfb/results/excluidas_3600.txt') if len(l.split())>1)
def parse(kind,r):
    a=r.split()
    if a[0]=='NA': return None
    t=float(a[1]); c=max(1.0,float(a[0]))
    ok=t<59.94 and not (kind=='sol' and a[2]=='3')
    return (max(t,1e-3),c,ok)
def leerdir(d):
    R={}
    for f in os.listdir(d):
        if not f.endswith('.csv'): continue
        kind=f.split('.')[0]
        n,arm,s,r=open(d+'/'+f).read().strip().split(',',3)
        R[(kind,n,arm,s)]=parse(kind,r)
    return R
N=leerdir(M+'/f60')
O={}
for row in csv.DictReader(open('/home/iaraya/ibex2024/examples/dfb/results/final_2026-10-14.csv')):
    O[(row['tipo'],row['instancia'],row['brazo'],row['semilla'])]=parse(row['tipo'],row['resultado'])
O.update({k:v for k,v in leerdir(M+'/ordcaja').items() if k[2]=='cajacerca'})
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
arms=sorted(set((k[0],k[2]) for k in N))
for kind,arm in arms:
    ks=[k for k in N if k[0]==kind and k[2]==arm and k[1] not in excl and k in O]
    tN=sum(1 for k in ks if N[k] and N[k][2]); tO=sum(1 for k in ks if O[k] and O[k][2])
    both=[k for k in ks if N[k] and O[k] and N[k][2] and O[k][2]]
    iden=sum(1 for k in both if N[k][1]==O[k][1])
    s=f'{kind} {arm:10s} pares {len(ks):4d} terminan {tO}->{tN} celdas iguales {iden}/{len(both)}'
    for th in (0.1,1,5):
        v=[k for k in both if O[k][0]>th]
        s+=f' | >{th}: n={len(v)} t x{gm([N[k][0]/O[k][0] for k in v]):.3f} c x{gm([N[k][1]/O[k][1] for k in v]):.3f}'
    print(s)
    ch=[k for k in ks if bool(N[k] and N[k][2])!=bool(O[k] and O[k][2])]
    for k in ch: print('   cambia termino',k,O[k],N[k])
    big=sorted((N[k][1]/O[k][1],k) for k in both if O[k][0]>1 and N[k][1]!=O[k][1])
    for r,k in big[:3]+big[-3:]: print(f'   celdas x{r:.2f}',k,'t %.1f->%.1f'%(O[k][0],N[k][0]))
