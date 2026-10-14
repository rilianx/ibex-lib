import csv,os,math,collections
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
T=60.0
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
def P(kind,r):
    if r=='NA': return None
    a=r.split(); t=float(a[1]); return (max(t,1e-3),max(1.0,float(a[0])),t<T*0.999 and not (kind=='sol' and a[2]=='3'))
ref={}
for k,n,arm,s,r in list(csv.reader(open(os.path.expanduser('~/ibex2024/examples/dfb/results/final_2026-10-14.csv'))))[1:]:
    if arm=='caja' and s=='1': ref[(k,n)]=P(k,r)
alt=collections.defaultdict(dict)
for f in os.listdir(M+'/ordcaja'):
    if not f.endswith('.csv'): continue
    k=f.split('.')[0]; n,arm,s,r=open(M+'/ordcaja/'+f).read().strip().split(',',3)
    alt[arm][(k,n)]=P(k,r)
for arm in sorted(alt):
  for k in ('sol','opt'):
    v=[];vc=[];w=[];fin=[0,0]
    for (kk,n),x in alt[arm].items():
        if kk!=k: continue
        r=ref.get((k,n)); fin[0]+=bool(x and x[2]); fin[1]+=bool(r and r[2])
        if x and r and x[2] and r[2]: v.append(x[0]/r[0]); vc.append(x[1]/r[1]); w.append((x[0]/r[0],n))
    w.sort()
    print('%s %-10s n=%d  t %.3f  celdas %.3f  terminan %d/%d  gana>10%% %d pierde>10%% %d'%(k,arm,len(v),gm(v),gm(vc),fin[0],fin[1],sum(a<0.9 for a in v),sum(a>1.1 for a in v)))
    print('    mejores',', '.join('%s %.2f'%(n,a) for a,n in w[:5]),' | peores',', '.join('%s %.2f'%(n,a) for a,n in w[-4:]))
