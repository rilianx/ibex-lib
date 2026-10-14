import csv,os,math,collections
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
T=60.0
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
def P(kind,r):
    if r=='NA': return None
    a=r.split(); t=float(a[1]); return (max(t,1e-3),max(1.0,float(a[0])),t<T*0.999 and not (kind=='sol' and a[2]=='3'))
D=collections.defaultdict(dict)
for k,n,arm,s,r in list(csv.reader(open(os.path.expanduser('~/ibex2024/examples/dfb/results/final_2026-10-14.csv'))))[1:]:
    if k=='sol': D[n][(arm,s)]=P(k,r)
for f in os.listdir(M+'/ordcaja'):
    if f.startswith('sol.') and 'cajacerca' in f:
        n,arm,s,r=open(M+'/ordcaja/'+f).read().strip().split(',',3); D[n][(arm,s)]=P('sol',r)
def cmp(x,ref):
    pairs=[]
    for n,d in D.items():
        for s in ('1','2'):
            a=d.get((x,s)); b=d.get((ref,s))
            if a and b and a[2] and b[2]: pairs.append((n,b[0],a[0]/b[0],a[1]/b[1]))
    o=[]
    for th in (0.1,1,5):
        v=[p for p in pairs if p[1]>th]; o.append('n=%d c%.3f t%.3f'%(len(v),gm([p[3] for p in v]),gm([p[2] for p in v])))
    fin=lambda z: sum(1 for d in D.values() for s in '12' if d.get((z,s)) and d[(z,s)][2])
    print('%-9s / %-9s'%(x,ref),' | '.join(o),' terminan %d/%d'%(fin(x),fin(ref)))
    return pairs
p=cmp('cajacerca','caja')
fam=collections.defaultdict(list)
for q in p:
    if q[1]>1: fam[q[0]].append((q[2],q[3]))
z=sorted((gm([a for a,b in v]),gm([b for a,b in v]),n) for n,v in fam.items())
print('mejores',', '.join('%s t%.2f c%.2f'%(n,a,b) for a,b,n in z[:8]))
print('peores ',', '.join('%s t%.2f c%.2f'%(n,a,b) for a,b,n in z[-6:]))
for x in ('sc','lejos'):
    for r in ('caja','cajacerca'): cmp(x,r)
print('--- terminan solo con uno')
for n,d in sorted(D.items()):
    for s in '12':
        a=d.get(('cajacerca',s)); b=d.get(('caja',s))
        if bool(a and a[2])!=bool(b and b[2]): print('  ',n,s,'cerca',a and round(a[0],1),'lejos',b and round(b[0],1))
