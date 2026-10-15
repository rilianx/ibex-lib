# Perfiles de desempeno SC (orden propio) contra la caja. uso: perfil.py brazo_bb_sol brazo_bb_opt [outdir]
import csv,os,sys,math,collections
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
T=60.0; EXC={'alkyl','alkylbis','bearing'}
def P(kind,r):
    if r=='NA': return None
    a=r.split(); t=float(a[1]); return (max(t,1e-3), t<T*0.999 and not (kind=='sol' and a[2]=='3'))
D=collections.defaultdict(dict)
for k,n,arm,s,r in list(csv.reader(open(os.path.expanduser('~/ibex2024/examples/dfb/results/final_2026-10-15.csv'))))[1:]:
    D[(k,n)][(arm,s)]=P(k,r)
for f in []:
    if f.endswith('.csv'):
        k=f.split('.')[0]; n,arm,s,r=open(M+'/ordcaja/'+f).read().strip().split(',',3); D[(k,n)][(arm,s)]=P(k,r)
grid=[1+0.01*i for i in range(100)]+[2+0.05*i for i in range(61)]+[5+0.25*i for i in range(21)]
def perfil(kind,bb,mode):
    R=[]
    for (k,n),d in D.items():
        if k!=kind or (kind=='opt' and n in EXC): continue
        for (a,s),x in d.items():
            if a!='sc': continue
            y=d.get((bb,s))
            if (bb,s) not in d: continue
            tx=x[0] if x and x[1] else math.inf; ty=y[0] if y and y[1] else math.inf
            if tx==ty==math.inf: continue
            m=min(tx,ty)
            ref={'min':m,'bb':y[0] if y else 0}[mode]
            if ref<=0.1: continue
            R.append((tx/m,ty/m))
    return R
mode=sys.argv[4] if len(sys.argv)>4 else 'min'
for kind,bb in (('sol',sys.argv[1]),('opt',sys.argv[2])):
    R=perfil(kind,bb,mode); N=len(R)
    f=lambda i,t: sum(1 for r in R if r[i]<=t+1e-12)/N
    print(kind,bb,'runs',N,'tau1 sc %.4f bb %.4f  tau2 sc %.3f bb %.3f  bb no termina %.3f sc no termina %.3f'%(f(0,1),f(1,1),f(0,2),f(1,2),sum(r[1]==math.inf for r in R)/N,sum(r[0]==math.inf for r in R)/N))
    if len(sys.argv)>3 and sys.argv[3]!='-':
        with open(os.path.join(sys.argv[3],'perfil_%s.dat'%kind),'w') as o:
            o.write('tau sc bb\n')
            for t in grid: o.write('%.3f %.4f %.4f\n'%(t,f(0,t),f(1,t)))
