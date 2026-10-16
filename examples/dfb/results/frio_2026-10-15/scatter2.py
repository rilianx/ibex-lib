import csv,os,collections,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
S='/tmp/claude-1001/-home-iaraya/2bd51c7a-9daf-4c06-b6ef-b89011baa9b5/scratchpad'
EXC=set(tuple(l.split()[:2]) for l in open(os.path.expanduser('~/ibex2024/examples/dfb/results/excluidas_3600.txt')) if l.strip() and not l.startswith('#'))
def P(kind,r,T):
    if r=='NA': return (T,False)
    a=r.split(); t=float(a[1]); return (max(t,1e-3), t<T*0.999 and not (kind=='sol' and a[2]=='3'))
D=collections.defaultdict(dict)
for k,n,arm,s,r in list(csv.reader(open(os.path.expanduser('~/ibex2024/examples/dfb/results/final_2026-10-15.csv'))))[1:]:
    D[(k,n)][(arm,s)]=P(k,r,60)
L=collections.defaultdict(dict)
for f in os.listdir(M+'/l3600z'):
    b,k=f.split('.')[:2]; n,arm,sd,r=open(M+'/l3600z/'+f).read().strip().split(',',3)
    a2={'lejos':'x','sc':'x','caja':'caja'}[arm]; L[(b,k,n)][(a2,sd)]=P(k,r,3600)
def pts(src,filt,x,y):
    o=[]
    for key,d in src.items():
        if not filt(key): continue
        for (a,s),v in d.items():
            if a==x and (y,s) in d:
                w=d[(y,s)]
                if v[1] or w[1]: o.append((w[0],v[0],v[1],w[1]))
    return o
def panel(ax,p,T,title,xl,yl,lo=0.01):
    both=[q for q in p if q[2] and q[3]]
    ax.scatter([q[0] for q in both],[q[1] for q in both],s=9,alpha=.55,c='#1f6fb2',lw=0,label='terminan los dos (%d)'%len(both))
    ox=[q for q in p if q[2] and not q[3]]; oy=[q for q in p if q[3] and not q[2]]
    if ox: ax.scatter([T*1.25]*len(ox),[q[1] for q in ox],s=22,marker='>',c='#2a9d3a',label='solo %s termina (%d)'%(yl.split()[0],len(ox)))
    if oy: ax.scatter([q[0] for q in oy],[T*1.25]*len(oy),s=22,marker='^',c='#d1495b',label='solo %s termina (%d)'%(xl.split()[0],len(oy)))
    ax.plot([lo,T],[lo,T],c='k',lw=.8); ax.plot([lo,T],[lo/2,T/2],c='gray',lw=.6,ls='--'); ax.plot([lo,T],[lo*2,T*2],c='gray',lw=.6,ls='--')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lo,T*1.6); ax.set_ylim(lo,T*1.6)
    ax.axvline(T,c='gray',lw=.4); ax.axhline(T,c='gray',lw=.4)
    ax.set_title(title,fontsize=10); ax.set_xlabel(xl+' (s)'); ax.set_ylabel(yl+' (s)'); ax.legend(fontsize=7,loc='upper left'); ax.grid(alpha=.2)
sol=lambda k:k[0]=='sol'; opt=lambda k:k[0]=='opt' and k[1] not in {'alkyl','alkylbis','bearing'}
fig,ax=plt.subplots(2,2,figsize=(11,10.5))
panel(ax[0,0],pts(D,sol,'sc','cajacerca'),60,'Sistemas, 60 s: SC (orden propio) vs SoPlex «cerca»','SoPlex','SC')
panel(ax[0,1],pts(D,opt,'sc','caja'),60,'Opt, 60 s: SC (orden propio) vs SoPlex «lejos»','SoPlex','SC')
panel(ax[1,0],pts(D,sol,'cajacerca','caja'),60,'Sistemas, 60 s: SoPlex «cerca» vs SoPlex «lejos»','SoPlex-lejos','SoPlex-cerca')
f36=lambda k:(k[1],k[2]) not in EXC
panel(ax[1,1],pts(L,f36,'x','caja'),3600,'3600 s: sistemas mismo orden; opt SC propio vs SoPlex «lejos»','SoPlex','SC',lo=0.1)
fig.suptitle('CPU por corrida (diagonal: igual; punteadas: ×2 y ×½; bajo la diagonal gana el eje y)',fontsize=10)
fig.tight_layout(); fig.savefig(M+'/scatter2.png',dpi=130)
