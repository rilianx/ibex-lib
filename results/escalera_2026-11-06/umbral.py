# Umbral del hibrido: por instancia (filas medias tras fusionar), SC base (R3) si filas <= T, simplex directo (R2) si no.
# Puntaje contra R0 (A/max : B/max), corridas > 1 s, y tambien media geometrica de t/t_R0.
import glob,os,re,math,collections,sys
import os; S=os.environ.get('S',os.path.join(os.path.dirname(os.path.abspath(__file__)),'work/S'))
M=os.environ.get('M',os.path.join(os.path.dirname(os.path.abspath(__file__)),'work/M'))
T=60.0; excl={'alkyl','alkylbis','bearing'}
filas={}
for f in os.listdir(M+'/filas_m2'):
    if f.count('.')!=3: continue
    k,n,lin,_=f.split('.')
    mm=re.search(r'([0-9.]+) por linealizacion',open(M+'/filas_m2/'+f).read())
    if mm: filas[(k,n,lin)]=float(mm.group(1))
runs=[]
for d,prog,lin in [('escf','opt','art'),('escf','sol','compo'),('escfa','sol','art')]:
    D=collections.defaultdict(dict)
    for f in glob.glob('%s/%s/%s.*.csv'%(S,d,prog)):
        n,p,s,r=open(f).read().strip().split(',',3)
        if n in excl or p not in ('R0','R2','R3'): continue
        v=r.split(); t=float(v[1]) if r!='NA' and len(v)>1 else T
        D[(n,s)][p]=min(max(t,1e-3),T)
    for (n,s),v in D.items():
        if len(v)==3 and (prog,n,lin) in filas: runs.append((prog,lin,n,filas[(prog,n,lin)],v['R0'],v['R2'],v['R3']))
def puntaje(sel,Tu):
    pr=[];g=[]
    for prog,lin,n,m,r0,r2,r3 in sel:
        h=r3 if m<=Tu else r2
        if max(h,r0)>1 and min(h,r0)<T*0.999:
            mx=max(h,r0); pr.append((h/mx,r0/mx))
        if r0>1: g.append(math.log(h/r0))
    A=sum(a for a,b in pr)/len(pr); B=sum(b for a,b in pr)/len(pr)
    return A/B, math.exp(sum(g)/len(g)), len(pr)
grupos={'todo':runs,'opt':[r for r in runs if r[0]=='opt'],'sol':[r for r in runs if r[0]=='sol']}
cands=[0,5,10,15,20,25,30,35,40,45,50,60,70,80,100,120,150,1e9]
print('%-6s'%'T'+''.join('%-26s'%g for g in grupos))
for Tu in cands:
    print('%-6s'%('inf' if Tu>1e8 else Tu)+''.join('%.3f  gm %.3f (%d)       '%puntaje(sel,Tu) for sel in grupos.values()))
# oraculo por instancia y los cruces SC vs simplex directo por filas
print('\nSC base / simplex directo por tramo de filas (media geom., R0>1s):')
tr=[0,10,15,20,25,30,40,50,60,80,120,1e9]
for prog in ('opt','sol'):
    out=[]
    for lo,hi in zip(tr,tr[1:]):
        v=[math.log(r3/r2) for p,l,n,m,r0,r2,r3 in runs if p==prog and lo<m<=hi and r0>1]
        out.append('%s-%s:%s'%(lo,'inf' if hi>1e8 else hi, '%.2f(%d)'%(math.exp(sum(v)/len(v)),len(v)) if v else '-'))
    print(' ',prog,' '.join(out))
