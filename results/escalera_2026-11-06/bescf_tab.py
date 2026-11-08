# Escalera por pasada sobre el banco (71 instancias x 1000 llamadas grabadas, afin). Cada peldano contra b0 caja.
import glob,os,re,collections,math,sys
import os; S=os.environ.get('S',os.path.join(os.path.dirname(os.path.abspath(__file__)),'work/S'))
M=os.environ.get('M',os.path.join(os.path.dirname(os.path.abspath(__file__)),'work/M'))
filas={}
for f in os.listdir(M+'/filas_m2'):
    if f.count('.')!=3: continue
    k,n,lin,_=f.split('.')
    mm=re.search(r'([0-9.]+) por linealizacion',open(M+'/filas_m2/'+f).read())
    if mm and lin=='art': filas[(k,n)]=float(mm.group(1))
D=collections.defaultdict(dict); contr=[]; noopt=collections.Counter()
for f in glob.glob(S+'/bescf/*.txt'):
    b=os.path.basename(f)[:-4]; prog,rest=b.split('.',1); n,pel=rest.rsplit('.',1)
    t=open(f).read()
    d=re.search(r'costo DFB  : pasada ([0-9.]+) s',t); c=re.search(r'costo caja : pasada ([0-9.]+) s',t)
    if not d or not c: print('sin dato',b); continue
    D[(prog,n)][pel+'c']=float(c.group(1)); D[(prog,n)][pel+'d']=float(d.group(1))
    m=re.search(r'contraccion media DFB=([0-9.e-]+) caja=([0-9.e-]+)',t)
    if m: contr.append((b,float(m.group(1)),float(m.group(2))))
    for side in ('DFB ','caja'):
        mm=re.search(r'\[abcaja\] %s.*no optimo=(\d+)'%side,t)
        if mm: noopt[(pel,side)]+=int(mm.group(1))
# optimizacion: solo la seleccion una por familia y clase (bopt_sel.py)
SEL=set(open(S+'/bopt_sel.txt').read().split())
for i in list(D):
    if i[0]=='opt' and i[1] not in SEL: print('fuera de la seleccion',i); del D[i]
# variantes del mismo problema (bis/ter/-1/_1): queda la de mas filas
def base(n): return re.sub(r'(bis|ter|-1|_1)$','',n) if re.match(r'ex\d+_\d+_\d+',n) else n
best={}
for (prog,n) in list(D):
    k=(prog,base(n))
    if k not in best or filas.get((prog,n),0)>filas.get(best[k],0): best[k]=(prog,n)
for i in list(D):
    if best[(i[0],base(i[1]))]!=i: print('variante descartada',i); del D[i]
for i in D: D[i]['_m']=filas[i]
def get(v,k):
    if k=='SC':  return (v['b0fd']+v['b1d'])/2 if 'b0fd' in v and 'b1d' in v else None
    if k=='SC+': return v.get('b2d')
    if k=='HYB': return get(v,'SC') if v['_m']<=30 else v.get('b2c')
    return v.get(k)
nom=[('b0fc','+ merge twin rows'),('b1c','+ SoPlex tuning'),('b2c','bare simplex (SPxSolver)'),('SC','SC base'),('HYB','hybrid (SC <=30 rows, else bare)')]
for prog in ('opt','sol'):
  for cl in ('<=30','>30','todas'):
    sel=[k for k in D if k[0]==prog and k in filas and (cl=='todas' or (filas[k]<=30)==(cl=='<=30'))]
    if not sel: continue
    print('== %s %s filas: %d instancias'%(prog,cl,len(sel)))
    for k,nm in nom:
        pr=[];g=[]
        for i in sel:
            a=get(D[i],k); b=D[i].get('b0c')
            if a is None or b is None: continue
            m=max(a,b); pr.append((a/m,b/m)); g.append(math.log(a/b))
        A=sum(x for x,y in pr)/len(pr); B=sum(y for x,y in pr)/len(pr)
        print('  %-30s %.2f:%.2f  %.2f   gm %.2f (%d)'%(nm,A,B,A/B,math.exp(sum(g)/len(g)),len(pr)))
    # SC vs bare directo
    g=[math.log(get(D[i],'SC')/D[i]['b2c']) for i in sel if get(D[i],'SC') and 'b2c' in D[i]]
    print('  SC base / bare  gm %.2f   win %d/%d'%(math.exp(sum(g)/len(g)),sum(x<0 for x in g),len(g)))
sinf=[k for k in D if k not in filas]; print('sin filas:',sinf)
dif=[c for c in contr if abs(c[1]-c[2])>1e-4*max(1,abs(c[2]))]; print('contraccion distinta: %d/%d'%(len(dif),len(contr)), dif[:8])
print('no optimo:',dict(noopt))
if '-v' in sys.argv:
    for i in sorted(D): print(i, filas.get(i), {k:round(get(D[i],k),3) for k in ('b0c','b0fc','b1c','b2c','SC','HYB') if get(D[i],k)})
# efecto de la fusion en cada motor: SC (lado DFB b0 vs b0f) y SoPlex de Ibex (caja b0 vs b0f)
print('== efecto de la fusion (con/sin)')
for prog in ('opt','sol'):
  for cl in ('<=30','>30'):
    sel=[i for i in D if i[0]==prog and (filas[i]<=30)==(cl=='<=30')]
    for nm,a,b in (('SC','b0fd','b0d'),('SoPlex Ibex','b0fc','b0c')):
        pr=[(D[i][a],D[i][b]) for i in sel if a in D[i] and b in D[i]]
        A=sum(x/max(x,y) for x,y in pr)/len(pr); B=sum(y/max(x,y) for x,y in pr)/len(pr)
        g=math.exp(sum(math.log(x/y) for x,y in pr)/len(pr))
        print('  %s %s %-12s %.2f:%.2f %.2f gm %.2f (%d)'%(prog,cl,nm,A,B,A/B,g,len(pr)))
