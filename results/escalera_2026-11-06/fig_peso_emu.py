# Figura de la frontera (fig_*.dat), peso de la pasada y emulacion del trabajo de SoPlex, desde un directorio S
# (por defecto work/S; la rehecha con limite 100 es S=work/S100). Escribe los .dat en el directorio $OUT (por defecto .).
import glob,os,re,math,statistics as st
H=os.path.dirname(os.path.abspath(__file__))
S=os.environ.get('S',H+'/work/S'); OUT=os.environ.get('OUT','.')
gm=lambda v: math.exp(sum(map(math.log,v))/len(v))
def carga(script):
    g={'__file__':H+'/'+script}
    src=open(H+'/'+script).read().split("nom=[")[0].replace('print(','(lambda *a,**k:None)(')
    exec(src,g); return g
# --- figura: pasada de SC / pasada del simplex directo contra filas
SEL_AFF=None
for script,tag in (('bescf_tab_c.py','aff'),('bcmp_tab_c.py','cmp')):
    g=carga(script); D,get,filas=g['D'],g['get'],g['filas']
    out={}
    for i in D:
        sc=get(D[i],'SC'); b=D[i].get('b2c')
        if sc and b and i in filas: out.setdefault(i[0],[]).append((filas[i],sc/b))
    for prog,v in out.items():
        fn='%s/fig_%s_%s.dat'%(OUT,prog,tag); open(fn,'w').write('m r\n'+''.join('%.2f %.4f\n'%x for x in sorted(v)))
        print('fig',prog,tag,len(v),'gm %.3f'%gm([x[1] for x in v]),'gana %d'%sum(x[1]<1 for x in v))
    if tag=='aff': SEL_AFF={i for i in D if 'b0c' in D[i] and get(D[i],'SC')}
# --- peso de la pasada: llamadas x costo de una pasada de SC en el banco / tiempo de la busqueda
for prog in ('opt','sol'):
    ws=[]
    for f in glob.glob('%s/peso/%s.*.txt'%(S,prog)):
        n=os.path.basename(f)[:-4].split('.',1)[1]
        if (prog,n) not in SEL_AFF: continue
        v=open(f).read().split()
        try: calls=int(v[1]); tm=float(v[3])
        except (ValueError,IndexError): continue
        b='%s/bescf_c/%s.%s.b1.txt'%(S,prog,n)
        m=re.search(r'costo DFB\s*: pasada ([0-9.]+) s',open(b).read()) if os.path.exists(b) else None
        if not m or tm<=0: continue
        ws.append(calls*float(m.group(1))/1000/tm)
    ws.sort(); print('peso',prog,len(ws),'mediana %.3f q1 %.3f q3 %.3f'%(st.median(ws),ws[len(ws)//4],ws[3*len(ws)//4]))
# --- emulacion: pasada de SC con el trabajo de SoPlex reinsertado, sobre la pasada del SoPlex ajustado (R2)
def pasadas(f):
    t=open(f).read() if os.path.exists(f) else ''
    d=re.search(r'costo DFB  : pasada ([0-9.]+) s',t); c=re.search(r'costo caja : pasada ([0-9.]+) s',t)
    return (float(d.group(1)),float(c.group(1))) if d and c else None
g=carga('bescf_tab_c.py'); filas=g['filas']
for prog,cl in (('opt','todas'),('sol','<=30')):
    R={k:[] for k in ('e0','e1','e2','e3','e4')}; n_ok=0; n_tot=0
    for i in SEL_AFF:
        if i[0]!=prog or (cl=='<=30' and filas.get(i,1e9)>30): continue
        n_tot+=1
        p={e:pasadas('%s/emu/%s.%s.%s.txt'%(S,prog,i[1],e)) for e in R}
        if any(v is None for v in p.values()): continue
        n_ok+=1
        for e in R: R[e].append(p[e][0]/p[e][1])
    r={e:gm(v) for e,v in R.items()}; gap=-math.log(r['e0'])
    sh=lambda a,b: math.log(r[b]/r[a])/gap
    print('emu %s %s: %d/%d  SC/R2 %.3f -> %.3f | refact %.0f%%  init %.0f%%  verif %.0f%%  copia %.0f%%'%(prog,cl,n_ok,n_tot,r['e0'],r['e4'],
          100*sh('e0','e1'),100*sh('e1','e2'),100*sh('e2','e3'),100*sh('e3','e4')))
