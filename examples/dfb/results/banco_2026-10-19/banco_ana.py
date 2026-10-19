# Tabla del banco de pasadas: tiempo por pasada de SC y de la caja sobre las mismas pasadas. uso: banco_ana.py dir
import os,re,sys,math,statistics as st
D=sys.argv[1]
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
R={'sol':[],'opt':[]}
for f in sorted(os.listdir(D)):
    k,n=f[:-4].split('.',1); t=open(D+'/'+f).read()
    a=re.search(r'costo DFB\s*: pasada ([\d.]+) s\s+pivotes=(\d+)\s+us/LP=([\d.]+)\s+piv/LP=([\d.]+)',t)
    b=re.search(r'costo caja\s*: pasada ([\d.]+) s\s+pivotes=(\d+)\s+us/LP=([\d.]+)\s+piv/LP=([\d.]+)',t)
    c=re.search(r'nodos=(\d+) \| contraccion media DFB=([\d.]+) caja=([\d.]+) \| nodos vaciados DFB=(\d+) caja=(\d+)',t)
    if not (a and b and c): print('falta',f); continue
    nod=int(c[1])
    R[k].append(dict(n=n,nod=nod,tsc=float(a[1])/nod*1e6,tbb=float(b[1])/nod*1e6,psc=float(a[4]),pbb=float(b[4]),usc=float(a[3]),ubb=float(b[3]),csc=float(c[2]),cbb=float(c[3]),vsc=int(c[4]),vbb=int(c[5])))
for k,L in R.items():
    if not L: continue
    print(f'== {k}: {len(L)} instancias, {sum(x["nod"] for x in L)} pasadas')
    print(f'   pasada SC/caja x{gm([x["tsc"]/x["tbb"] for x in L]):.2f} (mediana {st.median(x["tsc"] for x in L):.0f} vs {st.median(x["tbb"] for x in L):.0f} us), mas barata en {sum(x["tsc"]<x["tbb"] for x in L)}/{len(L)}')
    print(f'   us/LP x{gm([x["usc"]/x["ubb"] for x in L]):.2f}  piv/LP x{gm([max(x["psc"],.01)/max(x["pbb"],.01) for x in L]):.2f}')
    print(f'   contraccion media SC {st.mean(x["csc"] for x in L):.4f} caja {st.mean(x["cbb"] for x in L):.4f}; vaciados SC {sum(x["vsc"] for x in L)} caja {sum(x["vbb"] for x in L)}')
    L.sort(key=lambda x:x['tsc']/x['tbb'])
    for x in L: print(f'   {x["n"]:28s} pasada {x["tsc"]:8.0f} vs {x["tbb"]:8.0f} us x{x["tsc"]/x["tbb"]:.2f} | piv/LP {x["psc"]:5.2f} vs {x["pbb"]:5.2f} | contr {x["csc"]:.3f} vs {x["cbb"]:.3f}')
