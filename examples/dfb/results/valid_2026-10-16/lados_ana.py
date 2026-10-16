# tab:sides desde una sonda de lados (dual/primal) y la de franjas (sc). uso: lados_ana.py dir
import glob,re,statistics as st,math,sys
M=sys.argv[1]
D={}
for f in glob.glob(M+'/*.txt'):
    b=f.split('/')[-1][:-4]
    if b.count('.')<2: continue
    k,rest=b.split('.',1); n,arm=rest.rsplit('.',1)
    t=open(f).read()
    m=re.search(r'LPs tibios=(\d+)\s+m=([\d.]+)\s+columnas basicas nuevas=([\d.]+) \(([\d.]+)%.*pivotes=([\d.]+).*siguen: ([\d.]+), cambian de lado: ([\d.]+) \(([\d.]+)%',t)
    if m: D[(k,n,arm)]=dict(lps=int(m[1]),m=float(m[2]),nb=float(m[3]),pnb=float(m[4]),piv=float(m[5]),nn=float(m[6]),lado=float(m[7]),plado=float(m[8]))
    b2=re.search(r'variables ([\d.]+) \(cambian ([\d.]+)\), holguras ([\d.]+) \(cambian ([\d.]+)\).*franja .* media ([\d.e-]+), en las que cambian ([\d.e-]+)',t)
    if b2 and m: D[(k,n,arm)].update(nv=float(b2[1]),fv=float(b2[2]),nh=float(b2[3]),fh=float(b2[4]),rel=float(b2[5]),relf=float(b2[6]))
for k in ('sol','opt'):
    I=[n for (kk,n,a) in D if kk==k and a=='dual' and (k,n,'primal') in D and D[(k,n,'dual')]['lps']>=1000 and D[(k,n,'primal')]['lps']>=1000]
    ref='dual' if k=='sol' else 'primal'
    X=[D[(k,n,ref)] for n in I]
    print(k,'instancias',len(I))
    for key in ('m','pnb','lado','plado'):
        v=[x[key] for x in X]; qq=st.quantiles(v,n=4); print('  %s mediana %.2f  cuartiles %.2f / %.2f'%(key,st.median(v),qq[0],qq[2]))
    pd=[D[(k,n,'dual')]['piv'] for n in I]; pp=[D[(k,n,'primal')]['piv'] for n in I]
    print('  piv/LP mediana dual %.2f primal %.2f ; primal-dual mediana %.2f ; lados mediana %.2f'%(st.median(pd),st.median(pp),st.median([b-a for a,b in zip(pd,pp)]),st.median(D[(k,n,ref)]['lado'] for n in I)))
    Y=[D[(k,n,'sc')] for n in I if (k,n,'sc') in D and 'nv' in D[(k,n,'sc')] and D[(k,n,'sc')]['lps']>=1000]
    print('  franjas (sc) n=%d: var cambian %.1f%%  holguras cambian %.1f%%  franja media %.3f  en las que cambian %.3f'%(len(Y),
        100*st.median(y['fv']/y['nv'] for y in Y if y['nv']>0),100*st.median(y['fh']/y['nh'] for y in Y if y['nh']>0),st.median(y['rel'] for y in Y),st.median(y['relf'] for y in Y)))
