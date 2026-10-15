import os,re,sys,math,statistics
gm=lambda v: math.exp(sum(map(math.log,v))/len(v))
def parse(f):
    t=open(f).read()
    if 'nodos=' not in t: return None
    nod=int(re.search(r'nodos=(\d+)',t).group(1))
    lp={e:int(re.search(r'\[abcaja\] %s\s*: LPs=(\d+)'%e,t).group(1)) for e in ('DFB','caja')}
    c={}
    for e in ('DFB','caja'):
        m=re.search(r'costo %s\s*: pasada ([\d.]+) s  pivotes=(\d+)  us/LP=([\d.]+)  piv/LP=([\d.]+)  us/piv=([\d.]+)'%e,t)
        c[e]=tuple(map(float,m.groups()))
    return nod,lp,c
def rows(d,kind):
    R=[]
    for f in sorted(os.listdir(d)):
        if not f.startswith(kind+'.'): continue
        p=parse(os.path.join(d,f))
        if not p: continue
        nod,lp,c=p
        if nod<20 or min(lp.values())<100: continue
        b,s=c['caja'],c['DFB']
        R.append(dict(n=f[len(kind)+1:-4],tbb=b[0]/nod*1e6,tsc=s[0]/nod*1e6,pbb=b[1]/nod,psc=s[1]/nod,upbb=b[4],upsc=s[4],plbb=b[3],plsc=s[3],ulbb=b[2],ulsc=s[2]))
    return R
d=sys.argv[1]; kind=sys.argv[2]; out=sys.argv[3] if len(sys.argv)>3 else None
R=rows(d,kind)
piv=[r for r in R if r['plbb']>=2 and r['plsc']>=2]
print(kind,'n=%d pasada x%.3f (mediana %.0f vs %.0f us) mas barato %.0f%% | pivotes x%.3f | us/piv x%.3f (n=%d) | piv/LP %.2f vs %.2f'%(len(R),gm([r['tsc']/r['tbb'] for r in R]),statistics.median(r['tsc'] for r in R),statistics.median(r['tbb'] for r in R),100*sum(r['tsc']<r['tbb'] for r in R)/len(R),gm([r['psc']/r['pbb'] for r in R]),gm([r['upsc']/r['upbb'] for r in piv]),len(piv),gm([r['plsc'] for r in R]),gm([r['plbb'] for r in R])))
if out:
    with open(out+'/costo_pasada_%s.dat'%kind,'w') as o:
        o.write('tbb tsc pbb psc upbb upsc plbb plsc instancia\n')
        for r in R: o.write('%.3f %.3f %.3f %.3f %.4f %.4f %.3f %.3f %s\n'%(r['tbb'],r['tsc'],r['pbb'],r['psc'],r['upbb'],r['upsc'],r['plbb'],r['plsc'],r['n']))
    with open(out+'/costo_pivote_%s.dat'%kind,'w') as o:
        o.write('tbb tsc pbb psc upbb upsc plbb plsc instancia\n')
        for r in piv: o.write('%.3f %.3f %.3f %.3f %.4f %.4f %.3f %.3f %s\n'%(r['tbb'],r['tsc'],r['pbb'],r['psc'],r['upbb'],r['upsc'],r['plbb'],r['plsc'],r['n']))
    with open(out+'/costo_lp_%s.dat'%kind,'w') as o:
        o.write('bb sc instancia\n')
        for r in R: o.write('%.3f %.3f %s\n'%(r['ulbb'],r['ulsc'],r['n']))
