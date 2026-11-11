import re,glob,os,math,statistics as st,sys
def rd(f):
    t=open(f).read(); d={}
    for side in ('DFB','caja'):
        c=re.search(r'costo %s\s*: pasada ([0-9.]+) s\s+pivotes=(\d+)'%side,t)
        l=re.search(r'\[abcaja\] %s\s*: LPs=(\d+) .*no optimo=(\d+)'%side,t)
        if not c or not l: return None
        d[side]=(float(c.group(1)),int(c.group(2)),int(l.group(1)),int(l.group(2)))
    m=re.search(r'contraccion media DFB=([0-9.e-]+) caja=([0-9.e-]+)',t); d['con']=(float(m.group(1)),float(m.group(2)))
    return d
g=lambda x: math.exp(sum(map(math.log,x))/len(x)) if x else float('nan')
for bank,vs in (('cmp',('c0','c60','c150')),('aff',('c0','c150','t6'))):
    D={}
    for f in glob.glob('run/%s.*.txt'%bank):
        b=os.path.basename(f)[len(bank)+1:-4]; n,v=b.rsplit('.',1); r=rd(f)
        if r: D.setdefault(n,{})[v]=r
    full=[n for n in D if all(v in D[n] for v in vs)]
    print('== %s: completas %d de %d'%(bank,len(full),len(D)))
    for v in vs:
        x=[D[n][v] for n in full]
        if not x: continue
        print('  %-5s pasada SC/bare gm %.2f | pivotes gm %.2f | t/pivote gm %.2f | contr SC %.4f bare %.4f | sin optimo SC %d bare %d | SC gana %d/%d'%(v,
          g([r['DFB'][0]/r['caja'][0] for r in x]), g([r['DFB'][1]/max(1,r['caja'][1]) for r in x]),
          g([(r['DFB'][0]/max(1,r['DFB'][1]))/(r['caja'][0]/max(1,r['caja'][1])) for r in x]),
          st.mean(r['con'][0] for r in x), st.mean(r['con'][1] for r in x), sum(r['DFB'][3] for r in x), sum(r['caja'][3] for r in x),
          sum(r['DFB'][0]<r['caja'][0] for r in x),len(x)))
    if '-v' in sys.argv:
        for n in sorted(full,key=lambda n:D[n][vs[0]]['DFB'][0]/D[n][vs[0]]['caja'][0]):
            print('   %-26s'%n,' | '.join('%s %.2f piv %.2f con %.4f'%(v,D[n][v]['DFB'][0]/D[n][v]['caja'][0],D[n][v]['DFB'][1]/max(1,D[n][v]['caja'][1]),D[n][v]['con'][0]) for v in vs))
