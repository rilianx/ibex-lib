import subprocess,re,sys,collections
def load(f):
    out=subprocess.run(['perf','report','-i',f,'--stdio','--no-children','-n','--sort','symbol'],capture_output=True,text=True).stdout
    d={}
    for l in out.splitlines():
        m=re.match(r'\s*([0-9.]+)%\s+(\d+)\s+\[(.)\]\s+(.*)',l)
        if m: d[m.group(4).strip()]=d.get(m.group(4).strip(),0)+int(m.group(2))
    return d
def grp(s):
    if 'DFBSimplex' in s or 'DFBBasis' in s or 'DFBFact' in s: return 'SC propio'
    if 'CLUFactor' in s or 'SLUFactor' in s: return 'LU (codigo de SoPlex)'
    if 'soplex::' in s: return 'SoPlex resto'
    if 'CajaSPx' in s: return 'SoPlex resto'
    return 'otros'
for n in sys.argv[1:]:
    sc=load('perf.%s.sc.data'%n); ab=load('perf.%s.ab.data'%n)
    G=collections.Counter(); H=collections.Counter()
    for s,c in sc.items(): G[grp(s)]+=c
    for s,c in ab.items(): H[grp(s)]+=c
    print('==',n,'muestras: solo SC',sum(G.values()),' SC+SoPlex',sum(H.values()))
    for k in ('SC propio','LU (codigo de SoPlex)','SoPlex resto','otros'): print('  %-24s solo SC %6d   SC+SoPlex %6d   diferencia %6d'%(k,G[k],H[k],H[k]-G[k]))
    print('  -> SC: propio %d + LU %d = %d | SoPlex: propio %d + LU %d = %d'%(G['SC propio'],G['LU (codigo de SoPlex)'],G['SC propio']+G['LU (codigo de SoPlex)'],H['SoPlex resto'],H['LU (codigo de SoPlex)']-G['LU (codigo de SoPlex)'],H['SoPlex resto']+H['LU (codigo de SoPlex)']-G['LU (codigo de SoPlex)']))
    print('  SC propio, detalle:'); 
    for s,c in sorted(sc.items(),key=lambda x:-x[1]):
        if grp(s)=='SC propio' and c>=0.03*G['SC propio']: print('     %5d %s'%(c,s[:95]))
    print('  LU en SC, detalle:')
    for s,c in sorted(sc.items(),key=lambda x:-x[1]):
        if grp(s).startswith('LU') and c>=0.04*G['LU (codigo de SoPlex)']: print('     %5d (en SC+SoPlex %5d) %s'%(c,ab.get(s,0),s[:70]))
    print('  SoPlex resto, detalle:')
    for s,c in sorted(ab.items(),key=lambda x:-x[1]):
        if grp(s)=='SoPlex resto' and c>=0.04*H['SoPlex resto']: print('     %5d %s'%(c,s[:95]))
