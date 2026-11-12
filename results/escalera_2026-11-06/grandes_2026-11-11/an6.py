# Experimentos laterales con el limite de 100 iteraciones (run4 l100 como base, run6 las variantes).
import re,glob,os,math,statistics as st
exec(open('an3.py').read().split('g=lambda')[0])
g=lambda x: math.exp(sum(map(math.log,x))/len(x)) if x else float('nan')
def carga(pat):
    D={}
    for f in glob.glob(pat):
        b=os.path.basename(f)[:-4]; bank,rest=b.split('.',1); n,v=rest.rsplit('.',1); r=rd(f)
        if r: D.setdefault((bank,n),{})[v]=r
    return D
D=carga('run6/*.txt')
for k,v in carga('run4/*.txt').items():
    if 'l100' in v: D.setdefault(k,{})['l100']=v['l100']
def linea(bank,vs):
    full=[k for k in D if k[0]==bank and all(v in D[k] for v in vs)]
    print('== %s (%d)'%(bank,len(full)))
    for v in vs:
        x=[D[k][v] for k in full]
        print('  %-6s SC/SoPlex gm %.3f gana %d/%d | pivotes %.2f | contr SC %.4f SoPlex %.4f'%(v,g([r['DFB'][0]/r['caja'][0] for r in x]),
              sum(r['DFB'][0]<r['caja'][0] for r in x),len(x),g([r['DFB'][1]/max(1,r['caja'][1]) for r in x]),
              st.mean(r['con'][0] for r in x),st.mean(r['con'][1] for r in x)))
    return full
print('### representacion de la base (SoPlex directo)')
f=linea('cmp',('l100','col'))
print('  SoPlex col / filas gm %.2f'%g([D[k]['col']['caja'][0]/D[k]['l100']['caja'][0] for k in f]))
f=linea('aff',('l100','fil'))
print('### tolerancia 1e-6 en SC (afin)'); linea('aff',('l100','t6'))
print('### SoPlex parchado (sin refactorizar al optimo ni verificar precision)')
for grp,banks in (('opt',('optM','optS')),('aff',('aff',))):
    full=[k for k in D if k[0] in banks and 'pbase' in D[k] and 'pnrp' in D[k]]
    x=[D[k] for k in full]
    print('  %s (%d): SoPlex parchado/base gm %.3f | SC/base gm %.3f | SC/parchado gm %.3f gana %d | contr base %.4f parchado %.4f'%(grp,len(full),
          g([r['pnrp']['caja'][0]/r['pbase']['caja'][0] for r in x]), g([r['pbase']['DFB'][0]/r['pbase']['caja'][0] for r in x]),
          g([r['pnrp']['DFB'][0]/r['pnrp']['caja'][0] for r in x]), sum(r['pnrp']['DFB'][0]<r['pnrp']['caja'][0] for r in x),
          st.mean(r['pbase']['con'][1] for r in x), st.mean(r['pnrp']['con'][1] for r in x)))
