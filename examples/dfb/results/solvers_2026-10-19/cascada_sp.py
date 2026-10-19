# Cascada unica SC -> SoPlex sin presolve sobre el banco de pasadas (participacion en log).
import os,re,math,sys
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
def t(path,who):
    if not os.path.exists(path): return None
    s=open(path).read()
    a=re.search(r'costo %s\s*: pasada ([\d.]+) s'%who,s); c=re.search(r'nodos=(\d+)',s)
    return float(a[1])/int(c[1]) if a and c and float(a[1])>0 else None
REF=sys.argv[1] if len(sys.argv)>1 else 'banco_sp'
pasos=['sc_fact','sc_init','sc_verif','sc_copia']
for k in ('sol','opt'):
    L=[]
    for f in sorted(os.listdir(M+'/'+REF)):
        if not f.startswith(k+'.'): continue
        n=f[len(k)+1:-4]
        sc=t(M+'/'+REF+'/'+f,'DFB'); bb=t(M+'/'+REF+'/'+f,'caja')
        v=[t(f'{M}/cascada/{k}.{n}.{p}.txt','DFB') for p in pasos]
        if sc and bb and all(v): L.append([sc]+v+[bb])
    N=len(L); g=lambda i: math.exp(sum(math.log(x[i]/x[-1]) for x in L)/N)
    tot=math.log(g(0))
    print(f'{k}: {N} instancias, SC/BB x{g(0):.3f}')
    for i,p in enumerate(pasos,1): print(f'   +{p:9s} -> x{g(i):.3f}  participacion {100*(math.log(g(i))-math.log(g(i-1)))/(-tot):5.1f}%')
    print(f'   resto      -> x1.000  participacion {100*(-math.log(g(4)))/(-tot):5.1f}%')
