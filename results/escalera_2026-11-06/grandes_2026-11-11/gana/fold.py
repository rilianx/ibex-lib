import sys,re,subprocess,collections
def cls_spx(st):
    # st: list leaf->root
    s=' '.join(st[:40])
    if 'SPxSolver::factorize' in s or 'SPxBasis::factorize' in s or 'SLUFactor::load' in s: return 'refactorizar'
    if 'SPxSolver::init(' in s: return 'init (vectores, lado derecho, pricer)'
    if any(k in s for k in ('getPrimal','getDual','getRedCost','getSlacks','qualConstraint','unShift','unscale')): return 'extraer y verificar solucion'
    if any(k in s for k in ('PR::','SPxDevex','SPxSteep','SPxAutoPR')): return 'pricing'
    if any(k in s for k in ('RT::','SPxBoundFlipping','SPxFastRT')): return 'ratio test'
    if any(k in s for k in ('solveRight','solveLeft','vSolve','solve2right','forestUpdate','SLUFactor::change','solve4update')): return 'FTRAN/BTRAN/actualizar LU'
    if any(k in s for k in ('changeObj','changeBounds','changeLower','changeUpper','reLoad','setTerminationIter','setOpttol','setFeastol')): return 'cambiar objetivo/cotas'
    return 'resto del simplex'
def cls_sc(st):
    s=' '.join(st[:40])
    if 'DFBBasis::factorize' in s or 'DFBBasis::set_basis' in s or 'asegurar_fact' in s: return 'refactorizar'
    if any(k in s for k in ('solveRight','solveLeft','vSolve','forestUpdate','SLUFactor::change','solve4update')): return 'FTRAN/BTRAN/actualizar LU'
    return 'resto del simplex'
def run(f):
    p=subprocess.Popen(['perf','script','-i',f,'-F','ip,sym'],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,text=True)
    A=collections.Counter(); B=collections.Counter(); st=[]
    tot=0
    def flush():
        nonlocal tot
        if not st: return
        tot+=1
        s=' '.join(st)
        if 'CajaSPx::resolver' in s: A[cls_spx(st)]+=1
        elif 'DFBSimplex::solve' in s: B[cls_sc(st)]+=1
    for l in p.stdout:
        l=l.strip()
        if not l: flush(); st=[]; continue
        st.append(l.split(' ',1)[1] if ' ' in l else l)
    flush()
    return A,B,tot
for n in sys.argv[1:]:
    A,B,tot=run(n+'.dw.data'); a=sum(A.values()); b=sum(B.values())
    print('== %s  muestras: SoPlex %d  SC %d  (SoPlex/SC %.2f)'%(n,a,b,a/max(b,1)))
    print('  SoPlex:'); [print('    %-38s %5.1f%%  (%.2f x SC total)'%(k,100*v/a,v/b)) for k,v in A.most_common()]
    print('  SC:');     [print('    %-38s %5.1f%%'%(k,100*v/b)) for k,v in B.most_common()]

def leaves(f,bucket_fn,want,top=12):
    p=subprocess.Popen(['perf','script','-i',f,'-F','ip,sym'],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,text=True)
    C=collections.Counter(); st=[]
    for l in p.stdout:
        l=l.strip()
        if not l:
            if st and 'CajaSPx::resolver' in ' '.join(st) and bucket_fn(st)==want:
                # primer marco de soplex desde la hoja, y su llamador dentro de SPxSolver
                sx=[x for x in st if 'soplex::' in x]
                C[(re.sub(r'\(.*','',sx[0]) if sx else st[0])+'  <- '+next((re.sub(r'\(.*','',x) for x in st if 'SPxSolver::' in x or 'SPxBasis::' in x),'-')]+=1
            st=[]; continue
        st.append(l.split(' ',1)[1] if ' ' in l else l)
    for k,v in C.most_common(top): print('      %5d  %s'%(v,k[:130]))
