import os,math,sys
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
excl=set(l.split()[1] for l in open('/home/iaraya/ibex2024/examples/dfb/results/excluidas_3600.txt') if len(l.split())>1)
def leer(d,pat):
    R={}
    for f in os.listdir(d):
        if not f.startswith('sol.') or '.cajacerca.' not in f: continue
        a=open(d+'/'+f).read().strip().split(',')
        n,s=a[0],a[2]
        if n in excl: continue
        r=a[3].split()
        R[(n,s)]=None if r[0]=='NA' else (float(r[0]),float(r[1]))
    return R
A=leer(M+'/ordcaja',''); B=leer(M+'/fcaja','')
k=[x for x in A if x in B]
ter=lambda r: r is not None and r[1]<59.9
print('pares',len(k),'terminan sin',sum(ter(A[x]) for x in k),'con',sum(ter(B[x]) for x in k))
for th in (0.1,1,5):
    L=[(B[x][1]/A[x][1],B[x][0]/A[x][0]) for x in k if ter(A[x]) and ter(B[x]) and A[x][1]>th]
    print(f'>{th}: n={len(L)} t x{math.exp(sum(math.log(a) for a,b in L)/len(L)):.3f} celdas x{math.exp(sum(math.log(b) for a,b in L)/len(L)):.3f}')
ig=sum(1 for x in k if ter(A[x]) and ter(B[x]) and A[x][0]==B[x][0]); print('celdas identicas',ig)
D=sorted(((B[x][1]/A[x][1],x,A[x],B[x]) for x in k if ter(A[x]) and ter(B[x]) and A[x][1]>1))
for d in D[:5]+D[-5:]: print(d)
for x in k:
    if ter(A[x])!=ter(B[x]): print('cambia termino',x,A[x],B[x])
