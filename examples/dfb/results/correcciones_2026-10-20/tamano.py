import os,re,csv,math,collections,sys
M='/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad'
EXC={'alkyl','alkylbis','bearing'}
F={}
for f in os.listdir(M+'/filas_m2'):
    k,rest=f.split('.',1); n,l=rest[:-4].rsplit('.',1)
    t=open(M+'/filas_m2/'+f).read(); m=re.search(r'filas de la relajacion: ([\d.]+)',t)
    if m: F[(k,n,l)]=float(m[1])
def P(k,r):
    if r=='NA': return None
    x=r.split(); t=float(x[1]); return (max(t,1e-3),max(1.,float(x[0])),t<59.94 and not(k=='sol' and x[2]=='3'))
D=collections.defaultdict(dict)
for k,n,a,s,r in list(csv.reader(open('/home/iaraya/ibex2024/examples/dfb/results/final_2026-10-21.csv')))[1:]: D[(k,n)][(a,s)]=P(k,r)
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
arms={('opt','art'):('y_sc_art','y_caja_art'),('opt','compo'):('y_sc_compo','y_caja_compo'),('sol','art'):('z_sc_art','y_caja_art'),('sol','compo'):('w_sc_compo','y_caja_compo')}
bins=[(0,15),(15,30),(30,60),(60,120),(120,10**9)]
print(f'{"programa":5s} {"relaj":5s} '+' '.join(f'{f"m {a}-{b}" if b<10**9 else f"m >{a}":>16s}' for a,b in bins)+'   (SC/caja cpu, n instancias con >1 s)')
for (k,l),(x,ref) in arms.items():
    B=collections.defaultdict(list)
    for (kk,n),d in D.items():
        if kk!=k or (k=='opt' and n in EXC): continue
        m=F.get((k,n,l))
        if m is None: continue
        v=[]
        for s in '1234':
            X=d.get((x,s)); R=d.get((ref,s))
            if X and R and X[2] and R[2] and R[0]>1: v.append(X[0]/R[0])
        if v:
            for a,b in bins:
                if a<=m<b: B[(a,b)].append(gm(v))
    print(f'{k:5s} {l:5s} '+' '.join(f'{gm(B[bb]):7.2f} (n={len(B[bb]):3d})' for bb in bins))
print('instancias con tamano:',len(F))
