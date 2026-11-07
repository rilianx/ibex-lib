import glob,os,re,collections
exec(open('esc_fam.py').read().split("bancos=")[0])
R0=collections.defaultdict(float)
for f in glob.glob(S+'/escf/sol.*.R0.*.csv'):
    n,p,s,r=open(f).read().strip().split(',',3)
    v=r.split(); t=float(v[1]) if r!='NA' and len(v)>1 else T
    R0[n]=max(R0[n],t)
P={}
for l in open(S+'/escf_jobs.txt'):
    v=l.split()
    if len(v)>=2 and v[1]=='sol': P[os.path.basename(v[0])[:-4]]=v[0]
elig={}
for n,t in R0.items():
    if t<=1 or ('sol',n,'compo') not in filas: continue
    m=filas[('sol',n,'compo')]
    if m<=30: continue
    k=familia('sol',n)
    if k not in elig or m>elig[k][0]: elig[k]=(m,n)
sel=sorted(n for m,n in elig.values())
print(len(sel),'familias >30 con compo; filas',sorted(int(m) for m,n in elig.values()))
open(S+'/bcmp_lista.txt','w').write(''.join(P[n]+'\n' for n in sel))
