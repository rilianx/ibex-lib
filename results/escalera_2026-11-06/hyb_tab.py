# Hibrido real por nodo (H15/H20/H30) vs R0, R2 (bare), R3 (SC base), sistemas compo, mismas semillas 1-2, una instancia por familia y clase.
import glob,re,collections,sys
exec(open('esc_fam.py').read().split('nom=[')[0].replace("bancos=[('escf','opt','art'),('escf','sol','compo'),('escfa','sol','art')]","bancos=[('escf','sol','compo')]"))
for f in glob.glob(S+'/hyb/sol.*.csv'):
    n,p,s,r=open(f).read().strip().split(',',3)
    if n in excl: continue
    v=r.split(); t=float(v[1]) if r!='NA' and len(v)>1 else T
    R[('sol','compo',n,s)][p]=(min(t,T),t<T*0.999)
def tab(base):
  for c in ('<=30','>30'):
    sel=[(v[1],v[2]) for k,v in elig.items() if k[0]=='sol' and k[1]==c]
    print('== vs %s, %s: %d familias'%(base,c,len(sel)))
    for p in ('R2','R3','R7','H15','H20','H30'):
        out=[]
        for th in (1,5):
            pr=[]
            for (pg,lin,n,s),v in R.items():
                if (lin,n) not in sel or s not in ('1','2') or p not in v or base not in v: continue
                if not all(h in v for h in ('H15','H20','H30','R2','R3')): continue
                (a,fa),(b,fb)=v[p],v[base]
                if max(a,b)>th and (fa or fb): m=max(a,b); pr.append((a/m,b/m))
            if not pr: out.append('-'); continue
            A=sum(x for x,y in pr)/len(pr); B=sum(y for x,y in pr)/len(pr)
            out.append('%.2f (%d)'%(A/B,len(pr)))
        print('  %-5s >1s %-12s >5s %s'%(p,*out))
tab('R0'); tab('R2')
