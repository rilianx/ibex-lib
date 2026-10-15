import os,sys,math,collections
D=sys.argv[1]; T=60.0
EXC={'alkyl','alkylbis','bearing'}
R=collections.defaultdict(dict)
for f in os.listdir(D):
    if not f.endswith('.csv'): continue
    kind=f.split('.')[0]
    n,arm,seed,r=open(os.path.join(D,f)).read().strip().split(',',3)
    if r=='NA': R[(kind,n)][(arm,seed)]=None; continue
    a=r.split(); t=float(a[1]); c=max(1.0,float(a[0]))
    ok=t<T*0.999 and not (kind=='sol' and a[2]=='3')
    R[(kind,n)][(arm,seed)]=(max(t,1e-3),c,ok)
gm=lambda v: math.exp(sum(map(math.log,v))/len(v)) if v else float('nan')
def fin(kind,arm):
    runs=0; solved=0
    for (k,n),d in R.items():
        if k!=kind: continue
        v=[x for (a,s),x in d.items() if a==arm]
        runs+=sum(1 for x in v if x and x[2])
        if v and all(x and x[2] for x in v): solved+=1
    return runs,solved
def cmp(kind,x,ref,show=False):
    pairs=[]
    for (k,n),d in R.items():
        if k!=kind or (kind=='opt' and n in EXC): continue
        for (a,s),r in d.items():
            if a!=ref or not r: continue
            y=d.get((x,s))
            if y and r[2] and y[2]: pairs.append((n,r[0],y[0]/r[0],y[1]/r[1]))
    out=[]
    for th in (0.1,1,5):
        v=[p for p in pairs if p[1]>th]
        out.append((len(v),gm([p[2] for p in v]),gm([p[3] for p in v])))
    if show:
        fam=collections.defaultdict(list)
        for p in pairs:
            if p[1]>1: fam[p[0]].append(p[2])
        z=sorted((gm(v),n) for n,v in fam.items())
        print('     mejores:',', '.join('%s x%.2f'%(n,g) for g,n in z[:6]))
        print('     peores :',', '.join('%s x%.2f'%(n,g) for g,n in z[-6:]))
    return out
def line(kind,x,ref,show=False):
    o=cmp(kind,x,ref,show); fx=fin(kind,x); fr=fin(kind,ref)
    print('%s %-10s / %-6s cells %.3f/%.3f/%.3f  cpu %.3f/%.3f/%.3f  (n=%d/%d/%d)  runs %d/%d  solved %d/%d'%(
        kind,x,ref,o[0][2],o[1][2],o[2][2],o[0][1],o[1][1],o[2][1],o[0][0],o[1][0],o[2][0],fx[0],fr[0],fx[1],fr[1]))
for kind in ('opt','sol'):
    print('== principal')
    line(kind,'sc','caja',True); line(kind,'lejos','caja',True)
    print('== ablacion')
    for a in 'sinpert tibio sinlargo sindevex sintau tableau revisado sinfusion lejos'.split(): line(kind,a,'sc')

if len(sys.argv)>2:
  import re
  def famg(kind,x,ref,pat):
    v=[];inst=set();nf=[0,0]
    for (k,n),d in R.items():
      if k!=kind or not re.fullmatch(pat,n): continue
      for (a,s),r in d.items():
        if a!=ref: continue
        y=d.get((x,s))
        if r and y and r[2] and y[2]: v.append(y[0]/r[0]); inst.add(n)
        nf[0]+= bool(y and y[2]); nf[1]+= bool(r and r[2])
    return gm(v),len(inst),nf
  for pat in ['Brown-(1[89]|2\\d|30)(-1)?','Geneig-icse','Bratu-0(1[2-9]\\d|200)','Katsurabis-14']:
    for x in ('sc','lejos'): print(pat,x,famg('sol',x,'caja',pat))
  # worst per instance own order
  for x in ('sc','lejos'):
    w=[]
    for (k,n),d in R.items():
      if k!='sol': continue
      v=[d[(x,s)][0]/d[('caja',s)][0] for (a,s) in d if a=='caja' and d.get((x,s)) and d[('caja',s)] and d[(x,s)][2] and d[('caja',s)][2] and d[('caja',s)][0]>1]
      if v: w.append((gm(v),n))
    w.sort(); print(x,'peores >1s',w[-4:])
