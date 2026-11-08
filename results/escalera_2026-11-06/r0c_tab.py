import re,os,math,sys
def run(script,dirn,pref):
    g={'__file__':os.path.abspath(script)}; exec(open(script).read().split("nom=[")[0].replace('print(','(lambda *a,**k:None)('),g)
    D,get,filas=g['D'],g['get'],g['filas']
    for prog in ('opt','sol'):
      for cl in ('<=30','>30'):
        sel=[i for i in D if i[0]==prog and (filas[i]<=30)==(cl=='<=30')]
        pr=[];pr2=[];pr3=[];c=[]
        for i in sel:
            f=g['S']+'/r0c/%s.%s.b0fc.txt'%(pref or prog,i[1])
            if not os.path.exists(f) or not os.path.getsize(f): continue
            t=open(f).read(); m=re.search(r'costo caja : pasada ([0-9.]+) s',t)
            if not m: continue
            x=float(m.group(1)); pr.append((x,D[i]['b0fc'])); pr2.append((get(D[i],'SC'),x)); pr3.append((D[i]['b2c'],x))
            mm=re.search(r'contraccion media DFB=([0-9.]+) caja=([0-9.]+)',t); 
            if mm: c.append(float(mm.group(2)))
        if not pr: continue
        nq=lambda P: (sum(a/max(a,b) for a,b in P)/len(P))/(sum(b/max(a,b) for a,b in P)/len(P))
        gm=lambda P: math.exp(sum(math.log(a/b) for a,b in P)/len(P))
        print(script,prog,cl,len(pr),'R1cerca/R1nat nq %.2f gm %.2f | SC/R1cerca nq %.2f | bare/R1cerca nq %.2f | contr caja %.4f'%(nq(pr),gm(pr),nq(pr2),nq(pr3),sum(c)/len(c)))
run('bescf_tab.py','bescf',None)
run('bcmp_tab.py','bcmp','cmp')
