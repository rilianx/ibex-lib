import glob,re,os,collections,statistics as st,math
exec(open('bescf_tab.py').read().split("def get(v,k):")[0].replace("print('variante descartada',i);",""))
SEL=set(open(S+'/bopt_sel.txt').read().split()); drop={('opt','ex8_5_2'),('opt','concon'),('opt','pentagon'),('opt','ex8_2_4')}
X=collections.defaultdict(lambda: collections.defaultdict(list))
for f in glob.glob(S+'/bescf/*.txt'):
    b=os.path.basename(f)[:-4]; prog,rest=b.split('.',1); n,pel=rest.rsplit('.',1)
    if (prog,n) in drop or (prog=='opt' and n not in SEL): continue
    t=open(f).read(); cl=prog+('' if prog=='opt' else ('<=30' if filas[(prog,n)]<=30 else '>30'))
    for side,tag in (('DFB ','d'),('caja','c')):
        c=re.search(r'costo %s\s*: pasada ([0-9.]+) s\s+pivotes=(\d+)\s+us/LP=([0-9.]+)\s+piv/LP=([0-9.]+)'%side.strip(),t)
        l=re.search(r'\[abcaja\] %s\s*: LPs=(\d+) .*no optimo=(\d+)'%side.strip(),t)
        k=pel+tag
        X[cl][k+'_pass'].append(float(c.group(1))*1000); X[cl][k+'_uslp'].append(float(c.group(3))); X[cl][k+'_piv'].append(float(c.group(4)))
        X[cl][k+'_lps'].append(int(l.group(1))); X[cl][k+'_no'].append(int(l.group(2)))
    m=re.search(r'contraccion media DFB=([0-9.e-]+) caja=([0-9.e-]+)',t)
    X[cl][pel+'d_con'].append(float(m.group(1))); X[cl][pel+'c_con'].append(float(m.group(2)))
for cl in X:
    print('==',cl,len(X[cl]['b0c_pass']))
    for k,nm in (('b0c','Ibex'),('b0fc','Ibex+merge'),('b1c','tuned'),('b2c','bare'),('b1d','SC'),('b0d','SC nomerge')):
        x=X[cl]
        print('  %-11s pass med %7.0f us  us/LP med %6.1f  piv/LP med %5.2f  noopt %5.2f%%  contr %.4f'%(nm,st.median(x[k+'_pass']),st.median(x[k+'_uslp']),st.median(x[k+'_piv']),100*sum(x[k+'_no'])/sum(x[k+'_lps']),st.mean(x[k+'_con'])))
