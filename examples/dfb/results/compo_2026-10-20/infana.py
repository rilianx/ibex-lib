import sys,re,statistics as st
R=[]
for l in sys.stdin:
    m=re.search(r'(estructural|holgura) \d+ \| x=(\S+) \[(\S+), (\S+)\] ancho (\S+) \| viol (\S+) tol (\S+)',l)
    if m: R.append((m[1],abs(float(m[2])),float(m[5]),float(m[6]),float(m[7])))
if not R: print('   sin declaraciones'); sys.exit()
q=lambda v,p: sorted(v)[int(p*(len(v)-1))]
u=[r[3]/(2.2e-16*max(r[1],1.0)) for r in R]
print('   declaraciones %d (holguras %d) | mediana |x| %.3g, ancho %.3g, viol %.3g, tol %.3g'%(len(R),sum(r[0]=='holgura' for r in R),st.median(r[1] for r in R),st.median(r[2] for r in R),st.median(r[3] for r in R),st.median(r[4] for r in R)))
print('   violacion en ulps de |x|: mediana %.0f, p90 %.0f, max %.0f | tol en el piso 1e-15(1+|x|): %.0f%% | ancho < 1e-10 relativo: %.0f%%'%(st.median(u),q(u,0.9),max(u),100*sum(abs(r[4]-1e-15*(1+r[1]))<=1e-3*r[4] for r in R)/len(R),100*sum(r[2]<1e-10*(1+r[1]) for r in R)/len(R)))
