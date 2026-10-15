# dump2bch.py original.bch linea_dump salida.bch : cambia el dominio de x[n] por la caja volcada
import sys,re
src,lin,out=sys.argv[1],sys.argv[2],sys.argv[3]
v=lin.split()[2:]; n=len(v)//2
t=open(src).read()
t=re.sub(r'x\((\d+)\)',r'x\1',t)
dom=''.join(f'  x{i+1} in [{v[2*i]}, {v[2*i+1]}];\n' for i in range(n))
t=re.sub(r'\s*x\[\d+\] in \[[^\]]*\];\n',  '\n'+dom, t, count=1)
open(out,'w').write(t)
