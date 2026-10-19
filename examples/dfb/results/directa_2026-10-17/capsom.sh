#!/bin/bash
# captura 1500 cajas (una de cada 3 llamadas) y las reproduce con la sombra. $1 ruta $2 sol|opt
cd ~/ibex2024/examples/dfb
n=$(basename $1 .bch); C=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/cap; o=$C/$2.$n.sombra.txt
[ -s $o ] && exit 0
if [ $2 = opt ]; then DFBH_CAPTURA=$C/$2.$n.txt DFBH_CAPTURA_N=1500 DFBH_CAPTURA_CADA=3 timeout 100 ./ibexopt_cap $1 --filtering=acidhc4 --lr=dfbhull --seed=1 -t60 > /dev/null 2>&1; tol=1e-9
else DFBH_CAPTURA=$C/$2.$n.txt DFBH_CAPTURA_N=1500 DFBH_CAPTURA_CADA=3 timeout 100 ./ibexsolve_dfb_cap $1 --lr=dfbhull --lin=art -t 60 --seed=1 > /dev/null 2>&1; tol=1e-6; fi
DFBH_SOMBRA=1 DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_TOL=$tol timeout 900 ./ibexreplay $1 $2 $C/$2.$n.txt 2>&1 | grep -E 'sombra|replay' > $o
