#!/bin/bash
# Banco de pasadas: captura 1000 cajas (una de cada $3 llamadas) y las reproduce con DFBH_ABCAJA. $1 ruta $2 sol|opt $3 cada
cd ~/ibex2024/examples/dfb
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
n=$(basename $1 .bch); c=$M/banco/$2.$n.txt; o=$M/banco_ab/$2.$n.txt
ulimit -v 6000000
if [ ! -s $c ]; then
  if [ $2 = opt ]; then DFBH_CAPTURA=$c DFBH_CAPTURA_N=1000 DFBH_CAPTURA_CADA=$3 timeout 100 ./ibexopt_cap $1 --filtering=acidhc4 --lr=dfbhull --seed=1 -t60 > /dev/null 2>&1
  else DFBH_CAPTURA=$c DFBH_CAPTURA_N=1000 DFBH_CAPTURA_CADA=$3 timeout 100 ./ibexsolve_dfb_cap $1 --lr=dfbhull --lin=art -t 60 --seed=1 > /dev/null 2>&1; fi
fi
[ -s $o ] && exit 0
if [ $2 = opt ]; then tol=1e-9; else tol=1e-6; fi
DFBH_ABCAJA=1 DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol timeout 1800 ./ibexreplay $1 $2 $c 2>&1 | grep -E 'abcaja|replay' > $o
