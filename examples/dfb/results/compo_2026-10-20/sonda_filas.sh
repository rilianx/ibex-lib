#!/bin/bash
# Sonda con fusion de filas paralelas. $1 ruta $2 paralelas|quita $3 sol|opt (relajacion compo)
cd ~/ibex2024/examples/dfb
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
n=$(basename $1 .bch); o=$M/sonda_compo/$3.$n.compo_$2.txt
[ -s $o ] && exit 0
ulimit -v 6000000
if [ $3 = opt ]; then tol=1e-9; else tol=1e-6; fi
e="DFBH_FILAS=$2 DFBH_ABCAJA=1 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0"
if [ $3 = sol ]; then
  env $e timeout 200 ./ibexsolve_dfb_x $1 --lr=dfbhull --lin=compo -t 90 --seed=1 2>&1 | grep -E 'abcaja|filas|\[dfb\] piv|^[0-9]+ [0-9.]+ ' > $o
else
  env $e DFBH_LR=compo timeout 200 ./ibexopt_x $1 --filtering=acidhc4 --lr=dfbhull --seed=1 -t90 2>&1 | grep -E 'abcaja|filas|\[dfb\] piv|^[0-9]+ [0-9.]+ ' > $o
fi
