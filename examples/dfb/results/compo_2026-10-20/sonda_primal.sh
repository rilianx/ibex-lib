#!/bin/bash
# Sonda en la busqueda (DFBH_ABCAJA): SC y caja (SoPlex sin presolve) sobre los mismos LPs, afin o compo. $1 ruta $2 art|compo $3 sol|opt
cd ~/ibex2024/examples/dfb
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
n=$(basename $1 .bch); o=$M/sonda_compo/$3.$n.$2_primal.txt
[ -s $o ] && exit 0
ulimit -v 6000000
if [ $3 = opt ]; then tol=1e-9; else tol=1e-6; fi
e="DFB_SX_PIVOTEO=primal DFB_SX_REVPRIMAL=1 DFBH_ABCAJA=1 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0"
if [ $3 = sol ]; then
  env $e timeout 200 ./ibexsolve_dfb_x $1 --lr=dfbhull --lin=$2 -t 90 --seed=1 2>&1 | grep -E 'abcaja|^[0-9]' | tail -6 > $o
else
  [ $2 = compo ] && e="$e DFBH_LR=compo"
  env $e timeout 200 ./ibexopt_x $1 --filtering=acidhc4 --lr=dfbhull --seed=1 -t90 2>&1 | grep -E 'abcaja|^[0-9]' | tail -6 > $o
fi
