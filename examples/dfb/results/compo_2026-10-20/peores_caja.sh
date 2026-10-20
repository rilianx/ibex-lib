#!/bin/bash
# Caja negra (SoPlex sin presolve, tope igual) con la fusion de filas paralelas, compo. $1 ruta $2 sol|opt $3 compo $4 semilla
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/peores_caja/$2.$n.$3.$4.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
if [ $2 = sol ]; then tc=1e-6; else tc=1e-9; fi
e="DFBH_FILAS_PARALELAS=paralelas DFBH_ORDEN=cerca DFBH_MOTOR=caja DFBH_CAJA_TOL=$tc DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_CAJA_TOPE=igual"
if [ $2 = sol ]; then
  r=$(env $e timeout 150 ./ibexsolve_dfb_x "$1" --lr=dfbhull --lin=compo -t 60 --seed=$4 2>/dev/null | tail -1)
else
  r=$(env $e DFBH_LR=compo timeout 200 ./ibexopt_x "$1" --filtering=acidhc4 --lr=dfbhull --seed=$4 -t60 2>/dev/null | tail -1)
fi
[ -z "$r" ] && r="NA"
echo "$n,$3,$4,$r" > $o
