#!/bin/bash
# Sonda de lados a 60 s sobre el banco unico (2026-10-16), binarios con el arreglo afin. $1 ruta $2 brazo(dual|primal|sc) $3 sol|opt
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/lados60/$3.$n.$2.txt
[ -s "$o" ] && exit 0
ulimit -v 6000000
case $2 in
  sc)     e="X=1";;
  dual)   e="DFB_SX_PIVOTEO=dual";;
  primal) e="DFB_SX_PIVOTEO=primal DFB_SX_REVPRIMAL=1 DFB_SX_MOTOR_M=15";;
esac
if [ $3 = sol ]; then
  env $e DFBH_BASECOMP=200 timeout 200 ./ibexsolve_dfb_fix "$1" --lr=dfbhull --lin=art -t 60 --seed=1 2>&1 | grep basecomp | tail -2 > $o
else
  env $e DFBH_BASECOMP=200 timeout 250 ./ibexopt_fix "$1" --filtering=acidhc4 --lr=dfbhull --seed=1 -t60 2>&1 | grep basecomp | tail -2 > $o
fi
[ -s "$o" ] || echo NA > $o
