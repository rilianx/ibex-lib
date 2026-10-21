#!/bin/bash
# Remedios 1 (tolerancia por motor) y 2 (cambios de cota incrementales). $1 ruta $2 art|compo $3 variante $4 semilla
cd ~/ibex2024/examples/dfb
n=$(basename $1 .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/r12/$n.$2.$3.$4.csv
[ -s $o ] && exit 0
ulimit -v 6000000
case $3 in
  base)    e="X=1";;
  tolrev)  e="DFB_SX_TOL_REVISADO=1e-6";;
  flipinc) e="DFB_SX_FLIPINC=1";;
  ambos)   e="DFB_SX_TOL_REVISADO=1e-6 DFB_SX_FLIPINC=1";;
esac
r=$(env $e timeout 150 ./ibexsolve_dfb_y $1 --lr=dfbhull --lin=$2 -t 60 --seed=$4 2>/dev/null | tail -1)
echo "$n,$3,$4,$r" > $o
