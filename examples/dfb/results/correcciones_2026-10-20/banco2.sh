#!/bin/bash
# Banco de pasadas con SC corregido (ibexreplay_h2, 2026-10-21). $1 ruta $2 sol|opt $3 caja (sp|spx|gurobi|highs)
cd ~/ibex2024/examples/dfb
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
n=$(basename $1 .bch); c=$M/banco/$2.$n.txt; o=$M/banco2/$2.$n.$3.txt
[ -s $o ] && exit 0
ulimit -v 6000000
if [ $2 = opt ]; then tol=1e-9; else tol=1e-6; fi
case $3 in
  sp)     e="DFBH_CAJA_SIMP=0";;
  spx)    e="DFBH_CAJA_SPX=1";;
  gurobi) e="DFBH_CAJA_GUROBI=1 DFBH_CAJA_SIMP=0";;
  highs)  e="DFBH_CAJA_HIGHS=1 DFBH_CAJA_SIMP=0";;
esac
env $e DFBH_ABCAJA=1 DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol timeout 3600 ./ibexreplay_h2 $1 $2 $c 2>&1 | grep -E 'abcaja|replay' > $o
