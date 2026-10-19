#!/bin/bash
# Banco de pasadas con HiGHS (sin presolve) como caja. $1 ruta $2 sol|opt
cd ~/ibex2024/examples/dfb
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
n=$(basename $1 .bch); c=$M/banco/$2.$n.txt; o=$M/banco_h/$2.$n.txt
[ -s $o ] && exit 0
ulimit -v 6000000
if [ $2 = opt ]; then tol=1e-9; else tol=1e-6; fi
DFBH_CAJA_HIGHS=1 DFBH_CAJA_SIMP=0 DFBH_ABCAJA=1 DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol timeout 3600 ./ibexreplay_h $1 $2 $c 2>&1 | grep -E 'abcaja|replay' > $o
