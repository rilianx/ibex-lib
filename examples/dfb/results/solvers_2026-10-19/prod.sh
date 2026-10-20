#!/bin/bash
# PHull de produccion de Ibex (CtcPolytopeHull, --lr=xn) a 60 s, binarios _dir2 (2026-10-19).
# Solver: SoPlex con columnas (IBEX_LP_REP=col, la produccion mas fuerte); ibexopt: de fabrica.
# $1 ruta $2 brazo $3 semilla $4 sol|opt
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/prod/$4.$n.$2.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
if [ $4 = sol ]; then
  r=$(IBEX_LP_REP=col timeout 150 ./ibexsolve_dfb_dir2 "$1" --lr=xn --lin=art -t 60 --seed=$3 2>/dev/null | tail -1)
else
  r=$(timeout 200 ./ibexopt_dir2 "$1" --filtering=acidhc4 --lr=xn --seed=$3 -t60 2>/dev/null | tail -1)
fi
[ -z "$r" ] && r="NA"
echo "$n,$2,$3,$r" > $o
