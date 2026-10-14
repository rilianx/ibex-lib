#!/bin/bash
# SoPlex (caja) con otros ordenes de cotas, 60 s. $1 ruta $2 orden $3 semilla $4 sol|opt
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/ordcaja/$4.$n.caja$2.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
if [ $4 = sol ]; then tc=1e-6; else tc=1e-9; fi
e="DFBH_ORDEN=$2 DFBH_MOTOR=caja DFBH_CAJA_TOL=$tc"
if [ $4 = sol ]; then
  r=$(env $e timeout 150 ./ibexsolve_dfb "$1" --lr=dfbhull --lin=art -t 60 --seed=$3 2>/dev/null | tail -1)
else
  r=$(env $e timeout 200 ./ibexopt "$1" --filtering=acidhc4 --lr=dfbhull --seed=$3 -t60 2>/dev/null | tail -1)
fi
[ -z "$r" ] && r="NA"
echo "$n,caja$2,$3,$r" > $o
