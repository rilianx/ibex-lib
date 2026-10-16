#!/bin/bash
# 3600 s con binarios finales, reintento en frio en los dos motores (2026-10-15).
# $1 ruta $2 banco(l|all) $3 semilla $4 sol|opt $5 brazo: sc (orden propio), lejos (SC lejos), caja (SoPlex lejos)
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/l3600z/$2.$4.$n.$5.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
if [ $4 = sol ]; then tc=1e-6; else tc=1e-9; fi
case $5 in
  sc)    e="X=1";;
  lejos) e="DFBH_ORDEN=lejos";;
  caja)  e="DFBH_ORDEN=lejos DFBH_MOTOR=caja DFBH_CAJA_TOL=$tc";;
esac
if [ $4 = sol ]; then
  r=$(env $e timeout 3700 ./ibexsolve_dfb "$1" --lr=dfbhull --lin=art -t 3600 --seed=$3 2>/dev/null | tail -1)
else
  r=$(env $e timeout 3750 ./ibexopt "$1" --filtering=acidhc4 --lr=dfbhull --seed=$3 -t3600 2>/dev/null | tail -1)
fi
[ -z "$r" ] && r="NA"
echo "$n,$5,$3,$r" > $o
