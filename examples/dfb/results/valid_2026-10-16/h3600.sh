#!/bin/bash
# 102 problemas de optim/all fuera del banco del paper y no excluidos, 3600 s, binarios finales (2026-10-16).
# Reusa la corrida de 60 s si termino y la de l3600z (mismo binario, misma semilla) si existe. $1 ruta $2 brazo $3 semilla $4 opt
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=$M/h3600/$4.$n.$2.$3.csv
[ -s "$o" ] && exit 0
c=$M/h60/$4.$n.$2.$3.csv
if [ -s "$c" ] && awk -F, '{split($4,a," "); exit !(a[1]!="NA" && a[2]<59.9)}' "$c"; then cp "$c" "$o"; exit 0; fi
z=$M/l3600z/all.$4.$n.$2.$3.csv
if [ -s "$z" ]; then cp "$z" "$o"; exit 0; fi
ulimit -v 6000000
case $2 in
  sc)    e="X=1";;
  lejos) e="DFBH_ORDEN=lejos";;
  caja)  e="DFBH_ORDEN=lejos DFBH_MOTOR=caja DFBH_CAJA_TOL=1e-9";;
esac
r=$(env $e timeout 3750 ./ibexopt "$1" --filtering=acidhc4 --lr=dfbhull --seed=$3 -t3600 2>/dev/null | tail -1)
[ -z "$r" ] && r="NA"
echo "$n,$2,$3,$r" > $o
