#!/bin/bash
# R0/R0f en el solver con el orden de Ibex (mas cercana) (2026-11-07). $1 ruta $2 sol|opt $3 semilla $4 escf|escfa
H=$(cd "$(dirname "$0")" && pwd)   # this folder; work/ holds datos_trabajo.tgz unpacked
cd "$H/../../examples/dfb"
S=${S:-$H/work/S}/$4_r0c
n=$(basename "$1" .bch); ulimit -v 6000000
if [ $4 = escfa ]; then lin=art; else lin=compo; fi
caja_ibex="DFBH_MOTOR=caja DFBH_CAJA_TOL=1e-9 DFBH_ORDEN_CAJA=cerca DFBH_CORTE_CERT=1"
nofus="DFBH_FUSIONA=0 DFBH_FILAS_PARALELAS=no"
for p in R0 R0f; do
  o=$S/$2.$n.$p.$3.csv; [ -s "$o" ] && continue
  case $p in R0) e="$caja_ibex $nofus";; R0f) e="$caja_ibex";; esac
  if [ $2 = sol ]; then
    r=$(env $e timeout 150 ./ibexsolve_dfb_esc "$1" --lr=dfbhull --lin=$lin -t 60 --seed=$3 2>/dev/null | tail -1)
  else
    r=$(env $e timeout 200 ./ibexopt_esc "$1" --filtering=acidhc4 --lr=dfbhull --seed=$3 -t60 2>/dev/null | tail -1)
  fi
  [ -z "$r" ] && r="NA"
  echo "$n,$p,$3,$r" > $o
done
