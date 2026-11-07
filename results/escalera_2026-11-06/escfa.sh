#!/bin/bash
# Escalera con la relajacion afin en sistemas (2026-11-06): R3 SC base, R4 +tope 60, R5 +devex, R6 +paso largo (=SC final). $1 ruta $2 sol|opt $3 semilla
cd ~/ibex2024/examples/dfb
S=/tmp/claude-1001/-home-iaraya/2bd51c7a-9daf-4c06-b6ef-b89011baa9b5/scratchpad/escfa
n=$(basename "$1" .bch); ulimit -v 6000000
if [ $2 = sol ]; then tc=1e-6; lin=art; pel="R0 R0f R1 R2 R3 R4 R5 R6"; else tc=1e-9; lin=art; pel="R0 R0f R1 R2 R3 R4 R5 R6"; fi
caja_ibex="DFBH_MOTOR=caja DFBH_CAJA_TOL=1e-9 DFBH_ORDEN_CAJA=nat DFBH_CORTE_CERT=1"
caja_dir="DFBH_MOTOR=caja DFBH_CAJA_TOL=$tc DFBH_ORDEN_CAJA=cerca DFBH_CORTE_CERT=1 DFBH_CAJA_SIMP=0 DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_TOPE=igual"
nofus="DFBH_FUSIONA=0 DFBH_FILAS_PARALELAS=no"
for p in $pel; do
  o=$S/$2.$n.$p.$3.csv; [ -s "$o" ] && continue
  case $p in
    R0)  e="$caja_ibex $nofus";;
    R0f) e="$caja_ibex";;
    R1)  e="$caja_dir";;
    R2)  e="$caja_dir DFBH_CAJA_SPX=1";;
    R3)  e="DFB_SX_BFRT=0 DFB_SX_DEVEX=0 DFBH_LPCAP=0";;
    R4)  e="DFB_SX_BFRT=0 DFB_SX_DEVEX=0";;
    R5)  e="DFB_SX_BFRT=0";;
    R6)  e="X=1";;
    R7)  e="DFBH_HIBRIDO=30";;
  esac
  if [ $2 = sol ]; then
    r=$(env $e timeout 150 ./ibexsolve_dfb_esc "$1" --lr=dfbhull --lin=$lin -t 60 --seed=$3 2>/dev/null | tail -1)
  else
    r=$(env $e timeout 200 ./ibexopt_esc "$1" --filtering=acidhc4 --lr=dfbhull --seed=$3 -t60 2>/dev/null | tail -1)
  fi
  [ -z "$r" ] && r="NA"
  echo "$n,$p,$3,$r" > $o
done
