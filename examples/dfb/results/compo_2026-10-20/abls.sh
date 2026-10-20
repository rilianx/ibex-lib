#!/bin/bash
# Ablacion de sistemas con el orden nuevo por omision (costo_o_cerca), 60 s, binario _x (2026-10-20). $1 ruta $2 brazo $3 semilla
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/abls/sol.$n.$2.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
case $2 in
  tibio)     e="DFB_SX_PIVOTEO=primal DFB_SX_REVPRIMAL=1 DFB_SX_MOTOR_M=15";;
  sinpert)   e="DFB_SX_PERTURB=0";;
  sinlargo)  e="DFB_SX_BFRT=0";;
  sintau)    e="DFB_SX_TOLANCHO=0";;
  sindevex)  e="DFB_SX_DEVEX=0";;
  tableau)   e="DFB_SX_MOTOR=tableau";;
  revisado)  e="DFB_SX_MOTOR=revisado";;
  sinfusion) e="DFBH_FUSIONA=0";;
  lejos)     e="DFBH_ORDEN=lejos";;
  cerca)     e="DFBH_ORDEN=cerca";;
  conteo)    e="DFBH_ORDEN=costo";;
esac
r=$(env $e timeout 150 ./ibexsolve_dfb_x "$1" --lr=dfbhull --lin=art -t 60 --seed=$3 2>/dev/null | tail -1)
[ -z "$r" ] && r="NA"
echo "$n,$2,$3,$r" > $o
