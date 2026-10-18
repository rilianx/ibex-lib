#!/bin/bash
# Ablacion de opt con el orden cerca por omision del optimizador (2026-10-18), ibexopt_dir3. $1 ruta $2 brazo $3 semilla
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/abl3/opt.$n.$2.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
case $2 in
  tibio)     e="DFB_SX_PIVOTEO=dual";;
  sinpert)   e="DFB_SX_PERTURB=0";;
  sinlargo)  e="DFB_SX_BFRT=0";;
  sintau)    e="DFB_SX_TOLANCHO=0";;
  sindevex)  e="DFB_SX_DEVEX=0";;
  tableau)   e="DFB_SX_MOTOR=tableau";;
  revisado)  e="DFB_SX_MOTOR=revisado";;
  sinfusion) e="DFBH_FUSIONA=0";;
esac
r=$(env $e timeout 200 ./ibexopt_dir3 "$1" --filtering=acidhc4 --lr=dfbhull --seed=$3 -t60 2>/dev/null | tail -1)
[ -z "$r" ] && r="NA"
echo "$n,$2,$3,$r" > $o
