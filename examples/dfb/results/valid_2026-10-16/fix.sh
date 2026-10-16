#!/bin/bash
# Reejecuta con ibexopt_fix (guarda is_actif en el linearizador afin) los casos que abortaban (2026-10-16).
# $1 ruta $2 brazo $3 semilla $4 limite(60|3600) $5 etiqueta
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/fix/$5.opt.$n.$2.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
tibio="DFB_SX_PIVOTEO=dual"
case $2 in
  sc)        e="X=1";;
  lejos)     e="DFBH_ORDEN=lejos";;
  caja)      e="DFBH_ORDEN=lejos DFBH_MOTOR=caja DFBH_CAJA_TOL=1e-9";;
  tibio)     e="$tibio";;
  sinpert)   e="DFB_SX_PERTURB=0";;
  sinlargo)  e="DFB_SX_BFRT=0";;
  sintau)    e="DFB_SX_TOLANCHO=0";;
  sindevex)  e="DFB_SX_DEVEX=0";;
  tableau)   e="DFB_SX_MOTOR=tableau";;
  revisado)  e="DFB_SX_MOTOR=revisado";;
  sinfusion) e="DFBH_FUSIONA=0";;
esac
if [ $4 = 60 ]; then to=200; else to=3750; fi
r=$(env $e timeout $to ./ibexopt_fix "$1" --filtering=acidhc4 --lr=dfbhull --seed=$3 -t$4 2>/dev/null | tail -1)
[ -z "$r" ] && r="NA"
echo "$n,$2,$3,$r" > $o
