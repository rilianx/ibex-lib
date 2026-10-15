#!/bin/bash
# Banco final del paper con reintento en frio, binarios finales (2026-10-15). $1 ruta $2 brazo $3 semilla $4 sol|opt
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/f60/$4.$n.$2.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
if [ $4 = sol ]; then tc=1e-6; tibio="DFB_SX_PIVOTEO=primal DFB_SX_REVPRIMAL=1 DFB_SX_MOTOR_M=15"
else tc=1e-9; tibio="DFB_SX_PIVOTEO=dual"; fi
case $2 in
  sc)        e="X=1";;
  reint)     e="DFBH_REINTENTA=1";;
  scoc)      e="X=1";;
  lejos)     e="DFBH_ORDEN=lejos";;
  caja)      e="DFBH_ORDEN=lejos DFBH_MOTOR=caja DFBH_CAJA_TOL=$tc";;
  cajacerca) e="DFBH_ORDEN=cerca DFBH_MOTOR=caja DFBH_CAJA_TOL=$tc";;
  tibio)     e="$tibio";;
  sinpert)   e="DFB_SX_PERTURB=0";;
  sinlargo)  e="DFB_SX_BFRT=0";;
  sintau)    e="DFB_SX_TOLANCHO=0";;
  sindevex)  e="DFB_SX_DEVEX=0";;
  tableau)   e="DFB_SX_MOTOR=tableau";;
  revisado)  e="DFB_SX_MOTOR=revisado";;
  sinfusion) e="DFBH_FUSIONA=0";;
esac
if [ $4 = sol ]; then
  r=$(env $e timeout 150 ./ibexsolve_dfb "$1" --lr=dfbhull --lin=art -t 60 --seed=$3 2>/dev/null | tail -1)
else
  r=$(env $e timeout 200 ./ibexopt "$1" --filtering=acidhc4 --lr=dfbhull --seed=$3 -t60 2>/dev/null | tail -1)
fi
[ -z "$r" ] && r="NA"
echo "$n,$2,$3,$r" > $o
