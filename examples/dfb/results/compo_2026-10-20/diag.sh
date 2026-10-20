#!/bin/bash
# Diagnosticos de compo en sistemas. $1 ruta $2 art|compo $3 variante
cd ~/ibex2024/examples/dfb
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
n=$(basename $1 .bch); o=$M/diag_compo/$n.$2.$3.txt
[ -s $o ] && exit 0
ulimit -v 6000000
caja="DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=1e-6 DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0"
case $3 in
  resumen)     e="DFBH_ABCAJA=1 $caja";;
  sombra)      e="DFBH_SOMBRA=1 $caja";;
  sombra_base) e="DFBH_SOMBRA=1 DFBH_SOMBRA_BASE=1 $caja";;
  sinlargo)    e="DFB_SX_BFRT=0 DFBH_ABCAJA=1 $caja";;
  sinpert)     e="DFB_SX_PERTURB=0 DFBH_ABCAJA=1 $caja";;
  sindevex)    e="DFB_SX_DEVEX=0 DFBH_ABCAJA=1 $caja";;
  revisado)    e="DFB_SX_MOTOR=revisado DFBH_ABCAJA=1 $caja";;
  tableau)     e="DFB_SX_MOTOR=tableau DFBH_ABCAJA=1 $caja";;
  sinfusion)   e="DFBH_FUSIONA=0 DFBH_ABCAJA=1 $caja";;
  sintau)      e="DFB_SX_TOLANCHO=0 DFBH_ABCAJA=1 $caja";;
esac
env $e timeout 120 ./ibexsolve_dfb_x $1 --lr=dfbhull --lin=$2 -t 40 --seed=1 2>&1 | grep -E 'abcaja|sombra|\[dfb\]|^[0-9]+ [0-9.]+ ' > $o
