#!/bin/bash
# Diagnosticos compo, etapa 2: escalado, orden, perturbacion. $1 ruta $2 variante (relajacion compo, sistemas)
cd ~/ibex2024/examples/dfb
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
n=$(basename $1 .bch); o=$M/diag2/$n.$2.txt
[ -s $o ] && exit 0
ulimit -v 6000000
caja="DFBH_ABCAJA=1 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=1e-6 DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_CAJA_TOPE=igual"
case $2 in
  base)        e="$caja";;
  caja_sinesc) e="DFBH_CAJA_SCALER=0 $caja";;
  sc_escala)   e="DFB_SX_SCALE=1 $caja";;
  sc_cerca)    e="DFBH_ORDEN=cerca $caja";;
  sc_lejos)    e="DFBH_ORDEN=lejos $caja";;
  pert1e5)     e="DFB_SX_PERTURB=1e-5 $caja";;
  pert1e3)     e="DFB_SX_PERTURB=1e-3 $caja";;
  dump)        e="DFBH_CAJA_DUMP=$M/diag2/dump $caja";;
esac
env $e timeout 120 ./ibexsolve_dfb_x $1 --lr=dfbhull --lin=compo -t 40 --seed=1 2>&1 | grep -E 'abcaja|\[dfb\] piv|^[0-9]+ [0-9.]+ ' > $o
