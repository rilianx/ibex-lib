#!/bin/bash
# $1 banco (cmp|aff|opt)  $2 ruta .bch  $3 variante (base|l100): limite de iteraciones de SC y del simplex directo
R=/home/iaraya/ibex2024/results/escalera_2026-11-06; O=$(cd "$(dirname "$0")" && pwd)/run4
cd /home/iaraya/ibex2024/examples/dfb; n=$(basename $2 .bch); o=$O/$1.$n.$3.txt; [ -s $o ] && exit 0
ulimit -v 6000000
if [ $1 = opt ]; then tol=1e-9; else tol=1e-6; fi
dir="DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol DFBH_CAJA_SPX=1"
case $3 in
  base) e="DFBH_LPCAP=0 DFBH_CAJA_TOPE=igual";;
  l100) e="DFBH_LPCAP=100 DFBH_CAJA_TOPE=100";;
esac
case $1 in
  cmp) c=$R/work/S/bancocmp/sol.$n.txt; x="DFBH_REPLAY_LIN=compo"; bin=./ibexreplay_cmp; pr=sol;;
  aff) c=$R/work/M/banco/sol.$n.txt; x=""; bin=./ibexreplay_esc; pr=sol;;
  opt) c=$R/work/M/banco/opt.$n.txt; x=""; bin=./ibexreplay_esc; pr=opt;;
esac
env $dir DFB_SX_BFRT=0 DFB_SX_DEVEX=0 $e $x DFBH_ABCAJA=1 DFBH_CORTE_CERT=1 timeout 3600 $bin $2 $pr $c 2>&1 | grep -E 'abcaja|replay|estados|corte_cert|Assert|terminate' > $o
