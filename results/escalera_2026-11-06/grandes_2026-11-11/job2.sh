#!/bin/bash
# $1 banco (cmp|aff)  $2 ruta .bch  $3 variante
R=/home/iaraya/ibex2024/results/escalera_2026-11-06; O=$(cd "$(dirname "$0")" && pwd)/run2
cd /home/iaraya/ibex2024/examples/dfb; n=$(basename $2 .bch); o=$O/$1.$n.$3.txt; [ -s $o ] && exit 0
ulimit -v 6000000
dir="DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=1e-6 DFBH_CAJA_TOPE=igual DFBH_CAJA_SPX=1"
case $3 in
  c0)   e="DFBH_LPCAP=0";;
  c60)  e="DFBH_LPCAP=60";;
  c150) e="DFBH_LPCAP=150";;
  t6)   e="DFBH_LPCAP=0 DFB_SX_FEASTOL=1e-6 DFB_SX_DUALTOL=1e-6";;
  ralo) e="DFBH_LPCAP=0 DFB_SX_BTRAN=ralo";;
esac
if [ $1 = cmp ]; then c=$R/work/S/bancocmp/sol.$n.txt; x="DFBH_REPLAY_LIN=compo"; bin=./ibexreplay_cmp; else c=$R/work/M/banco/sol.$n.txt; x=""; bin=./ibexreplay_esc; fi
env $dir DFB_SX_BFRT=0 DFB_SX_DEVEX=0 $e $x DFBH_ABCAJA=1 DFBH_CORTE_CERT=1 timeout 3600 $bin $2 sol $c 2>&1 | grep -E 'abcaja|replay|estados|corte_cert|degen|Assert|terminate' > $o
