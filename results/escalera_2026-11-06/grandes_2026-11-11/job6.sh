#!/bin/bash
# Repeticion con el limite de 100 iteraciones por LP (SC y SoPlex) de los experimentos laterales del paper (2026-11-11).
# $1 banco (cmp|aff|optM|optS)  $2 ruta .bch  $3 variante:
#   col (SoPlex directo forzado a columnas), fil (forzado a filas), t6 (SC con tolerancia 1e-6),
#   pbase / pnrp (SoPlex parchado sin y con SOPLEX_SIN_REFACT_OPT + SOPLEX_SIN_PRECISION), psc (SC, mismo banco)
R=/home/iaraya/ibex2024/results/escalera_2026-11-06; O=$(cd "$(dirname "$0")" && pwd)/run6
cd /home/iaraya/ibex2024/examples/dfb; n=$(basename $2 .bch); o=$O/$1.$n.$3.txt; [ -s $o ] && exit 0
ulimit -v 6000000
L="DFBH_LPCAP=100 DFBH_CAJA_TOPE=100"
case $1 in
  cmp)  c=$R/work/S/bancocmp/sol.$n.txt; x="DFBH_REPLAY_LIN=compo"; bin=./ibexreplay_cmp; pr=sol; tol=1e-6;;
  aff)  c=$R/work/M/banco/sol.$n.txt; x=""; bin=./ibexreplay_esc; pr=sol; tol=1e-6;;
  optM) c=$R/work/M/banco/opt.$n.txt; x=""; bin=./ibexreplay_esc; pr=opt; tol=1e-9;;
  optS) c=$R/work/S/bancoopt/opt.$n.txt; x=""; bin=./ibexreplay_esc; pr=opt; tol=1e-9;;
esac
dir="DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol DFBH_CAJA_SPX=1"
case $3 in
  col)  e="DFBH_SPX_REP=columnas";;
  fil)  e="DFBH_SPX_REP=filas";;
  t6)   e="DFB_SX_FEASTOL=1e-6 DFB_SX_DUALTOL=1e-6";;
  pbase|psc) e=""; bin=./ibexreplay_sxp;;
  pnrp) e="SOPLEX_SIN_REFACT_OPT=1 SOPLEX_SIN_PRECISION=1"; bin=./ibexreplay_sxp;;
esac
env $dir $e DFB_SX_BFRT=0 DFB_SX_DEVEX=0 $L $x DFBH_ABCAJA=1 DFBH_CORTE_CERT=1 timeout 3600 $bin $2 $pr $c 2>&1 | grep -E 'abcaja|replay|estados|corte_cert|Assert|terminate' > $o
