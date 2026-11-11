#!/bin/bash
# $1 banco (optM|optS|aff)  $2 ruta .bch  $3 variante del SoPlex parchado (base|nr|np|nrp)
R=/home/iaraya/ibex2024/results/escalera_2026-11-06; O=$(cd "$(dirname "$0")" && pwd)/run5
cd /home/iaraya/ibex2024/examples/dfb; n=$(basename $2 .bch); o=$O/$1.$n.$3.txt; [ -s $o ] && exit 0
ulimit -v 6000000
case $1 in optM) c=$R/work/M/banco/opt.$n.txt; pr=opt; tol=1e-9;; optS) c=$R/work/S/bancoopt/opt.$n.txt; pr=opt; tol=1e-9;; aff) c=$R/work/M/banco/sol.$n.txt; pr=sol; tol=1e-6;; esac
case $3 in base) e="";; nr) e="SOPLEX_SIN_REFACT_OPT=1";; np) e="SOPLEX_SIN_PRECISION=1";; nrp) e="SOPLEX_SIN_REFACT_OPT=1 SOPLEX_SIN_PRECISION=1";; esac
dir="DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol DFBH_CAJA_TOPE=igual DFBH_CAJA_SPX=1"
env $dir $e DFB_SX_BFRT=0 DFB_SX_DEVEX=0 DFBH_LPCAP=0 DFBH_ABCAJA=1 DFBH_CORTE_CERT=1 timeout 3600 ./ibexreplay_sxp $2 $pr $c 2>&1 | grep -E 'abcaja|replay|Assert|terminate' > $o
