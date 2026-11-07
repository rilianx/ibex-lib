#!/bin/bash
# Escalera de ablacion sobre el banco de pasadas (2026-11-06). $1 ruta $2 sol|opt $3 peldano
# Cada corrida reproduce las 1000 llamadas grabadas con SC (lado DFB) y una caja negra (lado caja) sobre la misma caja.
#  b0 : caja = SoPlex como lo usa Ibex (LPSolver de Ibex, cotas como filas, presolve, 1e-9, orden natural, tope 100), sin fusion; DFB = SC base sin fusion
#  b0f: lo mismo con fusion; DFB = SC base
#  b1 : caja = SoPlex ajustado (instancia propia, sin presolve, tolerancia por programa, mas cercana, 20m+200); DFB = SC base
#  b2 : caja = simplex de SoPlex directo (SPxSolver); DFB = SC con devex, paso largo y tope de 60
cd ~/ibex2024/examples/dfb
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
S=/tmp/claude-1001/-home-iaraya/2bd51c7a-9daf-4c06-b6ef-b89011baa9b5/scratchpad
n=$(basename $1 .bch); c=$M/banco/$2.$n.txt; o=$S/bescf/$2.$n.$3.txt
[ -s $o ] && exit 0
ulimit -v 6000000
if [ $2 = opt ]; then tol=1e-9; else tol=1e-6; fi
ibex="DFBH_CAJA_TOL=1e-9 DFBH_ORDEN_CAJA=nat"
dir="DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol DFBH_CAJA_TOPE=igual"
case $3 in
  b0)  e="$ibex DFBH_FUSIONA=0 DFBH_FILAS_PARALELAS=no";;
  b0f) e="$ibex";;
  b1)  e="$dir";;
  b2)  e="$dir DFBH_CAJA_SPX=1 DFB_SX_DEVEX=1 DFB_SX_BFRT=1 DFBH_LPCAP=60";;
esac
env $e DFBH_ABCAJA=1 DFBH_CORTE_CERT=1 timeout 3600 ./ibexreplay_esc $1 $2 $c 2>&1 | grep -E 'abcaja|replay' > $o
