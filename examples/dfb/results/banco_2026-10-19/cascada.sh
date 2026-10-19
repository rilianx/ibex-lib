#!/bin/bash
# Cascada de atribucion sobre el banco de pasadas. $1 ruta $2 sol|opt $3 variante
cd ~/ibex2024/examples/dfb
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
n=$(basename $1 .bch); c=$M/banco/$2.$n.txt; o=$M/cascada/$2.$n.$3.txt
[ -s $o ] && exit 0
ulimit -v 6000000
if [ $2 = opt ]; then tol=1e-9; else tol=1e-6; fi
b=./ibexreplay
case $3 in
  sx_verif)  b=./ibexreplay_sx; e="SOPLEX_SIN_VERIF=1";;
  sx_refact) b=./ibexreplay_sx; e="SOPLEX_SIN_VERIF=1 SOPLEX_REFACT_FIN=1000000";;
  sx_presol) b=./ibexreplay_sx; e="SOPLEX_SIN_VERIF=1 SOPLEX_REFACT_FIN=1000000 DFBH_CAJA_SIMP=0";;
  sx_escala) b=./ibexreplay_sx; e="SOPLEX_SIN_VERIF=1 SOPLEX_REFACT_FIN=1000000 DFBH_CAJA_SIMP=0 DFBH_CAJA_SCALER=0";;
  sc_fact)   e="DFB_EMU_FACT=1";;
  sc_init)   e="DFB_EMU_FACT=1 DFB_EMU_INIT=1";;
  sc_verif)  e="DFB_EMU_FACT=1 DFB_EMU_INIT=1 DFB_EMU_VERIF=1";;
  sc_copia)  e="DFB_EMU_FACT=1 DFB_EMU_INIT=1 DFB_EMU_VERIF=1 DFB_EMU_COPIA=1";;
esac
env $e DFBH_ABCAJA=1 DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol timeout 1800 $b $1 $2 $c 2>&1 | grep -E 'abcaja|replay' > $o
