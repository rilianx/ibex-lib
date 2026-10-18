#!/bin/bash
# $1 ruta $2 sol|opt $3 nivel
cd ~/ibex2024/examples/dfb
n=$(basename $1 .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/emu/$2.$n.$3.txt
[ -s $o ] && exit 0
case $3 in
  base)  e="X=1";;
  fact)  e="DFB_EMU_FACT=1";;
  init)  e="DFB_EMU_FACT=1 DFB_EMU_INIT=1";;
  verif) e="DFB_EMU_VERIF=1";;
  copia) e="DFB_EMU_COPIA=1";;
  todo)  e="DFB_EMU_FACT=1 DFB_EMU_INIT=1 DFB_EMU_VERIF=1 DFB_EMU_COPIA=1";;
esac
c="DFBH_ABCAJA=1 DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_ORDEN_CAJA=cerca"
if [ $2 = sol ]; then env $e $c timeout 100 ./ibexsolve_dfb_emu $1 --lr=dfbhull --lin=art -t 20 --seed=1 2>&1 | grep 'costo\|abcaja\] DFB' > $o
else env $e $c DFBH_CAJA_TOL=1e-9 timeout 120 ./ibexopt_emu $1 --filtering=acidhc4 --lr=dfbhull --seed=1 -t20 2>&1 | grep 'costo\|abcaja\] DFB' > $o; fi
