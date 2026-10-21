#!/bin/bash
cd ~/ibex2024/examples/dfb
n=$(basename $1 .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/wm/$n.$2.$3.$4.txt
[ -s $o ] && exit 0
ulimit -v 6000000
case $3 in
  dual)   e="X=1";;
  primal) e="DFB_SX_PIVOTEO=primal DFB_SX_REVPRIMAL=1 DFB_SX_DEVEXP=1 DFB_SX_DEVEXPT=1";;
  medir)  e="DFBH_WARM=medir";;
esac
env $e timeout 150 ./ibexsolve_dfb_y $1 --lr=dfbhull --lin=$2 -t 60 --seed=$4 2>&1 | grep -E "warm medir|^[0-9]+ [0-9.]+ " > $o
