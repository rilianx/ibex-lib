#!/bin/bash
cd ~/ibex2024/examples/dfb
n=$(basename $1 .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/optc/$n.$2.$3.csv
[ -s $o ] && exit 0
ulimit -v 6000000
[ $2 = t6 ] && e="DFB_SX_FEASTOL=1e-6 DFB_SX_DUALTOL=1e-6" || e="X=1"
r=$(env $e DFBH_LR=compo timeout 200 ./ibexopt_y $1 --filtering=acidhc4 --lr=dfbhull --seed=$3 -t60 2>/dev/null | tail -1)
echo "$n,$2,$3,$r" > $o
