#!/bin/bash
cd ~/ibex2024/examples/dfb
n=$(basename $1 .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/cap/$n.$2.$3.csv
[ -s $o ] && exit 0
ulimit -v 6000000
[ $2 = sin ] && e="X=1" || e="DFBH_LPCAP=$2"
r=$(env $e timeout 150 ./ibexsolve_dfb_h2 $1 --lr=dfbhull --lin=compo -t 60 --seed=$3 2>/dev/null | tail -1)
echo "$n,$2,$3,$r" > $o
