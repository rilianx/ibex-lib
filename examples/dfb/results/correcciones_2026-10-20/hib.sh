#!/bin/bash
# Hibrido con umbral T. $1 ruta $2 art|compo $3 T $4 semilla $5 sol|opt
cd ~/ibex2024/examples/dfb
n=$(basename $1 .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/hib/$5.$n.$2.T$3.$4.csv
[ -s $o ] && exit 0
ulimit -v 6000000
if [ $5 = sol ]; then
  r=$(DFBH_HIBRIDO=$3 timeout 150 ./ibexsolve_dfb_h2 $1 --lr=dfbhull --lin=$2 -t 60 --seed=$4 2>/dev/null | tail -1)
else
  [ $2 = compo ] && e="DFBH_LR=compo" || e="X=1"
  r=$(env $e DFBH_HIBRIDO=$3 timeout 200 ./ibexopt_h2 $1 --filtering=acidhc4 --lr=dfbhull --seed=$4 -t60 2>/dev/null | tail -1)
fi
echo "$n,T$3,$4,$r" > $o
