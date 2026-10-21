#!/bin/bash
# Tope de pivotes por LP en las peores instancias. $1 ruta $2 sol|opt $3 art|compo $4 semilla $5 tope
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/capv/$2.$n.$3.m$5.$4.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
e="DFBH_LPCAP_M=$5"
if [ $2 = sol ]; then
  r=$(env $e timeout 150 ./ibexsolve_dfb_h3 "$1" --lr=dfbhull --lin=$3 -t 60 --seed=$4 2>/dev/null | tail -1)
else
  [ $3 = compo ] && e="$e DFBH_LR=compo"
  r=$(env $e timeout 200 ./ibexopt_h3 "$1" --filtering=acidhc4 --lr=dfbhull --seed=$4 -t60 2>/dev/null | tail -1)
fi
[ -z "$r" ] && r="NA"
echo "$n,capm$5,$4,$r" > $o
