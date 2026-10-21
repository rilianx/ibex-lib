#!/bin/bash
# Contractor hibrido (DFBH_HIBRIDO=30) con las correcciones (binarios _h2, 2026-10-21). $1 ruta $2 brazo (hib_art|hib_compo) $3 semilla $4 sol|opt
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/d60h/$4.$n.$2.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
lin=${2#*_}
e="DFBH_HIBRIDO=30"
if [ $4 = sol ]; then
  r=$(env $e timeout 150 ./ibexsolve_dfb_h2 "$1" --lr=dfbhull --lin=$lin -t 60 --seed=$3 2>/dev/null | tail -1)
else
  [ $lin = compo ] && e="$e DFBH_LR=compo"
  r=$(env $e timeout 200 ./ibexopt_h2 "$1" --filtering=acidhc4 --lr=dfbhull --seed=$3 -t60 2>/dev/null | tail -1)
fi
[ -z "$r" ] && r="NA"
echo "$n,$2,$3,$r" > $o
