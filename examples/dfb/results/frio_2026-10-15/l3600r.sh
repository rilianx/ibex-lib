#!/bin/bash
# SC a 3600 s con el binario final con frio_si_no_certifica (2026-10-15). $1 ruta $2 banco(l|all) $3 semilla $4 sol|opt
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/l3600r/$2.$4.$n.propio.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
if [ $4 = sol ]; then
  r=$(timeout 3700 ./ibexsolve_dfb "$1" --lr=dfbhull --lin=art -t 3600 --seed=$3 2>/dev/null | tail -1)
else
  r=$(timeout 3750 ./ibexopt "$1" --filtering=acidhc4 --lr=dfbhull --seed=$3 -t3600 2>/dev/null | tail -1)
fi
[ -z "$r" ] && r="NA"
echo "$n,propio,$3,$r" > $o
