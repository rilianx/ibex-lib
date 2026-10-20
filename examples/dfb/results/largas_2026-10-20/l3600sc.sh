#!/bin/bash
# SC en el solver de sistemas con la cota mas cercana primero (DFBH_ORDEN=cerca) a 3600 s, binarios _dir2 (2026-10-20).
# $1 ruta $2 brazo $3 semilla $4 sol
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/l3600sc/$4.$n.$2.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
r=$(DFBH_ORDEN=cerca timeout 3700 ./ibexsolve_dfb_dir2 "$1" --lr=dfbhull --lin=art -t 3600 --seed=$3 2>/dev/null | tail -1)
[ -z "$r" ] && r="NA"
echo "$n,$2,$3,$r" > $o
