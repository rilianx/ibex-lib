#!/bin/bash
# SC (orden cerca por omision) a 3600 s en los 24 problemas de opt de tab:long (2026-10-18). $1 ruta $2 semilla
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/l3600c/opt.$n.sc.$2.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
r=$(timeout 3750 ./ibexopt_dir3 "$1" --filtering=acidhc4 --lr=dfbhull --seed=$2 -t3600 2>/dev/null | tail -1)
[ -z "$r" ] && r="NA"
echo "$n,sc,$2,$r" > $o
