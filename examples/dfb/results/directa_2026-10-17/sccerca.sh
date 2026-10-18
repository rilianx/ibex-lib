#!/bin/bash
# SC con el orden de cotas "cerca" (DFBH_ORDEN=cerca), 60 s, binarios _dir2 (2026-10-18). $1 ruta $2 semilla $3 sol|opt
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/sccerca/$3.$n.sccerca.$2.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
if [ $3 = sol ]; then r=$(DFBH_ORDEN=cerca timeout 150 ./ibexsolve_dfb_dir2 "$1" --lr=dfbhull --lin=art -t 60 --seed=$2 2>/dev/null | tail -1)
else r=$(DFBH_ORDEN=cerca timeout 200 ./ibexopt_dir2 "$1" --filtering=acidhc4 --lr=dfbhull --seed=$2 -t60 2>/dev/null | tail -1); fi
[ -z "$r" ] && r="NA"
echo "$n,sccerca,$2,$r" > $o
