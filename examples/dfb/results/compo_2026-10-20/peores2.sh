#!/bin/bash
# SC con DFB_SX_PISO=error DFBH_FILAS_PARALELAS=paralelas en las peores instancias. $1 ruta $2 sol|opt $3 art|compo $4 semilla
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/peores2/$2.$n.$3.$4.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
if [ $2 = sol ]; then
  r=$(DFB_SX_PISO=error DFBH_FILAS_PARALELAS=paralelas timeout 150 ./ibexsolve_dfb_x "$1" --lr=dfbhull --lin=$3 -t 60 --seed=$4 2>/dev/null | tail -1)
else
  [ $3 = compo ] && x="DFBH_LR=compo" || x="X=1"
  r=$(env DFB_SX_PISO=error DFBH_FILAS_PARALELAS=paralelas $x timeout 200 ./ibexopt_x "$1" --filtering=acidhc4 --lr=dfbhull --seed=$4 -t60 2>/dev/null | tail -1)
fi
[ -z "$r" ] && r="NA"
echo "$n,$3,$4,$r" > $o
