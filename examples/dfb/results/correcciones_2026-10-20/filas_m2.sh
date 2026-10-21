#!/bin/bash
# filas de la relajacion por instancia. $1 ruta $2 art|compo $3 sol|opt
cd ~/ibex2024/examples/dfb
n=$(basename $1 .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/filas_m2/$3.$n.$2.txt
[ -s $o ] && exit 0
ulimit -v 6000000
if [ $3 = sol ]; then
  timeout 40 ./ibexsolve_dfb_y $1 --lr=dfbhull --lin=$2 -t 8 --seed=1 2>&1 | grep -E "filas de la relajacion" > $o
else
  [ $2 = compo ] && e="DFBH_LR=compo" || e="X=1"
  env $e timeout 40 ./ibexopt_y $1 --filtering=acidhc4 --lr=dfbhull --seed=1 -t8 2>&1 | grep -E "filas de la relajacion" > $o
fi
[ -s $o ] || echo "sin dato" > $o
