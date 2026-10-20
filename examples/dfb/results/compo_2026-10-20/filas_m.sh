#!/bin/bash
# filas de la relajacion por instancia (1 s), afin y compo. $1 ruta $2 art|compo
cd ~/ibex2024/examples/dfb
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
n=$(basename $1 .bch); o=$M/filas_m/$n.$2.txt
[ -s $o ] && exit 0
timeout 20 ./ibexsolve_dfb_x $1 --lr=dfbhull --lin=$2 -t 1 --seed=1 2>&1 | grep -E "filas de la relajacion|nb_var" > $o
