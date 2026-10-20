#!/bin/bash
# SC con DFBH_ORDEN=costo_o_cerca (conteo con tableau, mas cercana con el revisado), sistemas, 60 s, binario _x. $1 ruta $2 brazo $3 semilla $4 art|compo
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/d60coc/sol.$n.$2.$4.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
r=$(DFBH_ORDEN=costo_o_cerca timeout 150 ./ibexsolve_dfb_x "$1" --lr=dfbhull --lin=$4 -t 60 --seed=$3 2>/dev/null | tail -1)
[ -z "$r" ] && r="NA"
echo "$n,$2,$3,$r" > $o
