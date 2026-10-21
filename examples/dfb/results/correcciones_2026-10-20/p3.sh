#!/bin/bash
# Punto 3: tope de tolerancia por ancho en el optimizador. $1 ruta $2 brazo (sinancho|conancho) $3 semilla
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/p3/opt.$n.$2.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
[ $2 = conancho ] && e="DFB_SX_TOLANCHO=1e-3" || e="X=1"
r=$(env $e timeout 200 ./ibexopt_y "$1" --filtering=acidhc4 --lr=dfbhull --seed=$3 -t60 2>/dev/null | tail -1)
[ -z "$r" ] && r="NA"
echo "$n,$2,$3,$r" > $o
