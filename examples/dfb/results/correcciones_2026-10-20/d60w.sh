#!/bin/bash
# Bancos a 60 s con las correcciones por omision (binarios _y, 2026-10-20).
# $1 ruta $2 brazo (sc_art|sc_compo|caja_art|caja_compo) $3 semilla $4 sol|opt
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/d60w/$4.$n.$2.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
if [ $4 = sol ]; then tc=1e-6; else tc=1e-9; fi
motor=${2%_*}; lin=${2#*_}
e="X=1"
[ $motor = caja ] && e="DFBH_ORDEN=cerca DFBH_MOTOR=caja DFBH_CAJA_TOL=$tc DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0"
if [ $4 = sol ]; then
  r=$(env $e timeout 150 ./ibexsolve_dfb_y "$1" --lr=dfbhull --lin=$lin -t 60 --seed=$3 2>/dev/null | tail -1)
else
  [ $lin = compo ] && e="$e DFBH_LR=compo"
  r=$(env $e timeout 200 ./ibexopt_y "$1" --filtering=acidhc4 --lr=dfbhull --seed=$3 -t60 2>/dev/null | tail -1)
fi
[ -z "$r" ] && r="NA"
echo "$n,$2,$3,$r" > $o
