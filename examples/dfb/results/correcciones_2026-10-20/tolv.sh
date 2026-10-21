#!/bin/bash
# Variantes de tolerancia de SC. $1 ruta $2 art|compo $3 variante
cd ~/ibex2024/examples/dfb
n=$(basename $1 .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/tolv/$n.$2.$3.csv
[ -s $o ] && exit 0
ulimit -v 6000000
case $3 in
  t6)        e="X=1";;
  t9_ancho)  e="DFB_SX_FEASTOL=1e-9 DFB_SX_DUALTOL=1e-9 DFB_SX_TOLANCHO=1e-3";;
  t9)        e="DFB_SX_FEASTOL=1e-9 DFB_SX_DUALTOL=1e-9";;
  t9p6d)     e="DFB_SX_FEASTOL=1e-9 DFB_SX_DUALTOL=1e-6";;
  t6p9d)     e="DFB_SX_FEASTOL=1e-6 DFB_SX_DUALTOL=1e-9";;
  t7)        e="DFB_SX_FEASTOL=1e-7 DFB_SX_DUALTOL=1e-7";;
  t8)        e="DFB_SX_FEASTOL=1e-8 DFB_SX_DUALTOL=1e-8";;
esac
r=$(env $e timeout 150 ./ibexsolve_dfb_y $1 --lr=dfbhull --lin=$2 -t 60 --seed=1 2>/dev/null | tail -1)
echo "$n,$3,$r" > $o
