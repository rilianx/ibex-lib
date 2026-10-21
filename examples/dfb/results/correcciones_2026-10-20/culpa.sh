#!/bin/bash
# Atribucion de las perdidas de SC con afin: cada correccion apagada por separado. $1 ruta $2 variante
cd ~/ibex2024/examples/dfb
n=$(basename $1 .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/culpa/$n.$2.csv
[ -s $o ] && exit 0
ulimit -v 6000000
case $2 in
  todas)     e="X=1";;
  sin_rayo)  e="DFB_SX_RAYO=k";;
  sin_piso)  e="DFB_SX_PISO=no";;
  sin_filas) e="DFBH_FILAS_PARALELAS=no";;
  sin_tol)   e="DFB_SX_FEASTOL=1e-9 DFB_SX_DUALTOL=1e-9 DFB_SX_TOLANCHO=1e-3";;
  sin_orden) e="DFBH_ORDEN=costo";;
  ninguna)   e="DFB_SX_RAYO=k DFB_SX_PISO=no DFBH_FILAS_PARALELAS=no DFB_SX_FEASTOL=1e-9 DFB_SX_DUALTOL=1e-9 DFB_SX_TOLANCHO=1e-3 DFBH_ORDEN=costo";;
esac
r=$(env $e timeout 150 ./ibexsolve_dfb_y $1 --lr=dfbhull --lin=art -t 60 --seed=1 2>/dev/null | tail -1)
echo "$n,$2,$r" > $o
