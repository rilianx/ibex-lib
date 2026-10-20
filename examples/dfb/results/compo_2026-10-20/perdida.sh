#!/bin/bash
# Descomposicion de la perdida con compo y orden costo_o_cerca. $1 ruta $2 ab|sombra
cd ~/ibex2024/examples/dfb
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
n=$(basename $1 .bch); o=$M/perdida/$n.$2.txt
[ -s $o ] && exit 0
ulimit -v 6000000
c="DFBH_ORDEN=costo_o_cerca DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=1e-6 DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_CAJA_TOPE=igual"
case $2 in ab) e="DFBH_ABCAJA=1 $c";; sombra) e="DFBH_SOMBRA=1 $c";; abcol) e="DFBH_CAJA_REP=col DFBH_ABCAJA=1 $c";; sombracol) e="DFBH_CAJA_REP=col DFBH_SOMBRA=1 $c";; ab6) e="DFB_SX_FEASTOL=1e-6 DFBH_ABCAJA=1 $c";; sombra6) e="DFB_SX_FEASTOL=1e-6 DFBH_SOMBRA=1 $c";; esac
env $e timeout 150 ./ibexsolve_dfb_x $1 --lr=dfbhull --lin=compo -t 60 --seed=1 2>&1 | grep -E 'abcaja|sombra|\[dfb\] piv|\[dfb\] nodos=|\[dfb\] filas de|^[0-9]+ [0-9.]+ ' > $o
