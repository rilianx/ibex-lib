#!/bin/bash
cd ~/ibex2024/examples/dfb
n=$(basename $1 .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/ws/$n.$2.txt
[ $2 = primal ] && x="DFB_SX_PIVOTEO=primal DFB_SX_REVPRIMAL=1 DFB_SX_DEVEXP=1 DFB_SX_DEVEXPT=1" || x="X=1"
env DFB_SX_FEASTOL=1e-6 DFB_SX_DUALTOL=1e-6 DFB_SX_TOLANCHO=0 DFB_SX_PISO=error DFBH_FILAS_PARALELAS=paralelas DFBH_ABCAJA=1 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=1e-6 DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_CAJA_TOPE=igual DFBH_BASECOMP=200 $x timeout 150 ./ibexsolve_dfb_x $1 --lr=dfbhull --lin=compo -t 40 --seed=1 2>&1 | grep -E 'abcaja\] (DFB|caja|nodos|costo)|basecomp|\[dfb\] (estados|piv/LP|filas de)|^[0-9]+ [0-9.]+ ' > $o
