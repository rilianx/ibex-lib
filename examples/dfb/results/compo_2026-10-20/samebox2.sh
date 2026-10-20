#!/bin/bash
# tab:samebox con y sin el tope igualado: ambos motores sobre las mismas cajas, 60 s, afin. $1 ruta $2 tope100|topeigual
cd ~/ibex2024/examples/dfb
M=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad
n=$(basename $1 .bch); o=$M/samebox2/$n.$2.txt
[ -s $o ] && exit 0
ulimit -v 6000000
e="DFBH_ABCAJA=1 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=1e-6 DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0"
[ $2 = topeigual ] && e="$e DFBH_CAJA_TOPE=igual"
env $e timeout 150 ./ibexsolve_dfb_x $1 --lr=dfbhull --lin=art -t 60 --seed=1 2>&1 | grep -E 'abcaja|\[dfb\] estados|^[0-9]+ [0-9.]+ ' > $o
