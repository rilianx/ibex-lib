#!/bin/bash
# Hibrido (R7, SC hasta 30 filas) en las clases que faltaban de la Tabla 2:
# optimizacion (escf opt) y sistemas con la relajacion afin (escfa sol), 2026-11-13.
#   r7.sh escf|escfa <ruta> <prog> <semilla>     (mismas lineas que l100_jobs.txt)
H=$(cd "$(dirname "$0")" && pwd)
cd "$H/../../examples/dfb"
S=$H/work/S100; ulimit -v 6000000
SC="DFB_SX_BFRT=0 DFB_SX_DEVEX=0 DFBH_LPCAP=100"; TOPE="DFBH_CAJA_TOPE=100"
d=$1; n=$(basename "$2" .bch); pr=$3; s=$4
o=$S/$d/$pr.$n.R7.$s.csv; [ -s "$o" ] && exit 0
e="$SC $TOPE DFBH_HIBRIDO=30"
if [ $pr = sol ]; then r=$(env $e timeout 150 ./ibexsolve_dfb_esc "$2" --lr=dfbhull --lin=art -t 60 --seed=$s 2>/dev/null | tail -1)
else r=$(env $e timeout 200 ./ibexopt_esc "$2" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1); fi
[ -z "$r" ] && r="NA"; echo "$n,R7,$s,$r" > $o
