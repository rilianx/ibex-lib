#!/bin/bash
# Caja directa sin filas libres y SIN PRESOLVE (DFBH_CAJA_SIMP=0), 60 s, binarios _dir2 (2026-10-20).
# $1 ruta $2 brazo $3 semilla $4 sol|opt
cd ~/ibex2024/examples/dfb
n=$(basename "$1" .bch); o=/tmp/claude-1001/-home-iaraya/5f64b562-0cdb-4c11-89ea-a54f0d24b4aa/scratchpad/d60p/$4.$n.$2.$3.csv
[ -s "$o" ] && exit 0
ulimit -v 6000000
if [ $4 = sol ]; then tc=1e-6; else tc=1e-9; fi
case $2 in
  cajad)      e="DFBH_ORDEN=lejos DFBH_MOTOR=caja DFBH_CAJA_TOL=$tc DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1";;
  cajacercap) e="DFBH_ORDEN=cerca DFBH_MOTOR=caja DFBH_CAJA_TOL=$tc DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0";;
  scnew)      e="X=1";;
esac
if [ $4 = sol ]; then
  r=$(env $e timeout 150 ./ibexsolve_dfb_dir2 "$1" --lr=dfbhull --lin=art -t 60 --seed=$3 2>/dev/null | tail -1)
else
  r=$(env $e timeout 200 ./ibexopt_dir2 "$1" --filtering=acidhc4 --lr=dfbhull --seed=$3 -t60 2>/dev/null | tail -1)
fi
[ -z "$r" ] && r="NA"
echo "$n,$2,$3,$r" > $o
