#!/bin/bash
# Banco de pasadas, optimizacion una por familia (2026-11-06). $1 ruta. Captura todas las llamadas de una busqueda de 60 s con SC base,
# guarda 1000 repartidas uniformemente y corre los peldanos b0 b0f b1 b2 (como bescf.sh).
H=$(cd "$(dirname "$0")" && pwd)   # this folder; work/ holds datos_trabajo.tgz unpacked
cd "$H/../../examples/dfb"
S=${S:-$H/work/S}
n=$(basename $1 .bch); c=$S/bancoopt/opt.$n.txt; ulimit -v 6000000
if [ ! -s $c ]; then
  raw=$S/bancoopt/raw.$n
  DFB_SX_BFRT=0 DFB_SX_DEVEX=0 DFBH_LPCAP=0 DFBH_CAPTURA=$raw DFBH_CAPTURA_N=100000000 timeout 100 ./ibexopt_esc $1 --filtering=acidhc4 --lr=dfbhull --seed=1 -t60 > /dev/null 2>&1
  L=$(wc -l < $raw); awk -v L=$L 'BEGIN{k=-1} {j=int((NR-1)*1000/L); if(j!=k){print; k=j}}' $raw > $c; rm -f $raw
fi
tol=1e-9
ibex="DFBH_CAJA_TOL=1e-9 DFBH_ORDEN_CAJA=nat"
dir="DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol DFBH_CAJA_TOPE=igual"
for p in b0 b0f b1 b2; do
  o=$S/bescf/opt.$n.$p.txt; [ -s $o ] && continue
  case $p in
    b0)  e="$ibex DFBH_FUSIONA=0 DFBH_FILAS_PARALELAS=no";;
    b0f) e="$ibex";;
    b1)  e="$dir";;
    b2)  e="$dir DFBH_CAJA_SPX=1 DFB_SX_DEVEX=1 DFB_SX_BFRT=1 DFBH_LPCAP=60";;
  esac
  env $e DFBH_ABCAJA=1 DFBH_CORTE_CERT=1 timeout 3600 ./ibexreplay_esc $1 opt $c 2>&1 | grep -E 'abcaja|replay' > $o
done
