#!/bin/bash
# Peso de la pasada (2026-11-06): llamadas al contractor y tiempo de una busqueda de 60 s con SC base, relajacion afin. $1 ruta $2 sol|opt
H=$(cd "$(dirname "$0")" && pwd)   # this folder; work/ holds datos_trabajo.tgz unpacked
cd "$H/../../examples/dfb"
S=${S:-$H/work/S}
n=$(basename $1 .bch); o=$S/peso/$2.$n.txt; [ -s $o ] && exit 0; ulimit -v 6000000
c=$S/peso/$2.$n.cnt
if [ $2 = opt ]; then
  r=$(DFB_SX_BFRT=0 DFB_SX_DEVEX=0 DFBH_LPCAP=0 DFBH_CAPTURA=/dev/fd/3 DFBH_CAPTURA_N=1000000000 timeout 100 ./ibexopt_esc $1 --filtering=acidhc4 --lr=dfbhull --seed=1 -t60 2>/dev/null 3> >(wc -l > $c) | tail -1)
else
  r=$(DFB_SX_BFRT=0 DFB_SX_DEVEX=0 DFBH_LPCAP=0 DFBH_CAPTURA=/dev/fd/3 DFBH_CAPTURA_N=1000000000 timeout 100 ./ibexsolve_dfb_esc $1 --lr=dfbhull --lin=art -t 60 --seed=1 2>/dev/null 3> >(wc -l > $c) | tail -1)
fi
sleep 1; echo "$n $(cat $c) $r" > $o; rm -f $c
