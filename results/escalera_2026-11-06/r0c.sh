#!/bin/bash
# R0 con el orden de Ibex (mas cercana) en vez de natural (2026-11-07). $1 ruta $2 sol|opt|cmp
H=$(cd "$(dirname "$0")" && pwd)   # this folder; work/ holds datos_trabajo.tgz unpacked
cd "$H/../../examples/dfb"
M=${M:-$H/work/M}
S=${S:-$H/work/S}
n=$(basename $1 .bch); ulimit -v 6000000
if [ $2 = cmp ]; then c=$S/bancocmp/sol.$n.txt; pr=sol; x="DFBH_REPLAY_LIN=compo"; bin=./ibexreplay_cmp
else c=$S/bancoopt/$2.$n.txt; [ -s $c ] || c=$M/banco/$2.$n.txt; pr=$2; x=""; bin=./ibexreplay_esc; fi
[ -s $c ] || exit 0
o=$S/r0c/$2.$n.b0fc.txt; [ -s $o ] && exit 0
env $x DFBH_CAJA_TOL=1e-9 DFBH_ORDEN_CAJA=cerca DFBH_ABCAJA=1 DFBH_CORTE_CERT=1 timeout 3600 $bin $1 $pr $c 2>&1 | grep -E 'abcaja|replay' > $o
