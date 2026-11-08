#!/bin/bash
# Trabajo emulado con SC base (2026-11-07): config b1 (SC base vs SoPlex ajustado) + cascada acumulada DFB_EMU_*. $1 ruta $2 sol|opt
H=$(cd "$(dirname "$0")" && pwd)   # this folder; work/ holds datos_trabajo.tgz unpacked
cd "$H/../../examples/dfb"
M=${M:-$H/work/M}
S=${S:-$H/work/S}
n=$(basename $1 .bch); ulimit -v 6000000
c=$S/bancoopt/$2.$n.txt; [ -s $c ] || c=$M/banco/$2.$n.txt; [ -s $c ] || exit 0
if [ $2 = opt ]; then tol=1e-9; else tol=1e-6; fi
dir="DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol DFBH_CAJA_TOPE=igual"
for p in e0 e1 e2 e3 e4; do
  o=$S/emu/$2.$n.$p.txt; [ -s $o ] && continue
  case $p in
    e0) x="";;
    e1) x="DFB_EMU_FACT=1";;
    e2) x="DFB_EMU_FACT=1 DFB_EMU_INIT=1";;
    e3) x="DFB_EMU_FACT=1 DFB_EMU_INIT=1 DFB_EMU_VERIF=1";;
    e4) x="DFB_EMU_FACT=1 DFB_EMU_INIT=1 DFB_EMU_VERIF=1 DFB_EMU_COPIA=1";;
  esac
  env $dir $x DFBH_ABCAJA=1 DFBH_CORTE_CERT=1 timeout 3600 ./ibexreplay_esc $1 $2 $c 2>&1 | grep -E 'abcaja|replay|costo' > $o
done
