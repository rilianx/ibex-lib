#!/bin/bash
# Rehecha con el limite de iteraciones de Ibex (100 por LP) en TODOS los motores (2026-11-11).
# Escribe en work/S100. Los peldanos que solo usan SoPlex como lo usa Ibex (R0, R0f del solver)
# no cambian y se toman de work/S (escf_r0c, escfa_r0c).
#   l100.sh banco opt|sol|cmp <ruta>     banco de pasadas: b0 b0f (orden mas cercana, como r0c), b1, b2
#   l100.sh solver escf|escfa <ruta> <prog> <semilla>   R1 R2 R3 R6 (+R7 en sistemas compo)
#   l100.sh hyb <ruta> <semilla> | emu <ruta> <prog> | peso <ruta> <prog>
H=$(cd "$(dirname "$0")" && pwd)
cd "$H/../../examples/dfb"
M=${M:-$H/work/M}; S=$H/work/S100; ulimit -v 6000000
L="DFBH_LPCAP=100"                 # limite de SC
TOPE="DFBH_CAJA_TOPE=100"          # limite de SoPlex (R1, R2, hibrido)
SC="DFB_SX_BFRT=0 DFB_SX_DEVEX=0 $L"
SCMAS="DFB_SX_BFRT=1 DFB_SX_DEVEX=1 $L"   # SC + devex y paso largo, con el mismo limite de 100
case $1 in
banco)
  k=$2; n=$(basename $3 .bch)
  case $k in
    cmp) c=$S/bancocmp/sol.$n.txt; pr=sol; x="DFBH_REPLAY_LIN=compo"; bin=./ibexreplay_cmp; tol=1e-6; out=$S/bcmp/sol.$n;;
    opt) c=$S/bancoopt/opt.$n.txt; [ -s $c ] || c=$M/banco/opt.$n.txt; pr=opt; x=""; bin=./ibexreplay_esc; tol=1e-9; out=$S/bescf/opt.$n;;
    sol) c=$M/banco/sol.$n.txt; pr=sol; x=""; bin=./ibexreplay_esc; tol=1e-6; out=$S/bescf/sol.$n;;
  esac
  [ -s $c ] || exit 0
  ibex="DFBH_CAJA_TOL=1e-9 DFBH_ORDEN_CAJA=cerca"
  dir="DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol $TOPE"
  for p in b0 b0f b1 b2; do
    o=$out.$p.txt; [ -s $o ] && continue
    case $p in
      b0)  e="$SC $ibex DFBH_FUSIONA=0 DFBH_FILAS_PARALELAS=no";;
      b0f) e="$SC $ibex";;
      b1)  e="$SC $dir";;
      b2)  e="$dir DFBH_CAJA_SPX=1 $SCMAS";;
    esac
    env $x $e DFBH_ABCAJA=1 DFBH_CORTE_CERT=1 timeout 3600 $bin $3 $pr $c 2>&1 | grep -E 'abcaja|replay|costo|Assert|terminate' > $o
  done;;
solver)
  d=$2; n=$(basename "$3" .bch); pr=$4; s=$5
  if [ $d = escf ] && [ $pr = sol ]; then tc=1e-6; lin=compo; pel="R1 R2 R3 R6 R7"
  elif [ $pr = sol ]; then tc=1e-6; lin=art; pel="R1 R2 R3 R6"
  else tc=1e-9; lin=art; pel="R1 R2 R3 R6"; fi
  dir="DFBH_MOTOR=caja DFBH_CAJA_TOL=$tc DFBH_ORDEN_CAJA=cerca DFBH_CORTE_CERT=1 DFBH_CAJA_SIMP=0 DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 $TOPE"
  for p in $pel; do
    o=$S/$d/$pr.$n.$p.$s.csv; [ -s "$o" ] && continue
    case $p in
      R1) e="$dir";;
      R2) e="$dir DFBH_CAJA_SPX=1";;
      R3) e="$SC";;
      R6) e="$SCMAS";;
      R7) e="$SC $TOPE DFBH_HIBRIDO=30";;
    esac
    if [ $pr = sol ]; then r=$(env $e timeout 150 ./ibexsolve_dfb_esc "$3" --lr=dfbhull --lin=$lin -t 60 --seed=$s 2>/dev/null | tail -1)
    else r=$(env $e timeout 200 ./ibexopt_esc "$3" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1); fi
    [ -z "$r" ] && r="NA"; echo "$n,$p,$s,$r" > $o
  done;;
hyb)
  n=$(basename "$2" .bch)
  for T in 30 20 15; do
    o=$S/hyb/sol.$n.H$T.$3.csv; [ -s "$o" ] && continue
    r=$(env $SC $TOPE DFBH_HIBRIDO=$T timeout 150 ./ibexsolve_dfb_esc "$2" --lr=dfbhull --lin=compo -t 60 --seed=$3 2>/dev/null | tail -1)
    [ -z "$r" ] && r="NA"; echo "$n,H$T,$3,$r" > $o
  done;;
emu)
  n=$(basename $2 .bch); pr=$3
  c=$S/bancoopt/$pr.$n.txt; [ -s $c ] || c=$M/banco/$pr.$n.txt; [ -s $c ] || exit 0
  if [ $pr = opt ]; then tol=1e-9; else tol=1e-6; fi
  dir="DFBH_CAJA_DIRECTA=1 DFBH_CAJA_SINLIBRES=1 DFBH_CAJA_SIMP=0 DFBH_ORDEN_CAJA=cerca DFBH_CAJA_TOL=$tol $TOPE"
  for p in e0 e1 e2 e3 e4; do
    o=$S/emu/$pr.$n.$p.txt; [ -s $o ] && continue
    case $p in
      e0) x="";; e1) x="DFB_EMU_FACT=1";; e2) x="DFB_EMU_FACT=1 DFB_EMU_INIT=1";;
      e3) x="DFB_EMU_FACT=1 DFB_EMU_INIT=1 DFB_EMU_VERIF=1";; e4) x="DFB_EMU_FACT=1 DFB_EMU_INIT=1 DFB_EMU_VERIF=1 DFB_EMU_COPIA=1";;
    esac
    env $SC $dir $x DFBH_ABCAJA=1 DFBH_CORTE_CERT=1 timeout 3600 ./ibexreplay_esc $2 $pr $c 2>&1 | grep -E 'abcaja|replay|costo' > $o
  done;;
peso)
  n=$(basename $2 .bch); pr=$3; o=$S/peso/$pr.$n.txt; [ -s $o ] && exit 0; c=$S/peso/$pr.$n.cnt
  if [ $pr = opt ]; then
    r=$(env $SC DFBH_CAPTURA=/dev/fd/3 DFBH_CAPTURA_N=1000000000 timeout 100 ./ibexopt_esc $2 --filtering=acidhc4 --lr=dfbhull --seed=1 -t60 2>/dev/null 3> >(wc -l > $c) | tail -1)
  else
    r=$(env $SC DFBH_CAPTURA=/dev/fd/3 DFBH_CAPTURA_N=1000000000 timeout 100 ./ibexsolve_dfb_esc $2 --lr=dfbhull --lin=art -t 60 --seed=1 2>/dev/null 3> >(wc -l > $c) | tail -1)
  fi
  sleep 1; echo "$n $(cat $c) $r" > $o; rm -f $c;;
esac
