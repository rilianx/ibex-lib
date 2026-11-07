#!/bin/bash
# Hibrido real por nodo con SC base, compo, umbrales 15 20 30 (2026-11-06). $1 ruta $2 semilla
cd ~/ibex2024/examples/dfb
S=/tmp/claude-1001/-home-iaraya/2bd51c7a-9daf-4c06-b6ef-b89011baa9b5/scratchpad/hyb
n=$(basename "$1" .bch); ulimit -v 6000000
for T in 30 20 15; do
  o=$S/sol.$n.H$T.$2.csv; [ -s "$o" ] && continue
  r=$(env DFB_SX_BFRT=0 DFB_SX_DEVEX=0 DFBH_LPCAP=0 DFBH_HIBRIDO=$T timeout 150 ./ibexsolve_dfb_esc "$1" --lr=dfbhull --lin=compo -t 60 --seed=$2 2>/dev/null | tail -1)
  [ -z "$r" ] && r="NA"
  echo "$n,H$T,$2,$r" > $o
done
