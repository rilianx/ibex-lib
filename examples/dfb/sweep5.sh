#!/bin/bash
ff="$1"; var="$2"; s="$3"; out="/tmp/claude-1001/-home-iaraya/92eab584-7c36-4feb-85c5-d59a58d01d42/scratchpad/s5"
n=$(basename "$ff" .bch)
o="$out/$n.$var.$s.csv"
[ -s "$o" ] && exit 0
case $var in
  xn)     r=$(timeout 200 ./ibexopt "$ff" --filtering=acidhc4 --seed=$s -t60 2>/dev/null | tail -1) ;;
  base)   r=$(timeout 200 ./ibexopt "$ff" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1) ;;
  pivfix) r=$(DFB_SX_PIVFIX=1 timeout 200 ./ibexopt "$ff" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1) ;;
  reub)   r=$(DFB_SX_REUBICAR=1 timeout 200 ./ibexopt "$ff" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1) ;;
  both)   r=$(DFB_SX_PIVFIX=1 DFB_SX_REUBICAR=1 timeout 200 ./ibexopt "$ff" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1) ;;
esac
echo "$n,$var,$s,$r" > "$o"
