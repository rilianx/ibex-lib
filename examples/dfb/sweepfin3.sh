#!/bin/bash
ff="$1"; v="$2"; s="$3"; out="/tmp/claude-1001/-home-iaraya/92eab584-7c36-4feb-85c5-d59a58d01d42/scratchpad/final"
n=$(basename "$ff" .bch); o="$out/$n.$v.$s.csv"
[ -s "$o" ] && exit 0
case $v in
  prod) e="DFBH_LR=art"; x="" ;;
  dfb)  e="DFBH_LR=art DFBH_HULL=1 DFBH_HULLTAU=1e-6"; x="--lr=dfbhull" ;;
esac
r=$(env $e timeout 200 ./ibexopt "$ff" --filtering=acidhc4 $x --seed=$s -t60 2>/dev/null | tail -1)
echo "$n,$v,$s,$r" > "$o"
