#!/bin/bash
ff="$1"; v="$2"; s="$3"; out="/tmp/claude-1001/-home-iaraya/92eab584-7c36-4feb-85c5-d59a58d01d42/scratchpad/hc4b"
n=$(basename "$ff" .bch); o="$out/$n.$v.$s.csv"
[ -s "$o" ] && exit 0
case $v in
  base) e="" ;;
  k4)   e="DFBH_HC4MED=4" ;;
  r01)  e="DFBH_HC4RED=0.01" ;;
  r05)  e="DFBH_HC4RED=0.05" ;;
esac
r=$(env $e DFBH_LR=art timeout 200 ./ibexopt "$ff" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1)
echo "$n,$v,$s,$r" > "$o"
