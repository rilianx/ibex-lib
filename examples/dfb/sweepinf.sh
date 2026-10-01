#!/bin/bash
ff="$1"; v="$2"; s="$3"; out="/tmp/claude-1001/-home-iaraya/92eab584-7c36-4feb-85c5-d59a58d01d42/scratchpad/cic"
n=$(basename "$ff" .bch); o="$out/$n.$v.$s.csv"
[ -s "$o" ] && exit 0
case $v in
  fix)  e="DFBH_LR=art" ;;
  inf)  e="DFBH_NOFIX=1 DFBH_LR=art DFBH_CIRC=100000 DFBH_HC4VAR=0.3 DFBH_RELINCIC=100000" ;;
esac
r=$(env $e timeout 200 ./ibexopt "$ff" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1)
echo "$n,$v,$s,$r" > "$o"
