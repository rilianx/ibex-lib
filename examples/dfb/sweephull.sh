#!/bin/bash
ff="$1"; v="$2"; s="$3"; out="/tmp/claude-1001/-home-iaraya/92eab584-7c36-4feb-85c5-d59a58d01d42/scratchpad/hull"
n=$(basename "$ff" .bch); o="$out/$n.$v.$s.csv"
[ -s "$o" ] && exit 0
case $v in
  base) e="DFBH_LR=art" ;;
  h6)   e="DFBH_LR=art DFBH_HULL=1 DFBH_HULLTAU=1e-6" ;;
  h5)   e="DFBH_LR=art DFBH_HULL=1 DFBH_HULLTAU=1e-5" ;;
  h4)   e="DFBH_LR=art DFBH_HULL=1 DFBH_HULLTAU=1e-4" ;;
  h3)   e="DFBH_LR=art DFBH_HULL=1 DFBH_HULLTAU=1e-3" ;;
  h2)   e="DFBH_LR=art DFBH_HULL=1 DFBH_HULLTAU=1e-2" ;;
esac
r=$(env $e timeout 200 ./ibexopt "$ff" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1)
echo "$n,$v,$s,$r" > "$o"
