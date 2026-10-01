#!/bin/bash
ff="$1"; v="$2"; s="$3"; out="/tmp/claude-1001/-home-iaraya/92eab584-7c36-4feb-85c5-d59a58d01d42/scratchpad/cic"
n=$(basename "$ff" .bch); o="$out/$n.$v.$s.csv"
[ -s "$o" ] && exit 0
if [ "$v" = fix ]; then e="DFBH_LR=art"
else
  u=${v%%_*}; kk=${v##*_}
  e="DFBH_NOFIX=1 DFBH_LR=art DFBH_CIRC=8 DFBH_HC4VAR=0.${u#c} DFBH_RELINCIC=${kk#k}"
fi
r=$(env $e timeout 200 ./ibexopt "$ff" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1)
echo "$n,$v,$s,$r" > "$o"
