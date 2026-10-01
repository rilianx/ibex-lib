#!/bin/bash
ff="$1"; v="$2"; s="$3"; out="/tmp/claude-1001/-home-iaraya/92eab584-7c36-4feb-85c5-d59a58d01d42/scratchpad/th"
n=$(basename "$ff" .bch); o="$out/$n.$v.$s.csv"
[ -s "$o" ] && exit 0
[ "$v" = nf ] && e="" || e="DFBH_HC4VAR=0.${v#v}"
r=$(env DFBH_NOFIX=1 DFBH_LR=art $e timeout 200 ./ibexopt "$ff" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1)
echo "$n,$v,$s,$r" > "$o"
