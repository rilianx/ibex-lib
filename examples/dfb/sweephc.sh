#!/bin/bash
ff="$1"; k="$2"; s="$3"; out="/tmp/claude-1001/-home-iaraya/92eab584-7c36-4feb-85c5-d59a58d01d42/scratchpad/hc4m"
n=$(basename "$ff" .bch); o="$out/$n.$k.$s.csv"
[ -s "$o" ] && exit 0
r=$(DFBH_LR=art DFBH_HC4MED=$k timeout 200 ./ibexopt "$ff" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1)
echo "$n,$k,$s,$r" > "$o"
