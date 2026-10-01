#!/bin/bash
ff="$1"; ctc="$2"; r="$3"; s="$4"; out="/tmp/claude-1001/-home-iaraya/92eab584-7c36-4feb-85c5-d59a58d01d42/scratchpad/ratio"
n=$(basename "$ff" .bch)
o="$out/$n.$ctc.$r.$s.csv"
[ -s "$o" ] && exit 0
[ "$ctc" = dfbhull ] && extra=--lr=dfbhull || extra=""
v=$(DFBH_LR=art DFBH_RATIO=$r timeout 200 ./ibexopt "$ff" --filtering=acidhc4 $extra --seed=$s -t60 2>/dev/null | tail -1)
echo "$n,$ctc,$r,$s,$v" > "$o"
