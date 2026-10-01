#!/bin/bash
ff="$1"; ctc="$2"; lin="$3"; s="$4"; out="/tmp/claude-1001/-home-iaraya/92eab584-7c36-4feb-85c5-d59a58d01d42/scratchpad/fin"
n=$(basename "$ff" .bch)
o="$out/$n.$ctc.$lin.$s.csv"
[ -s "$o" ] && exit 0
[ "$ctc" = dfbhull ] && extra=--lr=dfbhull || extra=""
r=$(DFBH_LR=$lin timeout 200 ./ibexopt "$ff" --filtering=acidhc4 $extra --seed=$s -t60 2>/dev/null | tail -1)
echo "$n,$ctc,$lin,$s,$r" > "$o"
