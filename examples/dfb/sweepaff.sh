#!/bin/bash
ff="$1"; ctc="$2"; lin="$3"; out="/tmp/claude-1001/-home-iaraya/92eab584-7c36-4feb-85c5-d59a58d01d42/scratchpad/af"
n=$(basename "$ff" .bch)
o="$out/$n.$ctc.$lin.csv"
[ -s "$o" ] && exit 0
if [ "$ctc" = dfbhull ]; then extra=--lr=dfbhull; else extra=""; fi
r=$(DFBH_LR=$lin timeout 200 ./ibexopt "$ff" --filtering=acidhc4 $extra --seed=1 -t60 2>/dev/null | tail -1)
echo "$n,$ctc,$lin,$r" > "$o"
