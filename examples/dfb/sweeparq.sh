#!/bin/bash
ff="$1"; lin="$2"; arq="$3"; s="$4"; out="/tmp/claude-1001/-home-iaraya/92eab584-7c36-4feb-85c5-d59a58d01d42/scratchpad/arq"
n=$(basename "$ff" .bch)
o="$out/$n.$lin.$arq.$s.csv"
[ -s "$o" ] && exit 0
[ "$arq" = comp ] && pre="DFBH_COMPARTIDA=1" || pre=""
r=$(env $pre DFBH_LR=$lin timeout 200 ./ibexopt "$ff" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1)
echo "$n,$lin,$arq,$s,$r" > "$o"
