#!/bin/bash
ff="$1"; v="$2"; s="$3"; out="/tmp/claude-1001/-home-iaraya/92eab584-7c36-4feb-85c5-d59a58d01d42/scratchpad/nf"
n=$(basename "$ff" .bch); o="$out/$n.$v.$s.csv"
[ -s "$o" ] && exit 0
e="DFBH_NOFIX=1 DFBH_LR=art"
case $v in
  nf)        : ;;
  nf_hc4)    e="$e DFBH_HC4MED=4" ;;
  nf_hc4v)   e="$e DFBH_HC4VAR=0.3" ;;
  nf_compo)  e="DFBH_NOFIX=1 DFBH_LR=compo" ;;
  nf_cpohc4) e="DFBH_NOFIX=1 DFBH_LR=compo DFBH_HC4MED=4" ;;
  nf_nat)    e="$e DFBH_ORDEN=nat" ;;
  nf_cerca)  e="$e DFBH_ORDEN=cerca" ;;
  nf_nosk)   e="$e DFBH_NOSKIP=1" ;;
esac
r=$(env $e timeout 200 ./ibexopt "$ff" --filtering=acidhc4 --lr=dfbhull --seed=$s -t60 2>/dev/null | tail -1)
echo "$n,$v,$s,$r" > "$o"
