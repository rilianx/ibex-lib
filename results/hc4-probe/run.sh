#!/bin/sh
# Bisector with the drastic HC4 probe (results/sb-trajectory/README.md) vs
# lsmear-lffix, on the instances of results/ipopt-compo-fixed.csv.
# compo + ipoptxn, TLIM CPU seconds per run (default 600). One line of JSON
# per (instance, arm) in OUT; resumable (skips the pairs already in OUT).
#
#   git pull origin ml-branching && make -C build -j16
#   BIN=build/bin/ibexopt-ml JOBS=16 results/hc4-probe/run.sh
#   python3 results/hc4-probe/analyze.py
#
# Arms (ARMS="lffix hc4 hc4r025 hc4need" by default):
#   lffix     lsmear-lffix
#   hc4       + HC4 probe, 1 dimension x 4 parts, conservative r=0.5
#   hc4r025   + HC4 probe, 1 x 4, conservative r=0.25
#   hc4need   + HC4 probe, 1 x 4, deviate only if it empties more pieces
#   guard10   lsmear-guard:10 (reference)
cd "$(dirname "$0")/../.."
BIN=${BIN:-build/bin/ibexopt-ml}
OUT=${OUT:-results/hc4-probe/runs.jsonl}
TLIM=${TLIM:-600}; export TLIM
$BIN --help 2>/dev/null | grep -q probe-parts || { echo "$BIN no conoce --probe-parts: recompila"; exit 1; }
touch $OUT
for arm in ${ARMS:-lffix hc4 hc4r025 hc4need}; do while read i; do echo "$arm $i"; done < results/hc4-probe/instances.txt; done |
xargs -P ${JOBS:-8} -n2 sh -c '
  arm=$0; i=$1; grep -q "\"instance\": \"$i\", \"arm\": \"$arm\"" '"$OUT"' && exit 0
  case $arm in
    lffix)   a="--bisector lsmear-lffix" ;;
    hc4)     a="--oracle --oracle-score hc4 --probe-dims 1 --probe-parts 4 --sb-ratio 0.5 --bisector lsmear-lffix" ;;
    hc4r025) a="--oracle --oracle-score hc4 --probe-dims 1 --probe-parts 4 --sb-ratio 0.25 --bisector lsmear-lffix" ;;
    hc4need) a="--oracle --oracle-score hc4 --probe-dims 1 --probe-parts 4 --sb-ratio 0.5 --sb-need-pruned --bisector lsmear-lffix" ;;
    guard10) a="--bisector lsmear-guard --guard-horizon 10" ;;
  esac
  r=$(timeout $(($TLIM+200)) '"$BIN"' benchs/optim/all/$i --solve $a --relax both --loup ipoptxn \
       --random-seed 1 --timeout $TLIM 2>/dev/null | tail -1)
  [ -z "$r" ] && r="{\"status\": \"killed\"}"
  echo "{\"instance\": \"$i\", \"arm\": \"$arm\", \"result\": $r}" >> '"$OUT"'
'
