#!/bin/sh
# Reliability branching with pseudocosts relative to the bisector
# (--oracle-score pc, see MLNodeServer::set_oracle_pc) vs lsmear-lffix, on the
# instances of results/ipopt-compo-fixed.csv. compo + ipoptxn, TLIM CPU
# seconds per run (default 600). One line of JSON per (instance, arm) in OUT;
# resumable (skips the pairs already in OUT).
#
#   git pull fork ml-branching && make -C build -j16
#   BIN=build/bin/ibexopt-ml JOBS=16 results/pseudocost/run.sh
#   python3 results/pseudocost/analyze.py
#   git add results/pseudocost/runs.jsonl && git commit -m "Benchmark pseudocostos" && git push fork ml-branching
#
# Arms (ARMS="lffix pc pcr05" by default):
#   lffix   lsmear-lffix
#   pc      + pseudocosts, reliability 4, deviate when the mean says better
#   pcr05   + same, deviate only when the mean says at most half the volume
cd "$(dirname "$0")/../.."
BIN=${BIN:-build/bin/ibexopt-ml}
OUT=${OUT:-results/pseudocost/runs.jsonl}
TLIM=${TLIM:-600}; export TLIM
$BIN --help 2>/dev/null | grep -q pc-rel || { echo "$BIN no conoce --pc-rel: recompila"; exit 1; }
touch $OUT
for arm in ${ARMS:-lffix pc pcr05}; do while read i; do echo "$arm $i"; done < results/pseudocost/instances.txt; done |
xargs -P ${JOBS:-8} -n2 sh -c '
  arm=$0; i=$1; grep -q "\"instance\": \"$i\", \"arm\": \"$arm\"" '"$OUT"' && exit 0
  case $arm in
    lffix) a="--bisector lsmear-lffix" ;;
    pc)    a="--oracle --oracle-score pc --pc-rel 4 --bisector lsmear-lffix" ;;
    pcr05) a="--oracle --oracle-score pc --pc-rel 4 --sb-ratio 0.5 --bisector lsmear-lffix" ;;
  esac
  r=$(timeout $(($TLIM+200)) '"$BIN"' benchs/optim/all/$i --solve $a --relax both --loup ipoptxn \
       --random-seed 1 --timeout $TLIM 2>/dev/null | tail -1)
  [ -z "$r" ] && r="{\"status\": \"killed\"}"
  echo "{\"instance\": \"$i\", \"arm\": \"$arm\", \"result\": $r}" >> '"$OUT"'
'
