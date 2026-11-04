#!/bin/sh
# How much can node selection gain, per instance? Two runs of lsmear-lffix:
#   base  as usual; reports when the upper bound reached its final value
#         (loup_near_nodes/time: within 1e-3 relative; loup_final_*: exactly)
#   opt   the same with --initial-loup set to base's final loup: the tree when
#         the optimum is known from the start (no node order can do better,
#         up to ACID: acid_nbcidvar flags the runs where ACID tuned itself
#         differently, see results/fixed-loup/README.md)
#   INSTANCES=list TLIM=120 OUT=... BIN=... JOBS=... results/node-selection/headroom.sh
cd "$(dirname "$0")/../.."
BIN=${BIN:-build/bin/ibexopt-ml}
OUT=${OUT:-results/node-selection/headroom.jsonl}
TLIM=${TLIM:-600}; export TLIM BIN OUT
touch $OUT
xargs -P ${JOBS:-8} -n1 sh -c '
  i=$0; grep -q "\"instance\": \"$i\"" "$OUT" && exit 0
  common="--bisector lsmear-lffix --relax both --loup ipoptxn --random-seed 1 --timeout $TLIM"
  b=$(timeout $((TLIM+200)) $BIN benchs/optim/all/$i --solve $common 2>/dev/null | tail -1)
  [ -z "$b" ] && b="{\"status\": \"killed\"}"
  lo=$(echo "$b" | python3 -c "import json,sys;r=json.load(sys.stdin);print(r[\"loup\"] if r.get(\"status\")==\"complete\" else \"\")" 2>/dev/null)
  o="null"
  if [ -n "$lo" ]; then
    o=$(timeout $((TLIM+200)) $BIN benchs/optim/all/$i --solve $common --initial-loup $lo 2>/dev/null | tail -1)
    [ -z "$o" ] && o="{\"status\": \"killed\"}"
  fi
  echo "{\"instance\": \"$i\", \"base\": $b, \"opt\": $o}" >> "$OUT"
' < ${INSTANCES:-results/pseudocost/instances.txt}
