#!/bin/sh
# How much can node selection gain, per instance? Two runs of lsmear-lffix:
#   base  as usual; reports when the upper bound reached its final value
#         (loup_near_nodes/time: within 1e-3 relative; loup_final_*: exactly)
#   opt   the same with --late-loup set to base's final loup: the optimum is
#         given right after the root (no node order can do better). Not
#         --initial-loup: a bound from the start changes the root's own
#         contraction, and on some instances the whole tree (dnieper: 196 to
#         >7000 nodes depending on the bound's value).
# Both with ACID non-adaptive on all variables (--acid-fixed), so that the
# two histories do not tune it differently.
# EXTRA: more flags for both runs (the Ipopt base since results/ipopt-diam:
# EXTRA="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1").
#   INSTANCES=list TLIM=120 OUT=... BIN=... JOBS=... results/node-selection/headroom.sh
cd "$(dirname "$0")/../.."
BIN=${BIN:-build/bin/ibexopt-ml}
OUT=${OUT:-results/node-selection/headroom.jsonl}
TLIM=${TLIM:-600}; EXTRA=${EXTRA:-}; export TLIM BIN OUT EXTRA
touch $OUT
xargs -P ${JOBS:-8} -n1 sh -c '
  i=$0; grep -q "\"instance\": \"$i\"" "$OUT" && exit 0
  common="--bisector lsmear-lffix --relax both --loup ipoptxn --random-seed 1 --timeout $TLIM --acid-fixed 1000000 $EXTRA"
  b=$(timeout $((TLIM+200)) $BIN benchs/optim/all/$i --solve $common 2>/dev/null | tail -1)
  [ -z "$b" ] && b="{\"status\": \"killed\"}"
  lo=$(echo "$b" | python3 -c "import json,sys;r=json.load(sys.stdin);print(r[\"loup\"] if r.get(\"status\")==\"complete\" else \"\")" 2>/dev/null)
  o="null"
  if [ -n "$lo" ]; then
    o=$(timeout $((TLIM+200)) $BIN benchs/optim/all/$i --solve $common --late-loup $lo 2>/dev/null | tail -1)
    [ -z "$o" ] && o="{\"status\": \"killed\"}"
  fi
  echo "{\"instance\": \"$i\", \"base\": $b, \"opt\": $o}" >> "$OUT"
' < ${INSTANCES:-results/pseudocost/instances.txt}
