#!/bin/sh
# How often is the Lagrangian (multipliers of the last usable Ipopt run)
# certifiably convex on a node's box (interval Hessian + Gershgorin), and how
# often would its tangent-plane bound close the box? --cvx-stats, measurement
# only (the search is not changed). All 298 instances, TLIM s each.
#   TLIM=60 BIN=build/bin/ibexopt-ml JOBS=3 results/convexity/run.sh
cd "$(dirname "$0")/../.."
BIN=${BIN:-build/bin/ibexopt-ml}
OUT=${OUT:-results/convexity/runs.jsonl}
TLIM=${TLIM:-60}; export TLIM
touch $OUT
ls benchs/optim/all | grep '\.bch$' |
xargs -P ${JOBS:-8} -n1 sh -c '
  i=$0; grep -q "\"instance\": \"$i\"" '"$OUT"' && exit 0
  r=$(timeout $(($TLIM+200)) '"$BIN"' benchs/optim/all/$i --solve --bisector lsmear-lffix --relax both --loup ipoptxn \
       --ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --cvx-stats --random-seed 1 --timeout $TLIM 2>/dev/null | tail -1)
  [ -z "$r" ] && r="{\"status\": \"killed\"}"
  echo "{\"instance\": \"$i\", \"result\": $r}" >> '"$OUT"'
'
