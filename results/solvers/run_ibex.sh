#!/bin/sh
# IbexOpt (afull configuration) on the 298 instances with missing variable bounds
# replaced by +-1e12, the same domains SCIP gets with SCIP_INFBOUND=1e12.
# One JSON line per instance in OUT; resumable.
#   BIN=build-fix/bin/ibexopt-ml JOBS=3 TLIM=310 sh results/solvers/run_ibex.sh
cd "$(dirname "$0")/../.."
BIN=${BIN:-build-fix/bin/ibexopt-ml}
OUT=${OUT:-results/solvers/ibex_inf12.jsonl}
TLIM=${TLIM:-310}; export TLIM
# ARM/FLAGS: another configuration on the same domains (default: afull)
ARM=${ARM:-afull12}
FLAGS=${FLAGS:-"--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.05 --ipopt-guided-patience=20 --ipopt-guided-nodes --cvx-close --cvx-backoff=64 --cvx-test 2"}
DIR=results/solvers/bch12
mkdir -p $DIR
for f in benchs/optim/all/*.bch; do
  sed -E 's/-[ ]*oo\b/-1e12/g; s/\+?oo\b/1e12/g' $f > $DIR/$(basename $f)
done
touch $OUT
cat results/node-selection/closed-instances.txt results/anytime/open-instances.txt |
xargs -P ${JOBS:-3} -n1 sh -c '
  i=$0; grep -q "\"instance\": \"$i\", \"arm\": \"'"$ARM"'\"" '"$OUT"' && exit 0
  r=$(timeout $(($TLIM+200)) '"$BIN"' '"$DIR"'/$i --solve --bisector lsmear-lffix --relax both --loup ipoptxn \
       '"$FLAGS"' --random-seed 1 --timeout $TLIM 2>/dev/null | tail -1)
  [ -z "$r" ] && r="{\"status\": \"killed\"}"
  echo "{\"instance\": \"$i\", \"arm\": \"'"$ARM"'\", \"result\": $r}" >> '"$OUT"'
'
