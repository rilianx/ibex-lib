#!/bin/sh
# SCIP (PySCIPOpt) on the 298 benchmark instances, same time limit and gap rule as
# IbexOpt. One JSON line per instance in OUT; resumable.
#   PYTHONPATH=<pyscipopt> JOBS=1 TLIM=310 sh results/solvers/run.sh
cd "$(dirname "$0")/../.."
OUT=${OUT:-results/solvers/scip.jsonl}
TLIM=${TLIM:-310}; export TLIM
touch $OUT
cat results/node-selection/closed-instances.txt results/anytime/open-instances.txt |
xargs -P ${JOBS:-1} -n1 sh -c '
  i=$0; grep -q "\"instance\": \"$i\"" '"$OUT"' && exit 0
  r=$(timeout $(($TLIM+120)) python3 results/solvers/scip_run.py benchs/optim/all/$i $TLIM 2>/dev/null | tail -1)
  case "$r" in "{"*) ;; *) r="{\"status\": \"error\"}" ;; esac
  echo "{\"instance\": \"$i\", \"arm\": \"scip\", \"result\": $r}" >> '"$OUT"'
'
