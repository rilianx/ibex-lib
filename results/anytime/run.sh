#!/bin/sh
# The instances lsmear-lffix does not close in the 223-instance runs
# (open-instances.txt = benchs/optim/all minus node-selection/closed-instances.txt),
# with the Ipopt base of results/ipopt-diam and --bound-trace: the evolution of
# loup and uplo in time, for anytime metrics (final gap, primal and
# primal-dual integrals). One line of JSON per (instance, arm) in OUT; resumable.
#   TLIM=310 BIN=build/bin/ibexopt-ml JOBS=3 ARMS="ipoptfix" results/anytime/run.sh
cd "$(dirname "$0")/../.."
BIN=${BIN:-build/bin/ibexopt-ml}
OUT=${OUT:-results/anytime/runs.jsonl}
TLIM=${TLIM:-310}; export TLIM
$BIN --help 2>/dev/null | grep -q bound-trace || { echo "$BIN no conoce --bound-trace: recompila"; exit 1; }
touch $OUT
for arm in ${ARMS:-ipoptfix}; do while read i; do echo "$arm $i"; done < ${INSTANCES:-results/anytime/open-instances.txt}; done |
xargs -P ${JOBS:-8} -n2 sh -c '
  arm=$0; i=$1; grep -q "\"instance\": \"$i\", \"arm\": \"$arm\"" '"$OUT"' && exit 0
  case $arm in
    base)     x="" ;;
    ipoptfix) x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1" ;;
    fix2)     x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1" ;;  # ipoptfix, affine fix (9e69f77)
    tfinal)   x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1" ;;  # root time counted
    tms20)    x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-multistart=20" ;;
    tbase)    x="" ;;
    tg10)     x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.1" ;;
    tg10l)    x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.1 --ipopt-guided-local" ;;
    ms20)     x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-multistart=20" ;;
  esac
  r=$(timeout $(($TLIM+200)) '"$BIN"' benchs/optim/all/$i --solve --bisector lsmear-lffix --relax both --loup ipoptxn \
       $x --bound-trace --random-seed 1 --timeout $TLIM 2>/dev/null | tail -1)
  [ -z "$r" ] && r="{\"status\": \"killed\"}"
  echo "{\"instance\": \"$i\", \"arm\": \"$arm\", \"result\": $r}" >> '"$OUT"'
'
