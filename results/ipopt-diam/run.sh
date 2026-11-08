#!/bin/sh
# Ipopt's width threshold (--ipopt-diam; IbexOpt calls Ipopt only on boxes
# whose widest variable is at most 1e8) on lsmear-lffix, compo + ipoptxn.
# One line of JSON per (instance, arm) in OUT; resumable.
#   INSTANCES=list TLIM=600 BIN=build/bin/ibexopt-ml JOBS=16 results/ipopt-diam/run.sh
#   python3 results/ipopt-diam/analyze.py
# Arms (ARMS="d1e8 d1e12 dinf"): the threshold 1e8 (default), 1e12, 1e300;
# dnone: no threshold at all, unbounded boxes included (Ipopt starts from the
# point of the box closest to 0 when the box has no finite center).
# x1e5, x1e4: threshold 1e8 and the box around Ipopt's point where the
# correction search looks for a certified point widened from 1e-6 to 1e-5,
# 1e-4 (--ipopt-expansion); nx1e5: no threshold and 1e-5.
# c1e5, nc1e5: same, but the 1e-6 box is tried first (0.1 s) and the 1e-5 one
# only when it certifies nothing (--ipopt-expansion-first).
cd "$(dirname "$0")/../.."
BIN=${BIN:-build/bin/ibexopt-ml}
OUT=${OUT:-results/ipopt-diam/runs.jsonl}
TLIM=${TLIM:-600}; export TLIM
$BIN --help 2>/dev/null | grep -q ipopt-diam || { echo "$BIN no conoce --ipopt-diam: recompila"; exit 1; }
touch $OUT
for arm in ${ARMS:-d1e8 d1e12 dinf}; do while read i; do echo "$arm $i"; done < ${INSTANCES:-results/node-selection/closed-instances.txt}; done |
xargs -P ${JOBS:-8} -n2 sh -c '
  arm=$0; i=$1; grep -q "\"instance\": \"$i\", \"arm\": \"$arm\"" '"$OUT"' && exit 0
  case $arm in
    d1e8)  d=1e8 ;;
    d1e12) d=1e12 ;;
    dinf)  d=1e300 ;;
    dnone) d=-1 ;;
    x1e5)  d=1e8; x=--ipopt-expansion=1e-5 ;;
    x1e4)  d=1e8; x=--ipopt-expansion=1e-4 ;;
    nx1e5) d=-1;  x=--ipopt-expansion=1e-5 ;;
    c1e5)  d=1e8; x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1" ;;
    nc1e5) d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1" ;;
    fixaff) d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1" ;;  # final, with the Affine2 operator= fix
    base8) d=1e8;  x="" ;;  # d1e8, with the Affine2 operator= fix
    ms20)  d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-multistart=20" ;;
    # t*: with the time of the root counted (before, the timer started after the root)
    tbase)  d=1e8; x="" ;;
    tfinal) d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1" ;;
    tms5)   d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-multistart=5" ;;
    tms20)  d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-multistart=20" ;;
    tg10)   d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.1" ;;
    tg10l)  d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.1 --ipopt-guided-local" ;;
    tg10p)  d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.1 --ipopt-guided-patience=20" ;;
    tcvx)   d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --cvx-close" ;;
    gcvx)   d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.05 --ipopt-guided-patience=20 --cvx-close" ;;
    full)   d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.05 --ipopt-guided-patience=20 --cvx-close --cvx-backoff=64" ;;
    tbase2) d=1e8; x="" ;;  # tbase again, same machine as full (timing)
    g05p20) d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.05 --ipopt-guided-patience=20" ;;
    g20p20) d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.2 --ipopt-guided-patience=20" ;;
    g10p10) d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.1 --ipopt-guided-patience=10" ;;
    g10p50) d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.1 --ipopt-guided-patience=50" ;;
    final) d=-1;  x="--ipopt-expansion=1e-5 --ipopt-expansion-first=0.1" ;;  # nc1e5 with --ipopt-start-max (binary of 4197d3e+)
  esac
  r=$(timeout $(($TLIM+200)) '"$BIN"' benchs/optim/all/$i --solve --bisector lsmear-lffix --relax both --loup ipoptxn \
       --ipopt-diam=$d $x --random-seed 1 --timeout $TLIM 2>/dev/null | tail -1)
  [ -z "$r" ] && r="{\"status\": \"killed\"}"
  echo "{\"instance\": \"$i\", \"arm\": \"$arm\", \"result\": $r}" >> '"$OUT"'
'
