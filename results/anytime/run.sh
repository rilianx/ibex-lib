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
  L=ipoptxn
  case $arm in
    base)     x="" ;;
    ipoptfix) x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1" ;;
    fix2)     x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1" ;;  # ipoptfix, affine fix (9e69f77)
    tfinal)   x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1" ;;  # root time counted
    tms20)    x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-multistart=20" ;;
    tbase)    x="" ;;
    tg10)     x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.1" ;;
    tg10l)    x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.1 --ipopt-guided-local" ;;
    tg10p)    x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.1 --ipopt-guided-patience=20" ;;
    tcvx)     x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --cvx-close" ;;
    gcvx)     x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.05 --ipopt-guided-patience=20 --cvx-close" ;;
    full)     x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.05 --ipopt-guided-patience=20 --cvx-close --cvx-backoff=64" ;;
    abase) x="" ;;
    # tuned Ipopt schedule: adaptive backoff (max 100) + 1 s per call
    aad) x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-adaptive=100 --ipopt-call-limit=1" ;;
    cad) x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-adaptive=100 --ipopt-call-limit=1 --cvx-close --cvx-backoff=64 --cvx-test 2" ;;
    fad) x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-adaptive=100 --ipopt-call-limit=1 --ipopt-guided=0.05 --ipopt-guided-patience=20 --ipopt-guided-nodes --cvx-close --cvx-backoff=64 --cvx-test 2" ;;
    # ibex-ipopt schedule (freq 20, maxiter 1000 then 20, correction only when check fails and Ipopt converged)
    ami) x="--ipopt-freq 20 --ipopt-maxiter-feasible 20 --ipopt-corr-policy 2" ;;
    cmi) x="--ipopt-freq 20 --ipopt-maxiter-feasible 20 --ipopt-corr-policy 2 --cvx-close --cvx-backoff=64 --cvx-test 2" ;;
    acert0) x="--cvx-close --cvx-backoff=64 --cvx-test 0" ;;  # E2: certify only (Gerschgorin)
    acert1) x="--cvx-close --cvx-backoff=64 --cvx-test 1" ;;  # E2: certify only (+ comparison matrix, interval Cholesky)
    ascaled) x="--cvx-close --cvx-backoff=64 --cvx-test 3" ;;  # E3: scaled Gerschgorin shift
    acomb) x="--cvx-close --cvx-backoff=64 --cvx-test 4" ;;  # certify (Cholesky) else scaled shift
    acombr) x="--cvx-close --cvx-backoff=64 --cvx-test 4 --cvx-reduce" ;;  # acomb + range reduction with the bound
    acombr1) x="--cvx-close --cvx-backoff=64 --cvx-test 4 --cvx-reduce --cvx-floor 1" ;;  # acombr + exact minimum of each term
    acombr2) x="--cvx-close --cvx-backoff=64 --cvx-test 4 --cvx-reduce --cvx-floor 2" ;;  # acombr1 + interval-gradient floor
    astatsc) x="--cvx-close --cvx-backoff=64 --cvx-test 4 --cvx-stats" ;;  # E4 with acomb
    anobo) x="--cvx-close --cvx-test 2" ;;  # E3: shift without backoff
    astats) x="--cvx-close --cvx-backoff=64 --cvx-test 2 --cvx-stats" ;;  # E4: acvx0 + statistics
    anoip) x=""; L=default ;;  # no Ipopt at all (LoupFinderDefault)
    aipopt) x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1" ;;
    acvx0) x="--cvx-close --cvx-backoff=64 --cvx-test 2" ;;
    acvx) x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --cvx-close --cvx-backoff=64 --cvx-test 2" ;;
    afull) x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.05 --ipopt-guided-patience=20 --ipopt-guided-nodes --cvx-close --cvx-backoff=64 --cvx-test 2" ;;
    full2)    x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.05 --ipopt-guided-patience=20 --cvx-close --cvx-backoff=64 --cvx-test 2" ;;
    g05p20)   x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.05 --ipopt-guided-patience=20" ;;
    g20p20)   x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.2 --ipopt-guided-patience=20" ;;
    g10p10)   x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.1 --ipopt-guided-patience=10" ;;
    g10p50)   x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-guided=0.1 --ipopt-guided-patience=50" ;;
    ms20)     x="--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1 --ipopt-multistart=20" ;;
  esac
  r=$(timeout $(($TLIM+200)) '"$BIN"' benchs/optim/all/$i --solve --bisector lsmear-lffix --relax both --loup $L \
       $x --bound-trace --random-seed 1 --timeout $TLIM 2>/dev/null | tail -1)
  [ -z "$r" ] && r="{\"status\": \"killed\"}"
  echo "{\"instance\": \"$i\", \"arm\": \"$arm\", \"result\": $r}" >> '"$OUT"'
'
