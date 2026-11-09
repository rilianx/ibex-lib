# Certified-convex Lagrangian bound (--cvx-stats, --cvx-close)

At a node with box B, take the multipliers lambda of the last usable Ipopt
run (sign-projected for weak duality: >= 0 for f <= 0, <= 0 for f >= 0).
If the interval Hessian of L(y) = f(y) + sum lambda_k g_k(y) on B is
certified PSD (Gershgorin: nonnegative diagonal, diagonal dominance), L is
convex on B and, for x^ = Ipopt's point projected onto B,

    min { f(y) : y in B feasible } >= min_{y in B} L(x^) + grad L(x^)^T (y - x^)
                                       - sum_{equalities} |lambda_k| eps_h

(weak duality + convexity; the right side is linear, minimised exactly with
interval arithmetic; the eps_h term because the search relaxes equalities to
|h| <= eps_h). --cvx-close applies it as y >= lb in contract_and_bound.

## How often (measurement, 298 instances, 60 s, runs.jsonl)

Certified on 2.2% of the nodes, closing on 0.9%, concentrated on 24
instances: abel, fccu, abelter-1, abelbis, abelbis-1 (open; ~100% of nodes),
polak3 (open, 72%), harker (100%), hs113, srcpm, avgasa, sambal, schwefel5,
dipigri, hs100 (50-99%), and a few with < 50%.

## First test (--cvx-close, 120 s; base = tfinal, 310 s)

    instance    base                                   --cvx-close
    abel        timeout, uplo -1.03e6 (loup 225.19)    0.03 s, 0 nodes   (loup 225.1946, uplo 224.97)
    abelter-1   timeout, uplo 22.5                     0.03 s, 0 nodes
    fccu        timeout, uplo 5.85 (loup 11.149)       0.01 s, 2 nodes
    polak3      timeout, uplo 5.61 (loup 5.933)        28.2 s, 4706 nodes
    harker      125.2 s, 35078 nodes                   0.01 s, 0 nodes
    sambal      28.1 s, 9416 nodes                     0.01 s, 0 nodes
    hs113       0.78 s, 312                            0.02 s, 0
    schwefel5   0.58 s, 812                            0.04 s, 36
    dipigri     0.96 s, 212                            1.11 s, 48
    hs100       1.42 s, 222                            1.65 s, 66

Optima agree with the proved values (results/dataset-fixedloup/optima.txt)
within the solver's precision; every uplo stays below the optimum.

## 223 closed instances (tcvx = tfinal + --cvx-close, gcvx = g05p20 + --cvx-close)

    A      vs B       time diff [IC95]          nodes  total time (closed)
    tcvx   vs tfinal  +0.019 [+0.000, +0.036]   0.763  1384 vs 1669 s
    gcvx   vs g05p20  +0.007 [-0.011, +0.024]   0.784  1186 vs 1539 s
    gcvx   vs tbase   +0.008 [-0.017, +0.030]   0.523  1186 vs 1662 s

Big wins: harker 125 -> 0.0 s, sambal 28 -> 0.0 s, srcpm 1.4 -> 0.0 s (and
gcvx keeps matrix2 0.3 s, dixchlng, eigencco, launch, makela3, polak5).
Every optimum agrees with optima.txt (ex8_2_4's reference is the known
wrong one). The cost: the interval Hessian is evaluated at every node, also
where nothing is ever certified; ex6_2_13 (closes at ~293 s) times out with
0 raised bounds, and several instances are 1.2-1.4x slower. Next: back off
the check where it keeps failing.

## 75 open instances (310 s)

    arm      closed  primal int.  primal-dual int.   no loup
    tfinal   3       0.156        0.599              5
    g05p20   4       0.125        0.587              4
    tcvx     11      0.140        0.503              5 (hs118 closes)
    gcvx     11      0.138        0.509              4

Newly closed by the convex bound: abel, abelbis, abelbis-1, abelter-1 (at
the root; before: timeouts with lower bounds -1e6, -8e8, -313, 22.5 for an
optimum of 225.19), fccu, polak3; plus hs118 (thick-equality fix, in the
same binary) and mistake (borderline, ~280-300 s). Lower bounds improve on
hs101, hs102, linear, water-abs. gcvx vs tbase: primal-dual integral 0.509
vs 0.601 (27 better / 2 worse), primal 0.138 vs 0.178 (18/3).

Total: 222 of the 223 (ex6_2_13 times out from the per-node Hessian cost) +
11 of the 75 = 233 of 298, from 224 at the start of the Ipopt work.

## The per-node cost and --cvx-backoff (ex6_2_13, machine otherwise idle)

    ex6_2_13 (351770 nodes in every run)      time
    tfinal                                    221.8 s
    + --cvx-close                             227.2 s   (+2.4%)
    + --cvx-close --cvx-backoff=64            218.5 s

The timeout of tcvx/gcvx on ex6_2_13 came from running at ~293 s, next to
the limit, on a loaded machine; the Hessian checks cost ~2.4% there, and
the backoff (check 2x less often after each non-certified node, at most
every 64 nodes, reset by a certified one) removes that cost. On harker,
sambal, abel the backoff keeps the root closing; polak3 30.1 -> 32.6 s.

Recommended configuration:

    --ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1
    --ipopt-guided=0.05 --ipopt-guided-patience=20 --cvx-close --cvx-backoff=64

## Scaled Gershgorin shift (--cvx-test 3) vs plain (--cvx-test 2)

Full configuration, 310 s, same host:

    instance   test 2 (plain)          test 3 (plain or scaled, smaller penalty)
    robot      29.9 s, 36006 nodes     29.8 s, 35684
    ex6_2_13   193.1 s, 306264         183.6 s, 288602
    ex6_2_6    4.9 s, 13778            5.0 s, 13044
    ex8_4_1    27.9 s, 10016           27.0 s, 9878
    catena     0.0 s, 0                0.0 s, 0
    chakra     timeout                 timeout (uplo -188.21 vs -188.23)
    polak3     24.3 s, 5580            32.9 s, 6512
    haldmads   2.2 s, 1490             3.7 s, 1374
    dnieper    22.7 s, 2322            1.5 s, 78

The scaled shift trims nodes a little (2-6%) but brings no clear gain; the
differences in time are within the run-to-run variation, which is large
here: dnieper took 2.1 s with test 2 in the full2 run and 22.7 s now. The
guided-restart budget and the convexity backoff are measured in CPU time,
so the search is not deterministic across runs; path-sensitive instances
(dnieper, polak3, haldmads) swing. Plain Gershgorin (test 2) stays the
default.
