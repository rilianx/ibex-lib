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
