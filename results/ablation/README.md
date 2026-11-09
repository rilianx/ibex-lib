# Ablation (deterministic budgets, same host, 298 instances, 310 s)

Arms (all with the affine fix, the thick-equality fix and root time counted;
lsmear-lffix, `--relax both`, `--loup ipoptxn`):

| arm | adds | flags |
|---|---|---|
| abase | original Ipopt schedule (boxes <= 1e8 wide) | `--ipopt-diam=1e8` |
| aipopt | Ipopt everywhere + cascade certification | `--ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1` |
| acvx0 | Lagrangian bound only (old schedule) | `--ipopt-diam=1e8 --cvx-close --cvx-backoff=64 --cvx-test 2` |
| acvx | Ipopt + Lagrangian bound | aipopt + `--cvx-close --cvx-backoff=64 --cvx-test 2` |
| afull | + guided restarts (node budget) | acvx + `--ipopt-guided=0.05 --ipopt-guided-patience=20 --ipopt-guided-nodes` |

Data: `results/ipopt-diam/runs.jsonl` (223 instances closed by every arm) and
`results/anytime/runs.jsonl` (75 open instances, `--bound-trace`).

## Closed instances (of 298)

| arm | 223 set | 75 open set | total |
|---|---|---|---|
| abase | 223 | 5 (chenery-1, hs093, hs118, mistake, wall) | 228 |
| aipopt | 223 | 5 | 228 |
| acvx0 | 223 | 13 (+ abel, abelbis, abelbis-1, abelter-1, catena, ex8_4_1, fccu, polak3) | 236 |
| acvx | 223 | 13 | 236 |
| afull | 223 | 13 | 236 |

The 8 newly closed instances are all due to the Lagrangian bound.

## 223 closed instances

Time diff = mean of log10(max(t,1 s)) − log10(max(t_base,1 s)), bootstrap IC95.
Nodes = geometric mean of (nodes+1) ratios vs abase. Time/node = geometric
mean of per-instance ratios vs abase over instances with >= 100 nodes in both.
Breakdown = share of the total time.

| arm | total time | time diff vs abase [IC95] | nodes | time/node | Ipopt | correction | cvx bound |
|---|---|---|---|---|---|---|---|
| abase | 1184 s | — | 1.000 | 1.00 | 2.7 % | 0.5 % | 0 |
| aipopt | 1188 s | −0.000 [−0.007, +0.006] | 0.763 | 1.02 | 3.1 % | 0.4 % | 0 |
| acvx0 | 1030 s | −0.021 [−0.045, −0.004] | 0.798 | 1.02 | 3.0 % | 0.6 % | 0.4 % |
| acvx | 1048 s | −0.020 [−0.044, −0.002] | 0.567 | 1.05 | 3.5 % | 0.4 % | 0.5 % |
| afull | 942 s | −0.022 [−0.054, +0.003] | 0.513 | 1.12 | 5.3 % | 0.5 % | 0.4 % |

## 75 open instances

Integrals over [0, 310 s] normalized by 310 s (`analyze.py`), means over the 62
instances with metrics in every arm (9 unbounded, 2 infeasible, 2 killed
excluded). "No feasible point" = timeouts without a loup.

| arm | primal integral | primal-dual integral | no feasible point | time/node vs abase (timeouts) | Ipopt | correction | cvx bound |
|---|---|---|---|---|---|---|---|
| abase | 0.162 | 0.576 | 8 | 1.00 | 9.1 % | 4.3 % | 0 |
| aipopt | 0.145 | 0.574 | 4 | 1.00 (max 1.7) | 10.0 % | 3.1 % | 0 |
| acvx0 | 0.162 | 0.493 | 8 | 0.99 (max 1.1) | 10.0 % | 5.1 % | 0.6 % |
| acvx | 0.139 | 0.491 | 4 | 1.00 (max 1.7) | 12.3 % | 3.6 % | 0.7 % |
| afull | 0.108 | 0.481 | 4 | 1.01 (max 1.9) | 12.4 % | 4.7 % | 0.8 % |

## Reading

- The components are complementary. Ipopt everywhere improves the upper
  bound (primal integral −10 %, feasible points on 4 more instances, −24 %
  nodes) but not the lower bound. The Lagrangian bound improves only the lower
  bound (primal-dual integral 0.576 → 0.493, +8 closed) and is what closes
  instances. Guided restarts improve the primal side further (primal integral
  0.139 → 0.108).
- Together: half the nodes and −20 % total time on the 223 set, 236/298
  closed, primal integral −33 %, primal-dual integral −16 %.
- Time per node does not blow up: +2–5 % with Ipopt and the bound, +12 % with
  guided restarts (closed set); unchanged on the open set (max 1.9x on a single
  instance). The bound itself costs < 1 % of the time; the Ipopt correction
  (certification box) 0.4–5 %.
- The pairwise time metric barely separates the arms on the 223 set because
  most of those instances take < 1 s (1 s floor).
