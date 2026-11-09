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

## Relative times t(A)/max(t(A), t(abase))

Per instance r_A = t_A / max(t_A, t_base) and r_base = t_base / max(t_A, t_base)
(1 = the slower of the two). Timeouts count as 310 s. Only instances with
max(t_A, t_base) > 1 s. Diff = mean(r_A − r_base), bootstrap IC95; W/L = A
faster / slower by more than 5 %.

| arm | set | n | mean r_A | mean r_base | diff [IC95] | W/L |
|---|---|---|---|---|---|---|
| aipopt | 223 | 59 | 0.930 | 0.950 | −0.021 [−0.087, +0.040] | 8/10 |
| acvx0 | 223 | 57 | 0.882 | 0.988 | −0.106 [−0.174, −0.046] | 16/3 |
| acvx | 223 | 59 | 0.836 | 0.954 | −0.119 [−0.204, −0.035] | 20/13 |
| afull | 223 | 61 | 0.825 | 0.918 | −0.094 [−0.198, +0.001] | 24/22 |
| aipopt | 75 open | 60 | 0.998 | 1.000 | −0.001 [−0.004, +0.000] | 1/0 |
| acvx0 | 75 open | 60 | 0.870 | 1.000 | −0.130 [−0.215, −0.049] | 8/0 |
| acvx | 75 open | 60 | 0.868 | 1.000 | −0.131 [−0.216, −0.050] | 9/0 |
| afull | 75 open | 60 | 0.870 | 0.999 | −0.129 [−0.215, −0.048] | 8/0 |
| aipopt | all | 119 | 0.964 | 0.975 | −0.011 [−0.045, +0.019] | 9/10 |
| acvx0 | all | 117 | 0.876 | 0.994 | −0.118 [−0.176, −0.068] | 24/3 |
| acvx | all | 119 | 0.852 | 0.977 | −0.125 [−0.185, −0.066] | 29/13 |
| afull | all | 121 | 0.847 | 0.958 | −0.111 [−0.178, −0.048] | 32/22 |

Worst afull slowdowns (223 set): ex8_5_6 1.7 → 6.5 s, haldmads 1.0 → 3.3,
ex8_2_4 2.5 → 5.8, ex8_5_1bis 1.2 → 2.1. On the open set the gains are the 8
instances that close (the rest time out in both arms).

## Time ratio t(A)/t(abase)

Geometric mean of the per-instance ratio (bootstrap IC95) and median, over
instances with max(t_A, t_base) > 1 s; timeouts = 310 s.

| arm | set | n | geomean [IC95] | median |
|---|---|---|---|---|
| aipopt | 223 | 59 | 0.88 [0.69, 1.04] | 1.00 |
| acvx0 | 223 | 57 | 0.69 [0.45, 0.93] | 0.99 |
| acvx | 223 | 59 | 0.58 [0.36, 0.86] | 0.99 |
| afull | 223 | 61 | 0.55 [0.32, 0.86] | 1.01 |
| aipopt | 75 open | 60 | 1.00 [1.00, 1.00] | 1.00 |
| acvx0 | 75 open | 60 | 0.42 [0.21, 0.74] | 1.00 |
| acvx | 75 open | 60 | 0.34 [0.16, 0.69] | 1.00 |
| afull | 75 open | 60 | 0.35 [0.16, 0.70] | 1.00 |
| aipopt | all | 119 | 0.94 [0.83, 1.02] | 1.00 |
| acvx0 | all | 117 | 0.53 [0.35, 0.75] | 1.00 |
| acvx | all | 119 | 0.45 [0.28, 0.67] | 1.00 |
| afull | all | 121 | 0.44 [0.27, 0.66] | 1.00 |

The median is ~1: the typical instance does not change; the geometric mean is
driven by a minority with large speedups (instances the bound closes in
seconds that took minutes or timed out).
