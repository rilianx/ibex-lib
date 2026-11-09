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

All 298 instances, restricted to those that **some arm closes** and where some
arm takes > 1 s (n_total = 73); in parentheses, how many of them the arm
closes. Times floored at 1 s, a timeout counts as 310 s (ratios against a
timed-out base are upper bounds). Geometric mean (bootstrap IC95), median,
W/L = faster/slower than abase by > 5 %.

| arm | n_total (closes) | geomean [IC95] | median | W/L |
|---|---|---|---|---|
| abase | 73 (65) | 1.00 | 1.00 | — |
| aipopt | 73 (65) | 1.00 [0.94, 1.04] | 1.00 | 10/10 |
| acvx0 | 73 (73) | 0.52 [0.36, 0.73] | 0.99 | 24/3 |
| acvx | 73 (73) | 0.51 [0.34, 0.72] | 0.99 | 30/14 |
| afull | 73 (73) | 0.50 [0.33, 0.73] | 1.00 | 32/23 |

## Against no Ipopt at all (anoip: `--loup default`, LoupFinderDefault)

Same common-set rule, now including anoip: instances some arm closes and some
arm takes > 1 s: n_total = 81 (in parentheses, how many the arm closes).
Ratio t(A)/t(anoip), times floored at 1 s, timeout = 310 s.

| arm | closed /298 | n_total (closes) | geomean t/t(anoip) [IC95] | median | W/L vs anoip |
|---|---|---|---|---|---|
| anoip | 223 | 81 (69) | 1.00 | 1.00 | — |
| abase | 228 | 81 (73) | 0.62 [0.48, 0.77] | 0.98 | 38/13 |
| aipopt | 228 | 81 (73) | 0.62 [0.47, 0.77] | 0.98 | 38/13 |
| acvx0 | 236 | 81 (81) | 0.35 [0.24, 0.51] | 0.86 | 51/11 |
| acvx | 236 | 81 (81) | 0.34 [0.23, 0.51] | 0.84 | 50/13 |
| afull | 236 | 81 (81) | 0.34 [0.22, 0.50] | 0.85 | 53/18 |

Nodes (geometric mean vs anoip, 223 instances closed by every arm): abase
0.34, aipopt 0.27, acvx0 0.27, acvx 0.20, afull 0.18.

Open set (62 instances with metrics in every arm):

| arm | primal integral | primal-dual integral | no feasible point (timeouts) |
|---|---|---|---|
| anoip | 0.522 | 0.702 | 22 |
| abase | 0.162 | 0.576 | 8 |
| aipopt | 0.145 | 0.574 | 4 |
| acvx0 | 0.162 | 0.493 | 8 |
| acvx | 0.139 | 0.491 | 4 |
| afull | 0.108 | 0.481 | 4 |

Ipopt as primal oracle (anoip -> abase) is the largest single effect: 0.62 of
the time, a third of the nodes, primal integral 0.52 -> 0.16, feasible points
on 14 more open instances, +5 closed. As dual oracle (abase -> acvx0, the
Lagrangian bound from its multipliers): a further 0.35/0.62 = 0.56 of the
time and +8 closed. All together: about a third of the time of the solver
without Ipopt, 18 % of its nodes, +13 closed instances.

## Tuned Ipopt schedule (aad, cad, fad)

`--ipopt-adaptive=100 --ipopt-call-limit=1`: the next scheduled Ipopt call comes
right after a call that improved the loup and twice the previous gap later
(up to 100 loup-finder calls) after one that did not (instead of every 100
calls), and one Ipopt call is cut after 1 s of CPU (instead of 5 s).
aad = aipopt + this, cad = acvx + this, fad = afull + this. The 223 set was
rerun (aad2, cad2, fad2) without SCIP sharing the CPUs; aipopt2 (aipopt rerun
after the container restart) matches aipopt (1.00 [0.99, 1.01]).

223 closed instances (63 with some arm > 1 s; all arms close all 223):

| pair | time ratio [IC95] | faster/slower | nodes |
|---|---|---|---|
| aad / aipopt | 1.00 [0.94, 1.05] | 10/12 | 0.94 |
| cad / acvx | 0.98 [0.92, 1.03] | 17/11 | 0.95 |
| fad / afull | 0.94 [0.88, 1.00] | 11/8 | 0.97 |
| fad / abase | 0.79 [0.61, 0.96] | 24/19 | 0.50 |

Open set (69 instances with metrics in every arm, now including the 7
`no_lower_bound` stops with gap 1):

| arm | closed (open set) | primal integral | primal-dual integral |
|---|---|---|---|
| anoip | 4 | 0.550 | 0.732 |
| abase | 5 | 0.196 | 0.619 |
| aipopt | 5 | 0.181 | 0.617 |
| aad | 5 | 0.130 | 0.617 |
| acvx0 | 13 | 0.196 | 0.545 |
| acvx | 13 | 0.176 | 0.543 |
| cad | 13 | 0.145 | 0.542 |
| afull | 13 | 0.133 | 0.533 |
| fad | 13 | 0.133 | 0.542 |

The tuned schedule does not change the time to close (at most −6 %, within
noise) but improves the primal side on open instances (primal integral −28 %
for aipopt -> aad, −18 % for acvx -> cad), about as much as the guided
restarts; combining both (fad) adds nothing over either.
