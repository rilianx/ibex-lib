# Local runs (EXPERIMENTS.md §9)

## Machine

- CPU: 2 x Intel Xeon E5-2620 v4 @ 2.10 GHz (Broadwell, turbo 3.0 GHz), 16 physical
  cores, no hyperthreading; 62 GB RAM. Governor `schedutil`.
- OS: Ubuntu 24.04.1 LTS, gcc 13.3.0.
- Ipopt 3.11.9 (`coinor-libipopt-dev`, libmumps-seq 5.6.2). The machine also has an
  Ipopt 3.14.19 in `/usr/local` that pkg-config finds first, so the build was
  configured with
  `PKG_CONFIG_LIBDIR=/usr/lib/x86_64-linux-gnu/pkgconfig:/usr/lib/pkgconfig:/usr/share/pkgconfig`
  (CMake reports `Ipopt 3.11.9 found`; `ldd` shows `libipopt.so.1`).
- Build: `cmake .. -DCMAKE_BUILD_TYPE=Release -DLP_LIB=soplex -DIBEX_WITH_IPOPT=ON`
  (-O3), commit c605c1afb (last `src/` change: 95afa9ad5).
- All runs: JOBS=14, TLIM=600. Shared machine (other users logged in), load < 1
  before each run.

## E0: validation of the build

```bash
env OUT=results/local/runs-closed.jsonl BIN=build-fix/bin/ibexopt-ml JOBS=14 TLIM=600 ARMS="abase" sh results/ipopt-diam/run.sh
python3 results/paper/compare.py --tlim 600 --closed-file results/local/runs-closed.jsonl \
    --open-file results/local/runs-open.jsonl --set closed --ref abase abase
```

Wall time: 10 min (09:06-09:16, 2026-11-10).

- 223/223 complete, no killed runs.
- **Same node counts as the cloud `abase` on 221/223 instances (99.1 %) -> passes
  (>= 98 %).** Differences: `launch` (16 -> 26 nodes), `robot` (75742 -> 75882).
- **Speed: this machine is about 2.7x slower than the cloud container per run.**
  local/cloud time ratio on the 55 instances with equal nodes and > 1 s in the cloud:
  median 2.73, p10 2.36, p90 2.99 (geometric mean 2.35 on 73 with > 1 s in either).
  It is not load: `ex6_2_11` and `ex6_2_6` alone give the same times as under JOBS=14
  (10.2 s, 25.0 s vs 3.4 s, 8.2 s in the cloud), and the ratio is the same inside
  Ipopt (ipopt_time 3x). Consequence: TLIM=600 here is worth about 220 s of the cloud
  container, less than the 310 s of the cloud runs, so the 75 "open" instances are
  not expected to close more often here than in the cloud.

```
instances with every arm: 223, common set: 75, closed by every arm: 223
| arm | closed | n (closes) | time vs abase [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| abase | 223 | 75 (75) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 2.9 | 0.0 |
```

## Time limit: 1200 s from E1 on

Given the 2.7x slowdown above, TLIM was raised from 600 to **1200 s** for E1 and
later (about 440 s of the cloud container, above its 310 s). The E0 `abase` runs are
kept: all 223 closed in < 600 s (max 578 s), and the search is deterministic, so a
1200 s limit gives the same result. compare.py is run with `--tlim 1200`.

## E1: final ablation, 223 set

```bash
env OUT=results/local/runs-closed.jsonl BIN=build-fix/bin/ibexopt-ml JOBS=14 TLIM=1200 ARMS="anoip abase acvx0" sh results/ipopt-diam/run.sh
python3 results/paper/compare.py --tlim 1200 --closed-file results/local/runs-closed.jsonl \
    --open-file results/local/runs-open.jsonl --set closed --ref abase anoip abase acvx0
```

Wall time: 26 min (09:21-09:47, 2026-11-10; `abase` skipped, from E0). Max time of a
closed run: abase 578 s, acvx0 526 s.

```
instances with every arm: 223, common set: 93, closed by every arm: 219
| arm | closed | n (closes) | time vs abase [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| abase | 223 | 93 (93) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 2.9 | 0.0 |
| anoip | 219 | 93 (90) | 1.68 [1.36, 2.14] | 14/56 | 3.00 | 0.85 | 0.0 | 0.0 |
| acvx0 | 223 | 93 (93) | 0.85 [0.73, 0.95] | 18/0 | 0.80 | 1.01 | 3.2 | 0.3 |
```

acvx0 passes the acceptance rule against abase (0.85, IC95 [0.73, 0.95] < 1; 18 wins,
0 losses).

anoip does not close `eigencco`, `dnieper`, `ex7_2_3` (timeout at 1200 s) and
**`robot` is `killed`** (no output; `timeout` kills it at TLIM+200 s). Same as in the
cloud runs (`results/ipopt-diam/runs.jsonl`: robot/anoip `killed`). It is a hang, not a
crash: with `--timeout 30` it is still running at 120 s, at 100 % CPU and 14 MB RSS.
Backtrace (gdb, interrupted at 40 s):

```
ibex::Interval::pi()
ibex::ibwd_trigo(...)
ibex::CompiledFunction::backward<ibex::InHC4Revise>(...)
ibex::InHC4Revise::iproj(...)
ibex::LoupFinderInHC4::find(...)
ibex::LoupFinderDefault::find(...)
ibex::Optimizer::update_loup(...)
ibex::Optimizer::contract_and_bound(...)
```

i.e. the inner projection of a trigonometric term in LoupFinderInHC4 (used only by
`--loup default`, i.e. anoip) does not return and the time limit is not checked inside
it. compare.py excludes `killed` runs from the common set.

## E2: certify only vs shift, 223 set

```bash
env OUT=results/local/runs-closed.jsonl BIN=build-fix/bin/ibexopt-ml JOBS=14 TLIM=1200 ARMS="acert0 acert1" sh results/ipopt-diam/run.sh
python3 results/paper/compare.py --tlim 1200 --closed-file results/local/runs-closed.jsonl \
    --open-file results/local/runs-open.jsonl --set closed --ref abase abase acvx0 acert0 acert1
python3 results/paper/compare.py --tlim 1200 --closed-file results/local/runs-closed.jsonl \
    --open-file results/local/runs-open.jsonl --set closed --ref acvx0 acvx0 acert0 acert1
```

Wall time: 12.5 min (10:42-10:55, 2026-11-10). 223/223 complete in both arms (max
585 s), nothing killed.

```
instances with every arm: 223, common set: 76, closed by every arm: 223
| arm | closed | n (closes) | time vs abase [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| abase | 223 | 76 (76) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 2.9 | 0.0 |
| acvx0 | 223 | 76 (76) | 0.81 [0.67, 0.94] | 19/0 | 0.80 | 1.01 | 3.3 | 0.3 |
| acert0 | 223 | 76 (76) | 0.87 [0.71, 1.01] | 6/1 | 0.83 | 1.00 | 3.2 | 0.0 |
| acert1 | 223 | 76 (76) | 0.78 [0.62, 0.92] | 11/0 | 0.71 | 1.00 | 3.1 | 0.0 |

instances with every arm: 223, common set: 71, closed by every arm: 223
| arm | closed | n (closes) | time vs acvx0 [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| acvx0 | 223 | 71 (71) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 3.3 | 0.3 |
| acert0 | 223 | 71 (71) | 1.08 [1.03, 1.13] | 0/15 | 1.04 | 0.99 | 3.2 | 0.0 |
| acert1 | 223 | 71 (71) | 0.95 [0.83, 1.06] | 5/14 | 0.90 | 0.99 | 3.1 | 0.0 |
```

- acert0 (Gerschgorin certification only) is slower than the shift: 1.08 [1.03, 1.13]
  vs acvx0, 0 wins / 15 losses.
- acert1 (+ comparison matrix and interval Cholesky) is not distinguishable from
  acvx0 in time (0.95 [0.83, 1.06]), with 5 wins / 14 losses but **fewer nodes (0.90
  vs acvx0, 0.71 vs abase)**. Unexpected: a certify-only bound should not prune more
  than the shifted one; to look at per instance (which instances drive the node ratio).
- The common set differs from E1 (76 vs 93) because anoip is not in these tables.

Per instance (acert1 vs acvx0): 27 instances differ in nodes, 8 with fewer nodes in
acert1. The geometric mean is driven by 4 convex problems that acert1 closes at the
root (0 nodes) because the interval Cholesky certifies convexity there, while the
Gerschgorin shift of acvx0 does not prune at the root: `dualc5` (23186 nodes, 47.6 s
-> 0), `dualc1` (522 -> 0), `dualc1__coconut` (534 -> 0), `meanvar` (534 -> 0). The
shift is better elsewhere: `ex6_2_6` (13914 vs 27000 nodes), `ex6_2_8` (8660 vs
13672), `dnieper` (150 vs 434). Suggests combining both (certify with Cholesky, shift
when it fails).

## E3: variants (scaled shift, no backoff), 223 set

```bash
env OUT=results/local/runs-closed.jsonl BIN=build-fix/bin/ibexopt-ml JOBS=14 TLIM=1200 ARMS="ascaled anobo" sh results/ipopt-diam/run.sh
python3 results/paper/compare.py --tlim 1200 --closed-file results/local/runs-closed.jsonl \
    --open-file results/local/runs-open.jsonl --set closed --ref acvx0 acvx0 ascaled anobo
python3 results/paper/compare.py --tlim 1200 --closed-file results/local/runs-closed.jsonl \
    --open-file results/local/runs-open.jsonl --set closed --ref abase abase acvx0 ascaled anobo
```

Wall time: 11.5 min (10:55-11:07, 2026-11-10). 223/223 complete in both arms (max
509 s), nothing killed.

```
instances with every arm: 223, common set: 71, closed by every arm: 223
| arm | closed | n (closes) | time vs acvx0 [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| acvx0 | 223 | 71 (71) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 3.3 | 0.3 |
| ascaled | 223 | 71 (71) | 0.99 [0.95, 1.04] | 6/1 | 1.00 | 1.00 | 3.4 | 0.3 |
| anobo | 223 | 71 (71) | 1.02 [1.00, 1.05] | 3/1 | 1.00 | 1.00 | 3.3 | 2.3 |

instances with every arm: 223, common set: 76, closed by every arm: 223
| arm | closed | n (closes) | time vs abase [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| abase | 223 | 76 (76) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 2.9 | 0.0 |
| acvx0 | 223 | 76 (76) | 0.81 [0.67, 0.94] | 19/0 | 0.80 | 1.01 | 3.3 | 0.3 |
| ascaled | 223 | 76 (76) | 0.81 [0.66, 0.94] | 19/1 | 0.79 | 1.00 | 3.4 | 0.3 |
| anobo | 223 | 76 (76) | 0.83 [0.67, 0.96] | 19/2 | 0.79 | 1.01 | 3.3 | 2.3 |
```

- ascaled (scaled Gerschgorin): same as acvx0 (0.99 [0.95, 1.04], same nodes).
- anobo (no backoff): same nodes, bound time 2.3 % of the total instead of 0.3 %, so
  slightly slower (1.02 [1.00, 1.05]); the backoff only saves bound evaluations.

## E1: final ablation, 75 set (and all 298)

```bash
env OUT=results/local/runs-open.jsonl BIN=build-fix/bin/ibexopt-ml JOBS=14 TLIM=1200 ARMS="anoip abase acvx0" sh results/anytime/run.sh
python3 results/paper/compare.py --tlim 1200 --closed-file results/local/runs-closed.jsonl \
    --open-file results/local/runs-open.jsonl --set open --ref abase anoip abase acvx0
python3 results/paper/compare.py --tlim 1200 --closed-file results/local/runs-closed.jsonl \
    --open-file results/local/runs-open.jsonl --set all --ref abase anoip abase acvx0
```

Wall time: 4 h 21 min (11:07-15:28, 2026-11-10).

```
instances with every arm: 75, common set: 12, closed by every arm: 4
| arm | closed | n (closes) | time vs abase [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| abase | 6 | 12 (4) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 11.0 | 0.0 |
| anoip | 4 | 12 (2) | 1.23 [0.99, 1.82] | 1/3 | 1.10 | 0.85 | 0.0 | 0.0 |
| acvx0 | 14 | 12 (12) | 0.03 [0.00, 0.14] | 8/0 | 1.00 | 1.00 | 11.1 | 0.0 |

instances with every arm: 298, common set: 105, closed by every arm: 223
| arm | closed | n (closes) | time vs abase [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| abase | 229 | 105 (97) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 5.2 | 0.0 |
| anoip | 223 | 105 (92) | 1.62 [1.33, 2.04] | 15/59 | 2.94 | 0.85 | 0.0 | 0.0 |
| acvx0 | 237 | 105 (105) | 0.57 [0.41, 0.75] | 26/0 | 0.80 | 1.01 | 5.6 | 0.2 |
```

Instances of the 75 "open" set closed in 1200 s (time in s):

- anoip (4): hs118 (0), wall (0), hs093 (554), mistake (733)
- abase (6): hs118 (0), wall (0), chenery-1 (119), hs093 (603), mistake (639), hs102 (1112)
- acvx0 (14): the 6 of abase (same times) + abel (0), abelbis-1 (0), abelter-1 (0),
  catena (0), abelbis (3), fccu (20), ex8_4_1 (69), polak3 (85)

Statuses, every arm: `infeasible` ex7_3_6, launch__coconut; `no_lower_bound` (9):
3pk, ex8_4inf-1, ex8_5_6-1, like, like-1, modele_opti, ship, test_infinity3,
test_infinity4. Same as in the cloud.

**Killed (hangs):** abase on hs090, hs091; acvx0 on hs089, hs090, hs091 (no output,
killed at TLIM+200 = 1400 s). hs090 and hs091 are also `killed` in the cloud runs
(abase, acvx0). **hs089/acvx0 is new**: in the cloud it ended by `timeout` at 310 s
(6238 nodes); here it hangs after more than 310 s. anoip ends by timeout on all three.
Backtraces (gdb, same flags as run.sh, --timeout 1200, interrupted at 1250 s; one
sample each, binary without debug symbols):

```
hs089 / acvx0                                   hs090 / abase
CompiledFunction::backward<Gradient>            gaol::interval::operator*=
Gradient::jacobian                              CompiledFunction::forward<Eval>
CtcAcid::compute_smearorder                     Eval::eval
CtcAcid::contract                               LinearizerXTaylor::linear_relax
CtcCompo::contract                              LinearizerXTaylor::linearize
Optimizer::contract_and_bound                   LinearizerCompo::linearize
                                                CtcPolytopeHull::contract
                                                CtcCompo::contract
                                                CtcFixPoint::contract
                                                CtcCompo::contract
                                                Optimizer::contract_and_bound
```

Both are inside the contraction of one node (ACID; fixpoint of the polytope hull),
not in Ipopt nor in the Lagrangian bound, and the time limit is only checked between
nodes. One sample cannot tell an infinite loop from a node contraction that never
converges. Runs stopped before E2 on the 75 set, waiting for a decision.

## Hang fix (2d36f15a1) and reruns

Binary rebuilt at 2d36f15a1 ("Fix three hangs that ignored the time limit":
ibwd_trigo period scan, LoupFinderIpopt symbolic differentiation of huge
expressions, deadline checked inside CtcFixPoint/CtcAcid/CtcPropag), same Ipopt 3.11.9.

Check that the fix does not change other runs: 8 instances of the 223 set (ex6_2_6,
ex6_2_11, hs117, dualc1, ex8_4_4bis, hs100, dnieper, ex6_2_8) x {abase, acvx0}:
same node counts on 16/16, times within 5 %. The other runs are kept.

Deleted every line of robot, hs089, hs090, hs091 from `results/local/` and reran them
(robot: the 7 arms of the 223 set; hs089-091: anoip, abase, acvx0, acert0) with
`INSTANCES=<list>`, JOBS=7+7, TLIM=1200. No `killed` run left:

- robot closes with every arm: anoip 86214 nodes / 193 s (was `killed`), abase 75724 /
  173 s, acvx0 49994 / 86 s, acert0 76026 / 171 s, acert1 75972 / 171 s, ascaled
  50278 / 88 s, anobo 47382 / 84 s.
- hs089, hs090, hs091: `timeout` at 1200 s with every arm (was `killed` for abase,
  acvx0, acert0). hs090 gets a loup (1.5823 with abase, acvx0, acert0); hs089 and
  hs091 none.

Tables after the reruns (the E1-E3 tables above are the earlier ones; only the common
sets and robot change):

```
### 223 set, E1
instances with every arm: 223, common set: 94, closed by every arm: 220
| arm | closed | n (closes) | time vs abase [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| abase | 223 | 94 (94) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 2.9 | 0.0 |
| anoip | 220 | 94 (91) | 1.67 [1.36, 2.17] | 14/57 | 2.98 | 0.85 | 0.0 | 0.0 |
| acvx0 | 223 | 94 (94) | 0.85 [0.72, 0.95] | 19/0 | 0.80 | 1.00 | 3.3 | 0.3 |
### 223 set, E2
| acvx0 | 223 | 71 (71) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 3.3 | 0.3 |
| acert0 | 223 | 71 (71) | 1.08 [1.03, 1.13] | 0/15 | 1.04 | 0.99 | 3.2 | 0.0 |
| acert1 | 223 | 71 (71) | 0.95 [0.83, 1.06] | 5/14 | 0.90 | 0.99 | 3.1 | 0.0 |
### 223 set, E3
| acvx0 | 223 | 71 (71) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 3.3 | 0.3 |
| ascaled | 223 | 71 (71) | 0.99 [0.95, 1.04] | 6/1 | 1.00 | 1.00 | 3.4 | 0.3 |
| anobo | 223 | 71 (71) | 1.02 [1.00, 1.05] | 2/1 | 1.00 | 1.01 | 3.3 | 2.4 |
### 75 set, E1 (unchanged)
| abase | 6 | 12 (4) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 11.0 | 0.0 |
| anoip | 4 | 12 (2) | 1.23 [0.99, 1.82] | 1/3 | 1.10 | 0.85 | 0.0 | 0.0 |
| acvx0 | 14 | 12 (12) | 0.03 [0.00, 0.14] | 8/0 | 1.00 | 1.00 | 11.1 | 0.0 |
### all 298, E1
instances with every arm: 298, common set: 106, closed by every arm: 224
| abase | 229 | 106 (98) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 5.2 | 0.0 |
| anoip | 224 | 106 (93) | 1.61 [1.33, 2.04] | 15/60 | 2.93 | 0.85 | 0.0 | 0.0 |
| acvx0 | 237 | 106 (106) | 0.57 [0.40, 0.75] | 27/0 | 0.80 | 1.00 | 5.7 | 0.2 |
```

E2 on the 75 set is partial (acert0: 44 of 75, run with the old binary before the fix;
the paused runs were killed without writing lines and are rerun with the new one).

## E2: certify only vs shift, 75 set (and all 298)

```bash
env OUT=results/local/runs-open.jsonl BIN=build-fix/bin/ibexopt-ml JOBS=14 TLIM=1200 ARMS="acert0 acert1" sh results/anytime/run.sh
C="python3 results/paper/compare.py --tlim 1200 --closed-file results/local/runs-closed.jsonl --open-file results/local/runs-open.jsonl"
$C --set open --ref acvx0 acvx0 acert0 acert1
$C --set open --ref abase abase acvx0 acert0 acert1
$C --set all --ref acvx0 acvx0 acert0 acert1
```

acert0: 44 instances run with the binary before the hang fix (11-10, 15:28-16:48),
hs089-091 and the other 31 after it (11-11 and 11-12, 07:40-09:41); acert1: all 75
after the fix (the fix does not change node counts, see above). Server idle before
the run (load 0.4). No killed run.

```
instances with every arm: 75, common set: 10, closed by every arm: 12
| arm | closed | n (closes) | time vs acvx0 [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| acvx0 | 14 | 10 (9) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 6.9 | 0.1 |
| acert0 | 12 | 10 (7) | 2.65 [0.96, 11.17] | 1/2 | 1.00 | 0.98 | 6.5 | 0.1 |
| acert1 | 14 | 10 (9) | 0.63 [0.12, 2.27] | 3/1 | 0.98 | 0.98 | 6.7 | 0.1 |

instances with every arm: 75, common set: 13, closed by every arm: 6
| arm | closed | n (closes) | time vs abase [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| abase | 6 | 13 (4) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 6.7 | 0.0 |
| acvx0 | 14 | 13 (12) | 0.03 [0.01, 0.17] | 8/0 | 1.00 | 1.00 | 6.7 | 0.1 |
| acert0 | 12 | 13 (10) | 0.07 [0.01, 0.34] | 6/0 | 1.00 | 1.01 | 6.6 | 0.1 |
| acert1 | 14 | 13 (12) | 0.02 [0.00, 0.14] | 9/0 | 0.96 | 1.01 | 6.9 | 0.1 |

instances with every arm: 298, common set: 81, closed by every arm: 235
| arm | closed | n (closes) | time vs acvx0 [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| acvx0 | 237 | 81 (80) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 5.0 | 0.2 |
| acert0 | 235 | 81 (78) | 1.20 [1.05, 1.50] | 1/17 | 1.04 | 0.99 | 4.7 | 0.1 |
| acert1 | 237 | 81 (80) | 0.91 [0.72, 1.09] | 8/15 | 0.90 | 0.99 | 4.8 | 0.1 |
```

Closed in the 75 set (time in s):

- acert0 (12): abel, abelbis-1, abelter-1, hs118, wall (0), abelbis (3), fccu (20),
  polak3 (69), chenery-1 (124), hs093 (600), mistake (643), hs102 (1095). Misses
  catena and ex8_4_1 (closed by acvx0).
- acert1 (14): the same + catena (0), **portfl1 (0, closed by no other arm)**; misses
  ex8_4_1 (acvx0: 69 s). mistake 511 s (acvx0 643 s).

On all 298: acert0 is slower than the shift (1.20 [1.05, 1.50]); acert1 is not
distinguishable from acvx0 (0.91 [0.72, 1.09], 8 wins / 15 losses) and closes the
same number of instances (237), not the same ones (portfl1 vs ex8_4_1).

## E3: variants, 75 set (and all 298)

```bash
env OUT=results/local/runs-open.jsonl BIN=build-fix/bin/ibexopt-ml JOBS=14 TLIM=1200 ARMS="ascaled anobo" sh results/anytime/run.sh
$C --set open --ref acvx0 acvx0 ascaled anobo
$C --set all --ref acvx0 acvx0 ascaled anobo
```

Wall time 2 h 40 min (2026-11-12, 09:42-12:22), server idle before, no killed run.

```
instances with every arm: 75, common set: 10, closed by every arm: 14
| arm | closed | n (closes) | time vs acvx0 [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| acvx0 | 14 | 10 (8) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 6.7 | 0.2 |
| ascaled | 16 | 10 (10) | 1.04 [0.92, 1.24] | 2/1 | 1.11 | 0.92 | 6.7 | 0.8 |
| anobo | 14 | 10 (8) | 0.94 [0.83, 1.03] | 2/2 | 0.96 | 0.99 | 6.3 | 3.4 |

instances with every arm: 298, common set: 81, closed by every arm: 237
| arm | closed | n (closes) | time vs acvx0 [IC95] | W/L | nodes | t/node | ipopt% | bound% |
|---|---|---|---|---|---|---|---|---|
| acvx0 | 237 | 81 (79) | 1.00 [1.00, 1.00] | 0/0 | 1.00 | 1.00 | 4.9 | 0.3 |
| ascaled | 239 | 81 (81) | 1.00 [0.96, 1.04] | 8/2 | 1.01 | 0.99 | 5.0 | 0.6 |
| anobo | 237 | 81 (79) | 1.01 [0.98, 1.05] | 4/3 | 0.99 | 1.00 | 4.8 | 2.8 |
```

- ascaled closes the 14 of acvx0 plus **ex8_4_6 (983 s) and ex8_4_6ter (990 s)**, so
  239 vs 237 on all 298, same time otherwise (1.00 [0.96, 1.04]). chenery-1 is slower
  with it (258 s vs 120 s).
- anobo: same instances as acvx0, time not distinguishable (1.01 [0.98, 1.05]); the
  bound takes 2.8 % of the time instead of 0.3 %.

## E4: where the bound helps (astats = acvx0 + --cvx-stats), both sets

```bash
env OUT=results/local/runs-closed.jsonl BIN=build-fix/bin/ibexopt-ml JOBS=14 TLIM=1200 ARMS="astats" sh results/ipopt-diam/run.sh
env OUT=results/local/runs-open.jsonl  BIN=build-fix/bin/ibexopt-ml JOBS=14 TLIM=1200 ARMS="astats" sh results/anytime/run.sh
python3 results/local/e4/summarize.py      # -> results/local/e4/summary.md (new file, analysis only)
```

Wall time: 223 set 9.5 min, 75 set 1 h 21 min (2026-11-12, 12:22-13:53). No killed run.
Statuses identical to acvx0 (223 complete; 75 set: 14 complete, 50 timeout, 9
no_lower_bound, 2 infeasible).

**Unexpected:** --cvx-stats should only count, but the node counts equal acvx0's on
235/237 instances closed by both: ex9_2_5 (10 -> 6 nodes) and robot (49994 -> 50172)
differ. To check in the code (does the measurement touch any state used by the search?).

Counters summed over the instances (percent of `cvx_nodes`, the nodes where the bound
is evaluated):

| counter | 223 set | 75 set | all 298 | instances > 0 (298) |
|---|---|---|---|---|
| cvx_nodes | 5 152 511 | 53 957 296 | 59 109 807 | |
| cvx_with_dual (Ipopt multipliers available) | 12.43 % | 11.95 % | 11.99 % | |
| cvx_certified (Lagrangian Hessian certified PSD) | 0.00 % (200) | 0.85 % | 0.77 % | 20 |
| cvx_better_lb (bound above the node's lb) | 1.94 % | 0.44 % | 0.57 % | 74 |
| cvx_close (bound would close the node) | 0.73 % | 0.04 % | 0.10 % | 47 |
| cvx_raised | 1.23 % | 0.36 % | 0.44 % | |
| cvx_pruned | 0.52 % | 0.05 % | 0.09 % | 57 |
| cvx_center_indef | 0 | 0 | 0 | |

The bound acts on few nodes (0.1 % pruned) but on many instances (57 of 298), and
those prunings are what close the open instances (E1). Instances with the most
certified nodes (full list in `e4/summary.md`): polak6 (98 % of 424 001 nodes
certified, but never closes a node: timeout), genhs28 56 %, schwefel5(-abs) 50 %,
polak3 44 % (2 493 pruned, closes in 85 s), hs113 42 %, abelbis 41 %, linear 32 %
(12 419 pruned, still timeout), sambal 32 %, dipigri 31 %. robot: 50 certified but
5 567 pruned (the shift prunes without certification).

## E5: SCIP (in progress)

SCIP 10.0 through PySCIPOpt 6.3.0 (`pip install --target $HOME/pyscip pyscipopt`, same
SCIP version as the cloud runs); matplotlib 3.11.2 in `$HOME/pympl`.

### E5a: SCIP on the original domains

```bash
PYTHONPATH=$HOME/pyscip OUT=results/local/scip.jsonl JOBS=14 TLIM=1200 sh results/solvers/run.sh
```

Wall time 1 h 42 min (2026-11-12, 13:54-15:36), server idle before.

| status | local (1200 s) | cloud (310 s) |
|---|---|---|
| complete (optimal 93 + gaplimit 129) | 222 | 221 |
| timeout | 62 | 64 |
| error | 9 | 8 |
| infeasible | 4 | 4 |
| unbounded | 1 | 1 |

- Same statuses as the cloud except: **hs103 closes** (timeout at 310 s in the cloud),
  and **ex8_4_6bis is `error`** (timeout in the cloud). The 8 other errors are the
  same as in the cloud (brownal, cresc4, hs088-hs091, minlphi1, modele_opti: not
  translated). ex8_4_6bis runs fine with a 30 s limit (timeout, loup 7.837); being
  reproduced with 1200 s.
- infeasible according to SCIP: ex5_3_3, ex7_3_6, launch__coconut, wall (same as the
  cloud; IbexOpt closes wall with a finite optimum: see the cloud README on SCIP's
  infeasible claims).
- ex8_4inf-1 is `complete` with loup < -1e19 (SCIP's infinity; perfprof counts it as
  not solved).
- SCIP is also ~2.2x slower here than in the cloud (geometric mean on the 54 instances
  both close in > 1 s), like IbexOpt (2.7x).
