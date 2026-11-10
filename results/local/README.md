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
