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
