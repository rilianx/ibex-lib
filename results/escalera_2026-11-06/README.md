# Reproducing the experiments of the SC paper

This folder holds the raw results of every experiment in
*Certificate-driven polytope-hull contraction for interval branch and bound*
(`examples/dfb/paper/main.tex`), the scripts that ran them, and the scripts
that turn them into the tables. The code is in `examples/dfb` of the same
commit.

Comments inside the scripts are in Spanish, as are some run labels.

## 1. Build

Ibex is built in place from the repository root. It uses GAOL 4.2.0 and
SoPlex 4.0.2, which are the defaults:

    cmake -DCMAKE_BUILD_TYPE=Release -DINTERVAL_LIB=gaol -DLP_LIB=soplex . && make
    examples/dfb/build_paper.sh

This produces four binaries in `examples/dfb`:

- `ibexopt_esc`: the optimizer;
- `ibexsolve_dfb_esc`: the system solver;
- `ibexreplay_esc` and `ibexreplay_cmp`: the replay of a bank of passes.

The core build covers every rung, R0 to R5. Three other engines appear only in
side remarks of Section 5: SoPlex 8.1.0, HiGHS 1.15.1 and Gurobi 12.0.3. To
link them, use `build_paper.sh --externals`; the header of the script gives
the details. Without them, selecting one of those engines stops with a
message.

We checked the rebuilt binaries against the archived outputs. Replaying the
affine bank of Brent-10 gives the same output, line for line, on the SC side,
the tuned-SoPlex side and the bare-simplex side. A 60 s search with the affine
relaxation gives the same tree.

Runs with `compo` on systems of equations are not reproducible bit for bit,
not even with the original binaries. With `compo`, SC chooses between the
primal and the dual warm start by measuring CPU time per LP (Section 4.1 of the
paper), so the tree varies by a few percent from run to run. The banks of
passes are not affected, because they replay fixed calls.

## 2. Data

    mkdir work && tar xzf datos_trabajo.tgz -C work

Every script reads `work/S` and `work/M`, or the folders named by the
environment variables `S` and `M`.

The paper uses the limit of 100 iterations per LP in every engine. The runs
with that limit are in `work/S100` (unpack `datos_l100.tgz` into `work/`); run
the scripts with `S=work/S100`. `l100.sh` ran them, 14 at a time. `r7.sh` later added the hybrid (R5) on optimization
problems and on systems with the affine relaxation, 8 at a time (jobs in
`S100/r7_jobs.txt`), so that Table 2 has the hybrid in every class. Rungs that
did not change come from `work/S`: the recorded banks (symbolic links), and R0 and
R1 inside the search (copied into `escf_c`, `escfa_c`), which already used the limit of Ibex. `work/S` holds the
earlier runs, where the tuned SoPlex, the bare simplex and SC used the limit
20m+200. Most ratios of the paper differ between the two by 0.01–0.04; the
largest change is SC over the bare simplex on large relaxations (1.48 against
1.42 in geometric mean on the `compo` bank, 1.23 against 1.15 beyond 120 rows
inside the search).

| Folder in `work/` | Contents |
|---|---|
| `M/banco`, `S/bancoopt` | bank of passes, affine relaxation: 1000 recorded calls (box and relaxation) per instance |
| `S/bancocmp` | bank of passes with `compo`: systems with more than 30 rows, one per family |
| `M/filas_m2` | mean rows of the merged relaxation, per instance and relaxation |
| `S/bescf`, `S/bcmp` | bank replays: rungs b0, b0f, b1, b2 (see below) |
| `S/r0c` | bank replays of R0 and R1 with the nearest bound first |
| `S/bescf_c` | `bescf` with b0 and b0f taken from `r0c` |
| `S/escf`, `S/escfa` | solver runs, 60 s, seeds 1–4 (opt) and 1–2 (systems); `escfa` uses the affine relaxation on systems |
| `S/escf_r0c`, `S/escfa_r0c` | solver runs of R0 and R1 with the nearest bound first |
| `S/escf_c`, `S/escfa_c` | `escf`/`escfa` with R0 and R0f taken from the `_r0c` folders |
| `S/hyb` | real per-node hybrid, thresholds 15, 20 and 30 |
| `S/emu` | emulated work of SoPlex inside SC, cascade e0–e4 |
| `S/peso` | weight of the pass in the CPU time of a search |
| `S/kat2` | cost split on Katsurabis-14 |
| `S100/...` | the same folders, with the limit of 100 (`bescf_c` is a copy of `bescf`, whose b0 and b0f already use the nearest bound; `escf_c`, `escfa_c` take R0 and R0f from `S/escf_c`, `S/escfa_c`; `r0c` is empty) |
| `S/*_jobs.txt`, `S/*_lista.txt`, `S/bopt_sel.txt` | job lists and the selection of one instance per family |

The folders `banco/`, `banco_compo/`, `r0_cercana/`, `hibrido/`, `emulado/`,
`peso/`, `katsurabis_scbase/` and `solver_escf_escfa.tgz` hold the same results
unpacked or by topic.

### SoPlex 8.1.0 (footnote of R2)

`soplex8/job.sh` replays the banks with the tuned SoPlex 4.0.2 (`s4`) and with
SoPlex 8.1.0 through `ibex_Caja7.cpp`, with and without scaling (`s8`, `s8n`),
presolving off in all three; the binary `ibexreplay_s8` links SoPlex 8.1.0
built with `-Dsoplex=soplex8 -DBOOST=OFF -DGMP=OFF -DMPFR=OFF -DPAPILO=OFF
-DZLIB=OFF`. The run was stopped after 135 instances because other
experiments were loading the machine; the three variants of an instance ran
one after the other in the same job, so their ratio is comparable, but their
absolute times are not. The outputs are in `soplex8/out`.

## 3. Rungs of the ladder

| Paper | Solver run (`escf.sh`) | Bank side (`bescf.sh`, `bcmp.sh`) |
|---|---|---|
| R0: SoPlex as Ibex uses it | `R0` | box side of `b0` |
| R1: + merge twin rows | `R0f` | box side of `b0f` |
| R2: tuned SoPlex | `R1` | box side of `b1` |
| R3: bare simplex (SPxSolver) | `R2` | box side of `b2` |
| R4: SC base | `R3` | SC side, mean of `b0f` and `b1` |
| R5: hybrid, SC up to 30 rows | `R7` | R4 or R3 by the size of the relaxation |

The environment variables of each rung are written in the scripts.

R0 and R1 use the nearest bound first, as `CtcPolytopeHull` does in Ibex
(`DFBH_ORDEN_CAJA=cerca`). The scripts `r0c.sh`, `r0b.sh` and `r0o.sh` ran
them. The `nat` order in `escf.sh`, `bescf.sh` and `bcmp.sh` belongs to an
earlier R0 that the paper no longer uses; the `_c` folders replace it.

The SC rung of `l100.sh` (`R3`) sets `DFB_SX_BFRT=0 DFB_SX_DEVEX=0 DFBH_LPCAP=100`,
and the SoPlex rungs `DFBH_CAJA_TOPE=100`. These are the defaults of the
binaries (devex and the long step since commit f1b3717e5, the limit since the
commit that adds `work/S100`).

## 4. From results to the paper

| Paper | Script |
|---|---|
| Table 1, ladder per pass: optimization and systems up to 30 rows | `python3 bescf_tab_c.py` |
| Table 1, systems over 30 rows with `compo` | `python3 bcmp_tab_c.py` |
| Median pass times and contraction in Section 5.1 | `python3 bescf_extra_c.py` |
| Effect of the order of R0 | `python3 r0c_tab.py` |
| Table 2, solver | `python3 esc_fam_c.py` |
| Hybrid inside the search (offline threshold runs) | `python3 hyb_tab.py` (needs `escf` with R0: point `S` to a folder where `escf` links to `escf_c`) |
| Offline threshold and SC over R3 by rows | `python3 umbral.py` (same remark) |
| Figure 1, emulated work, weight of the pass | `OUT=examples/dfb/paper python3 fig_peso_emu.py` |
| Side experiments of Section 5.1 (representation, patched SoPlex, tolerance) | `grandes_2026-11-11/an6.py` (runs: `job6.sh`, `run6/`, base in `run4/`) |

The Katsurabis-14 split (`S/kat2`) was computed inline from the archived
outputs.

Each script prints its numbers for every class of instance, not only the ones
quoted in the paper.

All the runs were made on two Intel Xeon E5-2620 v4 (2.1 GHz, 16 cores), with
at most 6 GB each: `work/S` eight runs at a time, `work/S100` fourteen.
