# Experiments for the paper "Certified Lagrangian lower bounds from local-solver multipliers"

Everything needed to rerun the paper's experiments on a local machine, by hand or
by an independent agent. Branch: `ml-branching`. Draft paper: `results/paper/`.

## 1. Environment

Reference environment (cloud container, where the published numbers come from):
Ubuntu 24.04, gcc 13.3, Intel Xeon @ 2.30 GHz, 4 cores.

```bash
sudo apt-get install -y build-essential cmake pkg-config python3 \
     coinor-libipopt-dev          # Ipopt 3.11.9 with MUMPS (libmumps-seq)
```

Ipopt must be found by pkg-config (`pkg-config --modversion ipopt` -> 3.11.9). Other
Ipopt versions work but change the results (different iterates, different loups).

## 2. Build

```bash
git clone -b ml-branching https://github.com/rilianx/ibex-lib.git && cd ibex-lib
mkdir -p build-fix && cd build-fix
cmake .. -DCMAKE_BUILD_TYPE=Release -DLP_LIB=soplex -DIBEX_WITH_IPOPT=ON
make -j4 ibexopt-ml
cd ..
build-fix/bin/ibexopt-ml --help | grep -c cvx-close    # 1: the bound is compiled in
```

Gaol (intervals) and SoPlex (LP) are bundled. The scripts default to
`BIN=build-fix/bin/ibexopt-ml`.

## 3. Instances

298 instances of the COCONUT benchmark in Minibex form, `benchs/optim/all/`:

- `results/node-selection/closed-instances.txt`: 223 that IbexOpt closes in 310 s;
- `results/anytime/open-instances.txt`: 75 that it does not.

## 4. Configurations ("arms")

Common flags (set by the scripts):
`--solve --bisector lsmear-lffix --relax both --loup ipoptxn --random-seed 1 --timeout 310`.

| arm | meaning | extra flags |
|---|---|---|
| `anoip` | no local solver | `--loup default` instead of `ipoptxn` |
| `abase` | Ipopt for upper bounding (IbexOpt default) | `--ipopt-diam=1e8` |
| `acvx0` | + **Lagrangian bound** (paper's method) | `--ipopt-diam=1e8 --cvx-close --cvx-backoff=64 --cvx-test 2` |
| `afull` | acvx0 + Ipopt on any box, cascade certification, guided restarts | see `results/ipopt-diam/run.sh` |
| `acvx`, `aipopt`, `aad`, `cad`, `fad`, `ami`, `cmi` | rejected variants of the Ipopt call strategy | see the scripts and `results/ablation/README.md` |

Bound variants (section 5.4 of the paper): `--cvx-test 0|1` (certify only: Gerschgorin
/ + comparison matrix and interval Cholesky), `--cvx-test 2` (Gerschgorin shift,
default), `--cvx-test 3` (scaled Gerschgorin), `--cvx-backoff=N` (0: no backoff).
Other options: `--sqr-rewrite` (e*e -> sqr(e), off by default), `--cvx-stats`
(statistics on where the bound is certified / useful).

Arms are defined in two runners (same arm names): `results/ipopt-diam/run.sh` (223
set) and `results/anytime/run.sh` (75 set, adds `--bound-trace`). To add an arm, add a
`case` line to both.

## 5. Running

```bash
# 223 closed instances, then the 75 open ones; resumable (skips done lines)
env BIN=build-fix/bin/ibexopt-ml JOBS=3 TLIM=310 ARMS="anoip abase acvx0" sh results/ipopt-diam/run.sh
env BIN=build-fix/bin/ibexopt-ml JOBS=3 ARMS="anoip abase acvx0" sh results/anytime/run.sh
```

Output: one JSON line per (instance, arm) appended to `results/ipopt-diam/runs.jsonl`
and `results/anytime/runs.jsonl` (status, loup, uplo, nodes, time, ipopt_time,
ipopt_corr_time, cvx_time, ...). Cost: about 25 min per arm for the 223 set and
about 1.5-2 h per arm for the 75 set, with JOBS=3.

### Protocol (important)

1. **Same machine, no other load.** Never compare times of runs made on different
   machines or while something else runs (SCIP, other arms, builds). JOBS must be at
   most the number of physical cores minus one.
2. **Calibrate.** Every time the machine may have changed (new machine, container
   restart), rerun the reference arms there (e.g. as `abase2`, `acvx02`) and compare
   only against those. Observed here: the same arm 1.38x slower after a restart, and
   the slowdown differs between arms (1.27x vs 1.38x): a global factor is not enough.
3. **Determinism.** With the flags above the search is deterministic (budgets in
   nodes/calls, not CPU time): the same arm gives the same node counts. Check it on a
   rerun (same nodes on >= 98 % of the instances); a difference means a wrong binary or
   a different Ipopt.
4. **Acceptance rule** for a change: geometric-mean time ratio on the common set with
   the bootstrap 95 % interval below 1 (`compare.py`). Fewer nodes alone do not count.
5. When waiting for a run in a shell loop, do not use `pgrep -f <pattern>` with a
   pattern that also appears in the waiting command itself (it waits forever).

## 6. Analysis

```bash
python3 results/paper/compare.py --ref abase anoip abase acvx0        # ablation table, 298
python3 results/paper/compare.py --ref abase abase acvx0 --set closed # 223 set only
python3 results/anytime/analyze.py                                     # bounds over time (75 set)
```

`compare.py` prints the markdown table of the paper's ablation (definitions in its
docstring). Current numbers: `results/ablation/README.md`.

## 7. SCIP comparison

```bash
python3 -m pip install --target $HOME/pyscip pyscipopt      # SCIP 10.0 wheel
python3 -m pip install --target $HOME/pympl matplotlib
# SCIP on the original domains, and with missing bounds set to +-1e12
PYTHONPATH=$HOME/pyscip JOBS=3 sh results/solvers/run.sh
PYTHONPATH=$HOME/pyscip SCIP_INFBOUND=1e12 OUT=results/solvers/scip_inf12.jsonl JOBS=3 sh results/solvers/run.sh
# IbexOpt (afull) on the same +-1e12 domains
BIN=build-fix/bin/ibexopt-ml JOBS=3 TLIM=310 sh results/solvers/run_ibex.sh
# performance profiles
PYTHONPATH=$HOME/pympl python3 results/solvers/perfprof.py
PYTHONPATH=$HOME/pympl python3 results/solvers/perfprof12.py
```

`results/solvers/scip_run.py` translates Minibex to PySCIPOpt (7 instances it cannot
translate are reported as `error`). Run SCIP alone on the machine, like any arm.
Results and caveats: `results/solvers/README.md`.

## 8. Paper

`results/paper/intro-related.tex` + `refs.bib`; `sh results/paper/build.sh` builds a PDF
with pandoc + headless Chromium (no LaTeX needed); with LaTeX, `latexmk -pdf`.

## 9. Pending experiments (paper, section 5)

Write all local results to `results/local/` (never append to the cloud files: other
machine):

```bash
OUTC=results/local/runs-closed.jsonl; OUTO=results/local/runs-open.jsonl
env OUT=$OUTC BIN=build-fix/bin/ibexopt-ml JOBS=3 TLIM=310 ARMS="..." sh results/ipopt-diam/run.sh
env OUT=$OUTO BIN=build-fix/bin/ibexopt-ml JOBS=3 ARMS="..." sh results/anytime/run.sh
python3 results/paper/compare.py --closed-file $OUTC --open-file $OUTO --ref abase <arms>
```

| id | what | arms (both sets) |
|---|---|---|
| E0 | validation of the build | `abase` on the 223 set: node counts must match `results/ipopt-diam/runs.jsonl` (arm abase) on >= 98 % of the instances |
| E1 | final ablation | `anoip abase acvx0` |
| E2 | certify only vs shift | `acert0 acert1` (vs `acvx0`) |
| E3 | variants | `ascaled anobo` (vs `acvx0`) |
| E4 | where the bound helps | `astats` (counters `cvx_*` in the JSON) |
| E5 | SCIP, same domains | section 7, outputs under `results/local/` (`OUT=` for SCIP, `OUT=` for run_ibex.sh) |

Order: E0, then E1, E2, E3 on the 223 set (about 25 min per arm), then the same arms
on the 75 set (about 2 h per arm), then E4, then E5. Report each in
`results/local/README.md`: machine (CPU, cores, OS, Ipopt version), command lines,
the tables printed by `compare.py`, and anything unexpected. Commit and push to
`ml-branching` after each experiment (results only; ask before changing `src/`).
