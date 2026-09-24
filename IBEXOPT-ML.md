# ibexopt-ml

This tree is [ibex-lib](README.md) plus `ibexopt-ml`, a node-level
instrumentation of IbexOpt used to learn a branching (bisection) rule from
Python.

It adds, for every node of the search, the features a model can score variables
with, and — for each candidate variable — the size of the truncated branch &
bound ("dive") that bisecting on it would trigger. It also exposes contraction
and bisection as a line protocol, so the search itself can be driven from
Python with a learned rule in place of LSmear.

## Build

Same as plain Ibex. **An LP solver is mandatory** — `ibexopt-ml` needs it for
the same reason `ibexopt` does, and `LP_LIB` defaults to `none`, so it has to
be asked for explicitly. SoPlex is bundled and needs no download:

```bash
mkdir -p build && cd build
cmake .. -DLP_LIB=soplex          # -DLP_LIB=clp also works (bundled too)
make -j8
```

Leaving `-DLP_LIB` out builds cleanly but produces a binary that refuses to run
with `requires a LP Solver (use -DLP_LIB with cmake)`.

This produces `build/bin/ibexopt-ml` next to the usual `ibexopt` and
`ibexsolve`. The interval library defaults to Gaol, also bundled; nothing is
downloaded.

## Try it

```bash
# one JSON training sample per line
./build/bin/ibexopt-ml benchs/optim/easy/ex2_1_1.bch --collect -o data.jsonl \
    --budget 200 --sample-prob 0.2 --max-samples 20 --progress

# or talk to it as a server
printf '{"cmd":"info"}\n{"cmd":"quit"}\n' \
    | ./build/bin/ibexopt-ml benchs/optim/easy/ex2_1_1.bch

# from Python
python3 python/collect_dataset.py benchs/optim/easy/*.bch \
    -o data/train.jsonl --budget 200 --sample-prob 0.2 --max-samples 100
```

## Documentation

**[python/README.md](python/README.md)** — the full document: concepts, data
dictionary, protocol reference, Python API, a training path, DAgger, measured
costs, internals and known limitations.

## What was added

| | |
|---|---|
| `src/ml/ibex_MLNodeServer.{h,cpp}` | the instrumentation, on top of `Optimizer` |
| `src/ml/ibex_BscHijackGuard.{h,cpp}` | `--bisector lsmear-guard`: LSmear until one variable hijacks it, then RoundRobin |
| `src/ml/ibex_MLModel.{h,cpp}` | plain-text linear / GBDT models, scored inside the solver |
| `src/ml/ibex_Json.{h,cpp}` | dependency-free JSON reader / streaming writer |
| `src/bin/ibexopt-ml.cpp` | CLI, the four modes, command dispatch |
| `python/ibexml.py` | client, search drivers, tensor encoding, model export |
| `python/collect_dataset.py` | dataset collection over a set of instances |
| `python/make_package.py` | builds the standalone zip (binary + Python + benchmarks) |

Four pre-existing files were modified, all additively and all documented in
place; no public API changed:

| file | change | why |
|---|---|---|
| `src/optim/ibex_Optimizer.h` | search state `private` → `protected` | a subclass must be able to save and restore it around a dive |
| `src/tools/ibex_Random.h/.cpp` | `RNG::get_state()` / `set_state()` | `srand(s)` repositions the stream by *drawing* s numbers, so resuming an arbitrary position costs O(s); this saves and restores it in constant time |
| `src/contractor/ibex_CtcAcid.h` | `get_tuning()` / `set_tuning()` | ACID adapts across calls, so a speculative dive changes what the enclosing search does next unless its tuning is put back |
| `src/strategy/ibex_Sts.h` | `calls()` | read a call counter without parsing a report |

`src/CMakeLists.txt` and `src/bin/CMakeLists.txt` were extended to build the new
files.

Benchmarks are in `benchs/optim/` (`easy/`, `medium/`, `hard/`, …), in Minibex
format.
