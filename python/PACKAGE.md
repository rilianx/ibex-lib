# ibexopt-ml

Node-level instrumentation of IbexOpt, for learning a branching (bisection)
rule from Python.

For every node of the search it gives you the features a model can score
variables with, and — for each candidate variable — the size of the truncated
branch & bound ("dive") that bisecting on it would trigger. It also exposes
contraction and bisection as a line protocol, so the search itself can be
driven from Python with a learned rule in place of LSmear.

```
bin/ibexopt-ml            the executable, ready to run (Linux x86-64)
python/ibexml.py          client, a Python branch & bound loop, tensor encoding
python/collect_dataset.py dataset collection over a set of instances
python/test_admissible.py safety net: a no-information model must terminate
python/README.md          the full documentation
benchs/                   Minibex benchmarks: easy/, medium/, hard/, ...
```

Nothing to compile, nothing to install. The binary links only against
`libstdc++`, `libm`, `libgcc_s` and `libc`.

## Check it runs

```bash
./bin/ibexopt-ml --version
printf '{"cmd":"info"}\n{"cmd":"quit"}\n' | ./bin/ibexopt-ml benchs/easy/ex2_1_1.bch
```

## Collect a dataset

```bash
# the executable drives the search itself (fastest)
./bin/ibexopt-ml benchs/easy/ex2_1_1.bch --collect -o data.jsonl \
    --budget 200 --sample-prob 0.2 --max-samples 20 --progress

# run to completion under LSmear, or under a learned rule, and time it
./bin/ibexopt-ml benchs/easy/ex2_1_1.bch --solve
./bin/ibexopt-ml benchs/easy/ex2_1_1.bch --solve --model my.model

# or over a set of instances, from Python
python3 python/collect_dataset.py benchs/easy/*.bch \
    -o data/train.jsonl --budget 200 --sample-prob 0.2 --max-samples 100
```

Each line is one training sample: the node's features, plus — for every
candidate variable — the dive that bisecting on it triggers.

## Use it from Python

```python
import sys; sys.path.insert(0, "python")
from ibexml import IbexOptML, run, solve, encode, labels, pairs

with IbexOptML("benchs/easy/ex2_1_1.bch") as srv:
    root = srv.reset()["root"]
    s = srv.sample(root, budget=200)
    Xv, Xc, E   = encode(s["node"])       # variables, constraints, Jacobian edges
    y, mask, cz = labels(s, Xv.shape[0])  # dive size, candidate mask, censoring
    print(s["node"]["bisect_var"], y[mask], pairs(s))
```

## Put your own rule in the loop

`run()` lets the solver drive the real search and asks you for every branching
variable, so node counts stay comparable with `ibexopt`:

```python
def my_rule(node, candidates, state):
    return int(model(node).argmax())      # or None to let LSmear decide

run("benchs/easy/ex2_1_1.bch", my_rule, sample_prob=0.1, on_sample=save)
```

And once the model is a linear model or a tree ensemble, export it and let the
solver score it itself — the only way to measure wall-clock time:

```python
from ibexml import export_sklearn, solve
export_sklearn("my.model", gbr)
print(solve("benchs/easy/ex2_1_1.bch"))                   # LSmear
print(solve("benchs/easy/ex2_1_1.bch", model="my.model")) # the learned rule
```

`ibexml.py` locates the executable on its own; set `$IBEXOPT_ML` if you move it.
numpy is needed only for `encode()` and `labels()`.

## Read this next

**[python/README.md](python/README.md)** — concepts, the complete data
dictionary, the protocol reference, the Python API, a training path, DAgger,
measured costs, guarantees and known limitations.
