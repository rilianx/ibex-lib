#!/usr/bin/env python3
"""Collect a branching dataset from one or more Minibex files.

Each line of the output is one sample: the features of a visited node, plus --
for every admissible candidate variable -- the size of the dive that bisecting
on it triggers.

    {"format": 2, "node": {...}, "budget": 200, "labels": [...], "instance": "..."}

Two drivers, and the difference matters once a model is in the loop:

  --driver cpp   the search runs inside ibexopt-ml. Fastest, and the trajectory
                 is exactly the one ibexopt follows. With --model, the model
                 branches, so the data is on-policy for it.

  --driver run   the search still runs inside ibexopt-ml -- same node selection,
                 same contraction -- but every branching decision comes back
                 here, so you can plug a Python rule with --rule. This is the
                 DAgger loop: your model drives, the dives still tell the truth.

Examples
--------
    python3 python/collect_dataset.py benchs/optim/easy/*.bch \\
        -o data/train.jsonl --budget 200 --sample-prob 0.2 --max-samples 200

    # round 2: your own model branches
    python3 python/collect_dataset.py benchs/optim/easy/*.bch \\
        -o data/dagger.jsonl --driver run --rule mymodel:choose --sample-prob 0.2

    # or a model the solver scores itself
    python3 python/collect_dataset.py benchs/optim/easy/*.bch \\
        -o data/dagger.jsonl --model my.model --sample-prob 0.2
"""

import argparse
import importlib
import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ibexml import run, sample_format, DEFAULT_BINARY, SAMPLE_FORMAT   # noqa: E402


def load_rule(spec):
    """Import a branching rule given as ``module:function`` or ``module.py:function``."""
    mod, _, fun = spec.partition(":")
    if not fun:
        raise SystemExit("--rule must look like module:function")
    if mod.endswith(".py"):
        sys.path.insert(0, os.path.dirname(os.path.abspath(mod)) or ".")
        mod = os.path.basename(mod)[:-3]
    try:
        m = importlib.import_module(mod)
    except ImportError as e:
        raise SystemExit("cannot import %r: %s" % (mod, e))
    try:
        return getattr(m, fun)
    except AttributeError:
        raise SystemExit("%r has no %r" % (mod, fun))


def dive_args(args):
    """The flags that control diving and sampling, shared by both drivers."""
    argv = ["--budget", str(args.budget),
            "--budget-start", str(args.budget_start),
            "--max-depth", str(args.max_depth),
            "--topk", str(args.topk),
            "--random-seed", str(args.random_seed)]
    if args.no_prune:
        argv += ["--no-prune"]
    if args.no_goal:
        argv += ["--no-goal"]
    return argv


def emit(out, rec, bch, counter):
    if sample_format(rec) != SAMPLE_FORMAT:
        raise SystemExit("ibexopt-ml wrote sample format %d, this script expects %d "
                         "-- the binary and python/ are out of step"
                         % (sample_format(rec), SAMPLE_FORMAT))
    rec["instance"] = os.path.basename(bch)
    out.write(json.dumps(rec) + "\n")
    counter[0] += 1


def collect_cpp(bch, out, args, counter):
    """Let ibexopt-ml run the search and write the samples itself."""
    argv = [args.binary, bch, "--collect",
            "--sample-prob", str(args.sample_prob),
            "--max-samples", str(args.max_samples)] + dive_args(args)
    if args.model:     argv += ["--model", args.model]
    if args.max_nodes: argv += ["--max-nodes", str(args.max_nodes)]
    if args.timeout:   argv += ["--timeout", str(args.timeout)]

    proc = subprocess.Popen(argv, stdout=subprocess.PIPE, text=True, bufsize=1)
    for line in proc.stdout:
        line = line.strip()
        if line:
            emit(out, json.loads(line), bch, counter)
    proc.wait()


def collect_run(bch, out, args, counter, rule):
    """Drive the real search from here, sampling as we go."""
    import random
    rng = random.Random(int(args.random_seed))
    quota = [0]

    def should_sample(state):
        if args.max_samples and quota[0] >= args.max_samples:
            return False       # stop sampling, but let the search finish
        if rng.random() >= args.sample_prob:
            return False
        quota[0] += 1
        return True

    def on_sample(rec):
        emit(out, rec, bch, counter)
        out.flush()

    extra = ["--model", args.model] if args.model else []
    res = run(bch, rule,
              binary=args.binary,
              sample_when=should_sample,
              on_sample=on_sample,
              max_nodes=args.max_nodes,
              timeout=args.timeout,
              features=not args.no_features,
              include_goal=not args.no_goal,
              topk=args.topk,
              budget=args.budget,
              budget_start=args.budget_start,
              prune=not args.no_prune,
              max_depth=args.max_depth,
              random_seed=args.random_seed,
              extra_args=extra)
    if args.verbose and res:
        print("    search: %s" % {k: res[k] for k in ("status", "nodes", "loup", "uplo")},
              file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bch", nargs="+", help="Minibex files")
    p.add_argument("-o", "--output", required=True, help="output .jsonl file")
    p.add_argument("--driver", choices=("cpp", "run"), default="cpp")
    p.add_argument("--binary", default=DEFAULT_BINARY)

    g = p.add_argument_group("branching rule")
    g.add_argument("--model", help="a model file the solver scores itself")
    g.add_argument("--rule", help="a Python rule, as module:function "
                                  "(implies --driver run)")

    g = p.add_argument_group("diving")
    g.add_argument("--budget", type=int, default=200,
                   help="ceiling on the dive node budget (default: 200)")
    g.add_argument("--budget-start", type=int, default=25,
                   help="first budget probed, doubling until the dive closes; "
                        "0 goes straight to --budget (default: 25)")
    g.add_argument("--no-prune", action="store_true",
                   help="give every candidate the full budget instead of the "
                        "best result so far")
    g.add_argument("--max-depth", type=int, default=0, help="dive depth limit (0: none)")
    g.add_argument("--topk", type=int, default=0,
                   help="evaluate only the k most promising candidates (0: all). "
                        "The cost of a sample is linear in this")
    g.add_argument("--no-goal", action="store_true",
                   help="do not treat the objective variable as a candidate")

    g = p.add_argument_group("sampling and search")
    g.add_argument("--sample-prob", type=float, default=0.2)
    g.add_argument("--max-samples", type=int, default=200,
                   help="per instance; sampling stops, the search does not")
    g.add_argument("--max-nodes", type=int, default=0, help="search node limit (0: none)")
    g.add_argument("--timeout", type=float, default=0, help="search time limit, seconds")
    g.add_argument("--no-features", action="store_true",
                   help="with --driver run, do not ship the node features to the rule")
    g.add_argument("--random-seed", type=float, default=1.0)
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    rule = None
    if args.rule:
        rule = load_rule(args.rule)
        args.driver = "run"
    elif args.driver == "run":
        rule = lambda node, cand, state: None      # the solver's own bisector

    d = os.path.dirname(os.path.abspath(args.output))
    if d:
        os.makedirs(d, exist_ok=True)

    total = 0
    with open(args.output, "w") as out:
        for bch in args.bch:
            t0 = time.time()
            counter = [0]
            try:
                if args.driver == "cpp":
                    collect_cpp(bch, out, args, counter)
                else:
                    collect_run(bch, out, args, counter, rule)
            except SystemExit:
                raise
            except Exception as e:      # a bad instance must not stop the run
                print("  %-40s FAILED: %s" % (os.path.basename(bch), e), file=sys.stderr)
                continue
            total += counter[0]
            print("  %-40s %4d samples  (%.1fs)"
                  % (os.path.basename(bch), counter[0], time.time() - t0), file=sys.stderr)

    print("%d samples -> %s" % (total, args.output), file=sys.stderr)


if __name__ == "__main__":
    main()
