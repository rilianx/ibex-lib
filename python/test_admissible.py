#!/usr/bin/env python3
"""Safety net: a model that carries no information must still terminate.

A learned rule scores the candidates of a node, but it must not be able to pick
something the hand-written bisector would have refused, and a total tie must
fall back on what that bisector would have chosen. The cheapest way to check
both at once is to run the search with a model whose score is constant: every
candidate ties, so the choice is entirely decided by the admissible set and by
its ordering.

If this test regresses, a model is choosing outside the guards -- which is the
difference between a rule that is merely bad and one that never terminates. On
`ex8_1_2.bch` the failure mode was the objective variable: chosen at the root,
it turns 14 nodes into 250 000.

    python3 python/test_admissible.py
    python3 python/test_admissible.py --dir benchs/optim/medium --max-nodes 5000
"""

import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ibexml import solve, export_linear, MODEL_FEATURES   # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default="benchs/optim/easy", help="directory of .bch files")
    p.add_argument("--max-nodes", type=int, default=1000,
                   help="only check instances the default rule closes below this")
    p.add_argument("--base-timeout", type=float, default=10.0)
    p.add_argument("--model-timeout", type=float, default=20.0)
    p.add_argument("--model", default=None,
                   help="model file to use (default: a freshly written constant model)")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    what = "a constant-score model" if args.model is None else args.model
    model = args.model
    if model is None:
        import tempfile
        fd, model = tempfile.mkstemp(suffix=".model", prefix="ibexml-constant-")
        os.close(fd)
        export_linear(model, [0.0] * len(MODEL_FEATURES), bias=0.0)

    files = sorted(glob.glob(os.path.join(args.dir, "*.bch")))
    if not files:
        raise SystemExit("no .bch under %s" % args.dir)

    reference, failures, checked = {}, [], 0

    for f in files:
        try:
            r = solve(f, timeout=args.base_timeout)
        except Exception:
            continue
        if r["status"] != "complete" or r["nodes"] >= args.max_nodes:
            continue
        reference[f] = r

    for f, r in sorted(reference.items()):
        m = solve(f, model=model, timeout=args.model_timeout)
        checked += 1
        ok = m["status"] == "complete"
        if not ok:
            failures.append((f, m["status"], m["nodes"]))
        if args.verbose or not ok:
            print("  %-40s default %6d   model %8d  %s"
                  % (os.path.basename(f), r["nodes"], m["nodes"],
                     "ok" if ok else "FAILED (%s)" % m["status"]))

    print("checked %d instances the default rule closes in < %d nodes, with %s"
          % (checked, args.max_nodes, what))
    if failures:
        print("FAILED: %d instance(s) did not terminate with %s:" % (len(failures), what))
        for f, st, n in failures:
            print("   %-40s %s after %d nodes" % (os.path.basename(f), st, n))
        return 1
    print("PASSED: %s terminates on all of them" % what)
    if args.model is None:
        os.unlink(model)
    return 0


if __name__ == "__main__":
    sys.exit(main())
