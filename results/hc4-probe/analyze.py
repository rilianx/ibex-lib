#!/usr/bin/env python3
"""Pairwise comparison of the arms in results/hc4-probe/runs.jsonl.

Per pair of arms, over the instances both ran and at least one closed:
t(A)/max(t(A),t(B)) with a 1 s floor and 2 for a run that does not close
(lower is better), its mean with a 95% bootstrap interval of the difference,
and the instances each closes.

    python3 results/hc4-probe/analyze.py [runs.jsonl]
"""
import json, os, random, sys
from itertools import combinations
f = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "runs.jsonl")
R = {}
for l in open(f):
    r = json.loads(l); R.setdefault(r["arm"], {})[r["instance"]] = r["result"]
ok = lambda x: x.get("status") == "complete"
t = lambda x: max(x.get("time", 0), 1.0)
for arm, d in sorted(R.items()):
    print("%-8s %3d corridas, cierra %3d" % (arm, len(d), sum(ok(x) for x in d.values())))
for a, b in combinations(sorted(R), 2):
    common = [i for i in R[a] if i in R[b] and (ok(R[a][i]) or ok(R[b][i]))]
    if not common:
        continue
    sa, sb = [], []
    for i in common:
        x, y = R[a][i], R[b][i]
        tx, ty = t(x), t(y)
        sa.append(tx / max(tx, ty) if ok(x) else 2.0)
        sb.append(ty / max(tx, ty) if ok(y) else 2.0)
    d = [p - q for p, q in zip(sa, sb)]
    random.seed(0)
    bs = sorted(sum(random.choice(d) for _ in d) / len(d) for _ in range(2000))
    print("%s vs %s: n=%d  %s %.3f  %s %.3f  diff %.3f IC95 [%.3f, %.3f]  (menor es mejor)" % (
        a, b, len(common), a, sum(sa) / len(sa), b, sum(sb) / len(sb),
        sum(d) / len(d), bs[50], bs[1949]))
