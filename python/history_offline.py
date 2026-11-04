#!/usr/bin/env python3
"""Offline: does what bisecting variable j did at earlier nodes predict what it
does now?  Pseudocost / reliability-branching rules simulated on the samples
of a dataset, in search order.

Per sample, the strong-branching key of each candidate (pruned children, then
least open volume), relative to the sample's best: q_j = key_j - max key.
History rules (only previous samples of the same instance are used):
  pc-free  pseudocost from the nodes where j was actually bisected (the
           base's choice): the observation a search gets for free.
  pc-all   pseudocost from every probed candidate of previous samples (what
           reliability branching gathers while it still probes).
  rel-k    reliability branching: probe (true key) the candidates seen fewer
           than k times, pseudocost (pc-all history) for the others.
Unseen variables and ties fall back to the base's own order (LSmear rank).
Metric: per instance, sum of the chosen candidates' dive nodes over the
base's choice (geo over instances), instances above 1, and the share of
candidates probed.

    DATASET=dataset-fixedloup python3 python/history_offline.py
"""
import math, os, sys
from collections import defaultdict
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from oracle_sbfeat_offline import load, geo  # noqa


def main():
    XA, XP, XS, X2, Y, R = load()
    R = R.reset_index(drop=True)
    R["key"] = 100.0 * XP[:, 0] - XP[:, 3]          # pruned, then -log open volume
    R["var"] = -1
    R["q"] = R.key - R.groupby("sample").key.transform("max")
    R["lsrank"] = R.groupby("sample").cumcount()   # candidates come ranked by LSmear
    R["y"] = R.r
    # candidate variable index: recover from load()'s order -- labels are in
    # candidate order; we need the var id, so reload it alongside
    import glob, json
    DATA = os.path.join(HERE, "..", "results", os.environ.get("DATASET", "dataset-fixedloup"))
    vars_ = []
    for f in sorted(glob.glob(os.path.join(DATA, "samples", "*.jsonl"))):
        for line in open(f):
            s = json.loads(line)
            L = [l for l in s["labels"] if l["valid"]]
            done = [l["nodes"] for l in L if not l["censored"]]
            if len(L) < 2 or not done:
                continue
            vars_ += [l["var"] for l in L]
    assert len(vars_) == len(R), (len(vars_), len(R))
    R["var"] = vars_

    def run(rule, k=0):
        rows = []; probed = 0; total = 0
        for inst, G in R.groupby("instance", sort=False):
            free = defaultdict(list); full = defaultdict(list)
            for sid, S in G.groupby("sample", sort=True):
                base = S[S.guard]
                if base.empty:
                    for _, x in S.iterrows(): full[x["var"]].append(x.q)
                    continue
                b = base.iloc[0]
                score = {}
                for _, x in S.iterrows():
                    total += 1
                    h = free if rule == "pc-free" else full
                    if rule == "rel" and len(full[x["var"]]) < k:
                        score[x.name] = x.q; probed += 1
                    elif h[x["var"]]:
                        score[x.name] = float(np.mean(h[x["var"]]))
                    else:
                        score[x.name] = None
                known = {i: v for i, v in score.items() if v is not None}
                if rule == "sb":
                    ch = S.q.idxmax(); probed += len(S)
                elif known:
                    best = max(known.values())
                    # ties and unseen: the base's order (LSmear rank)
                    cand = [i for i, v in known.items() if v >= best - 1e-9]
                    ch = min(cand, key=lambda i: R.lsrank[i])
                    # never prefer a seen variable over the base unless strictly better
                    if b.name in known and known[b.name] >= best - 1e-9:
                        ch = b.name
                else:
                    ch = b.name
                rows.append((inst, R.y[ch], b.y))
                # observations for later samples
                for _, x in S.iterrows():
                    full[x["var"]].append(x.q)
                free[b["var"]].append(b.q)
        T = pd.DataFrame(rows, columns=["inst", "y", "yb"])
        ratio = T.groupby("inst").apply(lambda g: g.y.sum() / g.yb.sum())
        return geo(ratio), (ratio > 1.0001).mean(), probed / max(total, 1)

    print("%d muestras, %d instancias" % (R["sample"].nunique(), R.instance.nunique()))
    for name, rule, k in [("SB (sondea todo, referencia)", "sb", 0), ("pseudocosto gratis (solo lo bisectado)", "pc-free", 0),
                          ("pseudocosto con historia completa", "pc-all", 0),
                          ("reliability k=1", "rel", 1), ("reliability k=2", "rel", 2), ("reliability k=4", "rel", 4)]:
        g, w, p = run(rule, k)
        print("  %-40s suma/base %.3f   instancias >1 %3.0f%%   candidatas sondeadas %3.0f%%" % (name, g, 100 * w, 100 * p))


if __name__ == "__main__":
    main()
