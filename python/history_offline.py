#!/usr/bin/env python3
"""Offline: does what bisecting variable j did at earlier nodes predict what it
does now?  Pseudocost / reliability-branching rules simulated on the samples
of a dataset, in search order.

Per sample, the strong-branching key of each candidate (pruned children, then
least open volume), relative to the sample's best: q_j = key_j - max key.
History rules (only previous samples of the same instance are used):
  pc       pseudocost: the mean q of the variable's earlier observations.
  state    state-dependent pseudocost: a ridge regression fitted online, on
           this instance's earlier observations, from the candidate's state
           (LSmear rank and relative score, smear_sum_rel and rel_diam
           relative) and its history (mean q, number of observations).
  k        reliability: probe (true key) the candidates whose variable has
           fewer than k observations.
Observations are the probed candidates and the bisected one, as in a real
search; "ideal" counts every candidate of every earlier sample (upper bound).
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

    # state features per candidate (comparative, from the 30 model features):
    # lsmear rank / relmax, smear_sum_rel relmax, rel_diam relmax
    ST = XA[:, [15, 16, 22, 28]]

    def run(rule, k=0, ideal=False):
        """rule: sb | pc | state.  k: reliability threshold (probe a candidate
        whose variable has fewer than k observations).  ideal: every
        candidate of every earlier sample counts as observed (an upper bound:
        a real search only observes what it probes or bisects)."""
        rows = []; probed = 0; total = 0
        for inst, G in R.groupby("instance", sort=False):
            hist = defaultdict(list)          # var -> observed q
            Xtr, ytr = [], []                 # (state + history) -> q
            for sid, S in G.groupby("sample", sort=True):
                base = S[S.guard]
                if base.empty:
                    continue
                b = base.iloc[0]
                def feats(x):
                    h = hist[x["var"]]
                    return list(ST[x.name]) + [np.mean(h) if h else 0.0, math.log1p(len(h))]
                model = None
                if rule == "state" and len(ytr) >= 20:
                    A = np.array(Xtr); A1 = np.hstack([A, np.ones((len(A), 1))])
                    w = np.linalg.lstsq(A1.T @ A1 + 1e-3 * np.eye(A1.shape[1]), A1.T @ np.array(ytr), rcond=None)[0]
                    model = w
                score = {}; obs = []
                for _, x in S.iterrows():
                    total += 1
                    if rule == "sb" or len(hist[x["var"]]) < k:
                        score[x.name] = x.q; probed += 1; obs.append(x)
                    elif rule == "pc" and hist[x["var"]]:
                        score[x.name] = float(np.mean(hist[x["var"]]))
                    elif rule == "state" and model is not None:
                        score[x.name] = float(np.dot(model, feats(x) + [1.0]))
                    else:
                        score[x.name] = None
                known = {i: v for i, v in score.items() if v is not None}
                ch = b.name
                if known:
                    best = max(known.values())
                    cand = [i for i, v in known.items() if v >= best - 1e-9]
                    ch = min(cand, key=lambda i: R.lsrank[i])
                    if b.name in known and known[b.name] >= best - 1e-9:
                        ch = b.name
                rows.append((inst, R.y[ch], b.y))
                # what the search observes: the probed candidates and the bisected one
                seen = list(S.itertuples()) if ideal else obs + ([R.loc[ch]] if all(o.name != ch for o in obs) else [])
                for x in seen:
                    xx = R.loc[x.Index] if hasattr(x, "Index") else x
                    Xtr.append(feats(xx)); ytr.append(xx.q)
                    hist[xx["var"]].append(xx.q)
        T = pd.DataFrame(rows, columns=["inst", "y", "yb"])
        ratio = T.groupby("inst").apply(lambda g: g.y.sum() / g.yb.sum())
        return geo(ratio), (ratio > 1.0001).mean(), probed / max(total, 1)

    print("%d muestras, %d instancias" % (R["sample"].nunique(), R.instance.nunique()))
    for name, rule, k, ideal in [
            ("SB (sondea todo, referencia)", "sb", 0, False),
            ("pseudocosto, historia ideal", "pc", 0, True),
            ("estado, historia ideal", "state", 0, True),
            ("pseudocosto + reliability k=2", "pc", 2, False),
            ("pseudocosto + reliability k=4", "pc", 4, False),
            ("estado + reliability k=1", "state", 1, False),
            ("estado + reliability k=2", "state", 2, False),
            ("estado + reliability k=4", "state", 4, False)]:
        g, w, p = run(rule, k, ideal)
        print("  %-36s suma/base %.3f   instancias >1 %3.0f%%   candidatas sondeadas %3.0f%%" % (name, g, 100 * w, 100 * p))

if __name__ == "__main__":
    main()
