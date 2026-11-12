"""Pairwise performance profiles (Dolan & More): SCIP vs each IbexOpt arm, one pair
per panel so that a third solver never changes a pair's ratios.
Instances: translated for SCIP, closed by at least one of the two, and taking
more than 1 s for at least one of the two (otherwise both ratios are 1); times
floored at 1 s; an instance a solver does not close has ratio = oo.
Usage: python3 results/solvers/perfprof.py [out.png]"""
import json, sys, collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else 'results/local/perfprof.png'
S = {}
for l in open('results/local/scip.jsonl'):
    d = json.loads(l); S[d['instance']] = d['result']
X = collections.defaultdict(dict)
for f in ('results/local/runs-closed.jsonl', 'results/local/runs-open.jsonl'):
    for l in open(f):
        d = json.loads(l)
        X[d['arm']][d['instance']] = d['result']

FAIL = set(sys.argv[2:])  # SCIP errors counted as 'not solved' instead of dropped
def solved(r, scip=False):
    if r.get('status') != 'complete': return False
    if scip and r.get('loup') is not None and r['loup'] < -1e19: return False  # SCIP's infinity
    return True

PAIRS = [('anoip', 'IbexOpt, no Ipopt'), ('abase', 'IbexOpt + Ipopt'),
         ('acvx0', 'IbexOpt + Ipopt + Lagrangian bound')]
BLUE, ORANGE, INK, MUTED, GRID = '#2a78d6', '#eb6834', '#1f1f1e', '#6b6a64', '#e6e5e0'
TAUMAX = 1024
fig, axes = plt.subplots(1, len(PAIRS), figsize=(4.2 * len(PAIRS), 3.6), sharey=True)
for ax, (arm, label) in zip(axes, PAIRS):
    I = [i for i in S if (S[i].get('status') != 'error' or i in FAIL) and i in X[arm]
         and (solved(S[i], True) or solved(X[arm][i]))
         and max(S[i].get('time', 1200), X[arm][i].get('time', 1200)) > 1]
    n = len(I)
    def t(r, ok): return max(r['time'], 1.0) if ok else float('inf')
    ts = {i: t(S[i], solved(S[i], True)) for i in I}
    ta = {i: t(X[arm][i], solved(X[arm][i])) for i in I}
    for name, T, col, ls in (('SCIP 10', ts, ORANGE, '--'), (label, ta, BLUE, '-')):
        r = sorted(T[i] / min(ts[i], ta[i]) for i in I)
        xs, ys = [1], [sum(1 for v in r if v <= 1) / n]
        for v in r:
            if v == float('inf') or v <= 1: continue
            xs.append(v); ys.append(sum(1 for w in r if w <= v) / n)
        xs.append(TAUMAX); ys.append(ys[-1])
        ax.step(xs, ys, where='post', color=col, lw=2, ls=ls, label=name)
        ax.annotate('%s  %.0f%%' % ('SCIP' if name.startswith('SCIP') else 'IbexOpt', 100 * ys[-1]),
                    (TAUMAX, ys[-1]), xytext=(-4, 4 if col == BLUE else -12),
                    textcoords='offset points', ha='right', fontsize=8.5, color=INK)
    ax.set_xscale('log', base=2); ax.set_xlim(1, TAUMAX); ax.set_ylim(0, 1.02)
    ax.set_title('SCIP vs %s\n(n = %d: closed by one of the two, > 1 s for one)' % (label, n),
                 fontsize=9, color=INK)
    ax.set_xlabel('τ (time ratio to the faster of the two)', fontsize=8.5, color=MUTED)
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'): ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.legend(loc='lower center', fontsize=8, frameon=False)
axes[0].set_ylabel('fraction of instances within τ of the best', fontsize=8.5, color=MUTED)
fig.tight_layout()
fig.savefig(OUT, dpi=200)
fig.savefig(OUT.replace('.png', '.pdf'))
print('wrote', OUT)
