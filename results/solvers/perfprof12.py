"""Performance profile SCIP vs IbexOpt (afull) on the same domains: missing variable
bounds set to +-1e12 for both (scip_inf12.jsonl, ibex_inf12.jsonl). Instances SCIP can
read, closed by one of the two, > 1 s for one of the two; times floored at 1 s."""
import json, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT = sys.argv[1] if len(sys.argv) > 1 else 'results/solvers/perfprof12.png'
def load(f):
    D = {}
    for l in open(f):
        d = json.loads(l); D[d['instance']] = d['result']
    return D
S = load('results/solvers/scip_inf12.jsonl'); A = load('results/solvers/ibex_inf12.jsonl')
ok = lambda r: r.get('status') == 'complete'
I = [i for i in A if i in S and S[i].get('status') != 'error' and A[i].get('status') != 'killed'
     and (ok(A[i]) or ok(S[i])) and max(S[i].get('time', 310), A[i].get('time', 310)) > 1]
n = len(I)
t = lambda r: max(r['time'], 1.0) if ok(r) else float('inf')
ts = {i: t(S[i]) for i in I}; ta = {i: t(A[i]) for i in I}
BLUE, ORANGE, INK, MUTED, GRID = '#2a78d6', '#eb6834', '#1f1f1e', '#6b6a64', '#e6e5e0'
TAUMAX = 1024
fig, ax = plt.subplots(figsize=(5.2, 3.8))
for name, T, col, ls, short in (('SCIP 10', ts, ORANGE, '--', 'SCIP'),
                                ('IbexOpt + Lagrangian bound', ta, BLUE, '-', 'IbexOpt')):
    r = sorted(T[i] / min(ts[i], ta[i]) for i in I)
    xs, ys = [1], [sum(1 for v in r if v <= 1) / n]
    for v in r:
        if v == float('inf') or v <= 1: continue
        xs.append(v); ys.append(sum(1 for w in r if w <= v) / n)
    xs.append(TAUMAX); ys.append(ys[-1])
    ax.step(xs, ys, where='post', color=col, lw=2, ls=ls, label=name)
    ax.annotate('%s  %.0f%%' % (short, 100 * ys[-1]), (TAUMAX, ys[-1]),
                xytext=(-4, 4 if col == BLUE else -12), textcoords='offset points',
                ha='right', fontsize=8.5, color=INK)
ax.set_xscale('log', base=2); ax.set_xlim(1, TAUMAX); ax.set_ylim(0, 1.02)
ax.set_title('SCIP vs IbexOpt, same domains (missing bounds = ±1e12)\n(n = %d: closed by one of the two, > 1 s for one)' % n,
             fontsize=9, color=INK)
ax.set_xlabel('τ (time ratio to the faster of the two)', fontsize=8.5, color=MUTED)
ax.set_ylabel('fraction of instances within τ of the best', fontsize=8.5, color=MUTED)
ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
for s in ('left', 'bottom'): ax.spines[s].set_color(MUTED)
ax.tick_params(colors=MUTED, labelsize=8)
ax.legend(loc='lower center', fontsize=8, frameon=False)
fig.tight_layout(); fig.savefig(OUT, dpi=200); fig.savefig(OUT.replace('.png', '.pdf'))
print('wrote', OUT, 'n =', n)
