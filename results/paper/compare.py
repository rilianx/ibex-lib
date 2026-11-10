"""Compare configurations ("arms") on the 298 COCONUT instances: the table used in
the paper's ablation (results/ablation/README.md).

Usage:
  python3 results/paper/compare.py --ref abase anoip abase acvx0 [more arms]
  python3 results/paper/compare.py --ref abase3 abase3 ami --set closed   # 223 set only

Reads results/ipopt-diam/runs.jsonl (223 instances IbexOpt closes) and
results/anytime/runs.jsonl (75 it does not), one JSON line per (instance, arm).

Columns:
  closed     instances the arm closes (status "complete") over the chosen set
  n (closes) common set: instances some listed arm closes and some listed arm
             takes > 1 s; in parentheses, how many of them this arm closes
  time       geometric mean of t(arm)/t(ref) on the common set, times floored at
             1 s, a non-closed run counts as the time limit (310 s); bootstrap 95 %
             interval (2000 resamples, seed 1)
  W/L        instances where the arm is > 5 % faster / slower than ref
  nodes      geometric mean of (nodes+1) ratios on the instances every arm closes
  t/node     geometric mean of per-instance time-per-node ratios (>= 100 nodes in both)
  ipopt%, bound%   share of the total time (instances every arm closes) spent in
             Ipopt + its correction search, and in the Lagrangian bound

Only compare arms run on the same machine without other load (see EXPERIMENTS.md).
"""
import argparse, collections, json, math, random

ap = argparse.ArgumentParser()
ap.add_argument('arms', nargs='+')
ap.add_argument('--ref', required=True)
ap.add_argument('--set', choices=['all', 'closed', 'open'], default='all')
ap.add_argument('--tlim', type=float, default=310.0)
ap.add_argument('--floor', type=float, default=1.0)
args = ap.parse_args()

files = {'closed': ['results/ipopt-diam/runs.jsonl'], 'open': ['results/anytime/runs.jsonl'],
         'all': ['results/ipopt-diam/runs.jsonl', 'results/anytime/runs.jsonl']}[args.set]
arms = list(dict.fromkeys([args.ref] + args.arms))
R = collections.defaultdict(dict)
for f in files:
    for l in open(f):
        d = json.loads(l)
        if d['arm'] in arms:
            R[d['instance']][d['arm']] = d['result']

def ok(r): return r.get('status') == 'complete'
def t(r):
    s = r.get('status')
    if s == 'complete': return max(min(r['time'], args.tlim), args.floor)
    if s in ('timeout', 'no_lower_bound'): return args.tlim
    if s == 'unbounded' and (r.get('loup') or 0) > -1e300: return args.tlim   # old runs
    return None   # infeasible, proved unbounded, killed: not comparable

V = [v for v in R.values() if all(a in v for a in arms)]
I = [v for v in V if None not in [t(v[a]) for a in arms]
     and any(ok(v[a]) for a in arms) and max(t(v[a]) for a in arms) > args.floor]
C = [v for v in V if all(ok(v[a]) for a in arms)]
print('instances with every arm: %d, common set: %d, closed by every arm: %d' % (len(V), len(I), len(C)))
print('| arm | closed | n (closes) | time vs %s [IC95] | W/L | nodes | t/node | ipopt%% | bound%% |' % args.ref)
print('|---|---|---|---|---|---|---|---|---|')
gm = lambda xs: math.exp(sum(xs) / len(xs)) if xs else float('nan')
for a in arms:
    L = [math.log(t(v[a]) / t(v[args.ref])) for v in I]
    random.seed(1)
    bs = sorted(sum(random.choice(L) for _ in L) / len(L) for _ in range(2000))
    W = sum(t(v[a]) < 0.95 * t(v[args.ref]) for v in I)
    Lo = sum(t(v[args.ref]) < 0.95 * t(v[a]) for v in I)
    N = gm([math.log((v[a]['nodes'] + 1) / (v[args.ref]['nodes'] + 1)) for v in C])
    P = gm([math.log((v[a]['time'] / v[a]['nodes']) / (v[args.ref]['time'] / v[args.ref]['nodes']))
            for v in C if v[a]['nodes'] >= 100 and v[args.ref]['nodes'] >= 100])
    T = sum(v[a]['time'] for v in C) or 1
    ip = sum(v[a].get('ipopt_time', 0) + v[a].get('ipopt_corr_time', 0) for v in C)
    cv = sum(v[a].get('cvx_time', 0) for v in C)
    closed = sum(1 for v in V if ok(v[a]))
    print('| %s | %d | %d (%d) | %.2f [%.2f, %.2f] | %d/%d | %.2f | %.2f | %.1f | %.1f |' % (
        a, closed, len(I), sum(ok(v[a]) for v in I), gm(L), math.exp(bs[50]), math.exp(bs[1950]),
        W, Lo, N, P, 100 * ip / T, 100 * cv / T))
