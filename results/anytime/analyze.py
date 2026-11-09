"""Anytime metrics on results/anytime/runs.jsonl (runs with --bound-trace).

Per (instance, arm), over [0, T] (T = the run's time limit):
  gap      final relative gap (loup-uplo)/max(1,|loup|)
  pint     primal integral / T (Berthold): primal gap
           |loup(t)-best|/max(|loup(t)|,|best|) (1 with no loup, capped at 1),
           best = the best loup of any arm
  pdint    primal-dual integral / T: min(1,(loup-uplo)/max(|loup|,|uplo|,1))
           (1 when either bound is infinite)
  t_ub     time of the last loup change; ub_early: the loup at T/10 is the
           final one within 1e-6 (relative): then the remaining gap is the
           lower bound's
Usage: python3 results/anytime/analyze.py [runs.jsonl] [T]
"""
import json, re, sys, math, collections

F = sys.argv[1] if len(sys.argv) > 1 else 'results/anytime/runs.jsonl'
T = float(sys.argv[2]) if len(sys.argv) > 2 else 310.0
INF = float('inf')

R = collections.defaultdict(dict)
for l in open(F):
    d = json.loads(l)
    r = d['result']
    B = [[float(x) if x is not None else INF for x in b] for b in (r.get('bounds') or [])]
    st = r.get('status')
    # before 'no_lower_bound' existed, both stops were reported as 'unbounded':
    # only a loup of -oo (a feasible point with objective -oo) is evidence
    if st == 'unbounded' and r.get('loup') is not None and r['loup'] > -1e300: st = 'no_lower_bound'
    R[d['instance']][d['arm']] = (st, B, r.get('nodes'))

def pg(lo, best):
    if lo == INF or lo is None: return 1.0
    den = max(abs(lo), abs(best))
    return 0.0 if den == 0 else min(1.0, abs(lo - best) / den)

def pdg(lo, up):
    if lo == INF or up == -INF: return 1.0
    return min(1.0, max(0.0, lo - up) / max(abs(lo), abs(up), 1.0))

def integral(B, f):
    # before the root is processed there are no bounds: gap 1
    s = min(B[0][0], T) if B else T
    for k in range(len(B)):
        t0 = B[k][0]; t1 = B[k+1][0] if k+1 < len(B) else T
        t1 = min(t1, T)
        if t1 > t0: s += f(B[k]) * (t1 - t0)
    return s / T

def metrics(B, best, status):
    if not B: return None
    if status == 'complete':   # closed: gap 0 from the end on
        B = B + [[B[-1][0], B[-1][1], B[-1][2], B[-1][2]]]
    lo, up = B[-1][2], B[-1][3]
    gap = (lo - up) / max(1, abs(lo)) if lo < INF and up > -INF else INF
    lastl = 0.0
    for k in range(1, len(B)):
        if B[k][2] != B[k-1][2]: lastl = B[k][0]
    pre = [b for b in B if b[0] <= T/10]
    at10 = pre[-1][2] if pre else INF
    early = lo < INF and at10 < INF and abs(at10 - lo) <= 1e-6 * max(1, abs(lo))
    return dict(gap=gap, pint=integral(B, lambda b: pg(b[2], best)),
                pdint=integral(B, lambda b: pdg(b[2], b[3])), t_ub=lastl,
                ub_early=early, loup=lo, uplo=up)

arms = sorted({a for v in R.values() for a in v})
M = collections.defaultdict(dict)
for i, v in R.items():
    los = [B[-1][2] for (s, B, n) in v.values() if B and s not in ('unbounded',)]
    best = min(los) if los else INF
    for a, (s, B, n) in v.items():
        m = metrics(B, best, s) if s not in ('unbounded', 'infeasible') else None
        if m: m['status'] = s; M[i][a] = m

for a in arms:
    S = [M[i][a] for i in M if a in M[i]]
    st = collections.Counter(R[i][a][0] for i in R if a in R[i])
    fin = [m for m in S if m['gap'] < INF]
    print("%-9s runs %d %s" % (a, sum(st.values()), dict(st)))
    if not S: continue
    print("  mean primal integral %.3f   mean primal-dual integral %.3f" % (
        sum(m['pint'] for m in S)/len(S), sum(m['pdint'] for m in S)/len(S)))
    print("  finite final gap: %d, ub settled by T/10: %d" % (len(fin), sum(m['ub_early'] for m in fin)))

if len(arms) >= 2:
    a, b = arms[0], arms[1]
    print("\nper instance (%s vs %s): final loup, uplo, pdint" % (a, b))
    for i in sorted(M):
        if a in M[i] and b in M[i]:
            x, y = M[i][a], M[i][b]
            flag = ''
            if x['loup'] < y['loup'] - 1e-6*max(1, abs(y['loup'])): flag += ' loup<'
            if x['loup'] > y['loup'] + 1e-6*max(1, abs(y['loup'])): flag += ' loup>'
            if x['uplo'] > y['uplo'] + 1e-6*max(1, abs(y['uplo'])): flag += ' uplo>'
            if x['uplo'] < y['uplo'] - 1e-6*max(1, abs(y['uplo'])): flag += ' uplo<'
            print("  %-16s %-8s %13.6g %13.6g %.3f | %-8s %13.6g %13.6g %.3f%s" % (
                i.replace('.bch', ''), x['status'][:8], x['loup'], x['uplo'], x['pdint'],
                y['status'][:8], y['loup'], y['uplo'], y['pdint'], flag))
