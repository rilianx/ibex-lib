"""Solve a Minibex (.bch) optimization benchmark with SCIP (PySCIPOpt) and print one
JSON line with status, primal/dual bounds, time and nodes.

Usage: python3 scip_run.py file.bch [time_limit] [rel_gap] [abs_gap]

Translation: constants are substituted, variables get their bounds (oo -> infinite),
the objective goes into an epigraph variable (SCIP needs a linear objective), Minibex
^ -> **, ln -> log, abs -> fabs. Same stopping rule as IbexOpt by default
(relative gap 1e-3, absolute gap 1e-7). SCIP is not rigorous (floating-point
relaxations, feasibility tolerance 1e-6), so a 'complete' status is not a proof in the
interval sense.
"""
import json, math, re, sys, time

import pyscipopt as ps

path = sys.argv[1]
tlim = float(sys.argv[2]) if len(sys.argv) > 2 else 310.0
relgap = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
absgap = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-7

src = open(path).read()
src = re.sub(r'//[^\n]*|#[^\n]*', '', src)
low = src.lower()


def section(name, nxt):
    m = re.search(r'\b%s\b' % name, low)
    if not m:
        return ''
    end = len(src)
    for n in nxt:
        k = re.search(r'\b%s\b' % n, low[m.end():])
        if k:
            end = min(end, m.end() + k.start())
    return src[m.end():end]


consts_s = section('constants', ['variables'])
vars_s = section('variables', ['minimize', 'minimise', 'constraints'])
obj_s = section('minimize', ['constraints', 'end']) or section('minimise', ['constraints', 'end'])
ctr_s = section('constraints', ['end'])

model = ps.Model()
model.hideOutput()
NS = {'exp': ps.exp, 'log': ps.log, 'sqrt': ps.sqrt, 'sin': ps.sin, 'cos': ps.cos,
      'fabs': abs, 'oo': math.inf, 'pi': math.pi}


THICK = [0]


def thick(m):
    # a thick constant [a,b] means "some value in [a,b]": a fresh bounded variable
    THICK[0] += 1
    name = '__thick%d' % THICK[0]
    NS[name] = model.addVar(name, lb=float(ev(m.group(1))), ub=float(ev(m.group(2))))
    return name


def py(e):
    e = ' '.join(e.split())
    e = re.sub(r'\[([^\[\],]*),([^\[\],]*)\]', thick, e)
    e = e.replace('^', '**')
    e = re.sub(r'\bln\s*\(', 'log(', e)
    e = re.sub(r'\babs\s*\(', 'fabs(', e)
    return e


def ev(e):
    return eval(py(e), {'__builtins__': {}}, NS)


for st in consts_s.split(';'):
    if '=' in st:
        k, v = st.split('=', 1)
        NS[k.strip()] = float(ev(v.strip()))

INF = model.infinity()
for st in vars_s.split(';'):
    st = st.strip()
    if not st:
        continue
    m = re.match(r'(\w+)\s*(?:in\s*\[(.*),(.*)\])?$', st, re.S)
    if not m:
        raise SystemExit('cannot parse variable: %r' % st)
    lb = -INF if m.group(2) is None else float(ev(m.group(2)))
    ub = INF if m.group(3) is None else float(ev(m.group(3)))
    lb = -INF if lb == -math.inf else lb
    ub = INF if ub == math.inf else ub
    NS[m.group(1)] = model.addVar(m.group(1), lb=lb, ub=ub)

obj = obj_s.strip().rstrip(';')
if obj:
    z = model.addVar('__obj', lb=-INF, ub=INF)
    model.addCons(ev(obj) - z <= 0)
    model.setObjective(z, 'minimize')

for st in ctr_s.split(';'):
    st = st.strip()
    if not st:
        continue
    m = re.match(r'(.*?)(<=|>=|=)(.*)$', st, re.S)
    if not m:
        raise SystemExit('cannot parse constraint: %r' % st)
    lhs, op, rhs = ev(m.group(1)), m.group(2), ev(m.group(3))
    d = lhs - rhs
    model.addCons(d <= 0 if op == '<=' else d >= 0 if op == '>=' else d == 0)

model.setParam('limits/time', tlim)
model.setParam('limits/gap', relgap)
model.setParam('limits/absgap', absgap)
model.setParam('parallel/maxnthreads', 1)
model.setParam('randomization/randomseedshift', 0)
t0 = time.time()
model.optimize()
st = model.getStatus()
res = {
    'status': {'optimal': 'complete', 'gaplimit': 'complete', 'timelimit': 'timeout',
               'infeasible': 'infeasible', 'unbounded': 'unbounded',
               'inforunbd': 'unbounded'}.get(st, st),
    'scip_status': st,
    'loup': model.getPrimalbound() if model.getNSols() > 0 else None,
    'uplo': model.getDualbound(),
    'nodes': model.getNNodes(),
    'time': model.getSolvingTime(),
    'wall': time.time() - t0,
}
print(json.dumps(res))
