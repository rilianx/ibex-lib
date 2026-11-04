#!/usr/bin/env python3
"""Why does the final upper bound appear so late? Follow the path of boxes
from the root to the optimum x* (bisecting with the search's own bisector,
contracting as the search does) and, at each level:
  - does the loup finder (Ipopt at every call) reach the optimum in this box?
  - the box's size (max width over the original variables);
  - the lower bound of the child containing x* vs its sibling's (what a
    best-first selection would compare).

    python3 python/ub_path.py instance [instance ...]
x* and the optimum come from a solve with IBEX_LOUP_POINT=1.
"""
import json, math, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ibexml  # noqa
ROOT = os.path.join(HERE, "..")
BIN = os.environ.get("BIN", os.path.join(ROOT, "build-fix", "bin", "ibexopt-ml"))
COMMON = ["--relax", "both", "--bisector", "lsmear-lffix", "--acid-fixed", "1000000"]


def solve(bch):
    env = dict(os.environ, IBEX_LOUP_POINT="1")
    out = subprocess.run([BIN, bch, "--solve", "--loup", "ipoptxn", "--random-seed", "1",
                          "--timeout", "300"] + COMMON, capture_output=True, text=True, env=env).stdout
    return json.loads(out.strip().splitlines()[-1])


def contains(box, x, gv, tol=1e-6):
    k = 0
    for j, (a, b) in enumerate(box):
        if j == gv:
            continue
        v = x[k]; k += 1
        t = tol * (1 + abs(v))
        if v < a - t or v > b + t:
            return False
    return True


def width(box, gv):
    return max((b - a) for j, (a, b) in enumerate(box) if j != gv)


def main():
    for inst in sys.argv[1:]:
        bch = os.path.join(ROOT, "benchs", "optim", "all", inst + ".bch")
        r = solve(bch)
        if "loup_point" not in r:
            print(inst, "sin punto óptimo:", r.get("status")); continue
        opt = r["loup"]; xs = [p[0] for p in r["loup_point"]]
        tol_f = 1e-3 * abs(opt) + 1e-7
        print("\n== %s: óptimo %.6g, %d nodos (%s); loup final en el nodo %s" % (
            inst, opt, r["nodes"], r["status"], r.get("loup_final_nodes")))
        print("  prof  ancho      loup del buscador     encuentra  lb(hijo con x*)  lb(hermano)  hermano")
        with ibexml.IbexOptML(bch, binary=BIN, loup="ipoptxn", random_seed=1, ipopt_freq=1,
                              extra_args=COMMON) as srv:
            gv = srv.goal_var
            box = srv.reset(loup=ibexml.INF)["root"]
            first = None
            for depth in range(400):
                if not contains(box, xs, gv):
                    print("  %4d  x* fuera de la caja contraída: fin" % depth); break
                srv.set_loup(ibexml.INF)
                c = srv.contract(box, keep_loup=False)
                lo = c["loup_after"]
                found = lo is not None and lo <= opt + tol_f
                if found and first is None:
                    first = depth
                try:
                    var, left, right = srv.bisect(box)
                except ibexml.IbexError:
                    print("  %4d  no se puede bisectar: fin" % depth); break
                kids = []
                for b in (left, right):
                    srv.set_loup(ibexml.INF)
                    k = srv.contract(b, keep_loup=False)
                    kids.append((k["status"], k["box"]))
                inside = [i for i in (0, 1) if kids[i][0] == "open" and contains(kids[i][1], xs, gv)]
                if not inside:
                    print("  %4d  %.3g  ningún hijo abierto contiene x*: fin" % (depth, width(box, gv))); break
                i = inside[0]; s = 1 - i
                lbx = kids[i][1][gv][0]
                lbs = kids[s][1][gv][0] if kids[s][0] == "open" else None
                if depth < 12 or found or depth % 10 == 0:
                    print("  %4d  %-9.3g  %-20s  %-9s  %-15.6g  %-11s  %s" % (
                        depth, width(box, gv), ("%.6g" % lo) if lo is not None and lo < 1e300 else "-",
                        "SI" if found else "no", lbx, ("%.6g" % lbs) if lbs is not None else "-", kids[s][0]))
                if found:
                    break
                box = kids[i][1]
            print("  -> el buscador encuentra el óptimo a profundidad %s del camino" % first)


if __name__ == "__main__":
    main()
