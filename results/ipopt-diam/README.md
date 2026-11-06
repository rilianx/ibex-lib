# Ipopt's width threshold (`--ipopt-diam`)

IbexOpt calls Ipopt only on boxes whose widest variable is at most 1e8
(`ipopt_diam`, src/ipopt/ibex_LoupFinderIpopt.cpp). Variables with
"infinite" bounds written as +-1e8 keep many boxes at width 2e8, so Ipopt
never runs there (results/node-selection/README.md).

`run.sh` on the 223 instances lsmear-lffix closes, 310 s here (~600 s on the
user's machine), `runs.jsonl`. Pairwise time metric (lower is better):

    1e12 vs 1e8    0.976 vs 0.990   diff -0.014  IC95 [-0.030, -0.001]
    none vs 1e8    0.978 vs 0.984   diff -0.006  IC95 [-0.027,  0.017]
    1e12 vs none   0.985 vs 0.993   diff -0.008  IC95 [-0.027,  0.007]

1e12: all 223 close; faster (>1.2x) on 7, slower on 4 (geo time 0.97).
Big gains where the 1e8 threshold blocked Ipopt: ex8_2_4 21.5 -> 2.5 s,
eigencco 4.5 -> 0.6, launch 4.9 -> 2.0, dixchlng 2.3 -> 0.3. Worst losses:
makela3 1.7 -> 2.6 s, ex6_2_12 1.0 -> 1.4, haldmads 1.6 -> 2.2, ex14_1_7
15.8 -> 19.7.

No threshold: 222 close (matrix2 times out); slower on 6, among them robot
79.7 -> 203 s and makela3 1.7 -> 7.5 s: Ipopt on truly huge boxes costs.

1e12 is a small, significant, nearly free gain: Ipopt on boxes with +-1e8
"infinite" bounds, not on the really unbounded ones.

## When does Ipopt pay off? (`IBEX_IPOPT_LOG=<file>`: one line per Ipopt call)

lsmear-lffix, --ipopt-diam 1e12, on the instances with node-selection
headroom. For each: Ipopt calls in the whole search, the first call whose
value is the optimum (within 1e-3), and the first call whose point is
accepted as the new loup:

    instance    nodes  calls  Ipopt returns the optimum      accepted
    launch         24     3   call 0, depth 0                never
    hs116         118    29   call 0, depth 0                never
    ex8_2_4        42     4   call 41, depth 21              never
    concon        610     1   call 300, depth 295            never
    eigmaxc       304     1   call 156, depth 151            never
    eigencco       12     2   call 10, depth 6               same call
    ex7_3_4        34     2   call 10, depth 5               same call
    dixchlng      110     4   call 50, depth 10              same call
    ex8_1_6bis    176     2   call 89, depth 88              same call
    polak5       3504     1   call 2847, depth 51            same call
    polak4       1770     1   call 1300, depth 883           same call
    polak2       1920     1   call 1914, depth 959           same call
    makela3      2070     2   call 1953, depth 1035          same call
    ship-1        680    28   never                          never

Two separate failures:
1. Ipopt is called almost never (1-4 times in thousands of nodes): the
   schedule fires at loup-finder calls 0, 10, 20, 50 and every 100, but the
   call is skipped if that node's box is wider than ipopt_diam, and it is not
   retried at the next eligible node. With unbounded variables (polak,
   makela3) the boxes at the scheduled calls are wide; Ipopt only runs when
   forced (probing found a point below ymax), and then it finds the optimum
   at once.
2. Ipopt finds the optimum but its point is not accepted (launch and hs116
   at the root): the point must be certified, which IbexOpt attempts with a
   correction search around it (1 s of CPU); that fails, and the loup comes
   much later from elsewhere.

## Why the point is rejected, and the fixes (arms dnone, x1e5, c1e5, nc1e5)

The IBEX_IPOPT_LOG line now also gives the goal at Ipopt's point, the inner
test, the max violation in the normalized system and the correction search's
status, nodes, time and loup.

- launch: the point violates the equalities by 1e-4 (Ipopt's default
  constr_viol_tol); the correction search in the box of +-1e-6 around it times
  out (1 s) without a certified point. Tightening constr_viol_tol does not
  change the returned point.
- hs116: violation 1e-8, but the +-1e-6 box is proved infeasible at the root.
- ex8_2_4, concon: the width threshold (fixed by dnone / 1e12). ex8_2_4 also
  paid 1 s per call certifying Ipopt *failures* (status -13, value 0): the
  correction now runs only on a usable Ipopt result (21 s -> 3.4 s already at
  1e8).
- eigmaxc: Ipopt only runs at the end; nothing to fix.

Fixes:
- `--ipopt-diam=-1`: no threshold; Ipopt starts from the point of the box
  closest to 0 when the box is unbounded.
- `--ipopt-expansion=1e-5`: the correction box is +-1e-5. launch 98 -> 2
  nodes, hs116 118 -> 2. But a wider box certifies a slightly different point
  even when the narrow one works (haldmads: 2.5927330 vs 2.5927322 at the first
  call), the path diverges and haldmads goes 1.5 -> 2.7 s (seeds 1, 2, 4).
- `--ipopt-expansion-first=0.1`: try the 1e-6 box first (0.1 s) and widen
  only if it certifies nothing. haldmads and dipigri as in the base, launch
  and hs116 at 2 nodes.

Results, 223 instances, 310 s. Time metric: t/max(t_A,t_B), 1 s floor
(less is better). Nodes: geometric mean of (nodes_A+1)/(nodes_B+1) over the
instances both close, and instances with > 20% fewer / more nodes.

    A      vs B       time diff [IC95]          nodes gmean  better  worse
    c1e5   vs d1e8    -0.006 [-0.020, +0.006]   0.967         4       2
    nc1e5  vs d1e8    -0.007 [-0.025, +0.012]   0.740        19       6
    nc1e5  vs d1e12   +0.007 [-0.006, +0.021]   0.815        10       3
    nc1e5  vs dnone   +0.014 [+0.004, +0.023]   0.966         5       3
    dnone  vs d1e8    -0.018 [-0.002, +0.037]   0.771        17       5

The time metric is dominated by noise here: most instances close in 1-5 s and
the arms ran on different days (c1e5/nc1e5 after a container restart). The
nc1e5 losses to dnone are instances with the same node count and ~20% more
time (ex8_5_6, hs099, srcpm, ex6_2_12, himmel16, ex8_5_2); run back to back
they are within 0-12%. Node counts are deterministic and say: no threshold
is the big effect (0.77), the cascade adds launch and hs116 on top (0.97).

Per instance (time s / nodes):

    instance   d1e8          d1e12        dnone        c1e5         nc1e5
    launch     4.93/104      1.96/24      3.29/54      0.12/2       0.08/2
    hs116      0.56/118      0.59/118     0.64/118     0.01/2       0.01/2
    ex8_2_4    21.51/1818    2.50/42      2.73/46      3.38/1872    2.72/42
    eigencco   4.53/1188     0.65/12      0.68/12      4.43/1188    0.67/12
    dixchlng   2.28/2038     0.27/110     0.35/110     2.65/2462    0.35/112
    polak5     1.42/3504     1.32/3504    0.00/4       1.22/3504    0.00/4
    makela3    1.73/2078     2.60/2070    0.00/0       2.00/2078    0.00/0
    concon     0.63/610      0.52/610     0.02/10      0.51/610     0.02/10
    dnieper    5.17/434      5.37/434     5.36/434     7.64/624     6.96/624
    matrix2    130.3/4.57M   139.3/4.57M  timeout      143.6/4.57M  killed

Remaining problems:
- matrix2 without threshold: Ipopt fails (status -13, invalid number at the
  starting point) on every wide box; no longer costs a correction search, but
  the run still does not close in 310 s. A perturbed starting point is the
  next thing to try.
- dnieper is path-sensitive (it was 196 to >7000 nodes depending on the
  initial loup); the cascade moves it 434 -> 624.

## matrix2 and the final arm

matrix2 without threshold had two problems:
1. Boxes "at infinity" (x in [9e307, 1.8e308]): the starting point closest
   to 0 is itself ~1e308, the evaluations overflow, Ipopt fails at once.
2. Forced calls start from the point another loup finder just found; in
   matrix2 those are at ~1e20, and there Ipopt hangs inside one line search
   (max_iter and max_cpu_time are only checked between iterations): the run
   ignored its 310 s timeout.

Fix: with `--ipopt-diam` negative, Ipopt is not called when its starting
point (forced or not) has a coordinate beyond `--ipopt-start-max` (1e8).
matrix2 closes in 125 s (base 130 s). `--ipopt-max-time` (Ipopt's
max_cpu_time) also exists but does not help with that hang.

Arm `final`: `--ipopt-diam=-1 --ipopt-expansion=1e-5
--ipopt-expansion-first=0.1` with that binary. Closes 223/223 (total
1632 s vs 1651 s for the base).

    final vs   time diff [IC95]          nodes gmean  better  worse
    d1e8       -0.014 [-0.031, -0.000]   0.763        19       6
    d1e12      -0.000 [-0.010, +0.008]   0.838        10       3
    dnone      +0.006 [-0.008, +0.018]   0.989         5       3

Worse than the base by > 25% in nodes: ex14_1_7 (3990 -> 5200 nodes,
15.8 -> 20.1 s), dnieper (434 -> 624), haldmads (398 -> 982), and three
instances under 0.05 s (ackley_5, ex8_1_4, ex8_1_5).

Recommended configuration: `--ipopt-diam=-1 --ipopt-expansion=1e-5
--ipopt-expansion-first=0.1`.

## After the affine relaxation fix (arms base8, fixaff)

`--relax both` used LinearizerAffine2, which could add invalid cuts (an
Affine2 assignment kept the previous form's coefficients when the new one
was inactive; fixed in 9e69f77). base8 = d1e8 and fixaff = final, rerun
with the fixed binary:

    base8  vs d1e8    time +0.001 [-0.010, +0.010]  nodes gmean 0.999
    fixaff vs final   time +0.010 [+0.001, +0.017]  nodes gmean 1.006
    fixaff vs base8   time -0.006 [-0.021, +0.008]  nodes gmean 0.768
                      (19 instances with > 20% fewer nodes, 6 with > 25% more)

The invalid cuts barely mattered on these 223 runs (haldmads is the one
instance whose node count moves: 982 -> more with the fix); every optimum
is the same. The conclusions above hold.
