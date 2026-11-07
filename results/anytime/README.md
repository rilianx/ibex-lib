# Anytime behaviour on the instances that do not close

`open-instances.txt`: the 75 instances of benchs/optim/all that lsmear-lffix
does not close in 310 s (all minus results/node-selection/closed-instances.txt).
`run.sh` runs them with `--bound-trace` (loup and uplo in time); `analyze.py`
computes, over [0, 310 s]:

- final relative gap,
- primal integral / T (Berthold): primal gap |loup(t)-best|/max(|loup|,|best|),
  1 with no loup, best = best loup of any arm,
- primal-dual integral / T: min(1, (loup-uplo)/max(|loup|,|uplo|,1)),
- time before the root is processed counts as gap 1.

Arms (binary with the affine fix 9e69f77 and the root's time counted ae66555;
the older arms `ipoptfix`/`base` predate both):

    tbase   lsmear-lffix, both, ipoptxn (Ipopt threshold 1e8)
    tfinal  + --ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1
    tms20   tfinal + --ipopt-multistart=20 (20 extra Ipopt runs from random
            points of the root box)

Statuses (each arm): 59-60 timeouts, 8-9 unbounded, 2 infeasible, 3 complete
(chenery-1, hs093, wall), 2 killed (hs090, hs091: Ipopt's symbolic Hessian
blows up at construction).

Means over the 62 instances with metrics in both arms (lower is better):

    A       vs B       primal integral          primal-dual integral     final loup better / worse
    tfinal  vs tbase   0.138 vs 0.160 (14/5)    0.599 vs 0.601 (6/4)     orthrege, prolog / cresc4, polak6
    tms20   vs tfinal  0.106 vs 0.138 (6/29)    0.597 vs 0.599 (6/18)    cresc4, linear, orthrege, polak6, rk23 / hs089, polygon25
    tms20   vs tbase   0.106 vs 0.160 (12/27)   0.597 vs 0.601 (6/20)

(n/m: instances where A is better / worse by more than 1e-3.)

Instances with no feasible point at 310 s:

    tbase   discs ex5_3_3 hs089 hs107 hs118 minlphi minlphi1 orthregb rk23
    tfinal  discs ex5_3_3 hs118 minlphi minlphi1
    tms20   ex5_3_3 hs089 hs118 minlphi minlphi1

Reading:
- The upper bound settles early: with tfinal, in 44 of 50 instances with a
  finite gap the loup at 31 s is the final one (1e-6). The primal-dual
  integral barely moves between arms: what remains is the lower bound's gap.
- The new Ipopt base finds feasible points the old one does not (hs089,
  hs107, orthregb, rk23) and lowers the primal integral.
- Random multistart has a few large wins (discs: a feasible point; rk23,
  orthregb, polak6) and many small losses (the root's extra Ipopt runs delay
  the first loup by seconds; hs089 loses its feasible point, polygon25 ends
  worse).

## Guided Ipopt restarts (tg10, tg10l, tg10p)

tg10 = tfinal + `--ipopt-guided=0.1` (extra Ipopt runs at nodes, < 10% of
the CPU time, start = farthest of 10 random points of the node's box from
previous Ipopt solutions); tg10l = Ipopt bounded by the node's box; tg10p =
tg10 + `--ipopt-guided-patience=20`. Primal integrals use the best loup over
all arms (so the tfinal values here differ slightly from above).

    tg10p vs   primal integral          primal-dual integral     final loup better / worse
    tfinal     0.111 vs 0.143 (10/4)    0.569 vs 0.599 (12/3)    cresc4 linear orthrege polak6 polygon25 rk23 / prolog wall
    tg10       0.111 vs 0.120 (6/1)     0.569 vs 0.586 (17/1)    polak6 prolog / linear polygon25
    tg10l      0.111 vs 0.128 (6/6)     0.569 vs 0.589 (15/2)
    tbase      0.111 vs 0.165 (17/6)    0.569 vs 0.601 (12/5)    cresc4 linear orthrege polak6 polygon25 / wall

(n/m: instances better / worse by more than 1e-3.) Still no feasible point:
discs, ex5_3_3, hs118, minlphi, minlphi1.

Patience removes the cost of the budget on the lower bound (the primal-dual
integral is now better than tfinal's, not worse) and keeps the primal
gains. Recommended configuration:

    --ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1
    --ipopt-guided=0.1 --ipopt-guided-patience=20

## hs118: a modelling artifact, fixed (ef82025)

hs118 writes its range constraints as thick equalities, e.g.
`((-x1)+x4)-[-7, 6])=0` (i.e. -7 <= x4-x1 <= 6). NormalizedSystem kept the
interval constant inside f, so the inner test needed x4-x1 <= -7 and
x4-x1 >= 6 at once: no point could ever be certified (no loup in any arm,
while the lower bound reached the known optimum 664.82). Ipopt received
mid(f)=0, i.e. x4-x1 = -0.5, and reported infeasibility on all 914 calls.
With thick equalities normalized as ranges (a-eps <= f0 <= b+eps) and Ipopt
bounds g in [-rad, rad], hs118 closes at the root: loup 664.82045, 0 nodes,
0.014 s (tg10p config). Only hs118 and launch__coconut have interval
constants in their constraints; launch__coconut is proved infeasible at the
root by upstream ibexopt as well (0 cells), unchanged.

The sweep arms (g05p20, g20p20, g10p10, g10p50) still use the binary
without this fix.

## Patience 10 vs 20 (g10p10 vs tg10p)

    g10p10 vs  primal integral          primal-dual integral     final loup better / worse
    tg10p      0.126 vs 0.111 (3/6)     0.575 vs 0.569 (8/7)     - / polak6 prolog rk23
    tfinal     0.126 vs 0.143 (11/3)    0.575 vs 0.599 (13/5)    cresc4 linear orthrege polak6 polygon25 rk23 / prolog wall
    tbase      0.126 vs 0.165 (17/4)    0.575 vs 0.601 (14/4)

g10p10 closes 4 of the 75 (mistake in 281 s, besides chenery-1, hs093,
wall). The patience-20 advantage comes from rk23 (0.083 vs 0.444) and
polak6, where both arms make a single / a few guided runs: the budget is
measured in CPU time, so which node gets the run differs between runs, and
these instances are noise-sensitive. Patience 10 is the better choice on the
223 closed instances (a tie in time with tfinal, no dnieper loss) and about
the same here.

## Sweep of the guided restarts (budget rho, patience p), 75 open instances

Means over the 62 instances with metrics in every arm (best loup over all
arms as reference, so values differ slightly from the tables above):

    arm      rho   p    primal int.  primal-dual int.  closed  no loup   vs tfinal: primal / primal-dual better-worse
    tfinal   -     -    0.156        0.599             3       5
    tg10p    0.1   20   0.124        0.569             3       5
    g10p10   0.1   10   0.139        0.575             4       5         11/3   13/5
    g05p20   0.05  20   0.125        0.587             4       4         12/2   19/1
    g20p20   0.2   20   0.133        0.587             4       5         14/1   15/3
    g10p50   0.1   50   0.140        0.578             4       5          7/5   20/1

Every setting improves on tfinal in both integrals, and they are close to
one another (single-instance noise: rk23, polak6). With the 223 closed
instances (time vs tfinal: g05p20 +0.005, g10p10 -0.000, tg10p +0.012,
g10p50 +0.012, g20p20 +0.014), rho 0.05 / p 20 and rho 0.1 / p 10 are the
two best compromises; rho 0.05 also found a feasible point on one more
instance.
