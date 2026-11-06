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
