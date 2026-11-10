# SCIP on the 298 benchmark instances

`scip_run.py` translates a Minibex file to PySCIPOpt (SCIP 10.0, pip wheel) and
solves it with the same stopping rule as IbexOpt (relative gap 1e-3, absolute
1e-7), 310 s, one thread; `run.sh` runs the 298 instances (223 closed + 75 open
by IbexOpt). Objective in an epigraph variable; ^ -> **, ln -> log; a thick
constant [a,b] becomes a variable bounded in [a,b].

SCIP is **not rigorous** (floating-point relaxations, feasibility tolerance 1e-6):
its "solved" is not a proof in the interval sense. In the SCIP 8 paper
(Bestuzheva et al., arXiv 2301.00587), SCIP returned a wrong optimal value or an
infeasible solution on 39 of 1000 MINLPLib instances, BARON on 27.

## Status

SCIP: 221 complete, 64 timeouts, 4 infeasible, 1 unbounded, 8 not translated
(hs088-hs091: nesting too deep for Python's parser; cresc4, modele_opti: acos,
tan; minlphi1: a constraint without variables; brownal). One of the 221,
ex8_4inf-1, stops at SCIP's infinity (-1e20): the problem is unbounded, so 220
real closes.

Closed /298: SCIP 220, IbexOpt without Ipopt (anoip) 223, abase 228, afull 236.

- Only SCIP closes (15): 3pk-borne, 3pk, chenery, ex8_4_3, linear, minlphi,
  optmass, orthregb, pollut, portfl1, prolog, qp1-1, qp1, qp2, test_infinity3
  (mostly quadratic / bilinear models, where its LP relaxations are strong;
  3pk and test_infinity3 are among the IbexOpt "no_lower_bound" stops).
- Only afull closes (30): aircraftb, avion2, dixchlng, eigencco, ex6_1_1,
  ex6_1_3, ex6_1_3bis, ex6_2_9 to ex6_2_14, ex7_2_1, ex7_2_1bis, ex7_2_3,
  ex8_4_5bis, ex8_5_6, hs099, hs103, m4wd, pentagon, rosenbrock10, schwefel5,
  schwefel5-abs, srcpm, catena, hs093, polak3, wall (unbounded or huge
  variable domains, log-heavy GlobalLib models, trigonometric models).
- Different optimum: bearing (SCIP 1.92512, IbexOpt 1.95173; 12 equalities,
  SCIP accepts 1e-6 violations).

## Times (common set)

Instances some solver closes and some solver takes > 1 s: n_total = 110 (in
parentheses: closed by the solver). Ratio t/t(SCIP), times floored at 1 s,
timeout = 310 s; relative = t/max(t, t_SCIP) for the arm / for SCIP.

| solver | n_total (closes) | geomean t/t(SCIP) [IC95] | median | faster/slower than SCIP | relative arm / SCIP |
|---|---|---|---|---|---|
| SCIP 10 | 110 (81) | 1.00 | 1.00 | — | — |
| IbexOpt, no Ipopt (anoip) | 110 (83) | 1.56 [0.77, 3.08] | 1.35 | 43/62 | 0.682 / 0.578 |
| IbexOpt + Ipopt (abase) | 110 (87) | 1.10 [0.54, 2.13] | 1.00 | 47/61 | 0.653 / 0.633 |
| IbexOpt + Ipopt + Lagrangian bound (afull) | 110 (95) | 0.70 [0.37, 1.31] | 1.00 | 55/53 | 0.636 / 0.704 |

Without a local solver the rigorous interval solver is slower than SCIP; with
Ipopt as primal oracle it ties; with the Lagrangian bound from Ipopt's
multipliers it closes 16 more instances than SCIP on the 298 and is faster on
the common set (not significant: each solver wins big on different model
classes).

## IbexOpt "unbounded"

IbexOpt stopped with "unbounded" on 9 open instances. Only ex8_4inf-1 and
test_infinity4 have a feasible point with objective -oo (SCIP agrees). The
other 7 stop because a box that can no longer be bisected (at infinity) has an
objective enclosure unbounded below, which is an overestimation, not a proof:
3pk is a sum of squares written as products e*e (SCIP: optimum 1.72),
test_infinity3 (SCIP: -1e-8), like (SCIP finds 1005.5 < IbexOpt's 1218.5),
like-1, ex8_5_6-1, ship, modele_opti. ibexopt-ml now reports these as
`no_lower_bound` (not closed), and analyze.py reclassifies older runs.

## Performance profiles (pairwise)

`perfprof.py` -> `perfprof.png` / `perfprof.pdf`: one panel per pair (SCIP vs
one IbexOpt arm), so a third solver never changes a pair's ratios. Instances
translated for SCIP, closed by one of the two and taking > 1 s for one of the
two; times floored at 1 s. Fraction within τ of the faster of the two:

| pair | n | fastest (τ=1): SCIP / IbexOpt | closed at τ=1024: SCIP / IbexOpt |
|---|---|---|---|
| SCIP vs no Ipopt | 102 | 59 % / 41 % | 75 % / 77 % |
| SCIP vs + Ipopt | 100 | 54 % / 46 % | 73 % / 80 % |
| SCIP vs + Ipopt + Lagrangian bound | 94 | 46 % / 54 % | 69 % / 85 % |

## Same domains for both: missing bounds = ±1e12

As in the SCIP 8 paper, every missing variable bound replaced by ±1e12 for both
solvers (`SCIP_INFBOUND=1e12`, `run_ibex.sh` on sed-bounded copies of the
Minibex files); IbexOpt = afull. 291 instances (7 SCIP cannot read).

| | closed |
|---|---|
| IbexOpt afull ±1e12 | 240 of the 291 (241 of 298; 236 on the original domains: + brownal, ex8_4inf-1, rk23, test_infinity3, test_infinity4) |
| SCIP 10 ±1e12 | 223 (220: + ex8_4inf-1, polak3, test_infinity4) |

- Only IbexOpt (31): aircraftb, avion2, brownal, catena, dixchlng, eigencco,
  ex6_1_1, ex6_1_3, ex6_1_3bis, ex6_2_9 to ex6_2_14, ex7_2_1, ex7_2_1bis,
  ex7_2_3, ex8_4_5bis, ex8_5_6, hs093, hs099, hs103, m4wd, pentagon, rk23,
  rosenbrock10, schwefel5, schwefel5-abs, srcpm, wall.
- Only SCIP (14): 3pk, 3pk-borne, chenery, ex8_4_3, linear, minlphi, optmass,
  orthregb, pollut, portfl1, prolog, qp1, qp1-1, qp2.
- Different optimum (3): bearing (SCIP 1.92512 vs 1.95173); ex8_4inf-1
  (SCIP -1e16 vs IbexOpt -1e24) and test_infinity4 (SCIP -1.7e7 vs -6.3e48):
  unbounded problems whose bounded versions have their optimum at the
  ±1e12 corner, which SCIP reports as solved with a far worse value.

Times (n = 98: closed by one of the two, > 1 s for one; floor 1 s, timeout 310 s):
geomean t(IbexOpt)/t(SCIP) = 0.59 [0.28, 1.23], median 0.81, IbexOpt faster on
52, slower on 44; relative t/max: IbexOpt 0.569, SCIP 0.667.
Performance profile: `perfprof12.py` -> `perfprof12.png/pdf` (fastest at τ=1:
IbexOpt 54 %, SCIP 47 %; closed: 86 % vs 69 %).
