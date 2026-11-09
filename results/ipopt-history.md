# Solved instances since the Ipopt work started (benchs/optim/all, 298 instances, 310 s)

"Closed" = status complete (optimum proved within rel 1e-3 / abs 1e-7). The
298 instances = 223 that lsmear-lffix closed at the start
(node-selection/closed-instances.txt) + 75 open (anytime/open-instances.txt).
Times: total over closed instances; nodes: geometric mean vs d1e8 over the
instances both close.

| arm | date | what changed | closed /223 | total time | nodes vs d1e8 | closed /75 open | total /298 |
|---|---|---|---|---|---|---|---|
| d1e8 (+ open `base`) | 11-04 | start: lsmear-lffix, both, ipoptxn, Ipopt only on boxes <= 1e8 wide | 223 | 1651 s | 1.000 | 1 (chenery-1) | 224 |
| d1e12 | 11-04 | threshold 1e12 | 223 | 1618 s | 0.910 | - | - |
| dinf | 11-04 | threshold 1e300 | 222 (matrix2 TO) | 1585 s | 0.910 | - | - |
| dnone | 11-04 | no threshold, finite start on unbounded boxes | 222 (matrix2 TO) | 1427 s | 0.770 | - | - |
| x1e5 / c1e5 | 11-05 | correction box 1e-5 / cascade 1e-6 then 1e-5 | 222-223 | 1662-1759 s | 0.97-0.98 | - | - |
| nc1e5 | 11-05 | no threshold + cascade | 222 (matrix2) | 1519 s | 0.744 | - | - |
| final (+ open `ipoptfix`) | 11-05 | + no Ipopt from starts beyond 1e8 (matrix2 fixed) | 223 | 1632 s | 0.763 | 3 (chenery-1, hs093, mistake) | 226 |
| base8 / fixaff | 11-06 | affine relaxation fix (Affine2 operator=) | 223 / 223 | 1682 / 1707 s | 0.999 / 0.767 | - | - |
| ms20 | 11-06 | + 20 random Ipopt starts at the root | 221 (Ipopt hangs) | 1537 s | 0.649 | - | - |
| tbase | 11-06 | old base, all fixes, root time counted | 223 | 1662 s | 0.996 | 3 (chenery-1, hs093, wall) | 226 |
| tfinal | 11-06 | new Ipopt base, all fixes, root counted | 223 | 1669 s | 0.765 | 3 | 226 |
| tms5 / tms20 | 11-06 | + 5 / 20 root restarts (Ipopt call limit 5 s) | 223 / 223 | 1579 / 1674 s | 0.666 / 0.637 | - / 3 | - / 226 |
| tg10 / tg10l | 11-06 | + guided restarts (10% budget) global / local | 222 (ex6_2_13 TO) | 1429 / 1419 s | 0.686 / 0.684 | 3 | 225 |
| tg10p | 11-06 | + patience 20 (recommended) | 223 | 1599 s | 0.696 | 3 | 226 |
| tg10p + thick-equality fix | 11-06 | hs118 ranges certifiable | (223) | | | 4 (+ hs118) | 227 |
| g05p20 | 11-07 | guided restarts rho 0.05, patience 20 | 223 | 1539 s | 0.665 | 4 (+ mistake) | 227 |
| tcvx | 11-08 | + certified-convex Lagrangian bound (--cvx-close), thick fix | 222 (ex6_2_13 TO) | 1384 s | 0.763 | 11 (+ abel, abelbis, abelbis-1, abelter-1, fccu, polak3, hs118, mistake) | 233 |
| gcvx | 11-08 | g05p20 + --cvx-close | 222 (ex6_2_13 TO) | 1186 s | 0.523 vs tbase | 11 | 233 |

Other status changes on the open 75:
- wall: proved *infeasible* by the binaries before the affine fix, closes
  (optimum found) after it -- the invalid affine cuts emptied its root.
- mistake: closes at 299 s with the pre-fix binary without the root's time;
  ~310 s afterwards (borderline, timeout in the t* arms).
- matrix2 (223 set): 130-140 s in every arm without restarts, 0.3 s with
  guided restarts (tg10*) or root multistart.
- Proved unbounded: 8-9, infeasible: 2 (ex7_3_6, launch__coconut) in all
  arms; no feasible point at 310 s: 9 instances with the old base, 4 with
  tg10p + the thick fix (discs, ex5_3_3, minlphi, minlphi1).

## All together (arm `full`), and a clean timing comparison

full = --ipopt-diam=-1 --ipopt-expansion=1e-5 --ipopt-expansion-first=0.1
--ipopt-guided=0.05 --ipopt-guided-patience=20 --cvx-close --cvx-backoff=64
(+ affine fix, thick-equality fix, root time counted).

Closed: 223/223 + 11/75 = 234/298 (start: 224). Every optimum agrees with
optima.txt (ex8_2_4's reference is the known wrong one).

The container restarted right before this run on a ~30% faster host: against
gcvx, full has the same node counts (ratio 1.000) and ~30% lower times on
dozens of instances. tbase2 reruns the original base on the same host:

    A      vs B       time diff [IC95]          nodes  total time (closed)
    tbase2 vs tbase   -0.067 [-0.085, -0.050]   1.016  1229 vs 1662 s   (host speed only)
    full   vs tbase2  -0.007 [-0.029, +0.013]   0.516  1032 vs 1229 s

Same host: half the nodes, -16% total time, a tie on the pairwise metric
(1 s floor: most instances take < 1 s). Faster: harker 93 -> 0.0 s, sambal
18.6 -> 0.0, matrix2 94 -> 0.2, launch 9.8 -> 0.6, eigencco 3.6 -> 0.2,
dixchlng 1.7 -> 0.4, makela3 1.3 -> 0.0. Slower: dnieper 4.2 -> 5.5,
ex8_2_4 2.7 -> 6.0, ex9_2_5 1.0 -> 2.2, haldmads 1.1 -> 2.2.

Open 75 (62 with metrics in all arms): primal integral 0.182 (tbase) ->
0.135 (full), primal-dual 0.601 -> 0.505; no feasible point 9 -> 4 (discs,
ex5_3_3, minlphi, minlphi1).

## full2: full + shifted (alphaBB-like) Lagrangian bound (--cvx-test 2)

Same host as tbase2 / full.

    A      vs B       time diff [IC95]          nodes  total time (223)
    full2  vs full    -0.020 [-0.031, -0.009]   0.960  908 vs 1032 s
    full2  vs tbase2  -0.029 [-0.052, -0.008]   0.495  908 vs 1229 s

Faster than full: robot 62.2 -> 29.4 s, dnieper 5.5 -> 2.1, ex6_2_6 8.5 ->
5.1, ex6_2_8, ex6_2_11, ex8_4_4-1/-2; slower: haldmads 2.2 -> 3.4.
Open 75: 13 closed (+ catena, ex8_4_1), primal integral 0.120, primal-dual
0.487 (full: 0.135 / 0.505; tbase: 0.182 / 0.601). Optima agree with
optima.txt; catena and ex8_4_1 close at the loup every arm finds.

**Total: 223 + 13 = 236 of 298** (start: 224). Against the original base on
the same host: half the nodes, total time 1229 -> 908 s (-26%), pairwise
time -0.029 [-0.052, -0.008].
