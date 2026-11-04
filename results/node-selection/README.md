# How much can node selection gain?

`headroom.sh`: per instance, lsmear-lffix twice, both with ACID
non-adaptive on all variables (`--acid-fixed`): as usual (base), and with
the optimum given right after the root (`--late-loup` = base's final loup).
No node order can do better than the second; the root's own processing is
the same in both (unlike `--initial-loup`, which on dnieper changes the
whole tree: 196 to >7000 nodes depending on the bound's value).

49 instances (the 40 random of results/pseudocost + the 10 usual), 120 s,
`headroom-prelim-late.jsonl` (`headroom-prelim-adaptive.jsonl`: a first
version with --initial-loup and adaptive ACID, not to be used).

41 close; 21 with at least 50 nodes. The headroom is all or nothing:

    instance   nodes  final loup at node   nodes with the optimum
    polak4     1770   1766                 2
    ship-1      580    580                 0
    launch      104    100                 2
    ex7_3_4     124    102                 28
    (mconcon 18 -> 12, ex14_1_2 40 -> 28, ex14_1_6 30 -> 2: small trees)

Every other instance with >= 50 nodes is at 0.96-1.08: the final upper
bound is there at the root or within the first nodes (dnieper, ex6_2_8,
ex8_4_5bis, ex2_1_9fact, dipigri, schwefel5 ...). On the 4 with headroom,
the good point appears at the very end of the search: the tree is spent
proving bounds on regions while the optimum is found last.

Not covered: 5 timeouts (abelter-1, fccu, hs102, orthrege, water) and 3
unbounded/infeasible; with no proven optimum their headroom is unknown.

## All the instances lsmear-lffix closes (223), 310 s here (~600 s on the user's machine)

`headroom-all.jsonl`, `closed-instances.txt`. 222 close with ACID fixed
(1 timeout); 94 have at least 50 nodes. Nodes with the optimum given right
after the root, over the usual run:

    ratio <= 0.9   25 of 94
    ratio <= 0.8   22
    ratio <= 0.5   19
    ratio <= 0.2   15
    time <= 0.8x   10 (1 s floor)

The instances with headroom (ratio <= 0.8), most of them all-or-nothing --
the final upper bound appears within the last nodes:

    instance     nodes     time   final loup at node   nodes with the optimum
    ship-1         580      1.3        580                 0
    eigencco      1640      8.4       1640                 0
    makela3       2078      2.5       2078                 0
    matrix2    4511092    207.9    4511092                 0
    polak2        1920      0.5       1918                 2
    ex8_2_4       1846     24.8       1834                 2
    polak4        1770      0.2       1766                 2
    polak5        3190      1.4       3184                 4
    polak1         706      0.0        704                 2
    ex8_1_6bis     176      0.0        176                 2
    concon         598      0.6        590                10
    launch         104      4.8        100                 2
    hs116           98      0.8         94                 2
    dixchlng       698      1.5        602                18
    eigmaxc        300      1.3        298                30
    ex7_3_4        124      0.4        102                28
    m4wd            84      0.1         54                24
    ex4_1_5         54      0.0         38                22
    haldmads       262      1.3        124               110
    haifas         102      0.1         12                68
    ex7_3_4bis     264      0.2        132               198
    genhs28         58      0.2         10                46

None gets worse with the optimum (> 1.2x). In time, only 10 gain more than
20% (1 s floor): most of these trees are fast anyway; the ones that matter
are matrix2 (208 s), ex8_2_4 (25 s), eigencco (8 s), launch (5 s).

## Why does the final upper bound appear so late? (`python/ub_path.py`)

For each instance with headroom, the path of boxes from the root to the
optimum x* (the search's bisector, the search's contraction, ACID fixed): at
each level, whether the loup finder (Ipopt at every call) reaches the optimum
in that box, and the lower bound of the child containing x* vs its sibling's.
`ub-path.txt` (+ launch and polak5, run first).

    type                                        instances            loup finder reaches x* at
    A  optimum ~0, lower bound ~0 everywhere    ship-1, eigencco,    never on the path (ship-1: still
       (flat lb; termination needs a point      makela3, polak4,     1.09 at width 4e-14; matrix2:
       within abs_eps_f=1e-7 of 0)              matrix2              none at width 9e-5)
    B  lower bound misleads: x*'s region has    polak1, polak2,      depth 4 (polak1, polak5), 33
       lb = optimum, siblings lower             polak5               (polak2), boxes still unbounded
    C  path forced or best-first would follow   concon, eigmaxc,     depth 292, 137, 83, 31, 24, 5,
       it (siblings pruned or with higher lb),  ex8_1_6bis, ex7_3_4, 149: the finder fails on huge
       but the finder only succeeds deep         hs116, launch,       (often unbounded) boxes
                                                dixchlng
    D  lb nearly equal everywhere, finder       ex8_2_4              not by depth 30
       fails high in the tree

Only type B is a node-selection problem in the strict sense: the point is a
few levels down in a region best-first leaves for last, and the finder gets
it there. In A, C and D the bottleneck is the upper bounding on large boxes
(or the absolute precision near 0); a node order can only reach the deep box
sooner, which in C is already the forced path (concon: 292 levels, siblings
pruned). Launch is C with a schedule issue: x* is 5 levels down on the
lowest-lb path, but Ipopt runs every 100 loup-finder calls.

Caveat: Ipopt at every call here; in the real search Ipopt runs at calls 0,
10, 20, 50 and every 100, so "finder reaches x* at depth d" is an optimistic
bound for the default schedule.

## The real cause in types C and D: Ipopt is never called on wide boxes

The loup finder runs at every node (probing, XTaylor); Ipopt only on its
schedule (loup-finder calls 0, 10, 20, 50, then every 100) AND only when the
box's widest variable is at most `ipopt_diam = 1e8`
(src/ipopt/ibex_LoupFinderIpopt.cpp). Variables with "infinite" bounds
written as +-1e8 keep the box at width 2e8: on eigencco and ex8_2_4 Ipopt is
never called along the path; on dixchlng the optimum appears at depth 149,
exactly when the width drops to 1e8.

`--ipopt-diam` makes the threshold configurable. lsmear-lffix, 300 s:

    instance    ipopt_diam 1e8                 1e300 (no threshold)
    ex8_2_4     1832 nodes  20.6 s             42 nodes   2.5 s
    eigencco    1188         3.9 s             12         0.5 s
    dixchlng    2038         2.1 s             110        0.3 s
    launch       260         8.6 s             16         0.6 s
    ex7_3_4      126         0.4 s             34         0.2 s
    concon       610         0.4 s             610        0.6 s
    eigmaxc      304         1.0 s             304        1.0 s
    ex8_1_6bis   176         0.0 s             176        0.0 s

Not yet checked on the instances where it could cost (Ipopt on huge boxes
is slower and may fail): needs the full benchmark.
