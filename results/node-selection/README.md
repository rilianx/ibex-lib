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
