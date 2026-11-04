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
