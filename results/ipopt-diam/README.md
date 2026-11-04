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
