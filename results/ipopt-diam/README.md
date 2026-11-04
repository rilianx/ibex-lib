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
