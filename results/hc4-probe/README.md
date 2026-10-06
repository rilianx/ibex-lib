# Drastic HC4 probe vs lsmear-lffix, full benchmark

`run.sh` on the user's machine, 298 instances, compo + ipoptxn, 600 s;
`runs.jsonl`, `analyze.py` (t(A)/max(t(A),t(B)), 1 s floor, 2 for a run that
does not close, instances neither closes left out; lower is better).

    arm       closes   vs lffix                         diff   IC95
    lffix     223
    hc4       215      1.034 vs 0.894  (n=224)          0.140  [0.092, 0.194]
    hc4r025   219      1.014 vs 0.919  (n=223)          0.094  [0.062, 0.131]
    hc4need   0        (every run "killed": the binary used did not know
                        --sb-need-pruned; not a result)

On the instances both close: hc4 faster by >1.2x on 3, slower on 41 (geo
time 1.11, nodes 1.04); r=0.25 faster on 2, slower on 35. The gains are
dnieper (0.36x), ex7_2_3, polak5, ship-1; the losses are the mconcon pattern
(conservative deviations that are not conservative: launch 8x, concon 4.8x,
eigencco 4.6x, mconcon 3x, ex6_2_6 2.7x) plus the probe's own cost.

The 10 usual instances had said 0.846 vs 0.939: they over-represent dnieper
and ship-1. The drastic probe is not a better bisector than lsmear-lffix.
