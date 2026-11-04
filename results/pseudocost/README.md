# Reliability branching with pseudocosts relative to the bisector

`--oracle-score pc` (MLNodeServer::set_oracle_pc). An observation of j is
log(volume left by j / volume left by lsmear-lffix's choice), both probed at
the same node (bisect + full processing of both children). Candidates with
fewer than 4 observations are probed; a deviation on a pseudocost alone is
verified by probing the pair at that node.

Design steps that failed, on ship-1 / dnieper / schwefel5 / mconcon:
- absolute key (100 x pruned - log volume): biased by where a variable was
  probed (early, near the root) vs observed for free (deep): dnieper times
  out;
- bounded absolute volume left: same bias, schwefel5 times out;
- relative, without verification: once the mean says "better", the variable
  is chosen at every node and the pseudocost goes stale (a new hijack):
  schwefel5 times out, mconcon 16564 nodes.

## Preliminary run: 40 random instances (`prelim-instances.txt`), 120 s

`prelim.jsonl`, cloud machine (~1.9x faster than the user's). Pairwise time
metric (lower is better), `analyze.py`:

    lffix vs pc     n=32  0.936 vs 0.979   diff -0.043  IC95 [-0.098, 0.019]
    lffix vs pcr05  n=32  0.934 vs 0.996   diff -0.062  IC95 [-0.114, -0.017]

On the 32 both close: pc uses 0.83x the nodes (geo) but 1.05x the time;
faster than lffix (>1.2x) on 1, slower on 7. pcr05: 0.91x nodes, 1.08x time.
Gains in nodes (launch 226 -> 36, hs113 384 -> 266) are paid for by the
probes (launch 7.6 -> 10.2 s, eigminc 0.1 -> 1.7 s). Not worth the full run.
