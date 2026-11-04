# Step 2 dataset, dive labels with the loup fixed at the optimum

Collected with `collect.sh` (trajectory lsmear-lffix as usual; every dive
prunes with the optimum, `--dive-loup`, so its size measures the branching).
221 instances, 11159 samples; 11067 usable, 197 instances, 122 families.

## Offline (grouped 5-fold CV by family)

`python/oracle_sbfeat_offline.py` (sum = chosen candidates' dive nodes over
the base's own choice, geo over instances; share of instances above 1):

                              regret  best   sum/base  inst>1
    oracle                    1.000   100%   0.814     0%
    strong branching rule     1.042    78%   0.860     4%
    model P (probe)           1.035    80%   0.851     2%
    model A (30 features)     1.176    65%   0.998    23%
    base (lsmear-lffix)       1.172    63%   1.000     -

(dataset-lffix, upper bounding mixed in: oracle 0.757, SB 0.846.) The pure
branching headroom is smaller, and SB takes ~75% of it.

`python/deviate_offline.py` (`deviate-offline.txt`): conservative oracle
thr=0.7 0.828 (deviates on 21%); model A thr=0.9 0.988; model A+P thr=0.9
0.862 with no instance above 1.

`python/hc4probe_offline.py` (`hc4probe-offline.txt`, goal capped at the
search's ymax): HC4 alone 1.01-1.02, HC4+LP 1.01, full contractor 0.95.
Caveat: this offline probe does not reproduce the dive's first step well
(same pruned count 89%, same volume only 60%: fresh ACID tuning, no impact
context), so it underestimates probes; end to end is the test that counts.

## Search history: pseudocosts and reliability branching

`python/history_offline.py`: in search order, per instance, only earlier
samples are used. pc: a variable's mean relative SB key so far. state: a
ridge regression fitted online on this instance's observations, from the
candidate's state (LSmear rank and relative score, smear_sum_rel, rel_diam)
and its history (mean, count). Reliability k: probe the candidates whose
variable has fewer than k observations. Observations are what a search sees
(probed candidates and the bisected one); "ideal" counts every candidate of
every earlier sample.

    rule                                sum/base  inst>1  probed
    SB, probe everything                0.861      4%     100%
    pseudocost, ideal history           0.980     28%       0%
    state, ideal history                0.989     30%       0%
    pseudocost + reliability k=2        0.981     34%       3%
    pseudocost + reliability k=4        0.971     31%       5%
    state + reliability k=1/2/4         0.986 / 0.978 / 0.985

(A first version counted unprobed candidates as observed and reported 0.946
for k=4; this is the corrected accounting.) History carries a weak signal,
and conditioning on the state does not add to the per-variable mean. The
samples are ~150 per instance; a real search observes every node, so these
pseudocosts are built from far less data than a solver's would be.
