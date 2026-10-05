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
