statuses {('223', 'complete'): 223, ('75', 'complete'): 17, ('75', 'no_lower_bound'): 9, ('75', 'infeasible'): 2, ('75', 'timeout'): 47}
astatsc vs acomb nodes, both closed: 238/240 same [('launch.bch', 210, 26), ('robot.bch', 49980, 51610)]

## set 223: 223 instances, 5134929 nodes
| counter | total | % of nodes | instances with > 0 |
|---|---|---|---|
| cvx_nodes | 5134929 | 100.00 |  |
| cvx_with_dual | 632612 | 12.32 |  |
| cvx_certified | 1765 | 0.03 | 16 |
| cvx_close | 42520 | 0.83 | 41 |
| cvx_better_lb | 114105 | 2.22 | 67 |
| cvx_center_indef | 0 | 0.00 |  |
| cvx_raised | 71637 | 1.40 |  |
| cvx_pruned | 30994 | 0.60 | 47 |

## set 75: 75 instances, 54407112 nodes
| counter | total | % of nodes | instances with > 0 |
|---|---|---|---|
| cvx_nodes | 54407112 | 100.00 |  |
| cvx_with_dual | 6718224 | 12.35 |  |
| cvx_certified | 421943 | 0.78 | 8 |
| cvx_close | 86369 | 0.16 | 25 |
| cvx_better_lb | 802567 | 1.48 | 33 |
| cvx_center_indef | 0 | 0.00 |  |
| cvx_raised | 636687 | 1.17 |  |
| cvx_pruned | 86050 | 0.16 | 29 |

## set all: 298 instances, 59542041 nodes
| counter | total | % of nodes | instances with > 0 |
|---|---|---|---|
| cvx_nodes | 59542041 | 100.00 |  |
| cvx_with_dual | 7350836 | 12.35 |  |
| cvx_certified | 423708 | 0.71 | 24 |
| cvx_close | 128889 | 0.22 | 66 |
| cvx_better_lb | 916672 | 1.54 | 100 |
| cvx_center_indef | 0 | 0.00 |  |
| cvx_raised | 708324 | 1.19 |  |
| cvx_pruned | 117044 | 0.20 | 76 |

| instance | set | nodes | certified | % | would close | pruned | status | time |
|---|---|---|---|---|---|---|---|---|
| polak6.bch | 75 | 388101 | 374035 | 96 | 0 | 0 | timeout | 1200.9 |
| genhs28.bch | 223 | 25 | 14 | 56 | 4 | 1 | complete | 0.2 |
| schwefel5-abs.bch | 223 | 40 | 21 | 52 | 20 | 10 | complete | 0.1 |
| schwefel5.bch | 223 | 40 | 20 | 50 | 20 | 10 | complete | 0.1 |
| polak3.bch | 75 | 8226 | 3635 | 44 | 3159 | 2493 | complete | 85.6 |
| hs117.bch | 223 | 3568 | 1527 | 43 | 489 | 257 | complete | 46.5 |
| hs113.bch | 223 | 74 | 31 | 42 | 31 | 12 | complete | 0.5 |
| abelbis.bch | 75 | 339 | 138 | 41 | 138 | 84 | complete | 2.7 |
| hs099.bch | 223 | 42 | 16 | 38 | 16 | 16 | complete | 1.0 |
| haifas.bch | 223 | 17 | 6 | 35 | 5 | 1 | complete | 0.1 |
| dipigri.bch | 223 | 39 | 13 | 33 | 16 | 14 | complete | 0.3 |
| linear.bch | 75 | 121914 | 39915 | 33 | 6 | 12662 | timeout | 1200.0 |
| sambal.bch | 223 | 75 | 24 | 32 | 24 | 12 | complete | 0.3 |
| ex8_1_5.bch | 223 | 16 | 3 | 19 | 2 | 1 | complete | 0.1 |
| hs100.bch | 223 | 42 | 6 | 14 | 9 | 9 | complete | 1.5 |
| hs119.bch | 223 | 302 | 26 | 9 | 20 | 2 | complete | 5.4 |
| ex8_1_4.bch | 223 | 12 | 1 | 8 | 1 | 1 | complete | 0.0 |
| discs.bch | 75 | 26607 | 1200 | 5 | 0 | 1 | timeout | 1200.0 |
| mistake.bch | 75 | 78858 | 3000 | 4 | 155 | 76 | complete | 544.0 |
| odfits.bch | 223 | 1643 | 5 | 0 | 38 | 30 | complete | 7.6 |
| dixchlng.bch | 223 | 1173 | 3 | 0 | 2 | 0 | complete | 4.4 |
| robot.bch | 223 | 38569 | 49 | 0 | 6633 | 5590 | complete | 91.0 |
| water-abs.bch | 75 | 15568 | 10 | 0 | 0 | 0 | timeout | 1200.0 |
| fccu.bch | 75 | 24580 | 10 | 0 | 10 | 6 | complete | 19.9 |
