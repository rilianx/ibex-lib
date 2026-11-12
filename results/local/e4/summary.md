statuses {('223', 'complete'): 223, ('75', 'no_lower_bound'): 9, ('75', 'complete'): 14, ('75', 'infeasible'): 2, ('75', 'timeout'): 50}
astats vs acvx0 nodes, both closed: 235/237 same [('ex9_2_5.bch', 10, 6), ('robot.bch', 49994, 50172)]

## set 223: 223 instances, 5152511 nodes
| counter | total | % of nodes | instances with > 0 |
|---|---|---|---|
| cvx_nodes | 5152511 | 100.00 |  |
| cvx_with_dual | 640225 | 12.43 |  |
| cvx_certified | 200 | 0.00 | 13 |
| cvx_close | 37527 | 0.73 | 29 |
| cvx_better_lb | 100095 | 1.94 | 51 |
| cvx_center_indef | 0 | 0.00 |  |
| cvx_raised | 63141 | 1.23 |  |
| cvx_pruned | 27029 | 0.52 | 36 |

## set 75: 75 instances, 53957296 nodes
| counter | total | % of nodes | instances with > 0 |
|---|---|---|---|
| cvx_nodes | 53957296 | 100.00 |  |
| cvx_with_dual | 6449420 | 11.95 |  |
| cvx_certified | 456942 | 0.85 | 7 |
| cvx_close | 20658 | 0.04 | 18 |
| cvx_better_lb | 235481 | 0.44 | 23 |
| cvx_center_indef | 0 | 0.00 |  |
| cvx_raised | 194802 | 0.36 |  |
| cvx_pruned | 28788 | 0.05 | 21 |

## set all: 298 instances, 59109807 nodes
| counter | total | % of nodes | instances with > 0 |
|---|---|---|---|
| cvx_nodes | 59109807 | 100.00 |  |
| cvx_with_dual | 7089645 | 11.99 |  |
| cvx_certified | 457142 | 0.77 | 20 |
| cvx_close | 58185 | 0.10 | 47 |
| cvx_better_lb | 335576 | 0.57 | 74 |
| cvx_center_indef | 0 | 0.00 |  |
| cvx_raised | 257943 | 0.44 |  |
| cvx_pruned | 55817 | 0.09 | 57 |

| instance | set | nodes | certified | % | would close | pruned | status | time |
|---|---|---|---|---|---|---|---|---|
| polak6.bch | 75 | 424001 | 413809 | 98 | 0 | 0 | timeout | 1201.6 |
| genhs28.bch | 223 | 25 | 14 | 56 | 4 | 1 | complete | 0.2 |
| schwefel5-abs.bch | 223 | 40 | 21 | 52 | 20 | 10 | complete | 0.1 |
| schwefel5.bch | 223 | 40 | 20 | 50 | 20 | 10 | complete | 0.1 |
| polak3.bch | 75 | 8226 | 3635 | 44 | 3159 | 2493 | complete | 84.6 |
| hs113.bch | 223 | 74 | 31 | 42 | 31 | 12 | complete | 0.5 |
| abelbis.bch | 75 | 339 | 138 | 41 | 138 | 84 | complete | 2.7 |
| linear.bch | 75 | 120889 | 39240 | 32 | 6 | 12419 | timeout | 1200.0 |
| sambal.bch | 223 | 75 | 24 | 32 | 24 | 12 | complete | 0.3 |
| dipigri.bch | 223 | 39 | 12 | 31 | 16 | 14 | complete | 0.3 |
| hs099.bch | 223 | 41 | 11 | 27 | 12 | 12 | complete | 1.3 |
| ex8_1_5.bch | 223 | 16 | 3 | 19 | 2 | 1 | complete | 0.1 |
| hs100.bch | 223 | 42 | 6 | 14 | 8 | 8 | complete | 1.5 |
| ex8_1_4.bch | 223 | 12 | 1 | 8 | 1 | 1 | complete | 0.0 |
| hs119.bch | 223 | 302 | 2 | 1 | 3 | 3 | complete | 5.4 |
| discs.bch | 75 | 25501 | 100 | 0 | 0 | 0 | timeout | 1201.8 |
| odfits.bch | 223 | 1643 | 5 | 0 | 38 | 30 | complete | 7.6 |
| robot.bch | 223 | 37827 | 50 | 0 | 6477 | 5567 | complete | 86.7 |
| water-abs.bch | 75 | 16646 | 10 | 0 | 0 | 0 | timeout | 1200.0 |
| fccu.bch | 75 | 24580 | 10 | 0 | 10 | 6 | complete | 19.8 |
