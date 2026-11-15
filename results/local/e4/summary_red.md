statuses {('223', 'complete'): 223, ('75', 'complete'): 17, ('75', 'no_lower_bound'): 9, ('75', 'infeasible'): 2, ('75', 'timeout'): 47}
astatsr vs acombr nodes, both closed: 237/240 same [('ex9_2_5.bch', 4, 6), ('launch.bch', 26, 210), ('robot.bch', 33506, 33292)]

## set 223: 223 instances, 5107906 nodes
| counter | total | % of nodes | instances with > 0 |
|---|---|---|---|
| cvx_reduced | 39030 | 0.76 |  |
| cvx_reduce_emptied | 20618 | 0.40 |  |
| cvx_nodes | 5107906 | 100.00 |  |
| cvx_with_dual | 582641 | 11.41 |  |
| cvx_certified | 1365 | 0.03 | 17 |
| cvx_close | 25269 | 0.49 | 42 |
| cvx_better_lb | 90024 | 1.76 | 68 |
| cvx_center_indef | 0 | 0.00 |  |
| cvx_raised | 29811 | 0.58 |  |
| cvx_pruned | 42 | 0.00 | 1 |

## set 75: 75 instances, 53414645 nodes
| counter | total | % of nodes | instances with > 0 |
|---|---|---|---|
| cvx_reduced | 364257 | 0.68 |  |
| cvx_reduce_emptied | 43955 | 0.08 |  |
| cvx_nodes | 53414645 | 100.00 |  |
| cvx_with_dual | 6030451 | 11.29 |  |
| cvx_certified | 154745 | 0.29 | 8 |
| cvx_close | 42564 | 0.08 | 24 |
| cvx_better_lb | 614238 | 1.15 | 32 |
| cvx_center_indef | 0 | 0.00 |  |
| cvx_raised | 428567 | 0.80 |  |
| cvx_pruned | 13922 | 0.03 | 13 |

## set all: 298 instances, 58522551 nodes
| counter | total | % of nodes | instances with > 0 |
|---|---|---|---|
| cvx_reduced | 403287 | 0.69 |  |
| cvx_reduce_emptied | 64573 | 0.11 |  |
| cvx_nodes | 58522551 | 100.00 |  |
| cvx_with_dual | 6613092 | 11.30 |  |
| cvx_certified | 156110 | 0.27 | 25 |
| cvx_close | 67833 | 0.12 | 66 |
| cvx_better_lb | 704262 | 1.20 | 100 |
| cvx_center_indef | 0 | 0.00 |  |
| cvx_raised | 458378 | 0.78 |  |
| cvx_pruned | 13964 | 0.02 | 14 |

| instance | set | nodes | certified | % | would close | pruned | status | time |
|---|---|---|---|---|---|---|---|---|
| polak6.bch | 75 | 138402 | 107819 | 78 | 0 | 0 | timeout | 1200.3 |
| schwefel5-abs.bch | 223 | 40 | 21 | 52 | 20 | 0 | complete | 0.1 |
| schwefel5.bch | 223 | 40 | 20 | 50 | 20 | 0 | complete | 0.1 |
| genhs28.bch | 223 | 22 | 11 | 50 | 4 | 0 | complete | 0.2 |
| hs113.bch | 223 | 74 | 31 | 42 | 31 | 0 | complete | 0.5 |
| polak3.bch | 75 | 7956 | 3323 | 42 | 3090 | 0 | complete | 75.8 |
| abelbis.bch | 75 | 339 | 138 | 41 | 138 | 0 | complete | 2.7 |
| hs099.bch | 223 | 42 | 16 | 38 | 16 | 0 | complete | 1.0 |
| hs117.bch | 223 | 3117 | 1139 | 37 | 361 | 0 | complete | 39.4 |
| haifas.bch | 223 | 17 | 6 | 35 | 5 | 0 | complete | 0.1 |
| linear.bch | 75 | 122932 | 40745 | 33 | 6 | 12430 | timeout | 1200.0 |
| sambal.bch | 223 | 75 | 24 | 32 | 24 | 0 | complete | 0.3 |
| dipigri.bch | 223 | 30 | 8 | 27 | 9 | 0 | complete | 0.3 |
| ex8_1_5.bch | 223 | 16 | 3 | 19 | 2 | 0 | complete | 0.1 |
| hs100.bch | 223 | 38 | 6 | 16 | 9 | 0 | complete | 1.5 |
| hs119.bch | 223 | 292 | 28 | 10 | 26 | 0 | complete | 5.3 |
| ex8_1_4.bch | 223 | 12 | 1 | 8 | 1 | 0 | complete | 0.0 |
| ex4_1_5.bch | 223 | 28 | 2 | 7 | 1 | 0 | complete | 0.1 |
| mistake.bch | 75 | 77067 | 2100 | 3 | 94 | 2 | complete | 528.2 |
| discs.bch | 75 | 27101 | 600 | 2 | 0 | 9 | timeout | 1200.2 |
| odfits.bch | 223 | 1631 | 5 | 0 | 29 | 0 | complete | 7.6 |
| dixchlng.bch | 223 | 1173 | 3 | 0 | 2 | 0 | complete | 4.4 |
| robot.bch | 223 | 26580 | 41 | 0 | 1621 | 0 | complete | 45.6 |
| water-abs.bch | 75 | 15542 | 10 | 0 | 0 | 0 | timeout | 1200.0 |
| fccu.bch | 75 | 24580 | 10 | 0 | 10 | 0 | complete | 19.8 |
