# Results: roads-all-60.jsonl

Process baselines (empty run, used for records without their own base_kb): lpsolve 9.0 MB, lpsolve-ort 15.8 MB, classic 4.1 MB, flow-hipr 2.6 MB, flow-pseudoflow 1.9 MB, flow-bk 4.1 MB. 'net' = peak RSS minus baseline. Cells: solve time (s) / net peak memory (MB); 'TL' time limit, 'x' failed or wrong, '-' skipped after failing a smaller instance.


## sp-roadnet

| method | m=40,000 | m=160,000 | m=640,000 | m=2,545,844 | m=57,708,624 |
|---|---|---|---|---|---|
| classic:dijkstra | 0.00168 / 0.8 | 0.00778 / 3.7 | 0.0352 / 15.1 | 0.152 / 60.9 | 4.77 / 1375.1 |
| classic:lemon-dijkstra | 0.00128 / 1.8 | 0.00584 / 7.9 | 0.0274 / 32.6 | 0.132 / 129.7 | 4.76 / 2841.2 |
| classic:lemon-ns | 0.0444 / 6.2 | 0.506 / 25.5 | 4.39 / 101.4 | 53.8 / 407.9 | x(crash(rc=-6)) |
| clp-barrier | 0.29 / 19.9 | 1.51 / 73.5 | 9.08 / 298.2 | 45.4 / 1159.6 |  |
| clp-barrier-chol | 0.167 / 21.6 | 0.905 / 78.7 | TL / 336 | - |  |
| clp-dual | 1.51 / 21.2 | 25.8 / 79.2 | TL / 296 | - |  |
| clp-primal | 0.391 / 22.5 | 3.32 / 84.2 | 45.6 / 333.2 | TL / 1075 |  |
| glpk-dual | 8.47 / 35.6 | TL / 141 | - | - |  |
| glpk-ipt | 0.135 / 33.4 | 0.736 / 132.0 | 4.26 / 526.5 | 21.2 / 2093.3 |  |
| glpk-primal | 7.68 / 34.9 | TL / 139 | - | - |  |
| gurobi-barrier | 0.146 / 29.9 | 0.8 / 105.2 | 4.68 / 416.4 | 22.1 / 1634.9 |  |
| gurobi-barrier-cross | 0.178 / 38.8 | 0.967 / 127.7 | 6.07 / 508.9 | 32.2 / 1760.3 |  |
| gurobi-dual | 1.11 / 27.9 | 31 / 95.1 | TL / 319 | - |  |
| gurobi-pdhg | 21.2 / 29.5 | TL / 105 | - | - |  |
| gurobi-primal | 0.265 / 32.8 | 2.01 / 117.9 | 15.7 / 442.8 | TL / 1659 |  |
| highs-dual | 1.68 / 43.2 | 35.9 / 159.4 | TL / 537 | - |  |
| highs-ipx | x(Unknown, err 8e-01) | - | - | - |  |
| highs-ipx-cross | 2.17 / 41.1 | 38.5 / 161.5 | TL / 533 | - |  |
| highs-pdlp | TL / 24 | - | - | - |  |
| highs-primal | 0.6 / 46.2 | 8.99 / 157.6 | TL / 533 | - |  |
| ortools-glop | 2.17 / 53.7 | 36.3 / 232.1 | TL / 929 | - |  |
| ortools-glop-primal | 1 / 52.1 | 12.3 / 223.6 | TL / 931 | - |  |
| ortools-pdlp | TL / 11 | - | - | - |  |
| own-pdhg | TL / 8 | - | - | - |  |
| own-prow | TL / 7 | - | - | - |  |
| own-skm | TL / 6 | - |  |  |  |

## Net memory per input arc at the largest solved instance (bytes/arc)

| family | classic:dijkstra | classic:lemon-dijkstra | classic:lemon-ns | clp-barrier | clp-barrier-chol | clp-dual | clp-primal | glpk-dual | glpk-ipt | glpk-primal | gurobi-barrier | gurobi-barrier-cross | gurobi-dual | gurobi-pdhg | gurobi-primal | highs-dual | highs-ipx | highs-ipx-cross | highs-pdlp | highs-primal | ortools-glop | ortools-glop-primal | ortools-pdlp | own-pdhg | own-prow | own-skm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sp-roadnet | 25 (m=57708k) | 52 (m=57708k) | 168 (m=2545k) | 478 (m=2545k) | 515 (m=160k) | 519 (m=160k) | 546 (m=640k) | 933 (m=40k) | 862 (m=2545k) | 916 (m=40k) | 673 (m=2545k) | 725 (m=2545k) | 623 (m=160k) | 774 (m=40k) | 725 (m=640k) | 1045 (m=160k) |  | 1059 (m=160k) |  | 1033 (m=160k) | 1521 (m=160k) | 1466 (m=160k) |  |  |  |  |

## Correctness

Reference: {'sp': 'classic:dijkstra', 'mf': 'classic:lemon-preflow', 'mcf': 'highs-ipx-cross'}. Runs reporting optimal but not counted as solved (relative objective error > 0.001, or solve time over the limit): 0


Largest relative objective error among solved runs, per method: classic:dijkstra 0e+00, classic:lemon-dijkstra 0e+00, classic:lemon-ns 0e+00, clp-barrier 9e-13, clp-barrier-chol 0e+00, clp-dual 0e+00, clp-primal 0e+00, glpk-dual 0e+00, glpk-ipt 7e-09, glpk-primal 0e+00, gurobi-barrier 6e-13, gurobi-barrier-cross 6e-13, gurobi-dual 0e+00, gurobi-pdhg 7e-07, gurobi-primal 0e+00, highs-dual 0e+00, highs-ipx-cross 0e+00, highs-primal 0e+00, ortools-glop 0e+00, ortools-glop-primal 0e+00
