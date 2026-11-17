# Results: roads-all-60.jsonl

Process baselines (empty run): lpsolve 7.5 MB, lpsolve-ort 15.7 MB, classic 4.1 MB. 'net' = peak RSS minus baseline. Cells: solve time (s) / net peak memory (MB); 'TL' time limit, 'x' failed or wrong, '-' skipped after failing a smaller instance.


## sp-roadnet

| method | m=40,000 | m=160,000 | m=640,000 | m=2,545,844 | m=57,708,624 |
|---|---|---|---|---|---|
| classic:dijkstra | 0.00168 / 0.8 | 0.00778 / 3.7 | 0.0352 / 15.1 | 0.152 / 60.9 | 4.77 / 1375.1 |
| classic:lemon-dijkstra | 0.00128 / 1.8 | 0.00584 / 7.9 | 0.0274 / 32.6 | 0.132 / 129.7 | 4.76 / 2841.2 |
| classic:lemon-ns | 0.0444 / 6.2 | 0.506 / 25.5 | 4.39 / 101.4 | 53.8 / 407.9 | x(crash(rc=-6)) |
| clp-barrier | 0.29 / 20.0 | 1.51 / 73.6 | 9.08 / 298.2 | 45.4 / 1159.6 |  |
| clp-dual | 1.51 / 21.2 | 25.8 / 79.2 | TL / 296 | - |  |
| clp-primal | 0.391 / 22.6 | 3.32 / 84.3 | 45.6 / 333.2 | TL / 1075 |  |
| glpk-dual | 8.7 / 35.9 | TL / 142 | - | - |  |
| glpk-ipt | 0.135 / 33.4 | 0.735 / 132.1 | 4.32 / 523.4 | 21.2 / 2091.5 |  |
| glpk-primal | 8.34 / 34.2 | TL / 136 | - | - |  |
| highs-dual | 1.68 / 43.3 | 35.9 / 159.5 | TL / 537 | - |  |
| highs-ipx | x(Unknown, err 8e-01) | - | - | - |  |
| highs-ipx-cross | 2.17 / 41.2 | 38.5 / 161.6 | TL / 533 | - |  |
| highs-pdlp | TL / 24 | - | - | - |  |
| highs-primal | 0.6 / 46.3 | 8.99 / 157.7 | TL / 533 | - |  |
| ortools-glop | 2.17 / 53.8 | 36.3 / 232.2 | TL / 929 | - |  |
| ortools-glop-primal | 1 / 52.1 | 12.3 / 223.7 | TL / 931 | - |  |
| ortools-pdlp | TL / 11 | - | - | - |  |
| own-pdhg | TL / 9 | - | - | - |  |
| own-prow | TL / 7 | - | - | - |  |
| own-skm | TL / 6 | - |  |  |  |

## Net memory per input arc at the largest solved instance (bytes/arc)

| family | classic:dijkstra | classic:lemon-dijkstra | classic:lemon-ns | clp-barrier | clp-dual | clp-primal | glpk-dual | glpk-ipt | glpk-primal | highs-dual | highs-ipx | highs-ipx-cross | highs-pdlp | highs-primal | ortools-glop | ortools-glop-primal | ortools-pdlp | own-pdhg | own-prow | own-skm |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sp-roadnet | 25 (m=57708k) | 52 (m=57708k) | 168 (m=2545k) | 478 (m=2545k) | 519 (m=160k) | 546 (m=640k) | 941 (m=40k) | 861 (m=2545k) | 895 (m=40k) | 1045 (m=160k) |  | 1059 (m=160k) |  | 1034 (m=160k) | 1521 (m=160k) | 1466 (m=160k) |  |  |  |  |

## Correctness

Reference: {'sp': 'classic:dijkstra', 'mf': 'classic:lemon-preflow', 'mcf': 'highs-ipx-cross'}. Runs reporting optimal but with relative objective error > 0.001: 0


Largest relative objective error among solved runs, per method: classic:dijkstra 0e+00, classic:lemon-dijkstra 0e+00, classic:lemon-ns 0e+00, clp-barrier 9e-13, clp-dual 0e+00, clp-primal 0e+00, glpk-dual 0e+00, glpk-ipt 7e-09, glpk-primal 0e+00, highs-dual 0e+00, highs-ipx-cross 0e+00, highs-primal 0e+00, ortools-glop 0e+00, ortools-glop-primal 0e+00
