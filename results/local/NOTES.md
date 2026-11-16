# Local pilot rerun (2026-11-16)

Machine: Intel Core Ultra 5 125H, 14 GiB, Linux 7.2.11, GCC 16.2.1. Same plan and settings as the
sandbox pilot (`run.py --plan pilot --tl 60`, tolerance 1e-6, one seed, single thread).
Raw: `../pilot-local.jsonl`, `../pilot-local.log`; tables `summary.md`; figures `sp/mf/mcf.png`;
presolve on/off `presolve_check.txt`. 791 runs, no crashes, 0 runs reporting optimal with a wrong
objective.

## Differences from `../FINDINGS.md` (sandbox, 2-core VM)

Memory per arc is essentially unchanged; times are somewhat lower. Updated ranges (sp + mf families,
bytes per arc at the largest solved instance):

| | sandbox | local |
|---|---|---|
| own Dijkstra | 23–28 | 23–27 |
| LEMON Dijkstra | ~55 | 47–50 |
| general simplex | 440–1,600 | 432–1,365 |
| interior point | 500–3,100 | 507–2,279 |
| simplex / Dijkstra | 16–60× | 18–52× |
| simplex / Preflow | 7–27× | 7–21× |
| own PDHG | 184–219 | 189–205 |
| OR-Tools PDLP | 213–337 | 218–297 |
| HiGHS PDLP | 810–1,060 | 555–759 |
| proximal Hildreth | 170–245 | 155–180 |
| Garg–Könemann (mcf) | 60–90 | 43–81 |

* **Correction:** shortest paths via LP are **two to three** orders of magnitude slower than Dijkstra
  (fastest LP code per instance: 100–1,200×), not three to four. The sandbox's own example
  (61 s vs 0.05 s) was also ~1,200×.
* Clp primal now solves sp-road 512² (903k arcs) in 55 s; Clp barrier solves it in 28 s.
* mf-vision 512²: own PDHG 23 s, GLOP dual 23 s, HiGHS IPX+crossover 30 s, Preflow 15 s, LEMON NS 2.3 s.
  HiGHS IPX without crossover returns "Unknown" on vision 256².
* mf-rmf 159k: HiGHS IPX 2.7 s (5.9 s with crossover) vs Preflow 0.96 s.
* Two classical runs exceeded the 120 s wall cap (classical codes have no internal limit):
  LEMON network simplex on mf-rlg 512² and Garg–Könemann ε=0.05 on mcf-rand 5000.
* Presolve: HiGHS no longer peaks at 1.0 GB "in every mode"; PDLP peaks at 0.72 GB (0.64 without
  presolve). Dual simplex 1.04 → 0.77 GB, IPX+crossover 0.93 → 0.73 GB.
* OR-Tools PDLP objective error up to 9e-4 (sandbox 2e-4); still under the 1e-3 threshold.
