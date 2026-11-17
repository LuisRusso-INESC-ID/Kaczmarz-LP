# lpgraph — classical graph algorithms vs. sparse LP solvers

Benchmark harness for the question *how much time and, above all, how much memory does a
general sparse LP solver need on problems that have dedicated combinatorial algorithms?*
Problems: single-source shortest paths, s–t maximum flow, maximum concurrent multicommodity
flow. Every solve runs in its own process; peak memory is the process's own `VmHWM`.

## Layout

```
src/instance.h      DIMACS .gr / .max readers, own .mcf format, CSR helpers
src/gen.cpp         instance generators (re-implementations of the standard families)
src/classic.cpp     classical baselines (Dijkstra, LEMON Preflow / network simplex, Garg–Könemann)
src/lp.h            sparse LP model (CSC) + the three LP formulations
src/lpsolve.cpp     one LP backend per process: HiGHS, Clp, GLPK, OR-Tools, own prototypes
src/pdhg.h          own restarted PDHG (PDLP recipe: Ruiz + Pock–Chambolle, adaptive restarts)
src/skm.h           own Sampling Kaczmarz–Motzkin on the primal–dual optimality system
src/prowact.h       own proximal-point + Hildreth row action (Kaczmarz-style dual coordinate ascent)
src/ortools_backend.h  GLOP and PDLP through their direct C++ APIs
scripts/run.py      driver (plans "small" and "pilot"; resumable; JSONL output)
scripts/analyze.py  correctness check, tables (summary.md) and plots
scripts/fetch_dimacs_roads.sh  real 9th-DIMACS road graphs (run on your machine)
```

## Build

Dependencies (all built from source in the pilot): HiGHS ≥ 1.15, Clp + CoinUtils, GLPK,
LEMON 1.3.1 headers, optionally OR-Tools (C++). `scripts/build_deps.sh [PREFIX] [SRCDIR]` clones and
builds all of them. Then

```
make DEPS=/path/to/install LEMON=/path/to/lemon/include
```

The defaults are `DEPS=$HOME/lpdeps/install` and `LEMON=$HOME/lpdeps/src/lemon/include`, matching
`scripts/build_deps.sh ~/lpdeps/install ~/lpdeps/src`, so a plain `make` works after that.

OR-Tools is detected automatically if `$(DEPS)/ortools/include/ortools/pdlp/...` exists.

## LP formulations

| problem | variables | rows | notes |
|---|---|---|---|
| SSSP | one per arc, `x ≥ 0` | one per reachable node except the source | supply `r−1` at source, `−1` elsewhere; optimum = Σ dist(v) (checked against Dijkstra) |
| max flow | one per arc, `0 ≤ x ≤ u` | conservation at all nodes except s, t | maximize net outflow of s |
| max concurrent flow | one per (source group, arc) + λ | conservation per group and node (source row dropped) + one capacity row per arc | commodities sharing a source are aggregated (exact for concurrent flow) |

## Methods

* **Classical**: own binary-heap Dijkstra; LEMON Dijkstra; LEMON Preflow (highest-label
  push–relabel, phase 1); LEMON network simplex (the graph-aware simplex, on the same LP);
  Garg–Könemann/Fleischer FPTAS for concurrent flow (ε = 0.1, 0.05; approximate).
* **Simplex (sparse LU)**: HiGHS dual/primal (hyper-sparse revised simplex), Clp dual/primal,
  GLPK primal/dual, OR-Tools GLOP dual/primal.
* **Interior point**: HiGHS IPX (with and without crossover), Clp barrier, GLPK interior.
* **First-order**: HiGHS PDLP (cuPDLP-C), OR-Tools PDLP, own PDHG.
* **Kaczmarz-style**: own SKM (De Loera–Haddock–Needell) on the optimality system; own
  proximal Hildreth row action.

All solvers run single-threaded with relative tolerance 1e-6 (each solver's own definition);
answers are then checked against the exact combinatorial optimum (relative error ≤ 1e-3 counts
as solved).

## Instance families

Shortest path: `sp-grid` (4-neighbour grid), `sp-road` (jittered grid with dropped edges and
Euclidean lengths), `sp-rand` (random digraph, out-degree 4 + Hamiltonian cycle).
Max flow: `mf-rmf` (Goldfarb–Grigoriadis RMF), `mf-rlg` (Washington random level graph),
`mf-vision` (grid with terminal arcs). Multicommodity: `mcf-grid`, `mcf-rand`.
Real road networks: the 9th DIMACS Challenge site is offline, so `scripts/fetch_zenodo_roads.sh`
fetches the challenge's USA travel-time graphs (largest connected component) from Zenodo record
10.5281/zenodo.18026780 and converts them with `scripts/convert_roads.py`: `USA-road-t.NW`
(2.5M arcs), BFS regions of it with 40k/160k/640k arcs, and with `ALL=1` the full USA graph (58M
arcs). The readers accept the DIMACS formats directly.

## Running

```
python3 scripts/run.py --plan pilot --tl 60 --out results/pilot.jsonl
python3 scripts/analyze.py results/pilot.jsonl results/
```

Backends `highs-*` and `clp-*` accept a `-nopre` suffix (presolve off), e.g. `highs-dual-nopre`.
OR-Tools backends are dispatched to `bin/lpsolve-ort` automatically.

## Pilot results (in `results/`)

`FINDINGS.md` (main observations and caveats), `summary.md` (all tables), `sp.png`, `mf.png`,
`mcf.png` (memory and time vs. input size), `pilot-all.jsonl` (raw), `presolve_check.txt`.
Pilot settings: 60 s per solve, tolerance 1e-6, one seed, single thread, 2-core VM; a method is
skipped on larger instances of a family once it fails on a smaller one; SKM only on the two
smallest sizes. Re-generate the tables with
`python3 scripts/analyze.py results/pilot-all.jsonl results 60`.

## Convergence traces

Set `LPTRACE=file.csv` to log (time, iteration, objective, primal residual, dual residual, gap)
for `own-pdhg`, `own-prow` and `own-skm`. SKM diagnostics: `SKM_GAP_EVERY=k` (project the dense
gap row every k steps), `SKM_SCALE=1` (run on the equilibrated LP), `SKM_BALANCE=1` (project the
gap row in a metric that balances its primal and dual halves). `scripts/plot_traces.py` turns the
traces in `results/traces/` into `results/kaczmarz_convergence.png`. Stored traces are thinned to
≤30k rows and gzipped (`*.csv.gz`); the plot scripts read either form.

## Exploratory: hybrid prototype (not in the benchmark)

`own-hybrid-mem`, `own-hybrid-prox` and `own-hybrid-cg` (`src/hybrid.h`) are kept for
reproducibility only. They are not in `scripts/run.py`, the tables or the main plots, because they
are dominated by `own-pdhg` and `own-prow` on every instance tried.

The method keeps the objective out of the constraints: block projected-gradient steps in the null
space of the equality rows (balls of `HYB_RMAX` equality rows, default 64, doubled on stall),
bounds handled by an active set and a bent step, then repair of violated rows from a FIFO queue
with in-queue flags and incremental residuals (only rows of changed columns are re-checked).
Repairs: memoryless Kaczmarz (`mem`), Hildreth with dual memory, i.e. the exact projection
P(x + αd) (`prox`), or a least-norm correction solved by CG, local first with a global fallback
(`cg`). Only `cg` with whole-matrix blocks (`HYB_RMAX=1000000`) converges, and then it is plain
projected gradient. Plot: `python3 scripts/exploratory/plot_hybrid.py` →
`results/exploratory/hybrid_convergence.png`. Summary in `results/FINDINGS.md`.

## Exploratory: primal–dual bisection with Kaczmarz (`own-bisect-kz`, not in the benchmark)

`src/bisect.h`. Keeps bounds L ≤ opt ≤ U and a level t in between; interleaves an SKM run on the
primal slice {rows, box, c'x ≤ t} with a dual step towards {y : LB(y) ≥ t}, where LB is the
Lagrangian bound, valid for any sign-feasible y because every column has a finite box (columns
with cup = ∞ use valid implied bounds `LP::uhint` set by the builders). By strong duality exactly
one side can succeed, so no infeasibility certificate is needed. Dual step: projection onto the
supporting halfspace of {y : LB(y) ≥ t} (`BISECT_DUAL=polyak`), SKM on the explicit dual with z
variables (`skm`), the proximal dual (`prox`, the default; see below), or smoothed coordinate ascent
(`hildreth`): Hildreth row steps on
LB_μ(y) = rup'y + Σ min_box[(c − A'y)_j x_j + μ/2 x_j²], with μ halved when LB_μ stalls, floored at
1e-4·μ₀ and reset to 1e-2·μ₀ when the true bound sticks at a kink; `both` runs hildreth and polyak
on the same y each round. Row selection (default): FIFO queues of violated rows with in-queue
flags, primal rows violated beyond half the current tolerance and dual rows whose smoothed gradient
exceeds it, fed by the incremental residual updates (only rows of changed columns are re-checked);
an empty dual queue also triggers the next μ reduction. `BISECT_QUEUE=0` reverts to random sampling
(SKM for the primal, uniform rows for the dual).

Proximal dual (`BISECT_DUAL=prox`, default): the proximal point method on the dual,
y⁺ = argmax LB(y) − ‖y − y_k‖²/(2σ), carried out as the method of multipliers: x_d minimizes the
augmented Lagrangian c'x − y'(Ax − b) + σ/2‖Ax − b‖² over the box (implied bounds `uhint` for
columns with cup = ∞) by random column coordinate descent (5 sweeps per update), then
y ← y − σ(Ax_d − b) (y ≤ 0 on "<=" rows). σ starts at |c|/(largest box size), is halved (floor
10⁻³σ₀) when an update lowers the Lagrangian bound by more than 0.1%, and doubled after 10 updates
without primal-residual progress (`BISECT_SIGMA` overrides the start). x_d also serves as a primal
candidate. If the valid dual bound ever crosses the approximate primal bound, the primal point is
discarded and the primal tolerance tightened tenfold. Warm starts, jumps, adaptive split, over-relaxed objective rows,
gap-tied primal tolerance (measured in original units), equilibrated LP. Level rule: default is
bisection (t moves only after a success); `BISECT_ALPHA=a BISECT_BETA=b` uses split targets
recomputed every round, the primal asked to improve its bound by a·G and the dual by (1−b)·G
(value terms of a maximization, G = gap; b > a makes the slices overlap). Plot:
`python3 scripts/exploratory/plot_bisect.py` (bisection dashed, split α=1/4, β=3/4 solid).
