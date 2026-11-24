# lpgraph — classical graph algorithms vs. sparse LP solvers

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24335608.svg)](https://doi.org/10.5281/zenodo.24335608)

Benchmark harness for the question *how much time and, above all, how much memory does a
general sparse LP solver need on problems that have dedicated combinatorial algorithms?*
Problems: single-source shortest paths, s–t maximum flow, maximum concurrent multicommodity
flow. Every solve runs in its own process; peak memory is the process's own `VmHWM`.

This is the code and data for the paper *The Cost of Generality:
Memory and Time of Open-Source LP Solvers on Network Problems*,
submitted to ACM Transactions on Mathematical Software as a research
paper, with a Replicated Computational Results (RCR) review to be
requested. Every table and figure in the paper is regenerated from
the raw results in `results/` by the scripts below; see
[Replicating the paper](#replicating-the-paper).
Archived on Zenodo: [10.5281/zenodo.24335608](https://doi.org/10.5281/zenodo.24335608) (all
versions; v1.2, used in the paper, is [10.5281/zenodo.24368724](https://doi.org/10.5281/zenodo.24368724); v1.1 is [10.5281/zenodo.24341292](https://doi.org/10.5281/zenodo.24341292); v1.0 is [10.5281/zenodo.24335609](https://doi.org/10.5281/zenodo.24335609)). Citation
metadata is in `CITATION.cff`.

## Layout

```
src/instance.h      DIMACS .gr / .max readers, own .mcf format, CSR helpers
src/gen.cpp         instance generators (re-implementations of the standard families)
src/classic.cpp     classical baselines (Dijkstra, LEMON Preflow / network simplex, Garg–Könemann)
src/lp.h            sparse LP model (CSC) + the three LP formulations
src/lpsolve.cpp     one LP backend per process: HiGHS, Clp, GLPK, OR-Tools, Gurobi, own prototypes
src/pdhg.h          own restarted PDHG (PDLP recipe: Ruiz + Pock–Chambolle, adaptive restarts)
src/skm.h           own Sampling Kaczmarz–Motzkin on the primal–dual optimality system
src/prowact.h       own proximal-point + Hildreth row action (Kaczmarz-style dual coordinate ascent)
src/ortools_backend.h  GLOP and PDLP through their direct C++ APIs
src/gurobi_backend.h   Gurobi (commercial, needs a license) through its C API
src/flowref/        drivers for the reference max-flow codes HIPR, pseudoflow, Boykov–Kolmogorov
scripts/build_deps.sh   clone and build HiGHS, Clp, GLPK, LEMON, OR-Tools from source
scripts/run.py      driver (plans "small", "pilot", "roads", "roads-usa"; --seed; resumable; JSONL output)
scripts/memwrap.py  per-solve wrapper: address-space cap, timeout, CPU pinning (LPBENCH_CPU)
scripts/analyze.py  correctness check, tables (summary.md) and plots
scripts/seed_stats.py   ranges and run-to-run spread over seeds (results/seed_stats.txt)
scripts/paper_tables.py bytes/arc LaTeX table for the paper
scripts/plot_roads.py   road-network figure
scripts/plot_traces.py  convergence figure from results/traces/
scripts/fetch_flow_codes.sh    download HIPR 4.2, pseudoflow 3.23, BK maxflow 3.04 (checksummed)
scripts/depth_experiment.sh    hop-depth experiment for first-order methods (results/depth/)
scripts/hop_depth.py           hop depth of the shortest-path tree of a .gr instance
scripts/fetch_zenodo_roads.sh  9th-DIMACS USA road graphs from Zenodo, converted (convert_roads.py)
scripts/fetch_dimacs_roads.sh  same from the original challenge site, which is offline (404);
                    kept for reference, set DIMACS_URL to use a mirror
scripts/exploratory/    plots for the exploratory prototypes
```

## Build

Dependencies (all built from source): HiGHS 1.16, Clp + CoinUtils (with SuiteSparse AMD and
CHOLMOD), GLPK 5.0, LEMON 1.3.1 headers, optionally OR-Tools 9.15 (C++). `scripts/build_deps.sh [PREFIX] [SRCDIR]` clones and
builds all of them. Then

```
make DEPS=/path/to/install LEMON=/path/to/lemon/include
```

The defaults are `DEPS=$HOME/lpdeps/install` and `LEMON=$HOME/lpdeps/src/lemon/include`, matching
`scripts/build_deps.sh ~/lpdeps/install ~/lpdeps/src`, so a plain `make` works after that.

OR-Tools is detected automatically if `$(DEPS)/ortools/include/ortools/pdlp/...` exists.

Gurobi (optional, commercial; free academic licenses) is detected if `gurobi_c.h` is in
`$(GUROBI_HOME)/include` (default `/usr/include`); `make` then also builds `bin/lpsolve-grb`, a
separate binary so that libgurobi does not change the memory baseline of `bin/lpsolve`. The paper
used Gurobi 13.0.3 (`-lgurobi130`). Gurobi needs a valid license at run time (e.g. `~/gurobi.lic`);
set `LPSOLVE_GRB_LOG=1` to see Gurobi's log. Everything else runs without Gurobi.

The reference max-flow codes are not redistributed here because of their licenses (HIPR and
pseudoflow: academic, non-commercial; Boykov–Kolmogorov: GPL). `scripts/fetch_flow_codes.sh`
downloads them (default `~/lpdeps/src/flow`, the Makefile's `FLOW`), and `make` then builds
`bin/flow-hipr`, `bin/flow-pseudoflow` and `bin/flow-bk` from our drivers in `src/flowref/`, which
include the original sources unchanged. Each is a separate binary (GPL and academic licenses do not
mix), with the authors' optimization flags.

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
  Max flow also: HIPR 4.2 (Cherkassky–Goldberg highest-label push–relabel), pseudoflow 3.23
  (Chandran–Hochbaum, lowest label, FIFO buckets) and Boykov–Kolmogorov maxflow 3.04, as
  `classic:hipr`, `classic:pseudoflow`, `classic:bk`. Their reported value comes from phase 1 (HIPR,
  pseudoflow; as LEMON Preflow) or the full run (BK); RLG's 2^40 "infinite" capacities exceed the
  32-bit capacities of HIPR and pseudoflow, so `run.py` hands them a copy with those capacities
  replaced by the capacity of a finite s–t cut, which leaves the max-flow value unchanged.
* **Simplex (sparse LU)**: HiGHS dual/primal (hyper-sparse revised simplex), Clp dual/primal,
  GLPK primal/dual, OR-Tools GLOP dual/primal.
* **Interior point**: HiGHS IPX (with and without crossover), Clp barrier (native Cholesky,
  `clp-barrier`; CHOLMOD with AMD ordering, `clp-barrier-chol`), GLPK interior.
* **First-order**: HiGHS PDLP (cuPDLP-C), OR-Tools PDLP, own PDHG.
* **Commercial reference** (Gurobi 13): `gurobi-primal`, `gurobi-dual`, `gurobi-barrier` (no
  crossover), `gurobi-barrier-cross`, `gurobi-pdhg` (CPU, no crossover).
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

Runs must not overlap with other heavy work on the machine. `LPBENCH_CPU=k` pins every solver
process to CPU `k` (recorded in each result as `cpu`); on a hybrid CPU pick a performance core,
ideally at a fixed frequency, and keep the driver itself on other cores. The paper's runs:

```
export LPBENCH_CPU=2                      # P-core held at 3.0 GHz
for s in 1 2 3; do
  python3 scripts/run.py --plan pilot --tl 60 --seed $s --out results/pilot-pin-s$s.jsonl
  python3 scripts/analyze.py results/pilot-pin-s$s.jsonl results/pin-s$s 60
done
python3 scripts/seed_stats.py results/pin-s1/summary.md results/pilot-pin-s1.jsonl \
  results/pin-s2/summary.md results/pilot-pin-s2.jsonl results/pin-s3/summary.md results/pilot-pin-s3.jsonl > results/seed_stats.txt
scripts/fetch_zenodo_roads.sh             # ALL=1 also converts the full USA graph
python3 scripts/run.py --plan roads --tl 60 --out results/roads-local.jsonl
python3 scripts/run.py --plan roads --tl 600 --out results/roads-600-local.jsonl \
  --methods classic:dijkstra,classic:lemon-dijkstra,classic:lemon-ns,own-pdhg,ortools-pdlp,highs-pdlp,clp-barrier,clp-primal
python3 scripts/run.py --plan roads-usa --tl 600 --out results/roads-usa-local.jsonl \
  --methods classic:dijkstra,classic:lemon-dijkstra,classic:lemon-ns
python3 scripts/plot_roads.py results/roads/roads.png \
  results/roads-local.jsonl results/roads-600-local.jsonl results/roads-usa-local.jsonl
python3 scripts/paper_tables.py results/pin-s1/summary.md > bpa.tex
```

`results/roads/summary.md` is `analyze.py` run on the concatenation of `roads-local.jsonl` and
`roads-usa-local.jsonl` (named `roads-all-60.jsonl` in its header).
A solve counts as solved only if it reports optimal, finishes within the time limit, and is within
relative error 1e-3 of the combinatorial optimum. Net memory subtracts the median of nine
empty-process runs of the same binary.

Backends `highs-*` and `clp-*` accept a `-nopre` suffix (presolve off), e.g. `highs-dual-nopre`.
OR-Tools backends are dispatched to `bin/lpsolve-ort`, Gurobi backends to `bin/lpsolve-grb`, automatically.

## Replicating the paper

For the RCR review. The paper uses release v1.2
(`git checkout v1.2`). Later commits add the timer checks below
and documentation; the benchmark code is unchanged.

1. `scripts/build_deps.sh ~/lpdeps/install ~/lpdeps/src`,
   `scripts/fetch_flow_codes.sh`, then `make`. Gurobi is optional;
   without it, every result except the Gurobi records reproduces.
2. Run the commands in [Running](#running). Solver time is about
   3.5 h per seed (0.5 h of it Gurobi) and about 2 h for the roads.
3. Compare with the stored files, or regenerate the paper's
   tables and figures from the stored raw results directly:

| paper | file | produced by |
|---|---|---|
| Figures 1–3 | `results/pin-s1/{sp,mf,mcf}.png` | `analyze.py` |
| Table 3 (bytes/arc) | `bpa.tex` | `paper_tables.py` |
| ranges over seeds | `results/seed_stats.txt` | `seed_stats.py` |
| Figure 4 (roads) | `results/roads/roads.png` | `plot_roads.py` |
| hop depth | `results/depth/` | `depth_experiment.sh` |
| Figure 5 (Kaczmarz) | `results/kaczmarz_convergence.png` | `plot_traces.py` |
| Figure 6 (bisection) | `results/exploratory/` | `exploratory/plot_bisect.py` |

Memory should reproduce within a few percent on any Linux machine;
absolute times depend on the hardware, but ratios between methods
should hold on a similar CPU.

### Timer validation

`make timer-test` runs two checks; neither needs Gurobi or runs a
solver.

* `bin/timer-test` (`src/timer_test.cpp`) checks the timer the
  solvers use (`std::chrono::steady_clock`, the same `secs()` as
  `lpsolve.cpp` and `classic.cpp`): resolution, monotonicity over
  10^7 reads, and agreement with `CLOCK_MONOTONIC_RAW` over 1 ms to
  5 s (tolerance 10 µs + 100 ppm).
* `scripts/check_timers.py` cross-checks the stored results: each
  record's `t_solve` ≤ `t_total`, measured inside the solver, must
  not exceed `wall`, measured independently by `memwrap.py` in the
  parent process.

On the stored results of the paper (three seeds and the roads),
no record fails. For LP solvers `t_total` includes reading the
instance, and `wall` exceeds it by at most 1.9% on runs of 1 s or
more (median 0.05%): process start-up, exit and the 5 ms polling
of `memwrap.py`. The classical codes start their clock after
reading the instance, so for them `wall` also includes the read,
about 10 s on the full USA graph. On a second Linux machine
(2026-11-24) the timer self-test passed with 10 ns resolution and
`steady_clock` within 15 ppm of the raw clock.

## Results (in `results/`)

Current results, used in the paper (Intel Core Ultra 5 125H, every solve pinned to one P-core at
a fixed 3.0 GHz, 60 s per solve, tolerance 1e-6, single thread):

* `pilot-pin-s{1,2,3}.jsonl` / `.log` (raw, three seeds) and `pin-s{1,2,3}/` (summary tables and
  plots; the paper's figures and table use seed 1); `seed_stats.txt` (ranges over seeds).
* `roads-local`, `roads-600-local` (600 s limit, selected methods), `roads-usa-local` (full USA
  graph, classical codes) and `roads/` (summary, `roads.png`).
* The GLPK records (GLPK 5.0) and `clp-barrier-chol` (Clp barrier with CHOLMOD/AMD) were added to
  these files by a later rerun (`rerun-glpk5-chol.log`), and the HIPR, pseudoflow and BK records by
  `flowref-pin.log`; the GLPK 4.60 records they replace are in `superseded/`. The Gurobi records
  were added by `rerun-gurobi.log` (same machine and pinning).
* Net memory subtracts `base_kb`, the peak of an empty run of the same binary with the same shared
  libraries, which `run.py` stores in every record. Records written before it did so were annotated
  with the baseline measured for their build: 7,744 kB for `lpsolve` with the native-ordering Clp
  and GLPK 4.60 (the libraries kept in `~/lpdeps/superseded`), 9,192 kB for the CHOLMOD-linked
  build with GLPK 5.0 (the GLPK and `clp-barrier-chol` records), each the median of nine runs.
  Linking CHOLMOD and the BLAS raises the empty-process peak by 1.4 MB.
* `depth/`: hop depth of the shortest-path trees and the aspect-ratio experiment
  (`scripts/hop_depth.py`, `scripts/depth_experiment.sh`).

Earlier runs, kept for comparison and superseded by the above:

* `pilot-local.jsonl` and `local/`: the same plan on the same machine, one seed, not pinned;
  `local/NOTES.md` lists the differences from the first pilot; `local/presolve_check.txt`.
* `FINDINGS.md`, `summary.md`, `sp.png`, `mf.png`, `mcf.png`, `pilot-all.jsonl`, `pilot.jsonl`,
  `pilot-ort.jsonl`, `small.jsonl`, `presolve_check.txt`: the first pilot, run on a 2-core VM with one seed.

In all runs a method is skipped on larger instances of a family once it fails on a smaller one,
and SKM runs only on the two smallest sizes.

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
