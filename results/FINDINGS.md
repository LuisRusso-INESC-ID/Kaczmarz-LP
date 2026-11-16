# Pilot findings (60 s per solve, one seed, synthetic families)

Full tables: `summary.md`; plots: `sp.png`, `mf.png`, `mcf.png`; raw data: `pilot-all.jsonl`.
Memory is net peak RSS (process VmHWM minus the empty-process baseline of the same binary).

## Memory

* **Classical codes** need 23–28 B per input arc (own Dijkstra), about 55 (LEMON Dijkstra),
  about 60 (LEMON Preflow) and 125–140 (LEMON network simplex).
* **General sparse LP solvers** need 440–1,600 B/arc for simplex and 500–3,100 B/arc for
  interior point on shortest path and max flow. For simplex that is 16–60× Dijkstra and
  7–27× Preflow, on the same graph.
* **First-order methods** are the cheapest exact-ish LP option: own PDHG 184–219 B/arc,
  OR-Tools PDLP 213–337 B/arc, i.e. roughly 8× Dijkstra and 3× Preflow. They only store A
  plus a few vectors. HiGHS PDLP does not show this advantage (810–1,060 B/arc) because
  HiGHS keeps its model and presolve copies around.
* **Presolve is a large share of HiGHS's footprint.** On mf-vision 512² (m = 1.31M) HiGHS
  peaks at 1.0 GB in every mode; with presolve off this drops to 0.75 GB (dual simplex),
  0.72 GB (IPX) and 0.63 GB (PDLP). Clp barrier drops from 0.67 to 0.58 GB
  (`presolve_check.txt`).
* **Fill-in decides interior-point memory.** On planar grids and road-like graphs Clp barrier
  stays near 500 B/arc up to 1M arcs. On random digraphs it reaches 2 KB/arc at 20k arcs
  and exceeds 450 MB at 80k arcs (time limit). GLPK's interior point is the most
  memory-hungry, at 0.8–3 KB/arc.

## Time

* **Shortest paths via LP are 3–4 orders of magnitude slower than Dijkstra.** The best LP
  code on 1M-arc grids, Clp barrier, takes 61 s against Dijkstra's 0.05 s. All dual
  simplex codes hit the 60 s limit by 260k arcs. The primal simplex codes (Clp, GLOP, HiGHS
  on road graphs) still solve 225–320k arcs in 6–42 s but fail at about 1M. LEMON's network
  simplex, the graph-aware simplex on the same LP, sits in between: 8–150× slower than
  Dijkstra but 7–200× faster than the general simplex codes.
* **Max flow is where LP comes closest.** On vision grids at 1.3M arcs, own PDHG takes
  17 s and HiGHS IPX 20 s, against 15 s for LEMON Preflow. LEMON network simplex takes
  1.9 s, though, and push–relabel is not the strongest code for vision instances (see
  next steps). On RMF at 159k arcs, HiGHS IPX takes 3.1 s against 1.1 s for Preflow.
* **Multicommodity flow** (max concurrent, commodities aggregated by source) has no exact
  combinatorial competitor. Garg–Könemann uses 60–90 B/arc and stays within 0.02–0.9% of
  the optimum. The LP grows with K·m (K = 16 source groups at the largest size). There,
  first-order methods use the least memory (3.5–4.5 KB per input arc, about a third of
  simplex) and solve in 21–36 s. GLOP primal is faster (3.6–10 s) but uses 3× more memory.

## Kaczmarz-style methods

* **SKM on the primal–dual optimality system did not converge on any instance**, even with
  2.5k arcs and 60 s. The dense duality-gap row and the conditioning of the stacked system
  make Motzkin relaxation hopeless here. This is a clean negative result.
* **The proximal-point + Hildreth row-action prototype does converge** on max flow,
  multicommodity flow and small shortest-path instances. It has the smallest footprint of
  all LP methods (170–245 B/arc) but is 10–100× slower than PDHG.

## Exploratory: a hybrid feasibility/optimization method (not in the benchmark)

We also tried keeping the objective out of the constraints while still relying on row action:
gradient steps restricted to the null space of the equality rows (local circulations for network
matrices), with violated rows repaired afterwards. With Kaczmarz or Hildreth repairs, the method
stalls far from the optimum: clamping at bounds breaks flow conservation, and restoring it by row
projections is a slow Laplacian smoother (about 10⁶ row projections per step on an RMF instance
with 2,240 arcs). Replacing the repair by a CG solve makes it converge, but only with whole-matrix
blocks, at which point it is plain projected gradient (exact max-flow optimum in 0.4 s; concurrent
flow in 0.2 s; shortest paths still 13% above the optimum after 60 s), and it remains 3–7× slower
than PDHG. Local blocks fail because improving moves need long paths. The experiment isolates the
reason row-action methods struggle on these LPs, conservation repair, and was therefore excluded
from the benchmark. Code: `src/hybrid.h`; plot: `exploratory/hybrid_convergence.png`.

## Exploratory: primal–dual bisection with Kaczmarz (not in the benchmark)

A level t between a primal bound U and a dual bound L is tested by Kaczmarz on the primal slice
{c'x ≤ t} and, interleaved, by Polyak/Kaczmarz steps on the dual level set {LB(y) ≥ t}. Strong
duality makes this a race in which exactly one side can succeed, so no infeasibility certificate is
needed, and finite boxes make every dual iterate a valid bound. The framework is sound, but in 60 s
it closes the bracket only to 47% (shortest paths), 100% (max flow) and 44% (concurrent flow), against
PDHG's 1e-9 in under 2 s. Both inner solvers stall: the primal Kaczmarz cannot route flow to within
1e-3 feasibility while the objective row keeps pushing (its violation plateaus near 5e-3), and the
dual subgradient steps build node potentials slowly. Allowing a looser primal tolerance (1e-2) lets
the race start, but then the primal "bounds" come from 2.7%-infeasible points and overshoot the
optimum (0.79 against 0.655 on concurrent flow). The bisection only reorganizes the work; it does not
remove the slow inner problem, flow conservation by row action.

Overlapping split targets (primal asked to improve its bound by G/4, dual by G/4, recomputed every
round) fix the scheduling flaw of plain bisection, where the level moved only after a success
(2-3 levels in 60 s; now hundreds). In 60 s the primal bound reaches 99.8% of the optimum on max flow
(296,856 vs 297,584) and the concurrent-flow gap falls to 4.9% (from 44%). The dual bounds stay loose
(max flow: 3.4x the optimum), and on shortest paths the primal never improves on its first feasible
point (1.73x). The approximate primal bounds remain a risk: on concurrent flow the reported 0.6616
exceeds the optimum 0.6552 by 1% although the point violates the rows by only 7e-4 of the largest
right-hand side.

Replacing the Polyak dual step by smoothed coordinate ascent (Hildreth row steps on the dual with a
μ/2·‖x‖² smoothing, the true Lagrangian bound reported) tightens the dual bounds sharply in the first
seconds: in 60 s, max flow reaches [295,959, 318,757] (gap 3.7%; the dual bound 7.1% above the
optimum, from 3.4x) and concurrent flow's dual bound 0.705 (7.5% above, from 18%). Shortest paths
barely changes (dual bound at 71% of the optimum). After about 3 s the dual bounds stop improving:
as μ → 0 the coordinate ascent sticks at kinks of the nonsmooth dual, the known failure of single-
coordinate dual ascent on network problems; flooring and periodically re-raising μ helps only a little.
Escaping needs multi-coordinate moves (Bertsekas' relaxation and auction methods, i.e. network
algorithms again) or a proximal term on y (the proximal Hildreth method applied to the dual).

Selecting rows from FIFO queues of violated rows, instead of random sampling, reaches the same
bounds about twice as fast early on and ends somewhat tighter in 60 s: max flow [295,567, 315,310]
(gap 3.2%, from 3.7%), concurrent flow [0.6568, 0.6924] (1.5%; the primal overshoot drops from 1.7% to
0.24%), and on shortest paths the primal bound improves for the first time, 1.48x instead of 1.73x
(gap 34%, from 42%). The plateaus remain, so wasted sampling was a secondary cost; the primary limit is
the conditioning of row action on network matrices.

A proximal term on the dual removes the dual plateau. The proximal step
y⁺ = argmax LB(y) − ‖y − y_k‖²/(2σ) is the method of multipliers: minimize the augmented Lagrangian
over the box by column coordinate descent (sparse, 2–3 rows per column), then update y. Its x_d is
also a nearly feasible primal point, used as a primal candidate. Results in 60 s: max flow solved,
bracket [297,584.009, 297,584.016] (gap 1e-8) in 1.1 s; concurrent flow bracket [0.655237, 0.655423]
(gap 8e-5) after 2 s, then a plateau (an earlier, less cautious σ rule closed it to 7e-7 in 8 s);
shortest paths' dual bound reaches 99.6% of the optimum (9.01M, from 6.6M), while the primal bound
stays at 1.54x, now the bottleneck there. Two safeguards were needed: σ must adapt in both
directions (it otherwise exploded, or collapsed and let x_d diverge on unbounded arcs), and an
approximate primal bound that the valid dual bound crosses must be discarded with a tighter primal
tolerance (otherwise the bracket "closes" on a slightly infeasible point: 0.6616 vs 0.6552). With
this dual, the method is an augmented-Lagrangian solver with coordinate-descent inner steps more than
a Kaczmarz method, and PDHG is still faster (1e-9 in under 2 s on all three).

## Accuracy

Every exact method reproduces the combinatorial optimum. First-order methods stop at a
relative KKT error of 1e-6, which leaves objective errors of up to 2e-4 (OR-Tools PDLP),
3e-5 (HiGHS PDLP) and 5e-6 (own PDHG). In the tables, a relative objective error up to
1e-3 counts as solved.

## Caveats

* The 60 s time limit truncates the large sizes. Results come from a single seed.
* All instance families are synthetic; the DIMACS road graphs could not be downloaded here.
* The own prototypes were compiled with `-march=native`; the libraries use their default
  release flags.
* The OR-Tools runs shared the 2-core machine with the main sweep, so their timings are
  noisier.
* HiGHS PDLP returned status "Unknown" on vision 512² even though its objective was exact.

## Next steps

* Longer time limits (≥ 1 h), 3–5 seeds, and the real road networks
  (`scripts/fetch_dimacs_roads.sh`).
* Stronger max-flow baselines: hi_pr, HPF pseudoflow, Boykov–Kolmogorov / EIBFS for vision
  instances; CS2 for min-cost flow.
* GPU first-order runs (cuPDLP / cuOpt) for the memory-vs-time trade-off.
* Separate solver memory into factorization fill and working storage using solver
  statistics.
