#!/usr/bin/env python3
"""Road-network plot: net peak memory and solve time vs. arcs, combining the 60 s and 600 s runs.

usage: plot_roads.py out.png results/roads-local.jsonl results/roads-600-local.jsonl results/roads-usa-local.jsonl
A run counts as solved if its objective is within 1e-3 of Dijkstra's. Points from the 600 s runs
that were not solved within 60 s are drawn hollow. Net memory subtracts the empty-process baseline
of each binary (measured as in analyze.py).
"""
import json, os, statistics, subprocess, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out, f60, f600, fusa = sys.argv[1:5]
METHODS = [("classic:dijkstra", "Dijkstra", "k", "o"), ("classic:lemon-dijkstra", "LEMON Dijkstra", "0.45", "o"),
           ("classic:lemon-ns", "LEMON network simplex", "0.6", "s"), ("glpk-ipt", "GLPK interior point", "C4", "D"),
           ("clp-barrier", "Clp barrier", "C3", ">"), ("clp-primal", "Clp primal", "C1", "<"),
           ("highs-primal", "HiGHS primal", "C0", "^"), ("ortools-glop-primal", "GLOP primal", "C9", "v"),
           ("own-pdhg", "own PDHG", "C8", "P"), ("ortools-pdlp", "OR-Tools PDLP", "C6", "*")]


def baseline(cmd, reps=9):   # median of several empty runs, as in analyze.py
    v = []
    for _ in range(reps):
        o = subprocess.run(cmd, capture_output=True, text=True).stdout
        v += [float(t.split("=")[1]) for t in o.split() if t.startswith("hwm_end_kb=")]
    return statistics.median(v) if v else 0.0


probe = os.path.join(ROOT, "data", "sp-grid-3x3-s1.gr")
B = {"lp": baseline([os.path.join(ROOT, "bin", "lpsolve"), "build-only", probe]),
     "ort": baseline([os.path.join(ROOT, "bin", "lpsolve-ort"), "build-only", probe]),
     "classic": baseline([os.path.join(ROOT, "bin", "classic"), "dijkstra", probe])}
load = lambda p: [json.loads(l) for l in open(p)]
R60, R600, RU = load(f60), load(f600), load(fusa)
ref = {r["instance"]: r["obj"] for r in R60 + RU if r["method"] == "classic:dijkstra"}
arcs = {}
for name in ref:
    with open(os.path.join(ROOT, "data", name + ".gr")) as f:
        arcs[name] = next(int(l.split()[3]) for l in f if l.startswith("p "))


def solved(r):
    return r.get("status") in ("ok", "optimal") and isinstance(r.get("obj"), float) \
        and abs(r["obj"] - ref[r["instance"]]) <= 1e-3 * abs(ref[r["instance"]])


def net_mb(r):
    b = B["classic"] if r["method"].startswith("classic:") else B["ort"] if r["method"].startswith("ortools-") else B["lp"]
    return (r["hwm_end_kb"] - b) / 1024


fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
for m, label, col, mk in METHODS:
    pts = {}
    for r in R60 + RU:
        if r["method"] == m and solved(r):
            pts[r["instance"]] = (r, False)
    for r in R600:
        if r["method"] == m and solved(r) and r["instance"] not in pts:
            pts[r["instance"]] = (r, True)
    if not pts:
        continue
    xs = sorted(pts, key=lambda i: arcs[i])
    x = [arcs[i] for i in xs]
    for k, (val, ylab) in enumerate([(net_mb, "net peak memory (MB)"), (lambda r: r["t_solve"], "solve time (s)")]):
        y = [val(pts[i][0]) for i in xs]
        ax[k].plot(x, y, color=col, lw=1.2, label=label if k == 0 else None)
        for xi, yi, i in zip(x, y, xs):
            ax[k].plot(xi, yi, mk, color=col, ms=5, mfc="none" if pts[i][1] else col)
        ax[k].set_xscale("log"); ax[k].set_yscale("log"); ax[k].set_xlabel("input arcs m")
        ax[k].set_ylabel(ylab); ax[k].grid(True, which="both", alpha=0.25)
fig.legend(loc="lower center", ncol=5, fontsize=8, frameon=False)
fig.tight_layout(rect=(0, 0.1, 1, 1))
fig.savefig(out, dpi=150)
print("wrote", out)
