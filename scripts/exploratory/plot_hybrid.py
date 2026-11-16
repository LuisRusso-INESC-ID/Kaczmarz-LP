#!/usr/bin/env python3
"""EXPLORATORY (not part of the benchmark): hybrid block-gradient + repair variants vs proximal
Hildreth and PDHG.  Writes results/exploratory/hybrid_convergence.png."""
import csv, json, math, os
import gzip
import matplotlib
matplotlib.use("Agg")

def opentrace(p):
    """Open a trace CSV, transparently using the gzipped copy (p + '.gz') if present."""
    if os.path.exists(p): return open(p)
    return gzip.open(p + ".gz", "rt")

def tracexists(p):
    return os.path.exists(p) or os.path.exists(p + ".gz")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TD = os.path.join(ROOT, "results", "traces")
INST = [("sp-grid-32x32-s1", "shortest paths, 32×32 grid", +1, "classic:dijkstra"),
        ("mf-rmf-8x8-s1", "max flow, RMF 8×8×8", -1, "classic:lemon-preflow"),
        ("mcf-grid-10x10x4x5-s1", "max concurrent flow, 10×10 grid", -1, "highs-ipx-cross")]
METH = [("own-pdhg", "PDHG (reference)", "#2b7bba", "-"),
        ("own-prow", "Proximal Hildreth row action", "#e08214", "-"),
        ("own-hybrid-mem", "Hybrid, memoryless Kaczmarz repair", "#c2185b", "-"),
        ("own-hybrid-prox", "Hybrid, Hildreth (proximal) repair", "#2e7d32", "-"),
        ("own-hybrid-cg", "Hybrid, CG repair (64-row balls, doubling)", "#6a1b9a", "-"),
        ("own-hybrid-cg-full", "Hybrid, CG repair (whole-matrix blocks)", "#6a1b9a", "--")]
opt = {}
for line in open(os.path.join(ROOT, "results", "pilot-all.jsonl")):
    r = json.loads(line)
    for name, _, _, ref in INST:
        if r["instance"] == name and r["method"] == ref and r.get("status") in ("ok", "optimal"):
            opt[name] = r["obj"]

def load(p, maxpts=3000):
    rows = list(csv.DictReader(opentrace(p)))
    if len(rows) > maxpts:
        step = len(rows) / maxpts
        rows = [rows[int(k * step)] for k in range(maxpts)] + [rows[-1]]
    return rows

fig, axes = plt.subplots(3, 3, figsize=(14, 10), sharex="col")
for c, (name, title, sense, _) in enumerate(INST):
    v = opt[name]
    for fname, label, color, ls in METH:
        p = os.path.join(TD, f"{name}__{fname}.csv")
        if not tracexists(p): continue
        R = load(p)
        t = [max(float(r["time"]), 1e-4) for r in R]
        obj = [sense * float(r["pobj"]) for r in R]
        axes[0][c].plot(t, [o / v for o in obj], color=color, ls=ls, lw=1.4, label=label)
        axes[1][c].plot(t, [max(float(r["rp"]), 1e-12) for r in R], color=color, ls=ls, lw=1.4)
        if "hybrid" in fname:
            blocks = [max(float(r["iter"]), 1) for r in R]
            axes[2][c].plot(t, [float(r["extra"]) / b for r, b in zip(R, blocks)], color=color, lw=1.4)
    axes[0][c].axhline(1.0, color="k", lw=0.8)
    axes[0][c].set_title(f"{title}\noptimum = {v:.6g}", fontsize=11)
    axes[0][c].set_ylim(-0.1, 2.0)
    for r in range(3):
        axes[r][c].set_xscale("log"); axes[r][c].grid(alpha=.3, which="both", lw=.4)
        if r > 0: axes[r][c].set_yscale("log")
    axes[2][c].set_xlabel("time (s)")
axes[0][0].set_ylabel("objective / optimum")
axes[1][0].set_ylabel("primal infeasibility\n(max violation / max |rhs|)")
axes[2][0].set_ylabel("repair work per block step\n(row-ops; CG: rows x iterations)")
h, l = axes[0][0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, fontsize=9, frameon=False)
fig.tight_layout(rect=(0, 0.06, 1, 1))
out = os.path.join(ROOT, "results", "exploratory", "hybrid_convergence.png")
fig.savefig(out, dpi=140); print(out)
