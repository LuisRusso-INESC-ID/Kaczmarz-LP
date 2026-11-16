#!/usr/bin/env python3
"""Plot convergence traces (LPTRACE CSVs) of the Kaczmarz-style prototypes against PDHG,
relative to the known optimum (Dijkstra / Preflow / HiGHS IPX+crossover)."""
import csv, glob, json, math, os, sys
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

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TD = os.path.join(ROOT, "results", "traces")
INST = [("sp-grid-32x32-s1", "shortest paths, 32×32 grid", +1, "classic:dijkstra"),
        ("mf-rmf-8x8-s1", "max flow, RMF 8×8×8", -1, "classic:lemon-preflow"),
        ("mcf-grid-10x10x4x5-s1", "max concurrent flow, 10×10 grid", -1, "highs-ipx-cross")]
METH = [("pdhg", "own-pdhg", "PDHG (reference)", "#2b7bba", "-"),
        ("prow", "own-prow", "Proximal Hildreth row action", "#e08214", "-"),
        ("skm", "own-skm", "SKM on optimality system", "#c2185b", "-"),
        ("skm-gap1", "own-skm-gap1", "SKM, gap row every step", "#c2185b", ":"),
        ("skm-scaled", "own-skm-scaled", "SKM, equilibrated system", "#7b1fa2", "--"),
        ("skm-balanced", "own-skm-balanced", "SKM, balanced gap-row metric", "#2e7d32", "-")]

opt = {}
for line in open(os.path.join(ROOT, "results", "pilot-all.jsonl")):
    r = json.loads(line)
    for name, _, _, ref in INST:
        if r["instance"] == name and r["method"] == ref and r.get("status") in ("ok", "optimal"):
            opt[name] = r["obj"]

def load(path, maxpts=3000):
    rows = list(csv.DictReader(opentrace(path)))
    if len(rows) > maxpts:   # log-spaced thinning in iteration index
        idx = sorted({min(len(rows) - 1, int(round(10 ** (k * math.log10(len(rows)) / maxpts)))) for k in range(maxpts + 1)} | {0})
        rows = [rows[i] for i in idx]
    return rows

fig, axes = plt.subplots(4, 3, figsize=(14, 13), sharex="col")
for c, (name, title, sense, _) in enumerate(INST):
    v = opt[name]
    for tag, fname, label, color, ls in METH:
        p = os.path.join(TD, f"{name}__{fname}.csv")
        if not tracexists(p):
            continue
        R = load(p)
        t = [max(float(r["time"]), 1e-4) for r in R]
        obj = [sense * float(r["pobj"]) for r in R]
        axes[0][c].plot(t, [o / v for o in obj], color=color, ls=ls, lw=1.4, label=label)
        axes[1][c].plot(t, [max(abs(o - v) / abs(v), 1e-12) for o in obj], color=color, ls=ls, lw=1.4)
        axes[2][c].plot(t, [max(float(r["rp"]), 1e-12) for r in R], color=color, ls=ls, lw=1.4)
        if tag.startswith("skm"):
            # SKM logs |g|/(1+|pobj|); renormalize like PDHG: |g|/(1+|pobj|+|dobj|), dobj = pobj - g
            gap = []
            for r in R:
                po = float(r["pobj"]); g = float(r["gap"]) * (1 + abs(po))
                gap.append(max(g / (1 + abs(po) + abs(po - g)), 1e-12))
        else:
            gap = [max(float(r["gap"]), 1e-12) for r in R]
        axes[3][c].plot(t, gap, color=color, ls=ls, lw=1.4)
    axes[0][c].axhline(1.0, color="k", lw=0.8)
    axes[0][c].set_title(f"{title}\noptimum = {v:.6g}", fontsize=11)
    axes[0][c].set_ylim(-0.1, 1.3)
    for r in range(4):
        axes[r][c].set_xscale("log"); axes[r][c].grid(alpha=.3, which="both", lw=.4)
        if r > 0: axes[r][c].set_yscale("log")
    axes[3][c].set_xlabel("time (s)")
axes[0][0].set_ylabel("objective / optimum")
axes[1][0].set_ylabel("relative objective error")
axes[2][0].set_ylabel("primal infeasibility\n(max violation / max |rhs|)")
axes[3][0].set_ylabel("relative duality gap")
h, l = [], []
for ax in axes.flat:
    for hh, ll in zip(*ax.get_legend_handles_labels()):
        if ll not in l: h.append(hh); l.append(ll)
fig.legend(h, l, loc="lower center", ncol=3, fontsize=9, frameon=False)
fig.tight_layout(rect=(0, 0.045, 1, 1))
out = os.path.join(ROOT, "results", "kaczmarz_convergence.png")
fig.savefig(out, dpi=140)
print(out)
