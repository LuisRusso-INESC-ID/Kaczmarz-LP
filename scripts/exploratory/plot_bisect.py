#!/usr/bin/env python3
"""EXPLORATORY: primal-dual bisection with Kaczmarz (own-bisect-kz) -- the bracket [L, U] around
the optimum over time, against PDHG.  Writes results/exploratory/bisect_convergence.png."""
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
opt = {}
for line in open(os.path.join(ROOT, "results", "pilot-all.jsonl")):
    r = json.loads(line)
    for name, _, _, ref in INST:
        if r["instance"] == name and r["method"] == ref and r.get("status") in ("ok", "optimal"):
            opt[name] = r["obj"]

fig, axes = plt.subplots(2, 3, figsize=(14, 7.5), sharex="col")
for c, (name, title, sense, _) in enumerate(INST):
    v = opt[name]
    for fname, tagname, alpha_line, ls in (("own-bisect-kz-polyak", "Polyak dual", 0.3, ":"),
                                            ("own-bisect-kz", "smoothed Hildreth dual", 0.55, "--"),
                                            ("own-bisect-kz-prox", "proximal dual (multipliers)", 1.0, "-")):
        p = os.path.join(TD, f"{name}__{fname}.csv")
        if not tracexists(p): continue
        B = list(csv.DictReader(opentrace(p)))
        t = [max(float(r["time"]), 1e-4) for r in B]
        U = [sense * float(r["pobj"]) / v for r in B]
        L = [sense * float(r["extra"]) / v for r in B]
        primal = [u if math.isfinite(u) else float("nan") for u in U]
        dual = [l if math.isfinite(l) else float("nan") for l in L]
        axes[0][c].plot(t, primal, color="#c2185b", lw=1.6, ls=ls, alpha=alpha_line, label=f"{tagname}: primal bound")
        axes[0][c].plot(t, dual, color="#2e7d32", lw=1.6, ls=ls, alpha=alpha_line, label=f"{tagname}: dual bound")
        if fname == "own-bisect-kz-prox":
            lo = [min(a, b) for a, b in zip(primal, dual)]; hi = [max(a, b) for a, b in zip(primal, dual)]
            axes[0][c].fill_between(t, lo, hi, color="#bbbbbb", alpha=.3, lw=0)
        g = [float(r["gap"]) for r in B]
        axes[1][c].plot(t, [x if math.isfinite(x) else float("nan") for x in g], color="#6a1b9a", lw=1.6, ls=ls, alpha=alpha_line, label=f"{tagname}: relative gap")
    P = list(csv.DictReader(opentrace(os.path.join(TD, f"{name}__own-pdhg.csv"))))
    tp = [max(float(r["time"]), 1e-4) for r in P]
    axes[0][c].plot(tp, [sense * float(r["pobj"]) / v for r in P], color="#2b7bba", lw=1.2, label="PDHG objective")
    axes[0][c].axhline(1.0, color="k", lw=0.8)
    axes[0][c].set_ylim(-0.1, 2.6)
    axes[0][c].set_title(f"{title}\noptimum = {v:.6g}  (split targets ¼/¾)", fontsize=10.5)
    axes[1][c].plot(tp, [max(float(r["gap"]), 1e-12) for r in P], color="#2b7bba", lw=1.2, label="PDHG relative duality gap")
    for r in range(2):
        axes[r][c].set_xscale("log"); axes[r][c].grid(alpha=.3, which="both", lw=.4)
    axes[1][c].set_yscale("log"); axes[1][c].set_xlabel("time (s)")
axes[0][0].set_ylabel("bound / optimum")
axes[1][0].set_ylabel("relative gap")
h, l = [], []
for ax in axes.flat:
    for hh, ll in zip(*ax.get_legend_handles_labels()):
        if ll not in l: h.append(hh); l.append(ll)
fig.legend(h, l, loc="lower center", ncol=4, fontsize=8.5, frameon=False)
fig.tight_layout(rect=(0, 0.1, 1, 1))
out = os.path.join(ROOT, "results", "exploratory", "bisect_convergence.png")
fig.savefig(out, dpi=140); print(out)
