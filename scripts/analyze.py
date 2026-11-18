#!/usr/bin/env python3
"""Summarize a results .jsonl: correctness vs reference, time and peak-memory tables, plots.

usage: analyze.py results/pilot.jsonl [outdir]
Peak memory = VmHWM reported by the solving process itself (exact per exec'd process).
'Net' memory subtracts the empty-process baseline of the same binary (libraries, runtime).
"""
import statistics, json, os, subprocess, sys, math, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
res_path = sys.argv[1]
outdir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "results")
TL = float(sys.argv[3]) if len(sys.argv) > 3 else 60.0
SOLVED_ERR = 1e-3   # relative objective error accepted as "solved" (first-order methods stop at 1e-6 *KKT* error)
os.makedirs(outdir, exist_ok=True)

R = [json.loads(l) for l in open(res_path)]
inst_info = {}
for name in {r["instance"] for r in R}:
    for ext in ("gr", "max", "mcf"):
        p = os.path.join(ROOT, "data", f"{name}.{ext}")
        if os.path.exists(p):
            with open(p) as f:
                for line in f:
                    if line.startswith("p "):
                        t = line.split(); inst_info[name] = {"n": int(t[2]), "m": int(t[3]), "k": int(t[4]) if len(t) > 4 else 0}
                        break

def baseline(cmd, reps=9):
    # median of several empty runs: a single run varies by about 0.1 MB, which is up to ~30 B/arc
    # on the smallest instances
    vals = []
    for _ in range(reps):
        out = subprocess.run(cmd, capture_output=True, text=True).stdout
        vals += [float(t.split("=")[1]) for t in out.split() if t.startswith("hwm_end_kb=")]
    return statistics.median(vals) if vals else 0.0
probe = os.path.join(ROOT, "data", "sp-grid-3x3-s1.gr")
if not os.path.exists(probe):
    with open(probe, "w") as f:
        subprocess.run([os.path.join(ROOT, "bin", "gen"), "sp-grid", "3", "3", "1"], stdout=f)
BASE = {"lp": baseline([os.path.join(ROOT, "bin", "lpsolve"), "build-only", probe]),
        "lp-ort": baseline([os.path.join(ROOT, "bin", "lpsolve-ort"), "build-only", probe]) if os.path.exists(os.path.join(ROOT, "bin", "lpsolve-ort")) else 0.0,
        "classic": baseline([os.path.join(ROOT, "bin", "classic"), "dijkstra", probe])}

# reference objective
REF_METHOD = {"sp": "classic:dijkstra", "mf": "classic:lemon-preflow", "mcf": "highs-ipx-cross"}
ref = {}
for r in R:
    kind = r["family"].split("-")[0]
    if r["method"] == REF_METHOD[kind] and r.get("status") in ("ok", "optimal"):
        ref[r["instance"]] = r["obj"]
for r in R:   # fall back for mcf: any exact simplex optimum
    if r["instance"] not in ref and r["method"] in ("clp-dual", "highs-dual", "highs-primal", "ortools-glop", "ortools-glop-primal") and r.get("status") == "optimal":
        ref[r["instance"]] = r["obj"]

def relerr(r):
    if r["instance"] not in ref or "obj" not in r or not isinstance(r["obj"], float) or math.isnan(r["obj"]):
        return float("nan")
    v = ref[r["instance"]]
    return abs(r["obj"] - v) / max(1.0, abs(v))

def solved(r):
    if r.get("status") not in ("optimal", "ok"):
        return False
    # a solve must finish within the time limit; GLPK's interior-point method has no time-limit
    # parameter and is stopped only by memwrap.py's timeout (1.5 TL + 30 s)
    if r.get("t_solve", 0) > TL:
        return False
    if r["method"].startswith("classic:gk"):
        return True
    e = relerr(r)
    return not math.isnan(e) and e <= SOLVED_ERR

rows = []
for r in R:
    if r["method"] == "build-only":
        continue
    info = inst_info.get(r["instance"], {})
    kind = r["family"].split("-")[0]
    binkind = "classic" if r["method"].startswith("classic:") else ("lp-ort" if r["method"].startswith("ortools-") else "lp")
    hwm = r.get("hwm_end_kb", float("nan"))
    rows.append(dict(instance=r["instance"], family=r["family"], kind=kind, method=r["method"],
                     n=info.get("n"), m=info.get("m"), status=r.get("status"), solved=solved(r),
                     relerr=relerr(r), t=r.get("t_solve", float("nan")), hwm_mb=hwm / 1024,
                     net_mb=(hwm - BASE[binkind]) / 1024, nnz=r.get("nnz")))

# ---------------- markdown summary
fams = sorted({x["family"] for x in rows}, key=lambda f: ("sp mf mcf".split().index(f.split("-")[0]), f))
methods_order = sorted({x["method"] for x in rows}, key=lambda m: (not m.startswith("classic"), m))
md = [f"# Results: {os.path.basename(res_path)}\n",
      f"Process baselines (empty run): lpsolve {BASE['lp']/1024:.1f} MB, lpsolve-ort {BASE['lp-ort']/1024:.1f} MB, classic {BASE['classic']/1024:.1f} MB. "
      "'net' = peak RSS minus baseline. Cells: solve time (s) / net peak memory (MB); "
      "'TL' time limit, 'x' failed or wrong, '-' skipped after failing a smaller instance.\n"]
for fam in fams:
    F = [x for x in rows if x["family"] == fam]
    insts = sorted({(x["m"] or 0, x["instance"]) for x in F})
    md.append(f"\n## {fam}\n")
    md.append("| method | " + " | ".join(f"m={m:,}" for m, _ in insts) + " |")
    md.append("|---|" + "---|" * len(insts))
    for meth in methods_order:
        cells = []; any_ = False
        for m, inst in insts:
            x = next((y for y in F if y["instance"] == inst and y["method"] == meth), None)
            if x is None: cells.append(""); continue
            any_ = True
            st = x["status"] or ""
            if st.startswith("skipped"): cells.append("-")
            elif x["solved"]: cells.append(f"{x['t']:.3g} / {x['net_mb']:.1f}")
            elif st in ("time_limit", "killed(timeout)") or (x["t"] == x["t"] and x["t"] >= 0.95 * TL):
                cells.append(f"TL / {x['net_mb']:.0f}")
            else:
                e = x["relerr"]
                cells.append(f"x({st[:12]}{'' if math.isnan(e) else f', err {e:.0e}'})")
        if any_:
            md.append(f"| {meth} | " + " | ".join(cells) + " |")
# memory per arc at the largest instance each method solved
md.append("\n## Net memory per input arc at the largest solved instance (bytes/arc)\n")
md.append("| family | " + " | ".join(methods_order) + " |")
md.append("|---|" + "---|" * len(methods_order))
for fam in fams:
    cells = []
    for meth in methods_order:
        S = [x for x in rows if x["family"] == fam and x["method"] == meth and x["solved"] and x["m"]]
        if not S: cells.append(""); continue
        x = max(S, key=lambda y: y["m"])
        cells.append(f"{x['net_mb']*1024*1024/x['m']:.0f} (m={x['m']//1000}k)")
    md.append(f"| {fam} | " + " | ".join(cells) + " |")
# correctness
bad = [x for x in rows if x["status"] in ("optimal", "ok") and not x["solved"]]
md.append(f"\n## Correctness\n\nReference: {REF_METHOD}. Runs reporting optimal but with relative "
          f"objective error > {SOLVED_ERR:g}: {len(bad)}\n")
for x in bad:
    md.append(f"- {x['instance']} {x['method']}: rel. error {x['relerr']:.2e}")
md.append("\nLargest relative objective error among solved runs, per method: " + ", ".join(
    f"{mm} {max(x['relerr'] for x in rows if x['method']==mm and x['solved'] and x['relerr']==x['relerr']):.0e}"
    for mm in methods_order if not mm.startswith("classic:gk") and any(x['method']==mm and x['solved'] and x['relerr']==x['relerr'] for x in rows)))
gk = [x for x in rows if x["method"].startswith("classic:gk")]
if gk:
    md.append("\nGarg–Könemann (approximate) relative gap to the LP optimum: " +
              ", ".join(f"{x['instance'].split('-s')[0]} {x['method'][8:]}: {x['relerr']:.1e}" for x in gk if not math.isnan(x['relerr'])))
open(os.path.join(outdir, "summary.md"), "w").write("\n".join(md) + "\n")

# ---------------- plots
PLOT_METHODS = ["classic:dijkstra", "classic:lemon-preflow", "classic:lemon-ns", "classic:gk:0.05",
                "highs-dual", "highs-ipx-cross", "highs-pdlp", "clp-dual", "clp-barrier", "glpk-primal",
                "ortools-glop", "ortools-pdlp", "own-pdhg", "own-prow"]
STYLE = {m: (c, mk) for m, c, mk in zip(PLOT_METHODS,
         ["#111111", "#111111", "#666666", "#111111", "#1f77b4", "#2ca02c", "#17becf", "#ff7f0e", "#d62728",
          "#8c564b", "#9467bd", "#e377c2", "#bcbd22", "#7f7f7f"],
         ["o", "o", "s", "o", "^", "v", "D", "<", ">", "p", "h", "*", "P", "X"])}
for kind in ("sp", "mf", "mcf"):
    K = [x for x in rows if x["kind"] == kind]
    if not K: continue
    kf = sorted({x["family"] for x in K})
    fig, axes = plt.subplots(2, len(kf), figsize=(4.2 * len(kf), 7), squeeze=False)
    for col, fam in enumerate(kf):
        for row, (key, lab) in enumerate((("net_mb", "net peak memory (MB)"), ("t", "solve time (s)"))):
            ax = axes[row][col]
            for meth in PLOT_METHODS:
                S = sorted([x for x in K if x["family"] == fam and x["method"] == meth and x["solved"] and x["m"]], key=lambda y: y["m"])
                if not S: continue
                c, mk = STYLE[meth]
                ax.plot([x["m"] for x in S], [max(x[key], 1e-3) for x in S], marker=mk, color=c, lw=1.2, ms=4, label=meth)
            ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(alpha=.3, which="both", lw=.4)
            if row == 0: ax.set_title(fam)
            if row == 1: ax.set_xlabel("input arcs m")
            if col == 0: ax.set_ylabel(lab)
    h, l = [], []
    for ax in axes.flat:
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in l: h.append(hh); l.append(ll)
    if l:
        fig.legend(h, l, loc="lower center", ncol=min(7, len(l)), fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(outdir, f"{kind}.png"), dpi=150)
    plt.close(fig)
print(open(os.path.join(outdir, "summary.md")).read())
