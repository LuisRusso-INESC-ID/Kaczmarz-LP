#!/usr/bin/env python3
"""Seed robustness of the pilot: the ranges quoted in the paper, per seed, and run-to-run spread.

usage: seed_stats.py results/pin-s1/summary.md results/pilot-pin-s1.jsonl [summary2 jsonl2 ...]
  (pairs of analyze.py summary + raw results, one pair per seed)
Prints, per seed, the bytes/arc ranges and ratios over the sp and mf families (from the summary's
"Net memory per input arc" table), the fastest-LP / Dijkstra time ratios, and then the spread
(max/min over seeds) of solve time and peak memory for (family, size, method) solved in all seeds.
"""
import collections, json, math, re, statistics, sys

SPMF = ["sp-grid", "sp-rand", "sp-road", "mf-rlg", "mf-rmf", "mf-vision"]
SIMPLEX = ["highs-dual", "highs-primal", "clp-dual", "clp-primal", "glpk-primal", "glpk-dual",
           "ortools-glop", "ortools-glop-primal"]
IPM = ["highs-ipx", "highs-ipx-cross", "clp-barrier", "glpk-ipt"]   # default configurations
GROUPS = [("own Dijkstra", ["classic:dijkstra"]), ("LEMON Dijkstra", ["classic:lemon-dijkstra"]),
          ("Preflow", ["classic:lemon-preflow"]), ("LEMON NS", ["classic:lemon-ns"]),
          ("HIPR", ["classic:hipr"]), ("pseudoflow", ["classic:pseudoflow"]), ("BK", ["classic:bk"]),
          ("simplex", SIMPLEX), ("interior point", IPM), ("Clp barrier CHOLMOD", ["clp-barrier-chol"]), ("own PDHG", ["own-pdhg"]),
          ("OR-Tools PDLP", ["ortools-pdlp"]), ("HiGHS PDLP", ["highs-pdlp"]),
          ("prox. Hildreth", ["own-prow"])]


def bpa(summary):
    L = open(summary).read().splitlines()
    i = next(k for k, l in enumerate(L) if l.startswith("## Net memory per input arc"))
    rows = [l for l in L[i + 1:] if l.startswith("|")]
    hdr = [c.strip() for c in rows[0].strip("|").split("|")]
    T = {}
    for r in rows[2:]:
        c = [x.strip() for x in r.strip("|").split("|")]
        for h, v in zip(hdr[1:], c[1:]):
            m = re.match(r"(\d+)", v)
            if m:
                T[(c[0], h)] = int(m.group(1))
    return T


def rng(T, methods, fams=SPMF):
    v = [T[(f, m)] for f in fams for m in methods if (f, m) in T]
    return (min(v), max(v)) if v else None


def runs(path):
    R = [json.loads(l) for l in open(path)]
    by = collections.defaultdict(dict)
    for r in R:
        by[r["instance"]][r["method"]] = r
    return by


def ok(r, ref):
    return r.get("status") in ("optimal", "ok") and "t_solve" in r and isinstance(r.get("obj"), float) \
        and abs(r["obj"] - ref) <= 1e-3 * max(1.0, abs(ref))


pairs = list(zip(sys.argv[1::2], sys.argv[2::2]))
per_seed = []
for si, (summ, raw) in enumerate(pairs, 1):
    T = bpa(summ)
    by = runs(raw)
    out = {name: rng(T, ms) for name, ms in GROUPS}
    sd = [T[(f, m)] / T[(f, "classic:dijkstra")] for f in SPMF[:3] for m in SIMPLEX if (f, m) in T]
    sp = [T[(f, m)] / T[(f, "classic:lemon-preflow")] for f in SPMF[3:] for m in SIMPLEX if (f, m) in T]
    fd = [T[(f, m)] / T[(f, "classic:dijkstra")] for f in SPMF[:3] for m in ("own-pdhg", "ortools-pdlp") if (f, m) in T]
    fp = [T[(f, m)] / T[(f, "classic:lemon-preflow")] for f in SPMF[3:] for m in ("own-pdhg", "ortools-pdlp") if (f, m) in T]
    fh = [T[(f, m)] / T[(f, "classic:hipr")] for f in SPMF[3:] for m in ("own-pdhg", "ortools-pdlp") if (f, m) in T and (f, "classic:hipr") in T]
    sh = [T[(f, m)] / T[(f, "classic:hipr")] for f in SPMF[3:] for m in SIMPLEX if (f, m) in T and (f, "classic:hipr") in T]
    tr = []
    for inst, d in by.items():
        if not inst.startswith("sp-") or "classic:dijkstra" not in d:
            continue
        dj = d["classic:dijkstra"]
        best = [r["t_solve"] for m, r in d.items() if not m.startswith("classic") and m != "build-only" and ok(r, dj["obj"])]
        if best:
            tr.append(min(best) / max(dj["t_solve"], 1e-6))
    per_seed.append(out)
    print(f"seed {si}: " + "; ".join(f"{k} {v[0]}-{v[1]}" for k, v in out.items() if v))
    print(f"   simplex/Dijkstra {min(sd):.0f}-{max(sd):.0f}x, simplex/Preflow {min(sp):.0f}-{max(sp):.0f}x, "
          f"first-order/Dijkstra {min(fd):.1f}-{max(fd):.1f}x, first-order/Preflow {min(fp):.1f}-{max(fp):.1f}x, "
          f"fastest LP/Dijkstra time {min(tr):.0f}-{max(tr):.0f}x"
          + (f", first-order/HIPR {min(fh):.1f}-{max(fh):.1f}x, simplex/HIPR {min(sh):.0f}-{max(sh):.0f}x" if fh else ""))

# spread across seeds for runs solved in every seed
allruns = [runs(raw) for _, raw in pairs]
key = lambda inst: inst.rsplit("-s", 1)[0]
tspread, mspread = collections.defaultdict(list), collections.defaultdict(list)
for inst in allruns[0]:
    k = key(inst)
    insts = [next((i for i in br if key(i) == k), None) for br in allruns]
    if None in insts:
        continue
    for m, r0 in allruns[0][inst].items():
        rs = [br[i].get(m, {}) for br, i in zip(allruns, insts)]
        if all(r.get("status") in ("optimal", "ok") and r.get("t_solve", 0) > 0.05 for r in rs):
            t = [r["t_solve"] for r in rs]
            tspread[m].append(max(t) / min(t))
            h = [r.get("hwm_end_kb", 0) for r in rs]
            if min(h) > 0:
                mspread[m].append(max(h) / min(h))
print("\nspread over seeds (max/min) for runs solved in all seeds with t_solve > 0.05 s; median [max]:")
for m in sorted(tspread):
    print(f"  {m:22s} time {statistics.median(tspread[m]):.2f} [{max(tspread[m]):.2f}]   "
          f"memory {statistics.median(mspread[m]):.2f} [{max(mspread[m]):.2f}]   (n={len(tspread[m])})")
