#!/usr/bin/env python3
"""Benchmark driver: generate instances, run classical codes and LP backends, one process
per solve, recording wall time and peak RSS (ru_maxrss from wait4).

usage: run.py [--plan pilot|small|roads|roads-usa] [--tl SECONDS] [--seed S] [--out results/pilot.jsonl]
              [--only FAMILY] [--methods M1,M2,...]
The road plans read data/USA-road-t.*.gr (scripts/fetch_zenodo_roads.sh); --seed applies to generated families.
Re-running skips (instance, method) pairs already present in the output file.
"""
import argparse, json, os, resource, subprocess, sys, time, signal, shlex

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
DATA = os.path.join(ROOT, "data")

# ---------------------------------------------------------------- instance plans
PLANS = {
    "small": {
        "sp-grid": [(32, 32), (64, 64)],
        "sp-road": [(32, 32), (64, 64)],
        "sp-rand": [(1000, 4), (4000, 4)],
        "mf-rmf": [(8, 8), (12, 12)],
        "mf-rlg": [(64, 16), (64, 64)],
        "mf-vision": [(32, 32), (64, 64)],
        "mcf-grid": [(10, 10, 4, 5), (15, 15, 6, 5)],
        "mcf-rand": [(200, 3, 4, 5), (500, 3, 6, 5)],
    },
    "pilot": {
        "sp-grid": [(32, 32), (64, 64), (128, 128), (256, 256), (512, 512)],
        "sp-road": [(32, 32), (64, 64), (128, 128), (256, 256), (512, 512)],
        "sp-rand": [(1000, 4), (4000, 4), (16000, 4), (64000, 4), (256000, 4)],
        "mf-rmf": [(8, 8), (12, 12), (16, 16), (24, 24), (32, 32)],
        "mf-rlg": [(64, 16), (64, 64), (128, 128), (256, 256), (512, 512)],
        "mf-vision": [(32, 32), (64, 64), (128, 128), (256, 256), (512, 512)],
        "mcf-grid": [(10, 10, 4, 5), (15, 15, 6, 5), (20, 20, 8, 5), (30, 30, 12, 5), (40, 40, 16, 5)],
        "mcf-rand": [(200, 3, 4, 5), (500, 3, 6, 5), (1000, 3, 8, 5), (2500, 3, 12, 5), (5000, 3, 16, 5)],
    },
    # 9th DIMACS USA travel-time road graphs (largest connected component), BFS regions of NW
    "roads": {"sp-roadnet": [("NW-40k",), ("NW-160k",), ("NW-640k",), ("NW",)]},
    "roads-usa": {"sp-roadnet": [("USA",)]},
}
EXT = {"sp": "gr", "mf": "max", "mcf": "mcf"}

CLASSIC = {
    "sp": ["dijkstra", "lemon-dijkstra", "lemon-ns"],
    "mf": ["lemon-preflow", "lemon-ns"],
    "mcf": ["gk:0.1", "gk:0.05"],
}
LP_BACKENDS = ["build-only",
               "highs-dual", "highs-primal", "highs-ipx-cross", "highs-ipx", "highs-pdlp",
               "clp-dual", "clp-primal", "clp-barrier",
               "glpk-primal", "glpk-dual", "glpk-ipt",
               "ortools-glop", "ortools-glop-primal", "ortools-pdlp",
               "own-pdhg", "own-prow", "own-skm"]

MEM_LIMIT = 6 * 1024**3   # address-space cap per solve


def lp_binary(backend):
    return os.path.join(BIN, "lpsolve-ort" if backend.startswith("ortools-") else "lpsolve")


def available_lp_backends():
    _, probe = instance_path("sp-grid", (3, 3), 1)
    out = []
    for b in LP_BACKENDS:
        if not os.path.exists(lp_binary(b)):
            continue
        r = subprocess.run([lp_binary(b), b, probe, "5"], capture_output=True, text=True)
        if "unknown or unavailable" not in r.stderr:
            out.append(b)
    return out


def instance_path(fam, params, seed=1):
    if fam == "sp-roadnet":   # fixed files, not generated
        name = "USA-road-t." + params[0]
        path = os.path.join(DATA, name + ".gr")
        if not os.path.exists(path):
            sys.exit(f"missing {path}: run scripts/fetch_zenodo_roads.sh")
        return name, path
    name = f"{fam}-" + "x".join(map(str, params)) + f"-s{seed}"
    ext = EXT[fam.split("-")[0]]
    path = os.path.join(DATA, f"{name}.{ext}")
    if not os.path.exists(path):
        with open(path, "w") as f:
            subprocess.run([os.path.join(BIN, "gen"), fam, *map(str, params), str(seed)], stdout=f, check=True)
    return name, path




def parse_result(out):
    d0 = {}
    for line in out.splitlines():
        if line.startswith("HWM "):
            d0["hwm_end_kb"] = float(line.split("=")[1])
    for line in out.splitlines():
        if line.startswith("RESULT"):
            d = {}
            for tok in line.split()[1:]:
                k, _, v = tok.partition("=")
                try:
                    d[k] = float(v) if k not in ("status",) else v
                except ValueError:
                    d[k] = v
            d.update(d0)
            return d
    return None


def measured(cmd, timeout):
    """Fork via a tiny wrapper so wait4 gives the per-child max RSS precisely."""
    wrapper = [sys.executable, os.path.join(ROOT, "scripts", "memwrap.py"), str(MEM_LIMIT), str(timeout), *cmd]
    r = subprocess.run(wrapper, capture_output=True, text=True)
    meta = None
    for line in r.stderr.splitlines():
        if line.startswith("MEMWRAP "):
            meta = json.loads(line[8:])
    return r.stdout, r.stderr, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", default="small")
    ap.add_argument("--tl", type=float, default=120)
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=None)
    ap.add_argument("--only", default=None, help="comma-separated families")
    ap.add_argument("--methods", default=None, help="comma-separated methods to restrict to")
    a = ap.parse_args()
    out_path = a.out or os.path.join(ROOT, "results", f"{a.plan}.jsonl")
    os.makedirs(DATA, exist_ok=True); os.makedirs(os.path.dirname(out_path), exist_ok=True)
    done = set()
    if os.path.exists(out_path):
        for line in open(out_path):
            r = json.loads(line); done.add((r["instance"], r["method"]))
    lp_backs = available_lp_backends()
    print("LP backends:", " ".join(lp_backs), flush=True)
    plan = PLANS[a.plan]
    fams = a.only.split(",") if a.only else list(plan)
    restrict = set(a.methods.split(",")) if a.methods else None
    with open(out_path, "a") as fo:
        for fam in fams:
            kind = fam.split("-")[0]
            failed = set()   # methods that failed/timed out on a smaller size of this family
            for si, params in enumerate(plan[fam]):
                name, path = instance_path(fam, params, a.seed)
                methods = [("classic", m) for m in CLASSIC[kind]] + [("lp", b) for b in lp_backs]
                for typ, m in methods:
                    mname = m if typ == "lp" else "classic:" + m
                    if restrict and mname not in restrict and m not in restrict:
                        continue
                    if (name, mname) in done:
                        continue
                    if mname == "own-skm" and si >= 2:
                        continue   # SKM on the optimality system is only run on the two smallest sizes
                    if mname in failed:
                        rec = {"instance": name, "family": fam, "params": params, "method": mname,
                               "status": "skipped(smaller failed)"}
                        fo.write(json.dumps(rec) + "\n"); fo.flush(); continue
                    if typ == "classic":
                        algo, _, opt = m.partition(":")
                        cmd = [os.path.join(BIN, "classic"), algo, path] + ([opt] if opt else [])
                    else:
                        cmd = [lp_binary(m), m, path, str(a.tl), str(a.tol)]
                    out, err, meta = measured(cmd, a.tl * 1.5 + 30)
                    res = parse_result(out) or {}
                    rec = {"instance": name, "family": fam, "params": params, "method": mname, **res}
                    if meta:
                        rec["wall"] = meta["wall"]; rec["maxrss_kb"] = meta["maxrss_kb"]
                        if meta["killed"]:
                            rec["status"] = "killed(timeout)"
                        elif meta["returncode"] != 0 and "status" not in res:
                            rec["status"] = f"crash(rc={meta['returncode']})"
                            rec["stderr"] = err[-300:]
                    st = rec.get("status", "none")
                    if st not in ("optimal", "ok"):
                        failed.add(mname)
                    fo.write(json.dumps(rec) + "\n"); fo.flush()
                    print(f"{name:28s} {mname:22s} {st:16s} t={rec.get('t_solve', float('nan')):9.3f} "
                          f"hwm={rec.get('hwm_end_kb', 0)/1024:8.1f}MB obj={rec.get('obj', '')}", flush=True)


if __name__ == "__main__":
    main()
