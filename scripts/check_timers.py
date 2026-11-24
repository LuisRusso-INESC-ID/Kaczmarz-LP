#!/usr/bin/env python3
"""check_timers.py [results.jsonl ...] -- cross-check the times in
stored results. Each record has t_solve <= t_total, measured inside
the solver process (steady_clock), and wall, measured from outside
by memwrap.py (time.monotonic in the parent). A record fails if the
inner times exceed the outer one. Exit status 1 on any failure."""
import glob, json, statistics, sys

TOL = 1e-3
DEFAULT = ["results/pilot-pin-s1.jsonl", "results/pilot-pin-s2.jsonl",
           "results/pilot-pin-s3.jsonl", "results/roads-local.jsonl",
           "results/roads-600-local.jsonl",
           "results/roads-usa-local.jsonl"]


def timed(path):
    for line in open(path):
        r = json.loads(line)
        if all(k in r for k in ("t_solve", "t_total", "wall")):
            yield r


def failure(r):
    if r["t_solve"] > r["t_total"] + TOL:
        return "t_solve > t_total"
    if r["t_total"] > r["wall"] + TOL:
        return "t_total > wall"
    return None


def pct(xs, q):
    return xs[min(len(xs) - 1, int(q * len(xs)))]


def report(path):
    recs = list(timed(path))
    bad = [(r, failure(r)) for r in recs if failure(r)]
    gap = sorted(r["wall"] - r["t_total"] for r in recs)
    long = [r for r in recs if r["wall"] >= 1]
    share = sorted((r["wall"] - r["t_total"]) / r["wall"]
                   for r in long)
    print(f"{path}: {len(recs)} timed records, {len(bad)} failures")
    if gap:
        med = statistics.median(gap)
        print(f"  wall - t_total: median {med:.4f} s,"
              f" p99 {pct(gap, 0.99):.4f} s, max {gap[-1]:.4f} s")
    if share:
        print(f"  runs >= 1 s ({len(long)}): (wall - t_total)/wall"
              f" median {100 * statistics.median(share):.2f}%,"
              f" max {100 * share[-1]:.2f}%")
    for r, why in bad[:10]:
        print(f"  FAIL {r['instance']} {r['method']}: {why}"
              f" (t_solve={r['t_solve']}, t_total={r['t_total']},"
              f" wall={r['wall']:.4f})")
    return len(bad)


def main():
    paths = sys.argv[1:] or [p for p in DEFAULT if glob.glob(p)]
    nbad = sum(report(p) for p in paths)
    print("PASS" if nbad == 0 else f"FAIL: {nbad} records")
    sys.exit(1 if nbad else 0)


if __name__ == "__main__":
    main()
