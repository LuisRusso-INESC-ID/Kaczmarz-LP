#!/usr/bin/env python3
"""Convert the Zenodo USA road graphs (record 10.5281/zenodo.18026780; 9th DIMACS Challenge
travel-time graphs, largest connected component) to DIMACS .gr/.co, and cut regional subgraphs.

Input format: header "NAME t n m", then m lines "u v w" (0-based); coordinates: header
"NAME t n", then "v lon lat" (micro-degrees).

usage:
  convert_roads.py convert NW-t.graph NW-t.coordinates data/USA-road-t.NW
      -> data/USA-road-t.NW.gr, data/USA-road-t.NW.co
  convert_roads.py cut data/USA-road-t.NW 40000 160000 640000
      -> data/USA-road-t.NW-40k.gr (and .co) ...: nodes taken in BFS order from the vertex
         closest to the coordinate centroid until the induced subgraph has the target number
         of arcs; node ids renumbered in BFS order.
"""
import sys
import numpy as np


def convert(graph, coords, out):
    with open(graph) as f:
        name, kind, n, m = f.readline().split()
        n, m = int(n), int(m)
        a = np.loadtxt(f, dtype=np.int64, ndmin=2)
    assert a.shape == (m, 3), (a.shape, m)
    assert a[:, :2].min() >= 0 and a[:, :2].max() < n and a[:, 2].min() >= 0
    with open(out + ".gr", "w") as f:
        f.write(f"c 9th DIMACS Challenge USA road graph {name} ({'travel time' if kind == 't' else kind}),\n"
                f"c largest connected component, via Zenodo 10.5281/zenodo.18026780\n"
                f"p sp {n} {m}\n")
        a[:, :2] += 1
        np.savetxt(f, a, fmt="a %d %d %d")
    with open(coords) as f:
        f.readline()
        c = np.loadtxt(f, dtype=np.int64, ndmin=2)
    assert c.shape == (n, 3)
    with open(out + ".co", "w") as f:
        f.write(f"p aux sp co {n}\n")
        c[:, 0] += 1
        np.savetxt(f, c, fmt="v %d %d %d")
    print(f"{out}.gr: n={n} m={m}")


def read_gr(path):
    with open(path) as f:
        for line in f:
            if line.startswith("p "):
                n, m = map(int, line.split()[2:4])
                break
        a = np.loadtxt(f, dtype=np.int64, usecols=(1, 2, 3), ndmin=2)
    assert a.shape == (m, 3)
    return n, a - np.array([1, 1, 0])


def cut(base, targets):
    n, a = read_gr(base + ".gr")
    with open(base + ".co") as f:
        f.readline()
        co = np.loadtxt(f, dtype=np.int64, usecols=(1, 2), ndmin=2)
    # undirected adjacency (CSR) for the BFS
    u = np.concatenate([a[:, 0], a[:, 1]]); v = np.concatenate([a[:, 1], a[:, 0]])
    order = np.argsort(u, kind="stable"); u, v = u[order], v[order]
    start = np.searchsorted(u, np.arange(n + 1))
    centre = int(np.argmin(((co - co.mean(axis=0)) ** 2).sum(axis=1)))
    rank = np.full(n, -1, dtype=np.int64); rank[centre] = 0
    bfs = [centre]; head = 0
    while head < len(bfs):
        x = bfs[head]; head += 1
        for y in v[start[x]:start[x + 1]]:
            if rank[y] < 0:
                rank[y] = len(bfs); bfs.append(int(y))
    bfs = np.array(bfs)
    # an arc is in the subgraph of the first k BFS nodes iff max(rank[u], rank[v]) < k
    arc_k = np.maximum(rank[a[:, 0]], rank[a[:, 1]]) + 1
    arc_k_sorted = np.sort(arc_k)
    for t in sorted(targets):
        if t >= len(a):
            print(f"skip {t}: graph has only {len(a)} arcs"); continue
        k = int(arc_k_sorted[t - 1])                 # smallest node count giving >= t arcs
        sel = arc_k <= k
        sub = a[sel].copy()
        sub[:, 0] = rank[sub[:, 0]] + 1; sub[:, 1] = rank[sub[:, 1]] + 1
        tag = f"{t // 1000}k" if t % 1000 == 0 else str(t)
        out = f"{base}-{tag}"
        with open(out + ".gr", "w") as f:
            f.write(f"c BFS region of {base.split('/')[-1]} around vertex {centre + 1}: first {k} nodes in BFS order\n"
                    f"p sp {k} {len(sub)}\n")
            np.savetxt(f, sub, fmt="a %d %d %d")
        with open(out + ".co", "w") as f:
            f.write(f"p aux sp co {k}\n")
            np.savetxt(f, np.column_stack([np.arange(1, k + 1), co[bfs[:k]]]), fmt="v %d %d %d")
        print(f"{out}.gr: n={k} m={len(sub)}")


if __name__ == "__main__":
    if sys.argv[1] == "convert":
        convert(*sys.argv[2:5])
    elif sys.argv[1] == "cut":
        cut(sys.argv[2], [int(t) for t in sys.argv[3:]])
    else:
        sys.exit(__doc__)
