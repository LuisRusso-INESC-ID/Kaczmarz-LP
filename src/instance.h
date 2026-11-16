// instance.h -- common instance formats and readers.
//
//  Shortest path (DIMACS 9th challenge .gr):   p sp n m / a u v w
//      single-source problem; source given on the command line (default 1).
//  Max flow (DIMACS .max):                     p max n m / n s s / n t t / a u v cap
//  Max concurrent multicommodity flow (.mcf, own format):
//      p mcf n m k / a u v cap / d s t demand
//
// Nodes are 1-based in files, 0-based in memory.
#pragma once
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <stdexcept>

struct Arc { int u, v; int64_t w; };          // w = length or capacity
struct Demand { int s, t; int64_t d; };

struct Instance {
    std::string kind;     // "sp", "max", "mcf"
    int n = 0;
    std::vector<Arc> arcs;
    int s = 0, t = 0;     // max-flow terminals
    std::vector<Demand> dem; // mcf commodities
};

inline Instance read_instance(const char* path) {
    FILE* f = std::fopen(path, "r");
    if (!f) throw std::runtime_error(std::string("cannot open ") + path);
    Instance I;
    char line[512];
    size_t m = 0;
    while (std::fgets(line, sizeof line, f)) {
        char c = line[0];
        if (c == 'p') {
            char kind[16]; long n, mm, k = 0;
            int r = std::sscanf(line, "p %15s %ld %ld %ld", kind, &n, &mm, &k);
            if (r < 3) throw std::runtime_error("bad p line");
            I.kind = kind; I.n = (int)n; m = (size_t)mm;
            I.arcs.reserve(m);
            if (k) I.dem.reserve(k);
        } else if (c == 'a') {
            long u, v; long long w;
            std::sscanf(line, "a %ld %ld %lld", &u, &v, &w);
            I.arcs.push_back({(int)u - 1, (int)v - 1, (int64_t)w});
        } else if (c == 'n') {
            long id; char ty;
            std::sscanf(line, "n %ld %c", &id, &ty);
            if (ty == 's') I.s = (int)id - 1; else I.t = (int)id - 1;
        } else if (c == 'd') {
            long s, t; long long d;
            std::sscanf(line, "d %ld %ld %lld", &s, &t, &d);
            I.dem.push_back({(int)s - 1, (int)t - 1, (int64_t)d});
        }
    }
    std::fclose(f);
    if (I.arcs.size() != m)
        std::fprintf(stderr, "warning: header says %zu arcs, read %zu\n", m, I.arcs.size());
    return I;
}

// Forward-star adjacency (CSR over arcs by tail).
struct CSR {
    std::vector<int> start, arc; // arc ids sorted by tail
    void build(int n, const std::vector<Arc>& A, bool by_head = false) {
        start.assign(n + 1, 0);
        for (auto& a : A) start[(by_head ? a.v : a.u) + 1]++;
        for (int i = 0; i < n; i++) start[i + 1] += start[i];
        arc.resize(A.size());
        std::vector<int> pos(start.begin(), start.end() - 1);
        for (int e = 0; e < (int)A.size(); e++) arc[pos[by_head ? A[e].v : A[e].u]++] = e;
    }
};

// Commodities grouped by source (valid aggregation for max concurrent flow:
// a single-source multi-sink flow decomposes into per-sink paths).
struct SourceGroup { int s; int64_t total; std::vector<std::pair<int,int64_t>> sinks; };
inline std::vector<SourceGroup> group_by_source(const Instance& I) {
    std::vector<int> idx(I.n, -1);
    std::vector<SourceGroup> G;
    for (auto& d : I.dem) {
        if (idx[d.s] < 0) { idx[d.s] = (int)G.size(); G.push_back({d.s, 0, {}}); }
        auto& g = G[idx[d.s]];
        g.total += d.d; g.sinks.push_back({d.t, d.d});
    }
    return G;
}
