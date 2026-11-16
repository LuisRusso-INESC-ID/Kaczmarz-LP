// lp.h -- sparse LP model (column-major) and the three graph formulations.
//
// Model:  min c'x   s.t.  rlo <= A x <= rup,   clo <= x <= cup
// Rows are either equalities (rlo == rup) or upper-bounded (rlo = -inf).
#pragma once
#include "instance.h"
#include <cmath>
#include <limits>
#include <queue>
#include <cstdlib>

// Optional convergence trace: set LPTRACE=<path> to get CSV rows
//   method,time,iter,pobj,rp,rd,gap,extra
inline FILE* trace_file() {
    static FILE* f = nullptr; static bool init = false;
    if (!init) { init = true; if (const char* p = std::getenv("LPTRACE")) { f = std::fopen(p, "w");
        if (f) std::fprintf(f, "method,time,iter,pobj,rp,rd,gap,extra\n"); } }
    return f;
}
inline void trace(const char* meth, double t, long it, double pobj, double rp, double rd, double gap, double extra = 0) {
    if (FILE* f = trace_file()) { std::fprintf(f, "%s,%.6f,%ld,%.12g,%.6e,%.6e,%.6e,%.6g\n", meth, t, it, pobj, rp, rd, gap, extra); std::fflush(f); }
}

constexpr double INF = std::numeric_limits<double>::infinity();
constexpr int64_t BIGCAP = (int64_t)1 << 40;   // generator convention for "infinite" capacity

struct LP {
    int nrow = 0, ncol = 0;
    std::vector<double> c, clo, cup, rlo, rup;
    std::vector<int> start, index;   // CSC
    std::vector<double> value;
    double sense = 1.0;              // +1 min, -1 if original problem was max (c already negated)
    // Optional valid upper bounds for columns with cup = +inf: bounds that keep at least one
    // optimal solution feasible (used only to make Lagrangian dual bounds finite).
    std::vector<double> uhint;
    size_t nnz() const { return value.size(); }
};

// Column builder helper: columns appended in order, each with a few entries.
struct ColBuilder {
    LP& L;
    explicit ColBuilder(LP& l) : L(l) { L.start.assign(1, 0); }
    void add(double cost, double lo, double up, std::initializer_list<std::pair<int,double>> ent) {
        for (auto& p : ent) if (p.first >= 0) { L.index.push_back(p.first); L.value.push_back(p.second); }
        L.start.push_back((int)L.index.size());
        L.c.push_back(cost); L.clo.push_back(lo); L.cup.push_back(up); L.ncol++;
    }
};

static std::vector<char> reachable(const Instance& I, int src) {
    CSR G; G.build(I.n, I.arcs);
    std::vector<char> r(I.n, 0); std::vector<int> st{src}; r[src] = 1;
    while (!st.empty()) { int u = st.back(); st.pop_back();
        for (int k = G.start[u]; k < G.start[u + 1]; k++) { int v = I.arcs[G.arc[k]].v; if (!r[v]) { r[v] = 1; st.push_back(v); } } }
    return r;
}

// Shortest-path tree LP: supply (R-1) at the source, -1 at every other reachable node.
// Row of the source is dropped (it is implied), so A has full row rank.
inline LP build_sp(const Instance& I, int src) {
    LP L; auto r = reachable(I, src);
    std::vector<int> row(I.n, -1); long R = 0;
    for (int v = 0; v < I.n; v++) if (r[v]) { R++; if (v != src) row[v] = L.nrow++; }
    L.rlo.assign(L.nrow, -1.0); L.rup.assign(L.nrow, -1.0);
    ColBuilder B(L);
    L.index.reserve(2 * I.arcs.size()); L.value.reserve(2 * I.arcs.size());
    for (auto& a : I.arcs) if (r[a.u] && r[a.v])
        B.add((double)a.w, 0.0, INF, {{row[a.u], 1.0}, {row[a.v], -1.0}});
    L.uhint.assign(L.ncol, (double)(R - 1));   // an optimal flow never exceeds the total supply on an arc
    return L;
}

// Max flow: max net outflow of s; conservation at all other nodes except t.
inline LP build_max(const Instance& I) {
    LP L; std::vector<int> row(I.n, -1);
    for (int v = 0; v < I.n; v++) if (v != I.s && v != I.t) row[v] = L.nrow++;
    L.rlo.assign(L.nrow, 0.0); L.rup.assign(L.nrow, 0.0);
    L.sense = -1.0;
    ColBuilder B(L);
    for (auto& a : I.arcs) {
        double c = (a.u == I.s ? -1.0 : 0.0) + (a.v == I.s ? 1.0 : 0.0);   // minimize -(outflow of s)
        B.add(c, 0.0, a.w >= BIGCAP ? INF : (double)a.w, {{row[a.u], 1.0}, {row[a.v], -1.0}});
    }
    return L;
}

// Max concurrent flow, commodities aggregated by source:
//   max lambda  s.t. for each group g, node v != s_g:  out_g(v) - in_g(v) + d_{g,v} lambda = 0
//                    sum_g x_{g,e} <= u_e
inline LP build_mcf(const Instance& I) {
    LP L; auto G = group_by_source(I);
    int n = I.n, m = (int)I.arcs.size(), K = (int)G.size();
    // node rows: group g, node v -> g*(n-1) + rank of v skipping s_g
    auto nrow = [&](int g, int v) { int s = G[g].s; return v == s ? -1 : g * (n - 1) + (v < s ? v : v - 1); };
    int caprow0 = K * (n - 1);
    L.nrow = caprow0 + m;
    L.rlo.assign(L.nrow, 0.0); L.rup.assign(L.nrow, 0.0);
    for (int e = 0; e < m; e++) { L.rlo[caprow0 + e] = -INF; L.rup[caprow0 + e] = (double)I.arcs[e].w; }
    L.sense = -1.0;
    ColBuilder B(L);
    L.index.reserve((size_t)3 * K * m + n * K); L.value.reserve((size_t)3 * K * m + n * K);
    for (int g = 0; g < K; g++)
        for (int e = 0; e < m; e++) {
            auto& a = I.arcs[e];
            B.add(0.0, 0.0, INF, {{nrow(g, a.u), 1.0}, {nrow(g, a.v), -1.0}, {caprow0 + e, 1.0}});
        }
    // lambda column: +d at each sink row (sum over duplicate sinks)
    std::vector<double> coef(L.nrow, 0.0);
    for (int g = 0; g < K; g++) for (auto& sd : G[g].sinks) coef[nrow(g, sd.first)] += (double)sd.second;
    for (int i = 0; i < L.nrow; i++) if (coef[i] != 0) { L.index.push_back(i); L.value.push_back(coef[i]); }
    L.start.push_back((int)L.index.size());
    L.c.push_back(-1.0); L.clo.push_back(0.0); L.cup.push_back(INF); L.ncol++;
    {   // implied bounds: x_{g,e} <= u_e (capacity row); lambda * d_{g,v} <= capacity into v
        std::vector<double> incap(n, 0.0);
        for (auto& a : I.arcs) incap[a.v] += (double)a.w;
        double lam = INF;
        for (int g = 0; g < K; g++) for (auto& sd : G[g].sinks) {
            double d = coef[nrow(g, sd.first)];
            if (d > 0) lam = std::min(lam, incap[sd.first] / d);
        }
        L.uhint.resize(L.ncol);
        for (int g = 0; g < K; g++) for (int e = 0; e < m; e++) L.uhint[(size_t)g * m + e] = (double)I.arcs[e].w;
        L.uhint[L.ncol - 1] = lam;
    }
    // CSC requires sorted row indices inside a column for some solvers: entries are added
    // with arbitrary order; sort each column.
    for (int j = 0; j < L.ncol; j++) {
        int b = L.start[j], e = L.start[j + 1];
        for (int p = b + 1; p < e; p++)       // insertion sort, columns are tiny
            for (int q = p; q > b && L.index[q - 1] > L.index[q]; q--) {
                std::swap(L.index[q - 1], L.index[q]); std::swap(L.value[q - 1], L.value[q]); }
    }
    return L;
}

inline LP build_lp(const Instance& I, int src = 0) {
    LP L = I.kind == "sp" ? build_sp(I, src) : I.kind == "max" ? build_max(I) : build_mcf(I);
    // sort SP/max columns too (two entries)
    for (int j = 0; j < L.ncol; j++) {
        int b = L.start[j], e = L.start[j + 1];
        if (e - b == 2 && L.index[b] > L.index[b + 1]) { std::swap(L.index[b], L.index[b + 1]); std::swap(L.value[b], L.value[b + 1]); }
    }
    return L;
}

// Residuals of a primal point (max abs violation, relative to 1+max|rhs|).
inline double primal_residual(const LP& L, const std::vector<double>& x) {
    std::vector<double> ax(L.nrow, 0.0);
    double viol = 0, scale = 1;
    for (int j = 0; j < L.ncol; j++) {
        viol = std::max(viol, std::max(L.clo[j] - x[j], x[j] - L.cup[j]));
        for (int p = L.start[j]; p < L.start[j + 1]; p++) ax[L.index[p]] += L.value[p] * x[j];
    }
    for (int i = 0; i < L.nrow; i++) {
        viol = std::max(viol, std::max(L.rlo[i] - ax[i], ax[i] - L.rup[i]));
        if (std::isfinite(L.rup[i])) scale = std::max(scale, std::fabs(L.rup[i]));
    }
    return viol / scale;
}
