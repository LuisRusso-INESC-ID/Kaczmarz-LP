// gen.cpp -- instance generators (re-implementations of the standard families).
//
// usage: gen <family> <params...> <seed> > file
//
// Shortest path (.gr):
//   sp-grid   R C            bidirected R x C grid, lengths U[1,1000]
//   sp-rand   n d            random digraph, out-degree d + Hamiltonian cycle, lengths U[1,1000]
//   sp-road   R C            road-like: jittered grid, ~25% edges dropped (spanning tree kept),
//                            lengths = 1000 * Euclidean distance (rounded, >= 1)
// Max flow (.max):
//   mf-rmf    a b            Goldfarb-Grigoriadis RMF: b frames of a x a grids, c1=1, c2=10000
//   mf-rlg    R C            Washington random level graph: C levels of R nodes, 3 arcs to next level
//   mf-vision R C            4-connected grid with terminal arcs (computer-vision style)
// Max concurrent multicommodity flow (.mcf):
//   mcf-grid  R C S T        bidirected grid, caps U[50,150]; S sources x T sinks, demands U[1,100]
//   mcf-rand  n d S T        random digraph as sp-rand, caps U[50,150]
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>
#include <numeric>
#include <algorithm>
#include <string>

static std::mt19937_64 rng;
static int64_t U(int64_t lo, int64_t hi) { return std::uniform_int_distribution<int64_t>(lo, hi)(rng); }
static double Ud() { return std::uniform_real_distribution<double>(0, 1)(rng); }

struct A { int u, v; int64_t w; };

static void emit(const char* kind, int n, const std::vector<A>& arcs, const char* comment,
                 int s = -1, int t = -1, const std::vector<A>* dem = nullptr) {
    std::printf("c %s\n", comment);
    if (dem) std::printf("p %s %d %zu %zu\n", kind, n, arcs.size(), dem->size());
    else std::printf("p %s %d %zu\n", kind, n, arcs.size());
    if (s >= 0) std::printf("n %d s\nn %d t\n", s + 1, t + 1);
    for (auto& a : arcs) std::printf("a %d %d %lld\n", a.u + 1, a.v + 1, (long long)a.w);
    if (dem) for (auto& d : *dem) std::printf("d %d %d %lld\n", d.u + 1, d.v + 1, (long long)d.w);
}

static std::vector<A> grid(int R, int C, int64_t lo, int64_t hi) {
    std::vector<A> arcs;
    auto id = [&](int r, int c) { return r * C + c; };
    for (int r = 0; r < R; r++)
        for (int c = 0; c < C; c++) {
            if (c + 1 < C) { arcs.push_back({id(r, c), id(r, c + 1), U(lo, hi)}); arcs.push_back({id(r, c + 1), id(r, c), U(lo, hi)}); }
            if (r + 1 < R) { arcs.push_back({id(r, c), id(r + 1, c), U(lo, hi)}); arcs.push_back({id(r + 1, c), id(r, c), U(lo, hi)}); }
        }
    return arcs;
}

static std::vector<A> randgraph(int n, int d, int64_t lo, int64_t hi) {
    std::vector<A> arcs;
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    for (int i = 0; i < n; i++) arcs.push_back({perm[i], perm[(i + 1) % n], U(lo, hi)}); // strong connectivity
    for (int u = 0; u < n; u++)
        for (int j = 0; j < d; j++) {
            int v = (int)U(0, n - 1);
            if (v != u) arcs.push_back({u, v, U(lo, hi)});
        }
    return arcs;
}

// union-find for spanning-tree-preserving edge deletion
struct DSU { std::vector<int> p; DSU(int n) : p(n) { std::iota(p.begin(), p.end(), 0); }
    int f(int x) { while (p[x] != x) x = p[x] = p[p[x]]; return x; }
    bool u(int a, int b) { a = f(a); b = f(b); if (a == b) return false; p[a] = b; return true; } };

static std::vector<A> road(int R, int C) {
    int n = R * C;
    std::vector<double> X(n), Y(n);
    for (int r = 0; r < R; r++) for (int c = 0; c < C; c++) {
        X[r * C + c] = c + 0.8 * (Ud() - 0.5); Y[r * C + c] = r + 0.8 * (Ud() - 0.5); }
    std::vector<std::pair<int,int>> E;
    for (int r = 0; r < R; r++) for (int c = 0; c < C; c++) {
        int v = r * C + c;
        if (c + 1 < C) E.push_back({v, v + 1});
        if (r + 1 < R) E.push_back({v, v + C});
        if (c + 1 < C && r + 1 < R && Ud() < 0.1) E.push_back({v, v + C + 1}); // occasional diagonal
    }
    std::shuffle(E.begin(), E.end(), rng);
    DSU d(n);
    std::vector<char> keep(E.size(), 0);
    for (size_t i = 0; i < E.size(); i++) if (d.u(E[i].first, E[i].second)) keep[i] = 1;
    for (size_t i = 0; i < E.size(); i++) if (!keep[i] && Ud() < 0.66) keep[i] = 1; // drop ~25% overall
    std::vector<A> arcs;
    for (size_t i = 0; i < E.size(); i++) if (keep[i]) {
        int a = E[i].first, b = E[i].second;
        int64_t w = std::max<int64_t>(1, std::llround(1000.0 * std::hypot(X[a] - X[b], Y[a] - Y[b])));
        arcs.push_back({a, b, w}); arcs.push_back({b, a, w});
    }
    std::sort(arcs.begin(), arcs.end(), [](const A& x, const A& y) { return x.u != y.u ? x.u < y.u : x.v < y.v; });
    return arcs;
}

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: gen family params... seed\n"); return 1; }
    std::string fam = argv[1];
    rng.seed(std::strtoull(argv[argc - 1], nullptr, 10));
    char comment[256];
    std::snprintf(comment, sizeof comment, "generated by lpgraph gen:");
    std::string cm = comment;
    for (int i = 1; i < argc; i++) { cm += " "; cm += argv[i]; }

    if (fam == "sp-grid") {
        int R = atoi(argv[2]), C = atoi(argv[3]);
        emit("sp", R * C, grid(R, C, 1, 1000), cm.c_str());
    } else if (fam == "sp-rand") {
        int n = atoi(argv[2]), d = atoi(argv[3]);
        emit("sp", n, randgraph(n, d, 1, 1000), cm.c_str());
    } else if (fam == "sp-road") {
        int R = atoi(argv[2]), C = atoi(argv[3]);
        emit("sp", R * C, road(R, C), cm.c_str());
    } else if (fam == "mf-rmf") {
        int a = atoi(argv[2]), b = atoi(argv[3]);
        const int64_t c1 = 1, c2 = 10000;
        int fs = a * a, n = fs * b;
        std::vector<A> arcs;
        for (int k = 0; k < b; k++) {
            int base = k * fs;
            for (int r = 0; r < a; r++) for (int c = 0; c < a; c++) {
                int v = base + r * a + c;
                if (c + 1 < a) { arcs.push_back({v, v + 1, c2 * fs}); arcs.push_back({v + 1, v, c2 * fs}); }
                if (r + 1 < a) { arcs.push_back({v, v + a, c2 * fs}); arcs.push_back({v + a, v, c2 * fs}); }
            }
            if (k + 1 < b) {
                std::vector<int> perm(fs);
                std::iota(perm.begin(), perm.end(), 0);
                std::shuffle(perm.begin(), perm.end(), rng);
                for (int i = 0; i < fs; i++) arcs.push_back({base + i, base + fs + perm[i], U(c1, c2)});
            }
        }
        emit("max", n, arcs, cm.c_str(), 0, n - 1);
    } else if (fam == "mf-rlg") {
        int R = atoi(argv[2]), C = atoi(argv[3]);
        int n = R * C + 2, s = R * C, t = R * C + 1;
        std::vector<A> arcs;
        for (int r = 0; r < R; r++) arcs.push_back({s, r, (int64_t)1 << 40});
        for (int c = 0; c + 1 < C; c++)
            for (int r = 0; r < R; r++) {
                int v = c * R + r;
                for (int j = 0; j < 3; j++) arcs.push_back({v, (c + 1) * R + (int)U(0, R - 1), U(1, 10000)});
            }
        for (int r = 0; r < R; r++) arcs.push_back({(C - 1) * R + r, t, (int64_t)1 << 40});
        emit("max", n, arcs, cm.c_str(), s, t);
    } else if (fam == "mf-vision") {
        int R = atoi(argv[2]), C = atoi(argv[3]);
        int n = R * C + 2, s = R * C, t = R * C + 1;
        std::vector<A> arcs = grid(R, C, 1, 100);
        for (int v = 0; v < R * C; v++) {
            int64_t a = U(0, 200), b = U(0, 200);          // data term
            if (a > b) arcs.push_back({s, v, a - b}); else if (b > a) arcs.push_back({v, t, b - a});
        }
        emit("max", n, arcs, cm.c_str(), s, t);
    } else if (fam == "mcf-grid" || fam == "mcf-rand") {
        bool g = fam == "mcf-grid";
        int p1 = atoi(argv[2]), p2 = atoi(argv[3]), S = atoi(argv[4]), T = atoi(argv[5]);
        int n = g ? p1 * p2 : p1;
        std::vector<A> arcs = g ? grid(p1, p2, 50, 150) : randgraph(p1, p2, 50, 150);
        std::vector<A> dem;
        std::vector<int> nodes(n);
        std::iota(nodes.begin(), nodes.end(), 0);
        std::shuffle(nodes.begin(), nodes.end(), rng);
        for (int i = 0; i < S; i++) {
            int s = nodes[i];
            for (int j = 0; j < T; j++) {
                int t;
                do t = (int)U(0, n - 1); while (t == s);
                dem.push_back({s, t, U(1, 100)});
            }
        }
        emit("mcf", n, arcs, cm.c_str(), -1, -1, &dem);
    } else { std::fprintf(stderr, "unknown family %s\n", fam.c_str()); return 1; }
    return 0;
}
