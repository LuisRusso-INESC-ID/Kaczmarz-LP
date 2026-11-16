// classic.cpp -- classical combinatorial baselines.
//
// usage: classic <algo> <instance> [options]
//   sp:  dijkstra       own binary-heap Dijkstra on CSR (source = node 1)
//        lemon-dijkstra LEMON Dijkstra
//        lemon-ns       LEMON network simplex on the shortest-path-tree LP (graph-aware simplex)
//   max: lemon-preflow  LEMON Preflow (Goldberg-Tarjan, highest label, global+gap relabel)
//        lemon-ns       LEMON network simplex on the max-flow circulation
//   mcf: gk [eps]       Garg-Koenemann / Fleischer FPTAS for max concurrent flow (tree routing)
//
// Output: one line "RESULT key=value ..." parsed by the driver.
#include "instance.h"
#include <chrono>
#include <queue>
#include <cmath>
#include <limits>
#include <algorithm>
#include <lemon/smart_graph.h>
#include <lemon/dijkstra.h>
#include <lemon/preflow.h>
#include <lemon/network_simplex.h>
#include <fstream>
static long vm_hwm_kb() {
    std::ifstream f("/proc/self/status"); std::string l;
    while (std::getline(f, l)) if (l.rfind("VmHWM:", 0) == 0) return std::stol(l.substr(6));
    return -1;
}

using Clock = std::chrono::steady_clock;
static double secs(Clock::time_point a) { return std::chrono::duration<double>(Clock::now() - a).count(); }

// ---------------- own Dijkstra ----------------
static std::vector<int64_t> dijkstra(const Instance& I, const CSR& G, int src) {
    const int64_t INF = std::numeric_limits<int64_t>::max();
    std::vector<int64_t> d(I.n, INF);
    using P = std::pair<int64_t, int>;
    std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
    d[src] = 0; pq.push({0, src});
    while (!pq.empty()) {
        auto [du, u] = pq.top(); pq.pop();
        if (du != d[u]) continue;
        for (int k = G.start[u]; k < G.start[u + 1]; k++) {
            const Arc& a = I.arcs[G.arc[k]];
            int64_t nd = du + a.w;
            if (nd < d[a.v]) { d[a.v] = nd; pq.push({nd, a.v}); }
        }
    }
    return d;
}

static void sp(const std::string& algo, const Instance& I, int src) {
    auto t0 = Clock::now();
    int64_t sum = 0; long reach = 0;
    double tsolve = 0;
    if (algo == "dijkstra") {
        CSR G; G.build(I.n, I.arcs);
        auto t1 = Clock::now();
        auto d = dijkstra(I, G, src);
        tsolve = secs(t1);
        for (auto x : d) if (x != std::numeric_limits<int64_t>::max()) { sum += x; reach++; }
    } else {
        lemon::SmartDigraph g;
        g.reserveNode(I.n); g.reserveArc(I.arcs.size());
        std::vector<lemon::SmartDigraph::Node> N(I.n);
        for (int i = 0; i < I.n; i++) N[i] = g.addNode();
        lemon::SmartDigraph::ArcMap<int64_t> len(g);
        for (auto& a : I.arcs) len[g.addArc(N[a.u], N[a.v])] = a.w;
        if (algo == "lemon-dijkstra") {
            auto t1 = Clock::now();
            lemon::Dijkstra<lemon::SmartDigraph, lemon::SmartDigraph::ArcMap<int64_t>> D(g, len);
            D.run(N[src]);
            tsolve = secs(t1);
            for (int i = 0; i < I.n; i++) if (D.reached(N[i])) { sum += D.dist(N[i]); reach++; }
        } else if (algo == "lemon-ns") {
            // shortest-path-tree LP: supply (r-1) at source, -1 at every reachable node
            CSR G; G.build(I.n, I.arcs);
            auto d = dijkstra(I, G, src); // only for reachability (cheap preprocessing)
            lemon::SmartDigraph::NodeMap<int64_t> sup(g, 0);
            for (int i = 0; i < I.n; i++) if (d[i] != std::numeric_limits<int64_t>::max() && i != src) { sup[N[i]] = -1; reach++; }
            sup[N[src]] = reach; reach++;
            auto t1 = Clock::now();
            lemon::NetworkSimplex<lemon::SmartDigraph, int64_t, int64_t> NS(g);
            NS.costMap(len).supplyMap(sup);
            auto st = NS.run();
            tsolve = secs(t1);
            if (st != NS.OPTIMAL) { std::printf("RESULT status=fail\n"); return; }
            sum = NS.totalCost();
        } else { std::fprintf(stderr, "unknown sp algo\n"); std::exit(1); }
    }
    std::printf("RESULT status=ok obj=%lld reach=%ld t_solve=%.6f t_total=%.6f\n", (long long)sum, reach, tsolve, secs(t0));
}

static void maxflow(const std::string& algo, const Instance& I) {
    auto t0 = Clock::now();
    lemon::SmartDigraph g;
    g.reserveNode(I.n); g.reserveArc(I.arcs.size() + 1);
    std::vector<lemon::SmartDigraph::Node> N(I.n);
    for (int i = 0; i < I.n; i++) N[i] = g.addNode();
    lemon::SmartDigraph::ArcMap<int64_t> cap(g);
    for (auto& a : I.arcs) cap[g.addArc(N[a.u], N[a.v])] = a.w;
    int64_t val = 0; double tsolve = 0;
    if (algo == "lemon-preflow") {
        auto t1 = Clock::now();
        lemon::Preflow<lemon::SmartDigraph, lemon::SmartDigraph::ArcMap<int64_t>> P(g, cap, N[I.s], N[I.t]);
        P.runMinCut();   // phase 1 only: max preflow value = max flow value
        tsolve = secs(t1);
        val = P.flowValue();
    } else if (algo == "lemon-ns") {
        // circulation: add t->s arc with cost -1, all other costs 0
        int64_t big = 0; for (auto& a : I.arcs) if (a.u == I.s) big += a.w;
        auto ts = g.addArc(N[I.t], N[I.s]);
        cap[ts] = big;
        lemon::SmartDigraph::ArcMap<int64_t> cost(g, 0);
        cost[ts] = -1;
        auto t1 = Clock::now();
        lemon::NetworkSimplex<lemon::SmartDigraph, int64_t, int64_t> NS(g);
        NS.upperMap(cap).costMap(cost);
        auto st = NS.run();
        tsolve = secs(t1);
        if (st != NS.OPTIMAL) { std::printf("RESULT status=fail\n"); return; }
        val = -NS.totalCost();
    } else { std::fprintf(stderr, "unknown max algo\n"); std::exit(1); }
    std::printf("RESULT status=ok obj=%lld t_solve=%.6f t_total=%.6f\n", (long long)val, tsolve, secs(t0));
}

// ---------------- Garg-Koenemann / Fleischer max concurrent flow ----------------
static void gk(const Instance& I, double eps) {
    auto t0 = Clock::now();
    int m = (int)I.arcs.size(), n = I.n;
    CSR G; G.build(n, I.arcs);
    auto groups = group_by_source(I);
    std::vector<double> len(m), x(m, 0.0);
    double delta = std::pow((double)m / (1.0 - eps), -1.0 / eps); // GK choice of delta
    for (int e = 0; e < m; e++) len[e] = delta / (double)I.arcs[e].w;
    auto Dl = [&]() { double s = 0; for (int e = 0; e < m; e++) s += len[e] * I.arcs[e].w; return s; };
    std::vector<double> dist(n);
    std::vector<int> pred(n), order; order.reserve(n);
    std::vector<double> through(n);
    long phases = 0, dijkstras = 0;
    auto t1 = Clock::now();
    double D = Dl();
    while (D < 1.0) {
        bool done = false;
        for (auto& gp : groups) {
            // remaining demand of each sink for this phase
            std::vector<double> rem(gp.sinks.size());
            for (size_t j = 0; j < rem.size(); j++) rem[j] = (double)gp.sinks[j].second;
            double remtot = (double)gp.total;
            while (remtot > 1e-9 * gp.total) {
                // Dijkstra with real lengths
                std::fill(dist.begin(), dist.end(), std::numeric_limits<double>::infinity());
                std::fill(pred.begin(), pred.end(), -1);
                order.clear();
                using P = std::pair<double, int>;
                std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
                dist[gp.s] = 0; pq.push({0, gp.s});
                while (!pq.empty()) {
                    auto [du, u] = pq.top(); pq.pop();
                    if (du != dist[u]) continue;
                    order.push_back(u);
                    for (int k = G.start[u]; k < G.start[u + 1]; k++) {
                        int e = G.arc[k]; int v = I.arcs[e].v;
                        double nd = du + len[e];
                        if (nd < dist[v]) { dist[v] = nd; pred[v] = e; pq.push({nd, v}); }
                    }
                }
                dijkstras++;
                // tree flow: push remaining sink demands up the tree
                std::fill(through.begin(), through.end(), 0.0);
                for (size_t j = 0; j < rem.size(); j++) through[gp.sinks[j].first] += rem[j];
                double sigma = 1.0;
                for (int i = (int)order.size() - 1; i > 0; i--) {
                    int v = order[i]; int e = pred[v];
                    if (through[v] > 0) {
                        sigma = std::max(sigma, through[v] / (double)I.arcs[e].w);
                        through[I.arcs[e].u] += through[v];
                    }
                }
                // route through/sigma along the tree (recompute flows top-down via the same accumulation)
                std::fill(through.begin(), through.end(), 0.0);
                for (size_t j = 0; j < rem.size(); j++) through[gp.sinks[j].first] += rem[j] / sigma;
                for (int i = (int)order.size() - 1; i > 0; i--) {
                    int v = order[i]; int e = pred[v];
                    double f = through[v];
                    if (f > 0) {
                        x[e] += f;
                        double old = len[e];
                        len[e] *= 1.0 + eps * f / (double)I.arcs[e].w;
                        D += (len[e] - old) * (double)I.arcs[e].w;
                        through[I.arcs[e].u] += f;
                    }
                }
                for (auto& r : rem) r -= r / sigma;
                remtot -= remtot / sigma;
                if (D >= 1.0) { done = true; break; }
            }
            if (done) break;
        }
        if (!done) phases++;
    }
    double tsolve = secs(t1);
    // feasible solution: the flow of complete phases scaled by max congestion.
    // (x includes a partial last phase; that only helps feasibility scaling, so use phases as demand multiple)
    double cong = 0;
    for (int e = 0; e < m; e++) cong = std::max(cong, x[e] / (double)I.arcs[e].w);
    double lambda = cong > 0 ? (double)phases / cong : 0;
    std::printf("RESULT status=ok obj=%.9g phases=%ld dijkstras=%ld eps=%g t_solve=%.6f t_total=%.6f\n",
                lambda, phases, dijkstras, eps, tsolve, secs(t0));
}

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: classic algo instance [opts]\n"); return 1; }
    std::string algo = argv[1];
    auto tr = Clock::now();
    Instance I = read_instance(argv[2]);
    std::fprintf(stderr, "read %s: n=%d m=%zu in %.3fs\n", I.kind.c_str(), I.n, I.arcs.size(), secs(tr));
    if (I.kind == "sp") sp(algo, I, argc > 3 ? atoi(argv[3]) - 1 : 0);
    else if (I.kind == "max") maxflow(algo, I);
    else if (I.kind == "mcf") gk(I, argc > 3 ? atof(argv[3]) : 0.1);
    std::printf("HWM hwm_end_kb=%ld\n", vm_hwm_kb());
    return 0;
}
