// bk_driver.cpp -- lpgraph driver for Boykov and Kolmogorov's maxflow 3.04 (GPL; fetched by
// scripts/fetch_flow_codes.sh and not distributed here; the bin/flow-bk binary is GPL).
//
// usage: flow-bk instance.max
// Streams the DIMACS file straight into the BK graph (no intermediate copy of the instance):
// an arc u->v between inner nodes becomes add_edge(u, v, cap, 0); arcs s->v and v->t become
// terminal capacities; arcs into s or out of t cannot carry s-t flow and are dropped; s->t arcs
// add their capacity directly. t_solve times maxflow(), which returns the max-flow value and a
// flow; reading and graph construction are excluded. Output: RESULT/HWM lines as bin/classic.
#include "graph.cpp"          // graph construction; maxflow() comes from bk_maxflow.cpp
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

typedef Graph<int, long long, long long> G;   // captype, tcaptype, flowtype
template class Graph<int, long long, long long>;  // members defined in graph.cpp

static long vm_hwm_kb() {
    std::ifstream f("/proc/self/status"); std::string l;
    while (std::getline(f, l)) if (l.rfind("VmHWM:", 0) == 0) return std::stol(l.substr(6));
    return -1;
}
using Clock = std::chrono::steady_clock;
static double secs(Clock::time_point a) { return std::chrono::duration<double>(Clock::now() - a).count(); }

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: flow-bk instance.max\n"); return 1; }
    FILE* f = std::fopen(argv[1], "r");
    if (!f) { std::perror(argv[1]); return 1; }
    auto t0 = Clock::now();
    char line[512];
    long n = 0, m = 0, s = -1, t = -1;
    long long direct = 0;
    G* g = nullptr;
    while (std::fgets(line, sizeof line, f)) {
        if (line[0] == 'p') {
            char kind[16];
            if (std::sscanf(line, "p %15s %ld %ld", kind, &n, &m) != 3 || std::strcmp(kind, "max")) {
                std::printf("RESULT status=fail(not a max instance)\n"); return 1; }
            g = new G((int)n, (int)m);
            g->add_node((int)n);
        } else if (line[0] == 'n') {
            long v; char w;
            if (std::sscanf(line, "n %ld %c", &v, &w) == 2) { if (w == 's') s = v - 1; else if (w == 't') t = v - 1; }
        } else if (line[0] == 'a') {
            long u, v; long long c;
            if (std::sscanf(line, "a %ld %ld %lld", &u, &v, &c) != 3 || !g || s < 0 || t < 0) continue;
            u--; v--;
            if (u == v || v == s || u == t || c <= 0) continue;
            if (u == s && v == t) direct += c;
            else if (u == s) g->add_tweights((int)v, c, 0);
            else if (v == t) g->add_tweights((int)u, 0, c);
            else g->add_edge((int)u, (int)v, (int)c, 0);
        }
    }
    std::fclose(f);
    if (!g) { std::printf("RESULT status=fail(no p line)\n"); return 1; }
    auto t1 = Clock::now();
    long long value = g->maxflow() + direct;
    double ts = secs(t1);
    std::printf("RESULT status=ok obj=%lld t_solve=%.6f t_flow=%.6f t_total=%.6f\n", value, ts, ts, secs(t0));
    std::printf("HWM hwm_end_kb=%ld\n", vm_hwm_kb());
    return 0;
}
