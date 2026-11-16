// hybrid.h -- EXPLORATORY, not part of the benchmark (dominated by PDHG and the proximal Hildreth
// method; see README and results/FINDINGS.md).  Own hybrid prototype: Kaczmarz-type feasibility repair + block gradient steps
// that keep the equality rows satisfied.
//
//   min c'x  s.t.  A_E x = b_E (equality rows),  A_I x <= b_I (inequality rows),  lo <= x <= up
//
// State: a primal point x, the residual vector r = Ax (kept incrementally), and a FIFO queue of
// rows that may be violated (with an in-queue flag).  Changing x_j re-checks only the rows of
// column j (CSC), so neither phase ever rescans all constraints.
//
// Repeat:
//   1. Optimization step on a "middle ground" block: grow a ball R of equality rows around a
//      random column (rows are adjacent when they share a column), take the block columns B whose
//      equality rows all lie in R, and move along the projected gradient
//          d_B = -(c_B - M' y),  (M M' + delta I) y = M c_B,  M = A[R, B]
//      so that M d_B = 0: every equality row is unchanged (up to round-off).  For network
//      matrices these directions are local circulations/augmenting paths.  Columns at a bound
//      that d would push outward are dropped and d recomputed.  Ratio test on the bounds gives
//      the step alpha_b (exact line search: the objective is linear along d).
//   2. Repair: the step can violate inequality rows (and equalities drift by round-off).
//        variant "mem":  memoryless Kaczmarz projections (+ clamp) on queued violated rows;
//                        the step is a diminishing fraction of alpha_b once steps start creating
//                        violations (superiorization-style).
//        variant "prox": Hildreth row action with dual memory, which returns the Euclidean
//                        projection of the stepped point onto the feasible set, i.e. one exact
//                        proximal-point (projected-gradient) step  x+ = P(x + alpha_b d).
//   When blocks stop producing progress, the ball size doubles (eventually the whole matrix,
//   which is the plain projected gradient).  There is no dual certificate: the run stops at the
//   time limit or when even the full-size block gives no progress.
#pragma once
#include "lp.h"
#include <chrono>
#include <random>
#include <algorithm>
#include <numeric>

struct HybridResult { std::string status; double pobj, rp; long blocks, repairs; int rmax; std::vector<double> x; };

struct Hybrid {
    const LP& L; int m, n; bool prox; bool cgm;
    std::vector<int> rstart, rcol; std::vector<double> rval, rn;   // CSR + squared row norms
    std::vector<char> eq;
    std::vector<double> x, r;
    double obj = 0;
    // violated-row queue
    std::vector<int> q; size_t qh = 0; std::vector<char> inq;
    // Hildreth state (prox variant): sparse duals, A'y, recorded centre z
    std::vector<double> y, aty, z; std::vector<char> ztouched, yon; std::vector<int> ztl, ytl;
    double tol_eq, tol_in;
    std::vector<int> lr, lc; std::vector<char> frz;   // CG repair: local row/column maps, frozen columns

    Hybrid(const LP& l, bool p, bool cg = false) : L(l), m(l.nrow), n(l.ncol), prox(p && !cg), cgm(cg) {
        rstart.assign(m + 1, 0);
        for (size_t k = 0; k < L.nnz(); k++) rstart[L.index[k] + 1]++;
        for (int i = 0; i < m; i++) rstart[i + 1] += rstart[i];
        rcol.resize(L.nnz()); rval.resize(L.nnz());
        std::vector<int> pos(rstart.begin(), rstart.end() - 1);
        for (int j = 0; j < n; j++) for (int k = L.start[j]; k < L.start[j + 1]; k++) { int t = pos[L.index[k]]++; rcol[t] = j; rval[t] = L.value[k]; }
        rn.assign(m, 0.0); eq.assign(m, 0);
        double bmax = 1;
        for (int i = 0; i < m; i++) {
            for (int t = rstart[i]; t < rstart[i + 1]; t++) rn[i] += rval[t] * rval[t];
            eq[i] = L.rlo[i] == L.rup[i];
            if (std::isfinite(L.rup[i])) bmax = std::max(bmax, std::fabs(L.rup[i]));
        }
        tol_eq = 1e-7 * bmax; tol_in = 1e-7 * bmax;
        x.assign(n, 0.0); r.assign(m, 0.0); inq.assign(m, 0); lr.assign(m, -1); lc.assign(n, -1); frz.assign(n, 0);
        y.assign(m, 0.0); aty.assign(n, 0.0); z.assign(n, 0.0); ztouched.assign(n, 0); yon.assign(m, 0);
        for (int j = 0; j < n; j++) x[j] = std::clamp(0.0, L.clo[j], L.cup[j]);
        for (int j = 0; j < n; j++) { obj += L.c[j] * x[j]; for (int k = L.start[j]; k < L.start[j + 1]; k++) r[L.index[k]] += L.value[k] * x[j]; }
    }

    double viol(int i) const {            // signed violation (0 if satisfied within tolerance)
        double v = r[i] - L.rup[i];
        if (eq[i]) return std::fabs(v) > tol_eq ? v : 0.0;
        return v > tol_in ? v : 0.0;
    }
    bool wants(int i) const {             // should row i be (re)processed?
        if (viol(i) != 0) return true;
        return prox && yon[i] && !eq[i] && r[i] < L.rup[i] - tol_in && y[i] < 0;   // Hildreth: relax a slack row
    }
    void push(int i) { if (!inq[i] && wants(i)) { inq[i] = 1; q.push_back(i); } }
    void setx(int j, double v) {
        double d = v - x[j]; if (d == 0) return;
        x[j] = v; obj += L.c[j] * d;
        for (int k = L.start[j]; k < L.start[j + 1]; k++) { r[L.index[k]] += L.value[k] * d; push(L.index[k]); }
    }
    int pop() { int i = q[qh++]; inq[i] = 0; if (qh > 4096 && qh * 2 > q.size()) { q.erase(q.begin(), q.begin() + qh); qh = 0; } return i; }
    bool qempty() const { return qh >= q.size(); }

    // ---- repair, memoryless Kaczmarz (relaxation omega, then clamp)
    long repair_mem(long cap, double omega = 1.0) {
        long it = 0;
        while (!qempty() && it < cap) {
            int i = pop(); double v = viol(i); if (v == 0) continue;
            double s = omega * v / rn[i];
            for (int t = rstart[i]; t < rstart[i + 1]; t++) { int j = rcol[t]; setx(j, std::clamp(x[j] - s * rval[t], L.clo[j], L.cup[j])); }
            it++;
        }
        return it;
    }
    // ---- repair, Hildreth: Euclidean projection of the current point (taken as centre z)
    long repair_prox(long cap) {
        long it = 0;
        auto touch = [&](int j) { if (!ztouched[j]) { ztouched[j] = 1; z[j] = x[j]; ztl.push_back(j); } };   // pre-seeded entries keep their z
        while (!qempty() && it < cap) {
            int i = pop();
            if (!wants(i)) continue;
            double ny = y[i] + (L.rup[i] - r[i]) / rn[i];
            if (!eq[i]) ny = std::min(0.0, ny);
            double dy = ny - y[i]; if (dy == 0) continue;
            if (!yon[i]) { yon[i] = 1; ytl.push_back(i); }
            y[i] = ny;
            for (int t = rstart[i]; t < rstart[i + 1]; t++) {
                int j = rcol[t]; touch(j);
                aty[j] += rval[t] * dy;
                setx(j, std::clamp(z[j] + aty[j], L.clo[j], L.cup[j]));
            }
            it++;
        }
        // forget the duals: the next repair projects a new point
        for (int j : ztl) { ztouched[j] = 0; aty[j] = 0; } ztl.clear();
        for (int i : ytl) { yon[i] = 0; y[i] = 0; } ytl.clear();
        return it;
    }

    // ---- repair by CG: least-norm correction  dx = M' w,  (M M' + delta I) w = e,  where M is A
    // restricted to the active rows (equality rows of `rows`, plus violated inequality rows met
    // on the way) and the free columns of `cols`.  After clamping to the bounds, columns that hit
    // a bound are frozen and the solve is repeated (a few rounds).  Returns work in row-ops.
    long repair_cg(const std::vector<int>& rows, const std::vector<int>& cols, int rounds = 12) {
        long work = 0;
        std::vector<int> A, C;
        for (int j : cols) frz[j] = 0;
        for (int rd = 0; rd < rounds; rd++) {
            A.clear(); C.clear();
            for (int i : rows) if (eq[i] && lr[i] < 0) { lr[i] = (int)A.size(); A.push_back(i); }
            for (int j : cols) if (!frz[j]) {
                lc[j] = (int)C.size(); C.push_back(j);
                for (int k = L.start[j]; k < L.start[j + 1]; k++) { int i = L.index[k];
                    if (!eq[i] && lr[i] < 0 && r[i] > L.rup[i] + tol_in) { lr[i] = (int)A.size(); A.push_back(i); } }
            }
            int na = (int)A.size(), nc = (int)C.size();
            std::vector<double> e(na); double emax = 0, diag = 1;
            for (int h = 0; h < na; h++) { e[h] = L.rup[A[h]] - r[A[h]]; emax = std::max(emax, std::fabs(e[h])); diag = std::max(diag, rn[A[h]]); }
            bool done = emax <= 0.5 * tol_eq || nc == 0;
            if (!done) {
                auto Mv = [&](const std::vector<double>& v, std::vector<double>& out) {
                    for (int h = 0; h < na; h++) { double t = 0; int i = A[h];
                        for (int u = rstart[i]; u < rstart[i + 1]; u++) { int j = rcol[u]; if (lc[j] >= 0) t += rval[u] * v[lc[j]]; }
                        out[h] = t; } };
                auto MTw = [&](const std::vector<double>& w, std::vector<double>& out) {
                    std::fill(out.begin(), out.end(), 0.0);
                    for (int h = 0; h < na; h++) { int i = A[h];
                        for (int u = rstart[i]; u < rstart[i + 1]; u++) { int j = rcol[u]; if (lc[j] >= 0) out[lc[j]] += rval[u] * w[h]; } } };
                double delta = 1e-12 * diag;
                std::vector<double> w(na, 0.0), rr = e, p = e, ap(na), t(nc);
                double g = std::inner_product(rr.begin(), rr.end(), rr.begin(), 0.0), g0 = g;
                for (int it = 0; it < 10 * na + 100 && g > 1e-28 * std::max(g0, 1e-300); it++) {
                    MTw(p, t); Mv(t, ap); work += na;
                    for (int h = 0; h < na; h++) ap[h] += delta * p[h];
                    double pap = std::inner_product(p.begin(), p.end(), ap.begin(), 0.0);
                    if (!(pap > 0)) break;
                    double a = g / pap;
                    for (int h = 0; h < na; h++) { w[h] += a * p[h]; rr[h] -= a * ap[h]; }
                    double g2 = std::inner_product(rr.begin(), rr.end(), rr.begin(), 0.0);
                    for (int h = 0; h < na; h++) p[h] = rr[h] + (g2 / g) * p[h];
                    g = g2;
                }
                MTw(w, t);
                for (int b = 0; b < nc; b++) {
                    int j = C[b]; double v = x[j] + t[b], cv = std::clamp(v, L.clo[j], L.cup[j]);
                    if (cv != v) frz[j] = 1;
                    setx(j, cv);
                }
            }
            for (int i : A) lr[i] = -1;
            for (int j : C) lc[j] = -1;
            if (done) break;
        }
        for (int j : cols) frz[j] = 0;
        return work;
    }
    // drop queued rows that are satisfied; true if nothing violated remains
    bool settle() { while (!qempty()) { int i = pop(); if (viol(i) != 0) { inq[i] = 1; q.push_back(i); return false; } } return true; }
    std::vector<int> allrows, allcols;
    long repair_cg_global() {
        if (allrows.empty()) { allrows.resize(m); std::iota(allrows.begin(), allrows.end(), 0); allcols.resize(n); std::iota(allcols.begin(), allcols.end(), 0); }
        long w = repair_cg(allrows, allcols, 30);
        settle();
        return w;
    }
    long repair(long cap) { return prox ? repair_prox(cap) : repair_mem(cap); }

    double max_violation() const {
        double v = 0, s = 1;
        for (int i = 0; i < m; i++) {
            double e = r[i] - L.rup[i];
            v = std::max(v, eq[i] ? std::fabs(e) : std::max(0.0, e));
            if (std::isfinite(L.rup[i])) s = std::max(s, std::fabs(L.rup[i]));
        }
        return v / s;
    }

    HybridResult solve(double tlimit, int rmax0 = 64, uint64_t seed = 1, const char* tag = "hybrid") {
        auto t0 = std::chrono::steady_clock::now();
        auto el = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
        std::mt19937_64 rng(seed);
        // initial feasibility
        for (int i = 0; i < m; i++) push(i);
        long repairs = 0;
        if (cgm) { repairs += repair_cg_global(); while (!qempty() && el() < tlimit) repairs += repair_mem(1L << 20); }
        while (!qempty() && el() < tlimit) repairs += repair(1L << 20);
        int neq = 0; for (int i = 0; i < m; i++) neq += eq[i];
        int rmax = std::min(rmax0, std::max(1, neq));
        // block work arrays (stamped marks avoid clearing)
        std::vector<int> rmark(m, -1), cmark(n, -1), cloc(n, -1), R, B, F;
        std::vector<double> d, cF, yv, rv, pv, ap;
        long blocks = 0, stall = 0, violsteps = 0; double sstep = 0;
        double last_trace = -1;
        HybridResult res; res.status = "time_limit";
        std::uniform_int_distribution<int> pickc(0, n - 1);
        std::vector<int> csupp; for (int j = 0; j < n; j++) if (L.c[j] != 0) csupp.push_back(j);
        trace(tag, el(), 0, obj, max_violation(), NAN, NAN, (double)repairs);
        while (el() < tlimit) {
            int stamp = (int)blocks;
            // ---- grow a ball of equality rows around a random column
            R.clear();
            int seedc = pickc(rng);
            if (!csupp.empty() && (blocks & 1)) seedc = csupp[rng() % csupp.size()];
            for (int k = L.start[seedc]; k < L.start[seedc + 1]; k++) { int i = L.index[k]; if (eq[i] && rmark[i] != stamp) { rmark[i] = stamp; R.push_back(i); } }
            for (size_t h = 0; h < R.size() && (int)R.size() < rmax; h++)
                for (int t = rstart[R[h]]; t < rstart[R[h] + 1] && (int)R.size() < rmax; t++) {
                    int j = rcol[t];
                    for (int k = L.start[j]; k < L.start[j + 1] && (int)R.size() < rmax; k++) { int i = L.index[k]; if (eq[i] && rmark[i] != stamp) { rmark[i] = stamp; R.push_back(i); } }
                }
            // ---- block columns: all equality rows inside R
            B.clear();
            for (int i : R) for (int t = rstart[i]; t < rstart[i + 1]; t++) {
                int j = rcol[t]; if (cmark[j] == stamp) continue; cmark[j] = stamp;
                bool in = true;
                for (int k = L.start[j]; k < L.start[j + 1]; k++) if (eq[L.index[k]] && rmark[L.index[k]] != stamp) { in = false; break; }
                if (in) B.push_back(j);
            }
            F = B;
            // ---- projected gradient with bound-active-set refinement
            double dd = 0;
            for (int pass = 0; pass < 8 && !F.empty(); pass++) {
                int nf = (int)F.size(), nr = (int)R.size();
                for (int a = 0; a < nf; a++) cloc[F[a]] = a;
                // M v: for row R[h], sum over its columns in F
                auto Mv = [&](const std::vector<double>& v, std::vector<double>& out) {
                    for (int h = 0; h < nr; h++) { double s = 0; int i = R[h];
                        for (int t = rstart[i]; t < rstart[i + 1]; t++) { int j = rcol[t]; if (cloc[j] >= 0) s += rval[t] * v[cloc[j]]; }
                        out[h] = s; } };
                auto MTw = [&](const std::vector<double>& w, std::vector<double>& out) {
                    std::fill(out.begin(), out.end(), 0.0);
                    for (int h = 0; h < nr; h++) { int i = R[h];
                        for (int t = rstart[i]; t < rstart[i + 1]; t++) { int j = rcol[t]; if (cloc[j] >= 0) out[cloc[j]] += rval[t] * w[h]; } } };
                cF.resize(nf); for (int a = 0; a < nf; a++) cF[a] = L.c[F[a]];
                // CG on (M M' + delta I) y = M c
                std::vector<double> rhs(nr), Mt(nf);
                Mv(cF, rhs);
                double diag = 0; for (int i : R) diag = std::max(diag, rn[i]);
                double delta = 1e-12 * std::max(1.0, diag);
                yv.assign(nr, 0.0); rv = rhs; pv = rv; ap.assign(nr, 0.0);
                double rr = std::inner_product(rv.begin(), rv.end(), rv.begin(), 0.0), rr0 = rr;
                for (int it = 0; it < 4 * nr + 50 && rr > 1e-24 * std::max(rr0, 1e-300); it++) {
                    MTw(pv, Mt); Mv(Mt, ap);
                    for (int h = 0; h < nr; h++) ap[h] += delta * pv[h];
                    double pap = std::inner_product(pv.begin(), pv.end(), ap.begin(), 0.0);
                    if (pap <= 0) break;
                    double a = rr / pap;
                    for (int h = 0; h < nr; h++) { yv[h] += a * pv[h]; rv[h] -= a * ap[h]; }
                    double rr2 = std::inner_product(rv.begin(), rv.end(), rv.begin(), 0.0);
                    for (int h = 0; h < nr; h++) pv[h] = rv[h] + (rr2 / rr) * pv[h];
                    rr = rr2;
                }
                MTw(yv, Mt);
                d.assign(nf, 0.0);
                for (int a = 0; a < nf; a++) d[a] = -(cF[a] - Mt[a]);
                // drop columns at a bound that d pushes outward
                std::vector<int> keep; keep.reserve(nf);
                double dn = 0; for (double v : d) dn = std::max(dn, std::fabs(v));
                double eps = 1e-12 * std::max(1.0, dn);
                for (int a = 0; a < nf; a++) {
                    int j = F[a];
                    bool blocked = (d[a] < -eps && x[j] <= L.clo[j] + 1e-12) || (d[a] > eps && x[j] >= L.cup[j] - 1e-12);
                    if (!blocked) keep.push_back(j);
                }
                for (int a = 0; a < nf; a++) cloc[F[a]] = -1;
                if ((int)keep.size() == nf) {
                    for (int a = 0; a < nf; a++) cloc[F[a]] = a;   // keep mapping for the step below
                    dd = 0; for (double v : d) dd += v * v;
                    break;
                }
                F.swap(keep); dd = 0; d.clear();
            }
            blocks++;
            double progress = 0;
            if (!d.empty() && dd > 1e-24) {
                // ratio test on the bounds
                double ab = INF;
                for (size_t a = 0; a < F.size(); a++) {
                    int j = F[a];
                    if (d[a] > 0 && std::isfinite(L.cup[j])) ab = std::min(ab, (L.cup[j] - x[j]) / d[a]);
                    else if (d[a] < 0 && std::isfinite(L.clo[j])) ab = std::min(ab, (L.clo[j] - x[j]) / d[a]);
                }
                if (!std::isfinite(ab)) ab = 1.0;   // unbounded direction within the block
                // step length (in units of d): never shorter than the ratio-test step; beyond it the
                // step bends at the bounds and the repair projects back (avoids jamming).
                if (sstep <= 0) sstep = 4.0 * ab;
                double alpha = std::max(ab, (prox || cgm) ? sstep : sstep / std::sqrt(1.0 + violsteps / 50.0));
                if (!prox && !cgm) sstep = std::max(sstep, ab);
                double o0 = obj;
                if (prox) {
                    // centre z = x + alpha d (unclamped) on the block; Hildreth then returns P(z)
                    for (size_t a = 0; a < F.size(); a++) { int j = F[a];
                        if (!ztouched[j]) { ztouched[j] = 1; ztl.push_back(j); }
                        z[j] = x[j] + alpha * d[a];
                        setx(j, std::clamp(z[j], L.clo[j], L.cup[j])); }
                } else {
                    for (size_t a = 0; a < F.size(); a++) { int j = F[a]; setx(j, std::clamp(x[j] + alpha * d[a], L.clo[j], L.cup[j])); }
                }
                bool created = !qempty();
                long work = 0;
                if (created && cgm) {
                    violsteps++;
                    work += repair_cg(R, F);                       // local least-norm correction on the block
                    if (!settle()) work += repair_cg_global();     // fall back to a global solve
                    if (!settle()) { long cap = 200L * m + 10000; while (!qempty() && el() < tlimit) work += repair_mem(cap); }
                    repairs += work;
                } else if (created) {
                    violsteps++;
                    long cap = 200L * m + 10000;        // repair to completion (generous cap)
                    while (!qempty() && work < 20 * cap && el() < tlimit) work += repair(cap);
                    repairs += work;
                } else if (prox) { for (int j : ztl) ztouched[j] = 0; ztl.clear(); }
                progress = o0 - obj;
                // adapt the step: grow while the step is cheap to repair and keeps most of the
                // predicted decrease; shrink when the repair is expensive or the decrease is lost
                bool cheap = cgm ? true : work <= 4L * (long)(F.size() + R.size()) + 64;
                if (progress >= 0.25 * alpha * dd && cheap) sstep = 2.0 * alpha;
                else sstep = std::max(0.5 * alpha, 1e-300);
            }
            for (int j : F) cloc[j] = -1;   // always clear the local column map
            if (progress > 1e-12 * (1.0 + std::fabs(obj))) stall = 0;
            else if (++stall > std::max(50L, 4L * neq / rmax)) {
                if (rmax >= neq) { res.status = "stalled"; break; }
                rmax = std::min(neq, 2 * rmax); stall = 0;
            }
            double now = el();
            if (now - last_trace > 0.002 || blocks < 100) { trace(tag, now, blocks, obj, max_violation(), NAN, NAN, (double)repairs); last_trace = now; }
        }
        trace(tag, el(), blocks, obj, max_violation(), NAN, NAN, (double)repairs);
        res.pobj = obj; res.rp = max_violation(); res.blocks = blocks; res.repairs = repairs; res.rmax = rmax; res.x = x;
        return res;
    }
};
