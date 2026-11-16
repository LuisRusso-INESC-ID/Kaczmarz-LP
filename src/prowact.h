// prowact.h -- own Kaczmarz-style LP prototype #2: proximal point + Hildreth row action.
//
// Outer loop (proximal point / method of multipliers on x):
//     x_{k+1} = argmin  c'x + (eps/2)||x - x_k||^2   s.t.  rlo <= Ax <= rup,  clo <= x <= cup
// Inner loop: dual coordinate ascent on the row multipliers (Hildreth 1957; Bregman row
// action).  With L(x,y) = c'x + eps/2||x-x_k||^2 - y'(Ax-b) the primal minimizer is
//     x(y) = clamp(x_k + (A'y - c)/eps, clo, cup)
// and one row step is exactly a Kaczmarz projection of x(y) onto row i when no bound is
// active:   y_i += eps (b_i - a_i'x(y)) / ||a_i||^2   (then y_i <= 0 for "<=" rows).
// Only A (CSC + CSR) and O(1) vectors are stored.  Problem is Ruiz-scaled first.
#pragma once
#include "lp.h"
#include <chrono>
#include <random>
#include <numeric>
#include <algorithm>

struct ProwResult { std::string status; double pobj, dobj, rp, gap; long sweeps; int outer; std::vector<double> x; };

struct Prow {
    const LP& L0; int m, n;
    std::vector<double> val, D, E, c, lo, up, rlo, rup;   // scaled copy
    std::vector<int> rstart, rcol, rpos;                  // CSR with pointer into CSC values
    explicit Prow(const LP& L) : L0(L), m(L.nrow), n(L.ncol) {}

    void scale() {
        val = L0.value; D.assign(n, 1.0); E.assign(m, 1.0);
        std::vector<double> rmax(m), cmax(n);
        for (int it = 0; it < 10; it++) {
            std::fill(rmax.begin(), rmax.end(), 0.0);
            for (int j = 0; j < n; j++) { double cm = 0;
                for (int p = L0.start[j]; p < L0.start[j + 1]; p++) { double a = std::fabs(val[p]); cm = std::max(cm, a); rmax[L0.index[p]] = std::max(rmax[L0.index[p]], a); }
                cmax[j] = cm > 0 ? 1.0 / std::sqrt(cm) : 1.0; D[j] *= cmax[j]; }
            for (int i = 0; i < m; i++) { rmax[i] = rmax[i] > 0 ? 1.0 / std::sqrt(rmax[i]) : 1.0; E[i] *= rmax[i]; }
            for (int j = 0; j < n; j++) for (int p = L0.start[j]; p < L0.start[j + 1]; p++) val[p] *= cmax[j] * rmax[L0.index[p]];
        }
        c.resize(n); lo.resize(n); up.resize(n); rlo.resize(m); rup.resize(m);
        for (int j = 0; j < n; j++) { c[j] = L0.c[j] * D[j]; lo[j] = L0.clo[j] / D[j]; up[j] = L0.cup[j] / D[j]; }
        for (int i = 0; i < m; i++) { rlo[i] = L0.rlo[i] * E[i]; rup[i] = L0.rup[i] * E[i]; }
        rstart.assign(m + 1, 0);
        for (size_t p = 0; p < L0.nnz(); p++) rstart[L0.index[p] + 1]++;
        for (int i = 0; i < m; i++) rstart[i + 1] += rstart[i];
        rcol.resize(L0.nnz()); rpos.resize(L0.nnz());
        std::vector<int> pos(rstart.begin(), rstart.end() - 1);
        for (int j = 0; j < n; j++) for (int p = L0.start[j]; p < L0.start[j + 1]; p++) { int q = pos[L0.index[p]]++; rcol[q] = j; rpos[q] = p; }
    }

    ProwResult solve(double tol, double tlimit, int inner_sweeps = 4) {
        auto t0 = std::chrono::steady_clock::now();
        scale();
        std::vector<double> xc(n), y(m, 0.0), aty(n, 0.0), rn(m), ax(m);
        for (int j = 0; j < n; j++) xc[j] = std::clamp(0.0, lo[j], up[j]);
        for (int i = 0; i < m; i++) { double s = 0; for (int q = rstart[i]; q < rstart[i + 1]; q++) s += val[rpos[q]] * val[rpos[q]]; rn[i] = s; }
        double cmax = 0, bmax = 1;
        for (int j = 0; j < n; j++) cmax = std::max(cmax, std::fabs(c[j]));
        for (int i = 0; i < m; i++) if (std::isfinite(rup[i])) bmax = std::max(bmax, std::fabs(rup[i]));
        double xs = bmax;   // primal scale: largest finite bound or rhs
        for (int j = 0; j < n; j++) if (std::isfinite(up[j])) xs = std::max(xs, std::fabs(up[j]));
        double eps = std::max(1e-12, cmax / xs);
        double gap_ref = INF; int since = 0;
        std::vector<int> order(m); std::iota(order.begin(), order.end(), 0);
        std::mt19937_64 rng(1);
        auto xof = [&](int j) { return std::clamp(xc[j] + (aty[j] - c[j]) / eps, lo[j], up[j]); };
        ProwResult R; R.status = "iteration_limit"; long sweeps = 0; int outer = 0;
        std::vector<double> x(n);
        while (true) {
            for (int s = 0; s < inner_sweeps; s++) {
                std::shuffle(order.begin(), order.end(), rng);
                for (int i : order) {
                    double r = 0;
                    for (int q = rstart[i]; q < rstart[i + 1]; q++) r += val[rpos[q]] * xof(rcol[q]);
                    double b = rup[i];   // equality rows and "<=" rows both target rup; sign of y enforced below
                    double ny = y[i] + eps * (b - r) / rn[i];
                    if (!std::isfinite(rlo[i])) ny = std::min(0.0, ny);         // "<=" row: y <= 0
                    double d = ny - y[i];
                    if (d != 0) { y[i] = ny; for (int q = rstart[i]; q < rstart[i + 1]; q++) aty[rcol[q]] += val[rpos[q]] * d; }
                }
                sweeps++;
            }
            for (int j = 0; j < n; j++) x[j] = xof(j);
            outer++;
            // KKT check of the LP at (x, y)
            std::fill(ax.begin(), ax.end(), 0.0);
            double pobj = 0, dobj = 0, rp = 0;
            for (int j = 0; j < n; j++) { pobj += c[j] * x[j]; for (int p = L0.start[j]; p < L0.start[j + 1]; p++) ax[L0.index[p]] += val[p] * x[j]; }
            for (int i = 0; i < m; i++) { rp = std::max(rp, std::max(rlo[i] - ax[i], ax[i] - rup[i])); dobj += rup[i] * y[i]; }
            double rd = 0;
            for (int j = 0; j < n; j++) {
                double lam = c[j] - aty[j];
                if (lam > 0) { if (std::isfinite(lo[j])) dobj += lam * lo[j]; else rd = std::max(rd, lam); }
                else if (lam < 0) { if (std::isfinite(up[j])) dobj += lam * up[j]; else rd = std::max(rd, -lam); }
            }
            R.rp = rp / bmax; R.gap = std::fabs(pobj - dobj) / (1 + std::fabs(pobj) + std::fabs(dobj));
            R.pobj = pobj; R.dobj = dobj;
            double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (outer < 200 || outer % 10 == 0) trace("prow", el, sweeps, pobj, R.rp, rd / (1 + cmax), R.gap, eps);
            if (R.rp <= tol && R.gap <= tol && rd <= tol * (1 + cmax)) { R.status = "optimal"; break; }
            if (el > tlimit) { R.status = "time_limit"; break; }
            if (getenv("PROW_TRACE") && outer % 2000 == 0) std::fprintf(stderr, "outer %d eps %.3g pobj %.6g dobj %.6g rp %.2e rd %.2e xlast %.4g\n", outer, eps, pobj, dobj, R.rp, rd, x[n-1]);
            xc = x;   // proximal center update
            // adaptive proximal weight: shrink eps when the gap stalls
            (void)since; (void)gap_ref;
        }
        R.sweeps = sweeps; R.outer = outer;
        R.x.resize(n); for (int j = 0; j < n; j++) R.x[j] = x[j] * D[j];
        return R;
    }
};
