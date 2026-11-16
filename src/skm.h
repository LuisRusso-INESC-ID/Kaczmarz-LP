// skm.h -- own Kaczmarz-type LP prototype: Sampling Kaczmarz-Motzkin (De Loera, Haddock,
// Needell 2017) applied to the primal-dual optimality system of the LP.
//
// Unknowns  (x, y, zl, zu):   x primal, y row duals, zl/zu bound multipliers (zu only for
// finite upper bounds; all our models have clo = 0).
//   P_i :  rlo_i <= a_i' x <= rup_i                       (one sparse row per LP row)
//   D_j :  c_j - A_j' y - zl_j + zu_j = 0                (one sparse row per LP column)
//   G   :  c'x - (rup'y + clo'zl - cup'zu) <= 0         (dense duality-gap row)
//   simple sets: clo <= x <= cup, y_i <= 0 for "<=" rows, zl, zu >= 0  (exact projection)
// Each step samples `beta` rows from P u D, projects onto the most violated one (Motzkin
// choice inside the sample), then projects the touched coordinates onto the simple sets.
// The dense row G is projected once every `gap_every` steps.  Memory: A in CSC and CSR
// (2 nnz) plus O(1) vectors -- no factorization, no fill.
#pragma once
#include "lp.h"
#include <chrono>
#include <random>
#include <algorithm>

struct SkmResult { std::string status; double pobj, rp, rd, gap; long iters; std::vector<double> x; };

struct Skm {
    const LP& L; int m, n;
    std::vector<int> rstart, rcol; std::vector<double> rval;   // CSR copy
    std::vector<int> uidx;                                     // index into zu, -1 if infinite ub
    int nu = 0;
    explicit Skm(const LP& l) : L(l), m(l.nrow), n(l.ncol) {
        rstart.assign(m + 1, 0);
        for (int p = 0; p < (int)L.nnz(); p++) rstart[L.index[p] + 1]++;
        for (int i = 0; i < m; i++) rstart[i + 1] += rstart[i];
        rcol.resize(L.nnz()); rval.resize(L.nnz());
        std::vector<int> pos(rstart.begin(), rstart.end() - 1);
        for (int j = 0; j < n; j++) for (int p = L.start[j]; p < L.start[j + 1]; p++) {
            int q = pos[L.index[p]]++; rcol[q] = j; rval[q] = L.value[p]; }
        uidx.assign(n, -1);
        for (int j = 0; j < n; j++) if (std::isfinite(L.cup[j])) uidx[j] = nu++;
    }

    SkmResult solve(double tol, double tlimit, int beta = 32, uint64_t seed = 1) {
        const char* trace_tag = std::getenv("SKM_TAG") ? std::getenv("SKM_TAG") : "skm";
        auto t0 = std::chrono::steady_clock::now();
        std::mt19937_64 rng(seed);
        std::vector<double> x(n, 0.0), y(m, 0.0), zl(n, 0.0), zu(nu, 0.0), ax(m, 0.0), aty(n, 0.0);
        std::vector<double> pn(m), dn(n);
        for (int i = 0; i < m; i++) { double s = 0; for (int q = rstart[i]; q < rstart[i + 1]; q++) s += rval[q] * rval[q]; pn[i] = s; }
        for (int j = 0; j < n; j++) { double s = 1.0 + (uidx[j] >= 0); for (int p = L.start[j]; p < L.start[j + 1]; p++) s += L.value[p] * L.value[p]; dn[j] = s; }
        double gn = 0;
        for (int j = 0; j < n; j++) { gn += L.c[j] * L.c[j] + L.clo[j] * L.clo[j]; if (uidx[j] >= 0) gn += L.cup[j] * L.cup[j]; }
        for (int i = 0; i < m; i++) gn += L.rup[i] * L.rup[i];
        double cn2 = 0; for (int j = 0; j < n; j++) cn2 += L.c[j] * L.c[j];
        double dn2 = gn - cn2;                          // squared norm of the dual part of the gap row
        double W = 1.0;
        if (std::getenv("SKM_BALANCE")) W = std::sqrt(dn2 / std::max(cn2, 1e-300));
        // gap value maintained incrementally: g = c'x - rup'y - clo'zl + cup'zu
        double g = 0;
        auto set_x = [&](int j, double nv) {           // update x_j keeping ax and g consistent
            double d = nv - x[j]; if (d == 0) return; x[j] = nv; g += L.c[j] * d;
            for (int p = L.start[j]; p < L.start[j + 1]; p++) ax[L.index[p]] += L.value[p] * d; };
        auto set_y = [&](int i, double nv) {
            double d = nv - y[i]; if (d == 0) return; y[i] = nv; g -= L.rup[i] * d;
            for (int q = rstart[i]; q < rstart[i + 1]; q++) aty[rcol[q]] += rval[q] * d; };
        auto prow = [&](int i) {                       // signed violation of P_i
            if (ax[i] > L.rup[i]) return ax[i] - L.rup[i];
            if (ax[i] < L.rlo[i]) return ax[i] - L.rlo[i];
            return 0.0; };
        auto drow = [&](int j) { return L.c[j] - aty[j] - zl[j] + (uidx[j] >= 0 ? zu[uidx[j]] : 0.0); };
        std::uniform_int_distribution<long> pick(0, (long)m + n - 1);
        int gap_every = std::max(1, (m + n) / 16);
        if (const char* ge = std::getenv("SKM_GAP_EVERY")) gap_every = std::max(1, atoi(ge));
        long it = 0; SkmResult R; R.status = "iteration_limit";
        double bscale = 1, cscale = 1;
        for (int i = 0; i < m; i++) if (std::isfinite(L.rup[i])) bscale = std::max(bscale, std::fabs(L.rup[i]));
        for (int j = 0; j < n; j++) cscale = std::max(cscale, std::fabs(L.c[j]));
        long epoch = (long)m + n;
        long check = std::max(1L, epoch / 8);
        while (true) {
            // --- one SKM step
            long best = -1; double bestd = 0;
            for (int k = 0; k < beta; k++) {
                long r = pick(rng); double d;
                if (r < m) { double v = prow((int)r); d = pn[r] > 0 ? v * v / pn[r] : 0; }
                else { double v = drow((int)(r - m)); d = v * v / dn[r - m]; }
                if (d > bestd) { bestd = d; best = r; }
            }
            if (best >= 0) {
                if (best < m) {
                    int i = (int)best; double v = prow(i), s = v / pn[i];
                    for (int q = rstart[i]; q < rstart[i + 1]; q++) {
                        int j = rcol[q]; set_x(j, std::clamp(x[j] - s * rval[q], L.clo[j], L.cup[j])); }
                } else {
                    int j = (int)(best - m); double v = drow(j), s = v / dn[j];
                    // constraint: c_j - A_j'y - zl_j + zu_j = 0 ; gradient wrt (y, zl, zu) = (-A_j, -1, +1)
                    for (int p = L.start[j]; p < L.start[j + 1]; p++) {
                        int i = L.index[p]; double nv = y[i] + s * L.value[p];
                        if (!std::isfinite(L.rlo[i])) nv = std::min(0.0, nv);
                        set_y(i, nv); }
                    double nzl = std::max(0.0, zl[j] + s); g -= L.clo[j] * (nzl - zl[j]); zl[j] = nzl;
                    if (uidx[j] >= 0) { double& u = zu[uidx[j]]; double nz = std::max(0.0, u - s); g += L.cup[j] * (nz - u); u = nz; }
                }
            }
            it++;
            // --- dense gap row
            if (it % gap_every == 0 && g > 0) {
                // projection in the metric diag(1/W on x, W on duals); W = 1 is the Euclidean
                // (plain Kaczmarz) projection, W = Wbal balances the primal and dual halves.
                double s = g / (W * cn2 + dn2 / W), sx = s * W, sd = s / W;
                for (int j = 0; j < n; j++) set_x(j, std::clamp(x[j] - sx * L.c[j], L.clo[j], L.cup[j]));
                for (int i = 0; i < m; i++) { double nv = y[i] + sd * L.rup[i]; if (!std::isfinite(L.rlo[i])) nv = std::min(0.0, nv); set_y(i, nv); }
                for (int j = 0; j < n; j++) { double nz = std::max(0.0, zl[j] + sd * L.clo[j]); g -= L.clo[j] * (nz - zl[j]); zl[j] = nz;
                    if (uidx[j] >= 0) { double& u = zu[uidx[j]]; double nz2 = std::max(0.0, u - sd * L.cup[j]); g += L.cup[j] * (nz2 - u); u = nz2; } }
            }
            // --- convergence check once per epoch
            if (it % check == 0) {
                double rp = 0, rd = 0;
                for (int i = 0; i < m; i++) rp = std::max(rp, std::fabs(prow(i)));
                for (int j = 0; j < n; j++) rd = std::max(rd, std::fabs(drow(j)));
                double pobj = 0; for (int j = 0; j < n; j++) pobj += L.c[j] * x[j];
                R.rp = rp / bscale; R.rd = rd / cscale; R.gap = std::fabs(g) / (1 + std::fabs(pobj)); R.pobj = pobj;
                double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                trace(trace_tag, el, it, pobj, R.rp, R.rd, R.gap);
                if (R.rp <= tol && R.rd <= tol && R.gap <= tol) { R.status = "optimal"; break; }
                if (el > tlimit) { R.status = "time_limit"; break; }
            }
        }
        R.iters = it; R.x = x;
        return R;
    }
};
