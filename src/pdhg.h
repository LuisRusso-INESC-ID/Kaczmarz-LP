// pdhg.h -- own restarted PDHG prototype (a compact re-implementation of the PDLP recipe).
//
// Saddle point:  min_{clo<=x<=cup} max_{y in Y}  c'x - y'(Ax - b)
//   equality rows: y free;  rows "Ax <= b": y <= 0.
// Ingredients: Ruiz (10 passes) + Pock-Chambolle diagonal preconditioning, step size
// 1/||A||_2 (power iteration), primal weight with smoothing at restarts, KKT-based adaptive
// restarts to the better of current/average iterate (cuPDLP rule).  Memory: the scaled copy
// of A (CSC only; A'y is a per-column dot product) plus O(1) vectors of length m and n.
#pragma once
#include "lp.h"
#include <chrono>
#include <algorithm>
#include <numeric>

struct PdhgResult { std::string status; double pobj, dobj, rp, rd, gap; long iters; int restarts; std::vector<double> x; };

struct Pdhg {
    const LP& L0;
    int m, n;
    std::vector<double> val, D, E;      // scaled values, column scale D, row scale E (A_s = E A D)
    std::vector<double> c, lo, up, b;   // scaled data
    std::vector<char> leq;              // row is "<=" (dual sign constraint y<=0)
    explicit Pdhg(const LP& L) : L0(L), m(L.nrow), n(L.ncol) {}

    void Ax(const std::vector<double>& x, std::vector<double>& out) const {
        std::fill(out.begin(), out.end(), 0.0);
        for (int j = 0; j < n; j++) { double xj = x[j]; if (xj == 0) continue;
            for (int p = L0.start[j]; p < L0.start[j + 1]; p++) out[L0.index[p]] += val[p] * xj; }
    }
    void ATy(const std::vector<double>& y, std::vector<double>& out) const {
        for (int j = 0; j < n; j++) { double s = 0;
            for (int p = L0.start[j]; p < L0.start[j + 1]; p++) s += val[p] * y[L0.index[p]];
            out[j] = s; }
    }

    void scale() {
        val = L0.value; D.assign(n, 1.0); E.assign(m, 1.0);
        std::vector<double> rmax(m), cmax(n);
        for (int it = 0; it < 10; it++) {            // Ruiz: equilibrate inf-norms
            std::fill(rmax.begin(), rmax.end(), 0.0);
            for (int j = 0; j < n; j++) { double cm = 0;
                for (int p = L0.start[j]; p < L0.start[j + 1]; p++) { double a = std::fabs(val[p]); cm = std::max(cm, a); rmax[L0.index[p]] = std::max(rmax[L0.index[p]], a); }
                cmax[j] = cm; }
            for (int j = 0; j < n; j++) { double s = cmax[j] > 0 ? 1.0 / std::sqrt(cmax[j]) : 1.0; D[j] *= s; cmax[j] = s; }
            for (int i = 0; i < m; i++) { double s = rmax[i] > 0 ? 1.0 / std::sqrt(rmax[i]) : 1.0; E[i] *= s; rmax[i] = s; }
            for (int j = 0; j < n; j++) for (int p = L0.start[j]; p < L0.start[j + 1]; p++) val[p] *= cmax[j] * rmax[L0.index[p]];
        }
        {   // Pock-Chambolle alpha=1: rows by 1/sqrt(||row||_1), cols by 1/sqrt(||col||_1)
            std::vector<double> r1(m, 0.0), c1(n, 0.0);
            for (int j = 0; j < n; j++) for (int p = L0.start[j]; p < L0.start[j + 1]; p++) { double a = std::fabs(val[p]); c1[j] += a; r1[L0.index[p]] += a; }
            for (int j = 0; j < n; j++) c1[j] = c1[j] > 0 ? 1.0 / std::sqrt(c1[j]) : 1.0;
            for (int i = 0; i < m; i++) r1[i] = r1[i] > 0 ? 1.0 / std::sqrt(r1[i]) : 1.0;
            for (int j = 0; j < n; j++) { D[j] *= c1[j]; for (int p = L0.start[j]; p < L0.start[j + 1]; p++) val[p] *= c1[j] * r1[L0.index[p]]; }
            for (int i = 0; i < m; i++) E[i] *= r1[i];
        }
        // scaled problem in x' = x / D:  c' = D c, bounds / D, b' = E b
        c.resize(n); lo.resize(n); up.resize(n); b.resize(m); leq.resize(m);
        for (int j = 0; j < n; j++) { c[j] = L0.c[j] * D[j]; lo[j] = L0.clo[j] / D[j]; up[j] = L0.cup[j] / D[j]; }
        for (int i = 0; i < m; i++) { leq[i] = !std::isfinite(L0.rlo[i]); b[i] = L0.rup[i] * E[i]; }
    }

    double opnorm() const {   // power iteration on A'A
        std::vector<double> x(n, 1.0), ax(m), atax(n);
        double nrm = 0;
        for (int it = 0; it < 50; it++) {
            double s = std::sqrt(std::inner_product(x.begin(), x.end(), x.begin(), 0.0));
            for (auto& v : x) v /= s;
            Ax(x, ax); ATy(ax, atax);
            double nn = std::sqrt(std::inner_product(atax.begin(), atax.end(), atax.begin(), 0.0));
            x = atax;
            if (std::fabs(nn - nrm) < 1e-4 * nn) { nrm = nn; break; }
            nrm = nn;
        }
        return std::sqrt(nrm);
    }

    // KKT quantities in the *scaled* space (relative measures, as in PDLP termination).
    struct Kkt { double rp, rd, gap, pobj, dobj, err; };
    Kkt kkt(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& ax, std::vector<double>& aty, double w) const {
        Ax(x, ax); ATy(y, aty);
        double rp2 = 0, rd2 = 0, pobj = 0, dobj = 0, bn = 0, cn = 0;
        for (int i = 0; i < m; i++) {
            double r = b[i] - ax[i];
            if (leq[i]) r = std::min(0.0, r);    // only violation Ax > b counts
            rp2 += r * r; bn = std::max(bn, std::fabs(b[i])); dobj += b[i] * y[i];
        }
        for (int j = 0; j < n; j++) {
            double lam = c[j] - aty[j];          // reduced cost
            pobj += c[j] * x[j]; cn = std::max(cn, std::fabs(c[j]));
            double rdj = 0;
            if (lam > 0) { if (std::isfinite(lo[j])) dobj += lo[j] * lam; else rdj = lam; }
            else if (lam < 0) { if (std::isfinite(up[j])) dobj += up[j] * lam; else rdj = lam; }
            rd2 += rdj * rdj;
        }
        Kkt k; k.rp = std::sqrt(rp2) / (1 + bn); k.rd = std::sqrt(rd2) / (1 + cn);
        k.pobj = pobj; k.dobj = dobj;
        k.gap = std::fabs(pobj - dobj) / (1 + std::fabs(pobj) + std::fabs(dobj));
        k.err = std::sqrt(w * w * rp2 + rd2 / (w * w) + (pobj - dobj) * (pobj - dobj));
        return k;
    }

    PdhgResult solve(double tol, double tlimit, long maxit = 100000000) {
        auto t0 = std::chrono::steady_clock::now();
        scale();
        double nrmA = opnorm();
        double eta = 0.998 / nrmA;
        double cn = std::sqrt(std::inner_product(c.begin(), c.end(), c.begin(), 0.0));
        double bnorm = std::sqrt(std::inner_product(b.begin(), b.end(), b.begin(), 0.0));
        double w = (cn > 1e-10 && bnorm > 1e-10) ? cn / bnorm : 1.0;   // primal weight
        std::vector<double> x(n, 0.0), y(m, 0.0), xn(n), yn(m), xa(n, 0.0), ya(m, 0.0), xr(n), yr(m);
        for (int j = 0; j < n; j++) x[j] = std::clamp(0.0, lo[j], up[j]);
        std::vector<double> ax(m), aty(n), xbar(n), axbar(m);
        xr = x; yr = y;
        double wsum = 0; long it = 0, it_restart = 0; int restarts = 0;
        Kkt kr = kkt(x, y, ax, aty, w), klast = kr;
        double klast_err = kr.err;
        ATy(y, aty);
        PdhgResult R; R.status = "iteration_limit";
        Kkt kbest = kr;
        while (it < maxit) {
            double tau = eta / w, sig = eta * w;
            // x step
            for (int j = 0; j < n; j++) xn[j] = std::clamp(x[j] - tau * (c[j] - aty[j]), lo[j], up[j]);
            for (int j = 0; j < n; j++) xbar[j] = 2 * xn[j] - x[j];
            Ax(xbar, axbar);
            for (int i = 0; i < m; i++) { double v = y[i] + sig * (b[i] - axbar[i]); yn[i] = leq[i] ? std::min(0.0, v) : v; }
            x.swap(xn); y.swap(yn);
            ATy(y, aty);
            it++;
            // running average since last restart
            wsum += 1.0; double a = 1.0 / wsum;
            for (int j = 0; j < n; j++) xa[j] += a * (x[j] - xa[j]);
            for (int i = 0; i < m; i++) ya[i] += a * (y[i] - ya[i]);
            if (it % 64 == 0) {
                std::vector<double> t1(m), t2(n);
                Kkt kc = kkt(x, y, t1, t2, w), ka = kkt(xa, ya, t1, t2, w);
                bool useAvg = ka.err < kc.err;
                Kkt kcand = useAvg ? ka : kc;
                kbest = kcand;
                double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                trace("pdhg", el, it, kcand.pobj, kcand.rp, kcand.rd, kcand.gap, w);
                if (std::max({kcand.rp, kcand.rd, kcand.gap}) <= tol) {
                    if (useAvg) { x = xa; y = ya; }
                    R.status = "optimal"; break;
                }
                if (el > tlimit) { if (useAvg) { x = xa; y = ya; } R.status = "time_limit"; break; }
                long since = it - it_restart;
                bool restart = kcand.err <= 0.2 * kr.err ||
                               (kcand.err <= 0.8 * kr.err && kcand.err > klast_err) ||
                               since >= 0.36 * it;
                klast_err = kcand.err;
                if (restart) {
                    if (useAvg) { x = xa; y = ya; }
                    // primal weight update (PDLP, smoothing 0.5)
                    double dx = 0, dy = 0;
                    for (int j = 0; j < n; j++) dx += (x[j] - xr[j]) * (x[j] - xr[j]);
                    for (int i = 0; i < m; i++) dy += (y[i] - yr[i]) * (y[i] - yr[i]);
                    dx = std::sqrt(dx); dy = std::sqrt(dy);
                    if (dx > 1e-10 && dy > 1e-10) w = std::exp(0.5 * std::log(dy / dx) + 0.5 * std::log(w));
                    xr = x; yr = y; xa = x; ya = y; wsum = 0; it_restart = it; restarts++;
                    kr = kkt(x, y, t1, t2, w); klast_err = kr.err;
                    ATy(y, aty);
                }
            }
        }
        std::vector<double> t1(m), t2(n);
        Kkt kf = kkt(x, y, t1, t2, w);
        R.rp = kf.rp; R.rd = kf.rd; R.gap = kf.gap; R.pobj = kf.pobj; R.dobj = kf.dobj;
        R.iters = it; R.restarts = restarts;
        R.x.resize(n);
        for (int j = 0; j < n; j++) R.x[j] = x[j] * D[j];      // unscale
        return R;
    }
};
