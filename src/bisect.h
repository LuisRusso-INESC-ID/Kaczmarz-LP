// bisect.h -- own prototype: primal-dual bisection on the objective level, with interleaved
// Kaczmarz (SKM) runs on a primal slice and a dual slice.
//
// LP (min form):  min c'x  s.t.  rlo <= Ax <= rup (equality or "<=" rows),  clo <= x <= cup
// Dual:           max rup'y + clo'zl - cup'zu  s.t.  A'y + zl - zu = c,  y_i <= 0 on "<=" rows,
//                 zl, zu >= 0 (zu only where cup is finite)
//
// Keep bounds L <= opt <= U and a level t in (L, U).  By strong duality exactly one of
//     primal slice  P(t) = { x : rows, box, c'x <= t }
//     dual slice    D(t) = { (y,zl,zu) : dual rows, signs, dual objective >= t }
// is nonempty (both only at t = opt), so running Kaczmarz on both, interleaved, always makes one
// side succeed; its success certifies that the other side is infeasible -- no Farkas certificate
// is needed.  Optimizations (as discussed):
//   1. warm starts: both points persist; when t moves, each objective row is re-projected at once
//   2. jumps: primal success sets U = c'x (may beat t); the dual bound below may beat t
//   3. Lagrangian bound: with finite boxes (uhint fills infinite ones with valid implied bounds),
//      ANY sign-feasible y gives  L >= rup'y + sum_j min_{x_j in box} (c - A'y)_j x_j,  so the dual
//      side succeeds as soon as that bound reaches t -- dual feasibility is not required
//   4. adaptive split: t = L + theta (U - L); theta moves toward the side that keeps failing
//   5. over-relaxed projections (factor 1.5) on the two dense objective rows
//   6. primal feasibility tolerance tied to the current relative gap
//   7. both systems run on the Ruiz/Pock-Chambolle equilibrated LP (objective values unchanged)
// Each side is an SKM run (Motzkin choice among `beta` sampled sparse rows), the dense objective
// row is projected every `gap_every` steps.  Primal upper bounds are approximate: they come from
// points whose residual is within the (gap-tied) tolerance.
#pragma once
#include "lp.h"
#include <chrono>
#include <random>
#include <algorithm>

struct BisectResult { std::string status; double U, L, rp; long iters; int levels; std::vector<double> x; };

struct Bisect {
    const LP& L; int m, n;
    std::vector<int> rstart, rcol; std::vector<double> rval, pn, dn;
    std::vector<int> uidx; int nu = 0;
    std::vector<double> ubnd;          // finite upper bound used in the Lagrangian bound
    std::vector<double> rowunscale;    // 1/E_i if the LP was row-scaled: violations are measured in original units
    explicit Bisect(const LP& l) : L(l), m(l.nrow), n(l.ncol) {
        rstart.assign(m + 1, 0);
        for (size_t p = 0; p < L.nnz(); p++) rstart[L.index[p] + 1]++;
        for (int i = 0; i < m; i++) rstart[i + 1] += rstart[i];
        rcol.resize(L.nnz()); rval.resize(L.nnz());
        std::vector<int> pos(rstart.begin(), rstart.end() - 1);
        for (int j = 0; j < n; j++) for (int p = L.start[j]; p < L.start[j + 1]; p++) { int q = pos[L.index[p]]++; rcol[q] = j; rval[q] = L.value[p]; }
        uidx.assign(n, -1);
        for (int j = 0; j < n; j++) if (std::isfinite(L.cup[j])) uidx[j] = nu++;
        pn.assign(m, 0.0); dn.assign(n, 0.0);
        for (int i = 0; i < m; i++) for (int q = rstart[i]; q < rstart[i + 1]; q++) pn[i] += rval[q] * rval[q];
        for (int j = 0; j < n; j++) { dn[j] = 1.0 + (uidx[j] >= 0); for (int p = L.start[j]; p < L.start[j + 1]; p++) dn[j] += L.value[p] * L.value[p]; }
        ubnd.resize(n);
        for (int j = 0; j < n; j++) ubnd[j] = std::isfinite(L.cup[j]) ? L.cup[j] : (L.uhint.empty() ? INF : L.uhint[j]);
    }

    BisectResult solve(double tol, double tlimit, int beta = 32, uint64_t seed = 1) {
        auto t0 = std::chrono::steady_clock::now();
        auto el = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
        std::mt19937_64 rng(seed);
        // ---------------- primal state
        std::vector<double> x(n, 0.0), ax(m, 0.0);
        double cx = 0;
        for (int j = 0; j < n; j++) x[j] = std::clamp(0.0, L.clo[j], L.cup[j]);
        for (int j = 0; j < n; j++) { cx += L.c[j] * x[j]; for (int p = L.start[j]; p < L.start[j + 1]; p++) ax[L.index[p]] += L.value[p] * x[j]; }
        double cn2 = 0; for (int j = 0; j < n; j++) cn2 += L.c[j] * L.c[j];
        auto prow = [&](int i) { if (ax[i] > L.rup[i]) return ax[i] - L.rup[i]; if (ax[i] < L.rlo[i]) return ax[i] - L.rlo[i]; return 0.0; };
        if (rowunscale.empty()) rowunscale.assign(m, 1.0);
        // Violated-row queues (BISECT_QUEUE=0 falls back to random sampling): FIFO + in-queue flags,
        // fed by the incremental residual updates -- only rows of changed columns are re-checked.
        const bool useq = !(std::getenv("BISECT_QUEUE") && std::string(std::getenv("BISECT_QUEUE")) == "0");
        std::vector<int> qP, qD; size_t hP = 0, hD = 0; std::vector<char> inqP(m, 0), inqD(m, 0);
        double thr = 0;   // violation threshold, in original units (set every round from the tolerance)
        auto pushP = [&](int i) { if (!inqP[i] && std::fabs(prow(i)) * rowunscale[i] > thr) { inqP[i] = 1; qP.push_back(i); } };
        auto popq = [](std::vector<int>& q, size_t& h, std::vector<char>& f) { int i = q[h++]; f[i] = 0; if (h > 4096 && 2 * h > q.size()) { q.erase(q.begin(), q.begin() + h); h = 0; } return i; };
        auto set_x = [&](int j, double nv) { double d = nv - x[j]; if (d == 0) return; x[j] = nv; cx += L.c[j] * d;
            for (int p = L.start[j]; p < L.start[j + 1]; p++) { ax[L.index[p]] += L.value[p] * d; if (useq) pushP(L.index[p]); } };
        double bscale = 1; for (int i = 0; i < m; i++) if (std::isfinite(L.rup[i])) bscale = std::max(bscale, std::fabs(L.rup[i] * rowunscale[i]));
        auto pviol = [&]() { double v = 0; for (int i = 0; i < m; i++) v = std::max(v, std::fabs(prow(i)) * rowunscale[i]); return v / bscale; };
        // ---------------- dual state
        std::vector<double> y(m, 0.0), zl(n, 0.0), zu(nu, 0.0), aty(n, 0.0);
        double dobj = 0;   // rup'y + clo'zl - cup'zu  (explicit-z dual objective)
        auto set_y = [&](int i, double nv) { if (!std::isfinite(L.rlo[i])) nv = std::min(0.0, nv); double d = nv - y[i]; if (d == 0) return;
            y[i] = nv; dobj += L.rup[i] * d; for (int q = rstart[i]; q < rstart[i + 1]; q++) aty[rcol[q]] += rval[q] * d; };
        auto set_zl = [&](int j, double nv) { nv = std::max(0.0, nv); dobj += L.clo[j] * (nv - zl[j]); zl[j] = nv; };
        auto set_zu = [&](int j, double nv) { nv = std::max(0.0, nv); double& u = zu[uidx[j]]; dobj -= L.cup[j] * (nv - u); u = nv; };
        auto drow = [&](int j) { return L.c[j] - aty[j] - zl[j] + (uidx[j] >= 0 ? zu[uidx[j]] : 0.0); };
        double gd2 = 0;   // squared norm of the dual objective row
        for (int i = 0; i < m; i++) gd2 += L.rup[i] * L.rup[i];
        for (int j = 0; j < n; j++) { gd2 += L.clo[j] * L.clo[j]; if (uidx[j] >= 0) gd2 += L.cup[j] * L.cup[j]; }
        auto lagrange = [&]() {   // valid lower bound from y alone (optimization 3)
            double v = 0;
            for (int i = 0; i < m; i++) v += L.rup[i] * y[i];
            for (int j = 0; j < n; j++) { double rc = L.c[j] - aty[j];
                if (rc > 0) v += rc * L.clo[j]; else if (rc < 0) { if (!std::isfinite(ubnd[j])) return -INF; v += rc * ubnd[j]; } }
            return v; };
        // ---------------- bounds and level
        double Lb = lagrange(), Ub = INF, theta = 0.5;
        double t = Lb;
        int levels = 0;
        // Split mode (BISECT_ALPHA / BISECT_BETA, value-terms of a maximization): the primal is asked
        // to improve its own bound by alpha*G and the dual by (1-beta)*G, G = U - L, recomputed every
        // round.  In min form: tp = U - alpha G (primal slice c'x <= tp), td = L + (1-beta) G (dual
        // slice LB >= td).  beta > alpha makes the two slices overlap, so one is always nonempty.
        bool split = std::getenv("BISECT_ALPHA") != nullptr;
        double alpha = split ? atof(std::getenv("BISECT_ALPHA")) : 0.5;
        double betas = std::getenv("BISECT_BETA") ? atof(std::getenv("BISECT_BETA")) : 1.0 - alpha;
        double tp = INF, td = INF;
        std::vector<double> xbest = x;
        std::uniform_int_distribution<long> pick_p(0, std::max(0, m - 1)), pick_d(0, n - 1);
        long chunk = std::max(1000L, (long)(m + n));
        int gap_every = std::max(1, (m + n) / 16);
        const double omega = 1.5;
        long it = 0; int pfail = 0, dfail = 0;
        std::vector<double> gi(m);
        std::string dm = std::getenv("BISECT_DUAL") ? std::getenv("BISECT_DUAL") : "prox";
        int dual_mode = dm == "skm" ? 0 : dm == "polyak" ? 1 : dm == "hildreth" ? 2 : dm == "prox" ? 4 : 3;   // 3 = "both"
        bool dual_polyak = dual_mode == 1;
        // smoothing parameter for the Hildreth dual: start where x(y) spans the box for reduced
        // costs of the size of c, halve whenever the true bound stops improving
        double cmaxa = 0, ubmax = 1;
        for (int j = 0; j < n; j++) { cmaxa = std::max(cmaxa, std::fabs(L.c[j])); if (std::isfinite(ubnd[j])) ubmax = std::max(ubmax, ubnd[j]); }
        double mu = std::max(1e-12, cmaxa / ubmax), lb_ref = -INF; int mu_stall = 0;
        const double mu_floor = mu * 1e-4, mu_reset = mu * 1e-2; double lb_true_ref = -INF;
        double mu_q = -1; bool dual_settled = false;
        // proximal dual (method of multipliers): x_d minimizes the augmented Lagrangian
        //   c'x - y'(Ax - b) + sigma/2 |Ax - b|^2   over the box  (inequality rows: max(0, .) form)
        // by column coordinate descent; then y <- y - sigma (Ax_d - b)  (y <= 0 kept on "<=" rows).
        // This is the proximal point method on the dual: y+ = argmax LB(y) - |y - y_k|^2 / (2 sigma).
        std::vector<double> xd, rd, cn; double sigma = 0, rprev = INF; long cd_steps = 0, outer_d = 0; int r_stall = 0; double lb_prev_d = -INF, sigma0 = 0;
        if (dual_mode == 4) {
            xd.assign(n, 0.0); for (int j = 0; j < n; j++) xd[j] = std::clamp(0.0, L.clo[j], L.cup[j]);
            rd.assign(m, 0.0); for (int i = 0; i < m; i++) rd[i] = -L.rup[i];
            for (int j = 0; j < n; j++) for (int p = L.start[j]; p < L.start[j + 1]; p++) rd[L.index[p]] += L.value[p] * xd[j];
            cn.assign(n, 0.0); for (int j = 0; j < n; j++) for (int p = L.start[j]; p < L.start[j + 1]; p++) cn[j] += L.value[p] * L.value[p];
            // sigma ~ |c| / (size of the box): the penalty and the objective are of the same order
            sigma = std::getenv("BISECT_SIGMA") ? atof(std::getenv("BISECT_SIGMA")) : std::max(1e-12, cmaxa / ubmax);
            sigma0 = sigma;
        }
        BisectResult R; R.status = "time_limit";
        double tol_scale = 1.0; int crossings = 0;
        auto relgap = [&]() { return std::isfinite(Ub) ? (Ub - Lb) / (1 + std::fabs(Ub) + std::fabs(Lb)) : INF; };
        auto newlevel = [&]() {
            if (!std::isfinite(Ub)) t = INF;                      // no primal point yet: pure feasibility
            else t = Lb + theta * (Ub - Lb);
            tp = td = t;
            levels++;
        };
        auto settargets = [&]() {                                  // split mode, every round
            double otp = tp, otd = td;
            if (!std::isfinite(Ub)) { tp = INF; td = Lb + std::max(1.0, std::fabs(Lb)); }
            else { double G = Ub - Lb; tp = Ub - alpha * G; td = Lb + (1.0 - betas) * G; }
            if (tp != otp || td != otd) levels++;
        };
        if (split) settargets(); else newlevel();
        trace("bisect", el(), 0, Ub, pviol(), NAN, relgap(), Lb);
        while (el() < tlimit) {
            double tolp = tol_scale * (std::isfinite(Ub) ? std::clamp(0.01 * relgap(), tol, 1e-3) : 1e-3);   // optimization 6
            // ---------- primal chunk on P(t)
            thr = 0.5 * tolp * bscale;
            if (useq) {
                for (int i = 0; i < m; i++) pushP(i);          // threshold may have tightened
                for (long s = 0; s < chunk && m > 0; s++, it++) {
                    bool objv = std::isfinite(tp) && cx > tp && cn2 > 0;
                    if (hP >= qP.size()) { if (!objv) break; }   // slice solved to tolerance
                    else {
                        int i = popq(qP, hP, inqP);
                        if (std::fabs(prow(i)) * rowunscale[i] > thr) { double sc = prow(i) / pn[i];
                            for (int q = rstart[i]; q < rstart[i + 1]; q++) { int j = rcol[q]; set_x(j, std::clamp(x[j] - sc * rval[q], L.clo[j], L.cup[j])); } }
                        if (s % gap_every != 0) continue;
                    }
                    if (objv) {   // objective row c'x <= tp, over-relaxed; pushes the rows it breaks
                        double sc = omega * (cx - tp) / cn2;
                        for (int j = 0; j < n; j++) if (L.c[j] != 0) set_x(j, std::clamp(x[j] - sc * L.c[j], L.clo[j], L.cup[j]));
                    }
                }
            } else for (long s = 0; s < chunk && m > 0; s++, it++) {
                long best = -1; double bd = 0;
                for (int k = 0; k < beta; k++) { long r = pick_p(rng); double v = prow((int)r); double d = pn[r] > 0 ? v * v / pn[r] : 0; if (d > bd) { bd = d; best = r; } }
                if (best >= 0) { int i = (int)best; double sc = prow(i) / pn[i];
                    for (int q = rstart[i]; q < rstart[i + 1]; q++) { int j = rcol[q]; set_x(j, std::clamp(x[j] - sc * rval[q], L.clo[j], L.cup[j])); } }
                if (s % gap_every == 0 && std::isfinite(tp) && cx > tp && cn2 > 0) {   // objective row c'x <= tp, over-relaxed
                    double sc = omega * (cx - tp) / cn2;
                    for (int j = 0; j < n; j++) if (L.c[j] != 0) set_x(j, std::clamp(x[j] - sc * L.c[j], L.clo[j], L.cup[j]));
                }
            }
            double pv = pviol();
            bool psucc = pv <= tolp && cx <= (std::isfinite(tp) ? tp : INF) + 1e-12 * (1 + std::fabs(tp));
            if (psucc && cx < Ub) { Ub = cx; xbest = x; }
            // ---------- dual chunk
            if (dual_mode == 4) {
                auto wof = [&](int i) { double w = -y[i] + sigma * rd[i]; return std::isfinite(L.rlo[i]) ? w : std::max(0.0, w); };
                for (long s = 0; s < chunk; s++, it++, cd_steps++) {
                    int j = (int)pick_d(rng);
                    if (cn[j] <= 0) continue;
                    double g = L.c[j];
                    for (int p = L.start[j]; p < L.start[j + 1]; p++) g += L.value[p] * wof(L.index[p]);
                    double nx = std::clamp(xd[j] - g / (sigma * cn[j]), L.clo[j], ubnd[j]), d = nx - xd[j];   // implied bounds keep x_d bounded
                    if (d == 0) continue;
                    xd[j] = nx;
                    for (int p = L.start[j]; p < L.start[j + 1]; p++) rd[L.index[p]] += L.value[p] * d;
                    // multiplier update after every `sweeps` coordinate sweeps (inexact proximal step)
                }
                if (cd_steps >= 5L * n) {
                    cd_steps = 0; outer_d++;
                    double rmax = 0;
                    for (int i = 0; i < m; i++) {
                        double ny = y[i] - sigma * rd[i];
                        set_y(i, ny);   // set_y clamps y <= 0 on "<=" rows and maintains A'y
                        double v = std::isfinite(L.rlo[i]) ? std::fabs(rd[i]) : std::max(0.0, rd[i]);
                        rmax = std::max(rmax, v * rowunscale[i]);
                    }
                    if (std::getenv("BISECT_DEBUG")) { double al = 0; for (int j = 0; j < n; j++) al += L.c[j] * xd[j];
                        std::fprintf(stderr, "outer %ld t=%.2f sigma=%.3g rmax=%.3g cx_d=%.6g LB=%.6g Lb=%.6g\n", outer_d, el(), sigma, rmax, al, lagrange(), Lb); }
                    // sigma in both directions: an update that lowers the Lagrangian bound overshot
                    // (halve sigma); a stalled primal residual asks for more penalty (double it)
                    double lbn = lagrange();
                    if (std::isfinite(lb_prev_d) && lbn < lb_prev_d - 1e-3 * (1 + std::fabs(lb_prev_d))) { sigma = std::max(sigma0 * 1e-3, 0.5 * sigma); r_stall = 0; }
                    else if (rmax < 0.9 * rprev) { rprev = rmax; r_stall = 0; }
                    else if (++r_stall >= 10) { sigma *= 2.0; r_stall = 0; rprev = rmax; }
                    lb_prev_d = lbn;
                    // x_d is itself a nearly feasible primal point: use it for the primal bound too
                    double cxd = 0; for (int j = 0; j < n; j++) cxd += L.c[j] * xd[j];
                    if (rmax <= tolp * bscale && cxd < Ub) { Ub = cxd; xbest = xd; }
                }
            }
            if (dual_mode == 2 || dual_mode == 3) {
                // smoothed coordinate ascent (Hildreth row action on the dual):
                //   LB_mu(y) = rup'y + sum_j min_{x_j in box} [ (c - A'y)_j x_j + mu/2 x_j^2 ]
                //   x_j(y)   = clamp((A'y - c)_j / mu, clo_j, ubnd_j)
                //   row step y_i += mu (rup_i - a_i'x(y)) / |a_i|^2   (y_i <= 0 on "<=" rows)
                // The bound reported is the unsmoothed LB(y), valid for every y.
                auto xof = [&](int j) { return std::clamp((aty[j] - L.c[j]) / mu, L.clo[j], ubnd[j]); };
                auto grad = [&](int i) {    // smoothed dual gradient, zero where the sign constraint blocks
                    double r = 0; for (int q = rstart[i]; q < rstart[i + 1]; q++) r += rval[q] * xof(rcol[q]);
                    double g = L.rup[i] - r;
                    if (!std::isfinite(L.rlo[i]) && y[i] >= 0 && g > 0) g = 0;
                    return g; };
                auto pushD = [&](int k) { if (!inqD[k] && pn[k] > 0 && std::fabs(grad(k)) * rowunscale[k] > thr) { inqD[k] = 1; qD.push_back(k); } };
                if (useq) {
                    if (mu != mu_q) { for (int i = 0; i < m; i++) pushD(i); mu_q = mu; }   // all gradients changed
                    long s = 0;
                    for (; s < chunk && hD < qD.size(); s++, it++) {
                        int i = popq(qD, hD, inqD);
                        double g = grad(i);
                        if (std::fabs(g) * rowunscale[i] <= thr) continue;
                        set_y(i, y[i] + mu * g / pn[i]);
                        for (int q = rstart[i]; q < rstart[i + 1]; q++) { int j = rcol[q];        // rows sharing a column
                            for (int p = L.start[j]; p < L.start[j + 1]; p++) pushD(L.index[p]); }
                    }
                    if (hD >= qD.size()) dual_settled = true;    // stationary for this mu (to tolerance)
                } else for (long s = 0; s < chunk && m > 0; s++, it++) {
                    int i = (int)pick_p(rng);
                    if (pn[i] <= 0) continue;
                    double r = 0;
                    for (int q = rstart[i]; q < rstart[i + 1]; q++) r += rval[q] * xof(rcol[q]);
                    set_y(i, y[i] + mu * (L.rup[i] - r) / pn[i]);
                }
            }
            if (dual_polyak || dual_mode == 3) {
                // projection onto the supporting halfspace of {y : LB(y) >= target} (Polyak step):
                //   x(y) = box minimizer of the Lagrangian, g = rup - A x(y), y += w (target - LB) g / |g|^2
                double target = std::isfinite(td) ? td : Lb + std::max(1.0, std::fabs(Lb));
                long steps = std::max(1L, (long)(4.0 * chunk / std::max<size_t>(1, L.nnz())));
                for (long s = 0; s < steps; s++, it++) {
                    double lbv = 0; bool fin = true;
                    for (int i = 0; i < m; i++) { gi[i] = L.rup[i]; lbv += L.rup[i] * y[i]; }
                    for (int j = 0; j < n; j++) { double rc = L.c[j] - aty[j], xj;
                        if (rc > 0) xj = L.clo[j]; else if (rc < 0) xj = ubnd[j]; else xj = L.clo[j];
                        if (!std::isfinite(xj)) { fin = false; break; }
                        lbv += rc * xj;
                        if (xj != 0) for (int p = L.start[j]; p < L.start[j + 1]; p++) gi[L.index[p]] -= L.value[p] * xj; }
                    if (!fin || lbv >= target) break;
                    double g2 = 0; for (int i = 0; i < m; i++) { if (!std::isfinite(L.rlo[i]) && y[i] >= 0 && gi[i] > 0) gi[i] = 0; g2 += gi[i] * gi[i]; }
                    if (g2 <= 0) break;
                    double sc = omega * (target - lbv) / g2;
                    for (int i = 0; i < m; i++) if (gi[i] != 0) set_y(i, y[i] + sc * gi[i]);
                }
            } else if (dual_mode == 0) for (long s = 0; s < chunk; s++, it++) {
                long best = -1; double bd = 0;
                for (int k = 0; k < beta; k++) { long r = pick_d(rng); double v = drow((int)r); double d = v * v / dn[r]; if (d > bd) { bd = d; best = r; } }
                if (best >= 0) { int j = (int)best; double sc = drow(j) / dn[j];
                    for (int p = L.start[j]; p < L.start[j + 1]; p++) set_y(L.index[p], y[L.index[p]] + sc * L.value[p]);
                    set_zl(j, zl[j] + sc);
                    if (uidx[j] >= 0) set_zu(j, zu[uidx[j]] - sc); }
                double target = std::isfinite(td) ? td : Lb + std::max(1.0, std::fabs(Lb));
                if (s % gap_every == 0 && dobj < target && gd2 > 0) {   // objective row dobj >= target, over-relaxed
                    double sc = omega * (target - dobj) / gd2;
                    for (int i = 0; i < m; i++) if (L.rup[i] != 0) set_y(i, y[i] + sc * L.rup[i]);
                    for (int j = 0; j < n; j++) { if (L.clo[j] != 0) set_zl(j, zl[j] + sc * L.clo[j]); if (uidx[j] >= 0) set_zu(j, zu[uidx[j]] - sc * L.cup[j]); }
                }
            }
            double lb = lagrange();
            if (dual_mode == 2 || dual_mode == 3) {     // continuation on mu, driven by the smoothed objective
                double lbs = 0;
                for (int i = 0; i < m; i++) lbs += L.rup[i] * y[i];
                for (int j = 0; j < n; j++) { double rc = L.c[j] - aty[j], xj = std::clamp(-rc / mu, L.clo[j], ubnd[j]); lbs += rc * xj + 0.5 * mu * xj * xj; }
                if (std::getenv("BISECT_DEBUG")) std::fprintf(stderr, "t=%.2f mu=%.3g LBmu=%.6g LB=%.6g Lb=%.6g ymax=%.3g\n", el(), mu, lbs, lb, Lb,
                    *std::max_element(y.begin(), y.end(), [](double a, double b) { return std::fabs(a) < std::fabs(b); }));
                if (!dual_settled && (!std::isfinite(lb_ref) || lbs > lb_ref + 1e-4 * (1 + std::fabs(lb_ref)))) { lb_ref = lbs; mu_stall = 0; }
                else if (dual_settled || ++mu_stall >= 2) {
                    // shrink mu; at the floor, the true bound is stuck at a kink: re-smooth to escape
                    if (mu > mu_floor) mu = std::max(mu_floor, 0.5 * mu);
                    else if (lb <= lb_true_ref + 1e-7 * (1 + std::fabs(lb_true_ref))) mu = mu_reset;
                    lb_true_ref = std::max(lb_true_ref, lb);
                    mu_stall = 0; lb_ref = -INF;
                }
                dual_settled = false;
            }
            bool dsucc = std::isfinite(td) && std::isfinite(Ub) && lb >= td - 1e-12 * (1 + std::fabs(td));
            if (lb > Lb) Lb = lb;
            if (std::isfinite(Ub) && Lb > Ub + 1e-12 * (1 + std::fabs(Ub))) {
                // the valid dual bound crossed the primal "bound": that primal point was too infeasible.
                // Discard it and demand tighter primal feasibility from now on.
                tol_scale *= 0.1; Ub = INF; crossings++;
            }
            // ---------- update the level (optimizations 2 and 4)
            if (split) settargets();
            else if (psucc || dsucc || (!std::isfinite(t) && std::isfinite(Ub))) {
                if (psucc) { pfail = 0; dfail++; } if (dsucc) { dfail = 0; pfail++; }
                newlevel();
            }
            // theta: if the primal keeps succeeding, t is too easy for it -> move t down (harder
            // primal, easier dual), and vice versa
            if (!split && psucc) theta = std::max(0.1, theta * 0.8);
            if (!split && dsucc) theta = std::min(0.9, 1 - (1 - theta) * 0.8);
            double g = relgap();
            trace("bisect", el(), it, Ub, pv, NAN, g, Lb);
            if (g <= tol) { R.status = "optimal"; break; }
        }
        R.U = Ub; R.L = Lb; R.rp = pviol(); R.iters = it; R.levels = levels; R.x = xbest;
        return R;
    }
};
