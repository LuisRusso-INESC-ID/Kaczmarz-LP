// lpsolve.cpp -- solve a graph instance through its LP formulation with one backend.
//
// usage: lpsolve <backend> <instance> [time_limit=600] [tol=1e-6]
//   backends: build-only | highs-dual | highs-primal | highs-ipx | highs-ipx-cross |
//             highs-pdlp | glpk-primal | glpk-dual | glpk-ipt | clp-dual | clp-primal |
//             clp-barrier | own-pdhg | own-skm | (with -DHAVE_ORTOOLS) ortools-glop | ortools-pdlp
//             (with -DHAVE_GUROBI) gurobi-primal | gurobi-dual | gurobi-barrier[-cross] | gurobi-pdhg
//
// The process prints one "RESULT ..." line; peak memory is measured by the driver (wait4
// rusage), and also reported here as VmHWM at three checkpoints: after reading the
// instance, after building our CSC model, and after handing the model to the solver (our
// copy is freed at that point so the solver's own footprint is what remains).
#include "lp.h"
#include "pdhg.h"
#include "skm.h"
#include "prowact.h"
#include "hybrid.h"
#include "bisect.h"
#include <chrono>
#include <fstream>
#include <sstream>
#ifdef HAVE_HIGHS
#include "Highs.h"
#endif
#ifdef HAVE_GLPK
#include <glpk.h>
#endif
#ifdef HAVE_CLP
#include "ClpSimplex.hpp"
#include "ClpSolve.hpp"
#include "CoinMessageHandler.hpp"
#endif

using Clock = std::chrono::steady_clock;
static double secs(Clock::time_point a) { return std::chrono::duration<double>(Clock::now() - a).count(); }

static long vm_hwm_kb() {
    std::ifstream f("/proc/self/status"); std::string l;
    while (std::getline(f, l)) if (l.rfind("VmHWM:", 0) == 0) return std::stol(l.substr(6));
    return -1;
}

struct Out { std::string status = "unknown"; double obj = NAN, rp = NAN; long iters = -1; double t_solve = 0; std::string extra; };

template <class T> static void release(std::vector<T>& v) { std::vector<T>().swap(v); }
static void release(LP& L) { release(L.c); release(L.clo); release(L.cup); release(L.rlo); release(L.rup); release(L.start); release(L.index); release(L.value); }
#ifdef HAVE_ORTOOLS
#include "ortools_backend.h"
#endif
#ifdef HAVE_GUROBI
#include "gurobi_backend.h"
#endif

#ifdef HAVE_HIGHS
static Out run_highs(LP& L, std::string mode, double tl, double tol, long& hwm_load) {
    Highs h;
    if (mode.size() > 6 && mode.compare(mode.size() - 6, 6, "-nopre") == 0) { h.setOptionValue("presolve", "off"); mode.resize(mode.size() - 6); }
    h.setOptionValue("output_flag", false);
    h.setOptionValue("time_limit", tl);
    h.setOptionValue("threads", 1);
    h.setOptionValue("parallel", "off");
    {
        HighsLp lp;
        lp.num_col_ = L.ncol; lp.num_row_ = L.nrow;
        lp.col_cost_ = std::move(L.c); lp.col_lower_ = std::move(L.clo); lp.col_upper_ = std::move(L.cup);
        lp.row_lower_ = std::move(L.rlo); lp.row_upper_ = std::move(L.rup);
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = L.ncol; lp.a_matrix_.num_row_ = L.nrow;
        lp.a_matrix_.start_ = std::move(L.start); lp.a_matrix_.index_ = std::move(L.index); lp.a_matrix_.value_ = std::move(L.value);
        h.passModel(std::move(lp));
    }
    release(L);
    hwm_load = vm_hwm_kb();
    if (mode == "highs-dual" || mode == "highs-primal") {
        h.setOptionValue("solver", "simplex");
        h.setOptionValue("simplex_strategy", mode == "highs-dual" ? 1 : 4);
        h.setOptionValue("primal_feasibility_tolerance", tol);
        h.setOptionValue("dual_feasibility_tolerance", tol);
    } else if (mode == "highs-ipx" || mode == "highs-ipx-cross") {
        h.setOptionValue("solver", "ipx");
        h.setOptionValue("run_crossover", mode == "highs-ipx" ? "off" : "on");
        h.setOptionValue("ipm_optimality_tolerance", tol);
    } else if (mode == "highs-pdlp") {
        h.setOptionValue("solver", "pdlp");
        h.setOptionValue("pdlp_optimality_tolerance", tol);
        h.setOptionValue("pdlp_iteration_limit", 2000000000);
    }
    auto t1 = Clock::now();
    h.run();
    Out o; o.t_solve = secs(t1);
    auto st = h.getModelStatus();
    o.status = st == HighsModelStatus::kOptimal ? "optimal" : st == HighsModelStatus::kTimeLimit ? "time_limit" : h.modelStatusToString(st);
    o.obj = h.getInfo().objective_function_value;
    const auto& info = h.getInfo();
    o.iters = info.simplex_iteration_count + info.ipm_iteration_count + info.pdlp_iteration_count;
    o.rp = info.max_primal_infeasibility;
    return o;
}
#endif

#ifdef HAVE_GLPK
static Out run_glpk(LP& L, const std::string& mode, double tl, double tol, long& hwm_load) {
    glp_term_out(GLP_OFF);
    glp_prob* P = glp_create_prob();
    glp_set_obj_dir(P, GLP_MIN);
    glp_add_rows(P, L.nrow); glp_add_cols(P, L.ncol);
    for (int i = 0; i < L.nrow; i++) {
        if (L.rlo[i] == L.rup[i]) glp_set_row_bnds(P, i + 1, GLP_FX, L.rlo[i], L.rup[i]);
        else glp_set_row_bnds(P, i + 1, GLP_UP, 0.0, L.rup[i]);
    }
    for (int j = 0; j < L.ncol; j++) {
        glp_set_obj_coef(P, j + 1, L.c[j]);
        if (std::isfinite(L.cup[j])) glp_set_col_bnds(P, j + 1, GLP_DB, L.clo[j], L.cup[j]);
        else glp_set_col_bnds(P, j + 1, GLP_LO, L.clo[j], 0.0);
        // column-wise load avoids building global triplet arrays
        int len = L.start[j + 1] - L.start[j];
        std::vector<int> ind(len + 1); std::vector<double> v(len + 1);
        for (int k = 0; k < len; k++) { ind[k + 1] = L.index[L.start[j] + k] + 1; v[k + 1] = L.value[L.start[j] + k]; }
        glp_set_mat_col(P, j + 1, len, ind.data(), v.data());
    }
    release(L);
    hwm_load = vm_hwm_kb();
    Out o; auto t1 = Clock::now();
    if (mode == "glpk-ipt") {
        glp_iptcp ip; glp_init_iptcp(&ip); ip.msg_lev = GLP_MSG_OFF;
        int r = glp_interior(P, &ip);
        o.t_solve = secs(t1);
        int st = glp_ipt_status(P);
        o.status = (r == 0 && st == GLP_OPT) ? "optimal" : "fail(" + std::to_string(r) + "," + std::to_string(st) + ")";
        o.obj = glp_ipt_obj_val(P);
    } else {
        glp_smcp sp; glp_init_smcp(&sp); sp.msg_lev = GLP_MSG_OFF;
        sp.meth = mode == "glpk-dual" ? GLP_DUALP : GLP_PRIMAL;
        sp.presolve = GLP_ON; sp.tm_lim = (int)std::min(2e9, tl * 1000);
        sp.tol_bnd = tol; sp.tol_dj = tol;
        int r = glp_simplex(P, &sp);
        o.t_solve = secs(t1);
        int st = glp_get_status(P);
        o.status = (r == 0 && st == GLP_OPT) ? "optimal" : r == GLP_ETMLIM ? "time_limit" : "fail(" + std::to_string(r) + "," + std::to_string(st) + ")";
        o.obj = glp_get_obj_val(P);
    }
    glp_delete_prob(P);
    return o;
}
#endif

#ifdef HAVE_CLP
static Out run_clp(LP& L, std::string mode, double tl, double tol, long& hwm_load) {
    ClpSimplex M;
    bool nopre = mode.size() > 6 && mode.compare(mode.size() - 6, 6, "-nopre") == 0;
    if (nopre) mode.resize(mode.size() - 6);
    M.setLogLevel(0);
    for (auto& v : L.cup) if (!std::isfinite(v)) v = COIN_DBL_MAX;
    for (auto& v : L.rlo) if (!std::isfinite(v)) v = -COIN_DBL_MAX;
    M.loadProblem(L.ncol, L.nrow, L.start.data(), L.index.data(), L.value.data(), L.clo.data(), L.cup.data(), L.c.data(), L.rlo.data(), L.rup.data());
    release(L);
    hwm_load = vm_hwm_kb();
    M.setMaximumSeconds(tl);
    M.setPrimalTolerance(tol); M.setDualTolerance(tol);
    ClpSolve opt;
    if (mode == "clp-dual") opt.setSolveType(ClpSolve::useDual);
    else if (mode == "clp-primal") opt.setSolveType(ClpSolve::usePrimal);
    else opt.setSolveType(ClpSolve::useBarrierNoCross);
    // clp-barrier-chol: Cholesky by CHOLMOD with AMD ordering (ClpCholeskyUfl; needs Clp built
    // with CHOLMOD) instead of Clp's native factorization with its own minimum-fill ordering
    if (mode == "clp-barrier-chol") opt.setSpecialOption(4, 4);
    opt.setPresolveType(nopre ? ClpSolve::presolveOff : ClpSolve::presolveOn);
    auto t1 = Clock::now();
    M.initialSolve(opt);
    Out o; o.t_solve = secs(t1);
    int st = M.status();
    o.status = st == 0 ? "optimal" : st == 3 ? "time_limit" : "fail(" + std::to_string(st) + "," + std::to_string(M.secondaryStatus()) + ")";
    o.obj = M.objectiveValue();
    o.iters = M.numberIterations();
    o.rp = M.sumPrimalInfeasibilities();
    return o;
}
#endif

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: lpsolve backend instance [tl] [tol]\n"); return 1; }
    std::string be = argv[1];
    double tl = argc > 3 ? atof(argv[3]) : 600, tol = argc > 4 ? atof(argv[4]) : 1e-6;
    auto t0 = Clock::now();
    LP L;
    long hwm_read, hwm_build, hwm_load = -1;
    {
        Instance I = read_instance(argv[2]);
        hwm_read = vm_hwm_kb();
        L = build_lp(I, 0);
    }   // instance freed
    hwm_build = vm_hwm_kb();
    double t_build = secs(t0);
    int nrow = L.nrow, ncol = L.ncol; size_t nnz = L.nnz(); double sense = L.sense;
    Out o;
    if (be == "build-only") { o.status = "optimal"; o.obj = 0; }
#ifdef HAVE_HIGHS
    else if (be.rfind("highs-", 0) == 0) o = run_highs(L, be, tl, tol, hwm_load);
#endif
#ifdef HAVE_GLPK
    else if (be.rfind("glpk-", 0) == 0) o = run_glpk(L, be, tl, tol, hwm_load);
#endif
#ifdef HAVE_CLP
    else if (be.rfind("clp-", 0) == 0) o = run_clp(L, be, tl, tol, hwm_load);
#endif
#ifdef HAVE_ORTOOLS
    else if (be.rfind("ortools-", 0) == 0) o = run_ortools(L, be, tl, tol, hwm_load);
#endif
#ifdef HAVE_GUROBI
    else if (be.rfind("gurobi-", 0) == 0) o = run_gurobi(L, be, tl, tol, hwm_load);
#endif
    else if (be == "own-pdhg") {
        Pdhg P(L); auto t1 = Clock::now();
        auto r = P.solve(tol, tl);
        o.t_solve = secs(t1); o.status = r.status; o.obj = r.pobj; o.iters = r.iters;
        o.rp = primal_residual(L, r.x);
        char b[160]; std::snprintf(b, sizeof b, " rel_rp=%.2e rel_rd=%.2e rel_gap=%.2e restarts=%d", r.rp, r.rd, r.gap, r.restarts); o.extra = b;
    } else if (be == "own-skm") {
        LP LS;   // optional Ruiz + Pock-Chambolle equilibrated copy (SKM_SCALE=1): same objective values
        bool scaled = std::getenv("SKM_SCALE") != nullptr;
        if (scaled) {
            Pdhg P(L); P.scale();
            LS = L; LS.value = P.val; LS.c = P.c; LS.clo = P.lo; LS.cup = P.up;
            for (int i = 0; i < L.nrow; i++) { LS.rup[i] = L.rup[i] * P.E[i]; LS.rlo[i] = std::isfinite(L.rlo[i]) ? L.rlo[i] * P.E[i] : L.rlo[i]; }
        }
        Skm S(scaled ? LS : L); auto t1 = Clock::now();
        auto r = S.solve(tol, tl);
        o.t_solve = secs(t1); o.status = r.status; o.obj = r.pobj; o.iters = r.iters;
        o.rp = primal_residual(L, r.x);
        char b[160]; std::snprintf(b, sizeof b, " rel_rp=%.2e rel_rd=%.2e rel_gap=%.2e", r.rp, r.rd, r.gap); o.extra = b;
    } else if (be == "own-hybrid-mem" || be == "own-hybrid-prox" || be == "own-hybrid-cg") {
        bool px = be == "own-hybrid-prox", cg = be == "own-hybrid-cg";
        Hybrid H(L, px, cg); auto t1 = Clock::now();
        int rm = std::getenv("HYB_RMAX") ? atoi(std::getenv("HYB_RMAX")) : 64;
        auto r = H.solve(tl, rm, 1, cg ? "hybrid-cg" : px ? "hybrid-prox" : "hybrid-mem");
        o.t_solve = secs(t1); o.status = r.status; o.obj = r.pobj; o.iters = r.blocks; o.rp = primal_residual(L, r.x);
        char b[160]; std::snprintf(b, sizeof b, " rel_viol=%.2e repairs=%ld rmax=%d", r.rp, r.repairs, r.rmax); o.extra = b;
    } else if (be == "own-bisect-kz") {
        // optimization 7: run on the equilibrated LP (same objective values); unscale x at the end
        Pdhg P(L); P.scale();
        LP LS = L; LS.value = P.val; LS.c = P.c; LS.clo = P.lo; LS.cup = P.up;
        for (int i = 0; i < L.nrow; i++) { LS.rup[i] = L.rup[i] * P.E[i]; LS.rlo[i] = std::isfinite(L.rlo[i]) ? L.rlo[i] * P.E[i] : L.rlo[i]; }
        if (!L.uhint.empty()) for (int j = 0; j < L.ncol; j++) LS.uhint[j] = L.uhint[j] / P.D[j];
        Bisect S(LS); S.rowunscale.resize(L.nrow); for (int i = 0; i < L.nrow; i++) S.rowunscale[i] = 1.0 / P.E[i];
        auto t1 = Clock::now();
        auto r = S.solve(tol, tl);
        for (int j = 0; j < L.ncol; j++) r.x[j] *= P.D[j];
        o.t_solve = secs(t1); o.status = r.status; o.obj = r.U; o.iters = r.iters; o.rp = primal_residual(L, r.x);
        char b[200]; std::snprintf(b, sizeof b, " lower_bound=%.12g rel_gap=%.2e levels=%d", sense * r.L, (r.U - r.L) / (1 + std::fabs(r.U) + std::fabs(r.L)), r.levels); o.extra = b;
    } else if (be == "own-prow") {
        Prow S(L); auto t1 = Clock::now();
        auto r = S.solve(tol, tl);
        o.t_solve = secs(t1); o.status = r.status; o.obj = r.pobj; o.iters = r.sweeps;
        o.rp = primal_residual(L, r.x);
        char b[160]; std::snprintf(b, sizeof b, " rel_rp=%.2e rel_gap=%.2e outer=%d", r.rp, r.gap, r.outer); o.extra = b;
    } else { std::fprintf(stderr, "unknown or unavailable backend %s\n", be.c_str()); return 2; }
    std::printf("RESULT status=%s obj=%.12g iters=%ld rp=%.3e nrow=%d ncol=%d nnz=%zu t_build=%.4f t_solve=%.4f t_total=%.4f hwm_read_kb=%ld hwm_build_kb=%ld hwm_load_kb=%ld hwm_end_kb=%ld%s\n",
                o.status.c_str(), sense * o.obj, o.iters, o.rp, nrow, ncol, nnz, t_build, o.t_solve, secs(t0),
                hwm_read, hwm_build, hwm_load, vm_hwm_kb(), o.extra.c_str());
    return 0;
}
