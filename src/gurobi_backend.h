// gurobi_backend.h -- Gurobi (commercial; needs a license) through its C API, built into a
// separate binary, bin/lpsolve-grb, so that libgurobi does not change the memory baseline of
// bin/lpsolve. Modes:
//   gurobi-primal         primal simplex (Method=0)
//   gurobi-dual           dual simplex (Method=1)
//   gurobi-barrier        barrier without crossover (Method=2, Crossover=0), as clp-barrier
//   gurobi-barrier-cross  barrier with Gurobi's default crossover, as highs-ipx-cross
//   gurobi-pdhg           Gurobi's first-order PDHG (Method=6) without crossover, on the CPU
// Common settings: Threads=1, TimeLimit, presolve at its default (automatic). Tolerances:
// (Set LPSOLVE_GRB_LOG=1 to see Gurobi's log on stdout.)
// FeasibilityTol/OptimalityTol (simplex) and PDHGRelTol/PDHGAbsTol default to 1e-6 already and
// are set to tol; BarConvTol (default 1e-8) is set to tol, as HiGHS IPX's optimality tolerance.
#pragma once
#include "gurobi_c.h"

static Out run_gurobi(LP& L, const std::string& mode, double tl, double tol, long& hwm_load) {
    Out o;
    GRBenv* env = nullptr;
    if (GRBemptyenv(&env) || GRBsetintparam(env, "OutputFlag", getenv("LPSOLVE_GRB_LOG") ? 1 : 0) || GRBstartenv(env)) {
        o.status = std::string("fail(env: ") + (env ? GRBgeterrormsg(env) : "no env") + ")";
        if (env) GRBfreeenv(env);
        return o;
    }
    // rows are equalities (rlo == rup) or <= rup (rlo = -inf), as in run_glpk
    std::vector<char> sense(L.nrow);
    for (int i = 0; i < L.nrow; i++) sense[i] = L.rlo[i] == L.rup[i] ? GRB_EQUAL : GRB_LESS_EQUAL;
    std::vector<int> cnt(L.ncol);
    for (int j = 0; j < L.ncol; j++) cnt[j] = L.start[j + 1] - L.start[j];
    for (auto& v : L.cup) if (!std::isfinite(v)) v = GRB_INFINITY;
    GRBmodel* m = nullptr;
    int err = GRBloadmodel(env, &m, "lp", L.ncol, L.nrow, GRB_MINIMIZE, 0.0, L.c.data(), sense.data(),
                           L.rup.data(), L.start.data(), cnt.data(), L.index.data(), L.value.data(),
                           L.clo.data(), L.cup.data(), nullptr, nullptr, nullptr);
    if (!err) err = GRBupdatemodel(m);   // force Gurobi to build its copy before ours is freed
    release(L); std::vector<char>().swap(sense); std::vector<int>().swap(cnt);
    if (err) { o.status = std::string("fail(load: ") + GRBgeterrormsg(env) + ")"; GRBfreeenv(env); return o; }
    hwm_load = vm_hwm_kb();
    GRBenv* me = GRBgetenv(m);
    GRBsetintparam(me, "Threads", 1);
    GRBsetdblparam(me, "TimeLimit", tl);
    GRBsetdblparam(me, "FeasibilityTol", std::max(tol, 1e-9));
    GRBsetdblparam(me, "OptimalityTol", std::max(tol, 1e-9));
    if (mode == "gurobi-primal") GRBsetintparam(me, "Method", 0);
    else if (mode == "gurobi-dual") GRBsetintparam(me, "Method", 1);
    else if (mode == "gurobi-barrier" || mode == "gurobi-barrier-cross") {
        GRBsetintparam(me, "Method", 2);
        GRBsetdblparam(me, "BarConvTol", tol);
        if (mode == "gurobi-barrier") GRBsetintparam(me, "Crossover", 0);
    } else if (mode == "gurobi-pdhg") {
        GRBsetintparam(me, "Method", 6);
        GRBsetintparam(me, "Crossover", 0);
        GRBsetintparam(me, "PDHGGPU", 0);
        GRBsetdblparam(me, "PDHGRelTol", tol);
        GRBsetdblparam(me, "PDHGAbsTol", std::max(tol, 1e-9));
    } else { o.status = "fail(unknown mode)"; GRBfreemodel(m); GRBfreeenv(env); return o; }
    auto t1 = Clock::now();
    err = GRBoptimize(m);
    o.t_solve = secs(t1);
    int st = 0; GRBgetintattr(m, GRB_INT_ATTR_STATUS, &st);
    if (err) o.status = std::string("fail(") + GRBgeterrormsg(env) + ")";
    else if (st == GRB_OPTIMAL) o.status = "optimal";
    else if (st == GRB_TIME_LIMIT) o.status = "time_limit";
    else o.status = "fail(status=" + std::to_string(st) + ")";
    double v;
    if (!GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, &v)) o.obj = v;
    const char* itattr = mode == "gurobi-pdhg" ? "PDHGIterCount" : (mode.rfind("gurobi-barrier", 0) == 0 ? "BarIterCount" : "IterCount");
    double it; int iti;
    if (!GRBgetdblattr(m, itattr, &it)) o.iters = (long)it;
    else if (!GRBgetintattr(m, itattr, &iti)) o.iters = iti;
    if (!GRBgetdblattr(m, "MaxVio", &v)) o.rp = v;
    GRBfreemodel(m); GRBfreeenv(env);
    return o;
}
