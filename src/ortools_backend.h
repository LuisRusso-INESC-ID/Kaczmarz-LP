// ortools_backend.h -- OR-Tools GLOP (primal/dual simplex) and PDLP via their direct C++ APIs
// (not MPSolver, to avoid its per-variable object overhead in the memory figures).
#pragma once
#include "ortools/glop/lp_solver.h"
#include "ortools/glop/parameters.pb.h"
#include "ortools/lp_data/lp_data.h"
#include "ortools/pdlp/primal_dual_hybrid_gradient.h"
#include "ortools/pdlp/quadratic_program.h"
#include "ortools/pdlp/solvers.pb.h"
#include "ortools/pdlp/solve_log.pb.h"

static Out run_ortools(LP& L, const std::string& mode, double tl, double tol, long& hwm_load) {
    namespace glop = operations_research::glop;
    namespace pdlp = operations_research::pdlp;
    Out o;
    if (mode == "ortools-glop" || mode == "ortools-glop-primal") {
        glop::LinearProgram lp;
        for (int i = 0; i < L.nrow; i++) { auto r = lp.CreateNewConstraint(); lp.SetConstraintBounds(r, L.rlo[i], L.rup[i]); }
        for (int j = 0; j < L.ncol; j++) {
            auto col = lp.CreateNewVariable();
            lp.SetVariableBounds(col, L.clo[j], L.cup[j]);
            lp.SetObjectiveCoefficient(col, L.c[j]);
            for (int p = L.start[j]; p < L.start[j + 1]; p++) lp.SetCoefficient(glop::RowIndex(L.index[p]), col, L.value[p]);
        }
        lp.CleanUp();
        release(L);
        hwm_load = vm_hwm_kb();
        glop::GlopParameters prm;
        prm.set_max_time_in_seconds(tl);
        prm.set_use_dual_simplex(mode == "ortools-glop");
        prm.set_primal_feasibility_tolerance(tol);
        prm.set_dual_feasibility_tolerance(tol);
        glop::LPSolver solver; solver.SetParameters(prm);
        auto t1 = Clock::now();
        auto st = solver.Solve(lp);
        o.t_solve = secs(t1);
        o.status = st == glop::ProblemStatus::OPTIMAL ? "optimal"
                 : o.t_solve >= 0.98 * tl ? "time_limit" : glop::GetProblemStatusString(st);
        o.obj = solver.GetObjectiveValue();
        o.iters = solver.GetNumberOfSimplexIterations();
        return o;
    }
    // PDLP
    pdlp::QuadraticProgram qp(L.ncol, L.nrow);
    for (int j = 0; j < L.ncol; j++) { qp.objective_vector[j] = L.c[j]; qp.variable_lower_bounds[j] = L.clo[j]; qp.variable_upper_bounds[j] = L.cup[j]; }
    for (int i = 0; i < L.nrow; i++) { qp.constraint_lower_bounds[i] = L.rlo[i]; qp.constraint_upper_bounds[i] = L.rup[i]; }
    {
        std::vector<int64_t> st(L.start.begin(), L.start.end()), ix(L.index.begin(), L.index.end());
        Eigen::Map<const Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>> A(L.nrow, L.ncol, (int64_t)L.nnz(), st.data(), ix.data(), L.value.data());
        qp.constraint_matrix = A;
    }
    release(L);
    hwm_load = vm_hwm_kb();
    pdlp::PrimalDualHybridGradientParams prm;
    prm.set_num_threads(1);
    auto* tc = prm.mutable_termination_criteria();
    tc->set_time_sec_limit(tl);
    tc->mutable_simple_optimality_criteria()->set_eps_optimal_relative(tol);
    tc->mutable_simple_optimality_criteria()->set_eps_optimal_absolute(tol);
    auto t1 = Clock::now();
    pdlp::SolverResult r = pdlp::PrimalDualHybridGradient(std::move(qp), prm);
    o.t_solve = secs(t1);
    auto tr = r.solve_log.termination_reason();
    o.status = tr == pdlp::TERMINATION_REASON_OPTIMAL ? "optimal" : tr == pdlp::TERMINATION_REASON_TIME_LIMIT ? "time_limit" : pdlp::TerminationReason_Name(tr);
    o.iters = r.solve_log.iteration_count();
    const auto& ss = r.solve_log.solution_stats();
    if (ss.convergence_information_size() > 0) {
        o.obj = ss.convergence_information(0).primal_objective();
        o.rp = ss.convergence_information(0).l_inf_primal_residual();
    }
    return o;
}
