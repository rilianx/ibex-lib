#include "Highs.h"
#include "ibex_CajaHighs.h"
#include <vector>
#include <cmath>
namespace ibex {
CajaHighs::CajaHighs() : h(new Highs()), n_(0), m_(0) {
    Highs& s = *(Highs*) h;
    s.setOptionValue("output_flag", false);
    s.setOptionValue("threads", 1);
}
CajaHighs::~CajaHighs() { delete (Highs*) h; }
static double inf_h(double v) { return v >= 1e100 ? kHighsInf : (v <= -1e100 ? -kHighsInf : v); }
void CajaHighs::cargar(int n, int m, const double* lb, const double* ub, const double* lhs,
                       const double* rhs, const double* A, double tol, bool presolve) {
    Highs& s = *(Highs*) h;
    s.clearModel();
    s.setOptionValue("primal_feasibility_tolerance", tol);
    s.setOptionValue("dual_feasibility_tolerance", tol);
    s.setOptionValue("presolve", presolve ? "choose" : "off");
    HighsLp lp;
    lp.num_col_ = n; lp.num_row_ = m;
    lp.col_cost_.assign(n, 0.0);
    lp.col_lower_.resize(n); lp.col_upper_.resize(n);
    for (int i = 0; i < n; ++i) { lp.col_lower_[i] = inf_h(lb[i]); lp.col_upper_[i] = inf_h(ub[i]); }
    lp.row_lower_.resize(m); lp.row_upper_.resize(m);
    lp.a_matrix_.format_ = MatrixFormat::kRowwise;
    lp.a_matrix_.num_col_ = n; lp.a_matrix_.num_row_ = m;
    lp.a_matrix_.start_.assign(1, 0);
    for (int j = 0; j < m; ++j) {
        const double lo = inf_h(lhs[j]), up = inf_h(rhs[j]);
        lp.row_lower_[j] = lo; lp.row_upper_[j] = up;
        const bool libre = (lo == -kHighsInf && up == kHighsInf);
        if (!libre)
            for (int i = 0; i < n; ++i) {
                const double a = A[(long) j * n + i];
                if (a != 0.0) { lp.a_matrix_.index_.push_back(i); lp.a_matrix_.value_.push_back(a); }
            }
        lp.a_matrix_.start_.push_back((HighsInt) lp.a_matrix_.index_.size());
    }
    lp.sense_ = ObjSense::kMinimize;
    s.passModel(lp);
    n_ = n; m_ = m;
}
void CajaHighs::cota(int i, double lb, double ub) { ((Highs*) h)->changeColBounds(i, inf_h(lb), inf_h(ub)); }
int CajaHighs::resolver(int k, bool maximizar, bool frio, double* lam, double* x, int* iters) {
    Highs& s = *(Highs*) h;
    s.changeColCost(k, maximizar ? -1.0 : 1.0);
    if (frio) s.clearSolver();
    s.run();
    const HighsModelStatus st = s.getModelStatus();
    *iters = s.getInfo().simplex_iteration_count;
    int res = 3;
    if (st == HighsModelStatus::kOptimal) {
        const HighsSolution& sol = s.getSolution();
        if ((int) sol.row_dual.size() == m_ && (int) sol.col_value.size() == n_) {
            for (int j = 0; j < m_; ++j) lam[j] = sol.row_dual[j];
            for (int i = 0; i < n_; ++i) x[i] = sol.col_value[i];
            res = 0;
        }
    } else if (st == HighsModelStatus::kInfeasible) {
        bool tiene = false;
        std::vector<double> ray(m_, 0.0);
        if (s.getDualRay(tiene, ray.data()) == HighsStatus::kOk && tiene) {
            for (int j = 0; j < m_; ++j) lam[j] = ray[j];
            res = 1;
        } else res = 2;
    }
    s.changeColCost(k, 0.0);
    return res;
}
}
