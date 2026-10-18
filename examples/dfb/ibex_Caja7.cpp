/* Se compila con -Dsoplex=soplex7 contra los encabezados de SoPlex 7 (ver ibex_Caja7.h). */
#include "soplex.h"
#include "ibex_Caja7.h"
#include <vector>
#include <cstdlib>
using namespace soplex;
namespace ibex {
Caja7::Caja7() : sx(new SoPlex()) {}
Caja7::~Caja7() { delete (SoPlex*) sx; }
void Caja7::cargar(int n, int m, const double* lb, const double* ub, const double* lhs,
                   const double* rhs, const double* A, double tol, int iterlimit) {
    SoPlex& s = *(SoPlex*) sx;
    s.clearLPReal();
    s.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
    s.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
    s.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_REAL);
    s.setIntParam(SoPlex::TIMER, SoPlex::TIMER_OFF);
    s.setIntParam(SoPlex::ITERLIMIT, iterlimit);
    s.setRealParam(SoPlex::FEASTOL, tol);
    s.setRealParam(SoPlex::OPTTOL, tol);
    s.setBoolParam(SoPlex::ENSURERAY, true);
    /* CAJA7_SIMP=0 / CAJA7_SCALER=0: sin presolve / sin escalado (para buscar la mejor configuracion) */
    if (getenv("CAJA7_SIMP") && atoi(getenv("CAJA7_SIMP")) == 0) s.setIntParam(SoPlex::SIMPLIFIER, SoPlex::SIMPLIFIER_OFF);
    if (getenv("CAJA7_SCALER") && atoi(getenv("CAJA7_SCALER")) == 0) s.setIntParam(SoPlex::SCALER, SoPlex::SCALER_OFF);
    const double inf = s.realParam(SoPlex::INFTY);
    LPColSetReal cols(n); DSVectorReal vacio(0);
    for (int i = 0; i < n; ++i) cols.add(0.0, lb[i], vacio, ub[i]);
    s.addColsReal(cols);
    LPRowSetReal rows(m);
    for (int j = 0; j < m; ++j) {
        DSVectorReal r(n);
        const bool libre = lhs[j] <= -inf && rhs[j] >= inf;
        if (!libre) for (int i = 0; i < n; ++i) if (A[(long) j * n + i] != 0.0) r.add(i, A[(long) j * n + i]);
        rows.add(lhs[j], r, rhs[j]);
    }
    s.addRowsReal(rows);
}
void Caja7::cota(int i, double lb, double ub) { ((SoPlex*) sx)->changeBoundsReal(i, lb, ub); }
int Caja7::resolver(int k, bool maximizar, bool frio, double* lam, double* x, int* iters) {
    SoPlex& s = *(SoPlex*) sx;
    const int m = s.numRowsReal(), n = s.numColsReal();
    s.changeObjReal(k, maximizar ? -1.0 : 1.0);
    if (frio) s.clearBasis();
    s.optimize();
    s.ignoreUnscaledViolations();
    const SPxSolver::Status st = s.status();
    *iters = s.numIterations();
    int res = 3;
    if (st == SPxSolver::OPTIMAL) {
        if (s.getDualReal(lam, m) && s.getPrimalReal(x, n)) res = 0;
    } else if (st == SPxSolver::INFEASIBLE) {
        if (s.getDualFarkasReal(lam, m)) res = 1;
        else res = 2;
    }
    s.changeObjReal(k, 0.0);
    return res;
}
}
