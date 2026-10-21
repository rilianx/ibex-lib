#include "soplex.h"
#include "ibex_CajaSPx.h"
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>
using namespace soplex;
namespace ibex {
struct SPxCosas {
    SPxOut out;
    SPxSolver sx;
    SLUFactor lu;
    SPxBoundFlippingRT rt;
    SPxAutoPR pr;
    SPxEquiliSC esc;
    bool escalado = false;
    long nlp = 0, nit = 0, nprim = 0, nlp_prim = 0, nlp_mixto = 0;
    SPxCosas() : sx(SPxSolver::LEAVE, SPxSolver::COLUMN, Timer::OFF), esc(true) {
        out.setVerbosity(SPxOut::ERROR);
        sx.setOutstream(out);
        esc.setOutstream(out);
        lu.setUtype(SLUFactor::FOREST_TOMLIN);
        lu.changeTimer(Timer::OFF);           /* como la clase SoPlex con timer=0 */
        sx.setBasisSolver(&lu);
        sx.setTester(&rt);
        sx.setPricer(&pr);
        sx.setTerminationIter(100);          /* el tope que Ibex da a SoPlex */
        /* los mismos ajustes que la clase SoPlex con los parametros de Ibex */
        sx.setSparsePricingFactor(0.6);
        sx.setNonzeroFactor(10.0);
        sx.setFillFactor(5.0);
        sx.setMemFactor(1.5);
    }
};
#define C (*(SPxCosas*) p)
CajaSPx::CajaSPx() : p(new SPxCosas()), n_(0), m_(0) {}
CajaSPx::~CajaSPx() {
    if (getenv("DFBH_CAJA_VERB"))
        fprintf(stderr, "[spx] LPs %ld factorizaciones %d (%.3f por LP)\n", C.nlp, C.lu.getFactorCount(),
                C.nlp ? (double) C.lu.getFactorCount() / C.nlp : 0.0);
    delete (SPxCosas*) p;
}
static Real inf_s(double v) { return v >= 1e100 ? infinity : (v <= -1e100 ? -infinity : v); }
void CajaSPx::cargar(int n, int m, const double* lb, const double* ub, const double* lhs,
                     const double* rhs, const double* A, double tol, bool escalar) {
    SPxLP lp;
    lp.changeSense(SPxLP::MINIMIZE);
    LPColSet cols(n, n * m + 1);
    DSVector vacio(0);
    for (int i = 0; i < n; ++i) cols.add(0.0, inf_s(lb[i]), vacio, inf_s(ub[i]));
    lp.addCols(cols);
    LPRowSet rows(m, n * m + 1);
    DSVector fila(n);
    for (int j = 0; j < m; ++j) {
        fila.clear();
        const Real lo = inf_s(lhs[j]), up = inf_s(rhs[j]);
        const bool libre = (lo <= -infinity && up >= infinity);
        if (!libre)
            for (int i = 0; i < n; ++i) {
                const double a = A[(long) j * n + i];
                if (a != 0.0) fila.add(i, a);
            }
        rows.add(lo, fila, up);
    }
    lp.addRows(rows);
    C.escalado = escalar;
    if (escalar) C.esc.scale(lp, true);
    C.sx.loadLP(lp);
    C.sx.setFeastol(tol < 1e-12 ? 1e-12 : tol);
    C.sx.setOpttol(tol < 1e-12 ? 1e-12 : tol);
    /* la regla de representacion de la clase SoPlex (REPRESENTATION_AUTO, 1.2) */
    if ((n + 1) * 1.2 >= (m + 1)) { C.sx.setRep(SPxSolver::COLUMN); C.sx.setType(SPxSolver::LEAVE); }
    else                          { C.sx.setRep(SPxSolver::ROW);    C.sx.setType(SPxSolver::ENTER); }
    /* `DFBH_CAJA_TOPE=igual`: el tope de iteraciones de SC (20m+200) en vez del de Ibex (100) */
    static const bool tope_igual = getenv("DFBH_CAJA_TOPE") && std::string(getenv("DFBH_CAJA_TOPE")) == "igual";
    C.sx.setTerminationIter(tope_igual ? 20*m + 200 : 100);
    n_ = n; m_ = m;
}
void CajaSPx::cota(int i, double lb, double ub) {
    C.sx.changeBounds(i, inf_s(lb), inf_s(ub), C.escalado);
}
long CajaSPx::n_fact() const { return 0; }
int CajaSPx::resolver(int k, bool maximizar, bool frio, double* lam, double* x, int* iters) {
    SPxSolver& s = C.sx;
    ++C.nlp;
    static const bool verb = getenv("DFBH_CAJA_VERB") != NULL;
    if (verb && C.nlp % 2000 == 0)
        fprintf(stderr, "[spx] LPs %ld factorizaciones %d (%.3f por LP) | pivotes %ld, primales %.1f%% | LPs solo primales %.1f%%, mixtos %.1f%%\n",
                C.nlp, C.lu.getFactorCount(), (double) C.lu.getFactorCount() / C.nlp,
                C.nit, C.nit ? 100.0*C.nprim/C.nit : 0.0, 100.0*C.nlp_prim/C.nlp, 100.0*C.nlp_mixto/C.nlp);
    s.changeObj(k, maximizar ? -1.0 : 1.0, C.escalado);
    if (frio) s.reLoad();
    /* `DFBH_SPX_DUAL=1`: como la clase SoPlex, vuelve al simplex dual en cada
     * llamada; si no, SPxSolver sigue con el tipo en que termino */
    static const bool fija_dual = getenv("DFBH_SPX_DUAL") != NULL;
    if (fija_dual) {
        if (s.rep() == SPxSolver::COLUMN && s.type() != SPxSolver::LEAVE) s.setType(SPxSolver::LEAVE);
        if (s.rep() == SPxSolver::ROW && s.type() != SPxSolver::ENTER) s.setType(SPxSolver::ENTER);
    }
    SPxSolver::Status st;
    try { st = s.solve(); } catch (const SPxException& e) {
        st = SPxSolver::ERROR;
        static int ne = 5; if (ne-- > 0) fprintf(stderr, "[spx] excepcion: %s\n", e.what().c_str());
    }
    *iters = s.iterations();
    {
        const int it = s.iterations(), pr = s.primalIterations();
        C.nit += it; C.nprim += pr;
        if (it > 0 && pr == it) ++C.nlp_prim; else if (pr > 0) ++C.nlp_mixto;
    }
    static int dbg = verb ? 20 : 0;
    if (dbg > 0) { --dbg; fprintf(stderr, "[spx] k %d st %d it %d rep %d tipo %d filas %d cols %d\n", k, (int) st, s.iterations(), (int) s.rep(), (int) s.type(), s.nRows(), s.nCols()); }
    int res = 3;
    static long hist[32] = {0}; static long nh = 0;
    if (verb) {
        ++hist[(int) st + 16]; ++nh;
        if (nh % 20000 == 0) { fprintf(stderr, "[spx] estados:"); for (int q = 0; q < 32; ++q) if (hist[q]) fprintf(stderr, " %d:%ld", q - 16, hist[q]); fprintf(stderr, "\n"); }
    }
    if (st == SPxSolver::OPTIMAL) {
        DVector y(m_), z(n_);
        if (s.getDual(y) == SPxSolver::OPTIMAL && s.getPrimal(z) == SPxSolver::OPTIMAL) {
            if (C.escalado) { C.esc.unscaleDual(s, y); C.esc.unscalePrimal(s, z); }
            for (int j = 0; j < m_; ++j) lam[j] = y[j];
            for (int i = 0; i < n_; ++i) x[i] = z[i];
            res = 0;
        }
    } else if (st == SPxSolver::INFEASIBLE) {
        DVector f(m_);
        bool ok = false;
        try { ok = (s.getDualfarkas(f) == SPxSolver::INFEASIBLE); } catch (const SPxException&) { ok = false; }
        if (ok) {
            if (C.escalado) C.esc.unscaleDual(s, f);
            for (int j = 0; j < m_; ++j) lam[j] = f[j];
            res = 1;
        } else res = 2;
    }
    s.changeObj(k, 0.0, C.escalado);
    return res;
}
}
