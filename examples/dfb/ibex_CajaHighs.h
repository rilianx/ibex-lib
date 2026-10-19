/* Caja negra con HiGHS (2026-10-20): mismo LP que la caja directa de SoPlex
 * (filas y lados de la relajacion, cotas de columna, filas libres vacias), con
 * la base de HiGHS conservada entre LPs. Interfaz sin tipos de HiGHS. */
#ifndef __IBEX_CAJA_HIGHS_H__
#define __IBEX_CAJA_HIGHS_H__
namespace ibex {
class CajaHighs {
public:
    CajaHighs();
    ~CajaHighs();
    void cargar(int n, int m, const double* lb, const double* ub, const double* lhs,
                const double* rhs, const double* A, double tol, bool presolve);
    void cota(int i, double lb, double ub);
    /* 0 optimo, 1 infactible con rayo, 2 infactible sin rayo, 3 sin optimo */
    int resolver(int k, bool maximizar, bool frio, double* lam, double* x, int* iters);
private:
    void* h;
    int n_, m_;
};
}
#endif
