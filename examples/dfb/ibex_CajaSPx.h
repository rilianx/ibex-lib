/* Caja negra con el nucleo de SoPlex 4.0.2 (SPxSolver) usado directamente
 * (2026-10-19): el LP queda cargado y entre cotas solo cambian objetivo y
 * cotas de columna, asi que la factorizacion de la base se conserva entre LPs
 * (la clase SoPlex la rehace en cada llamada). Mismo LP que la caja directa,
 * misma regla de representacion y algoritmo dual que la clase SoPlex, sin
 * presolve ni verificacion de la solucion. Interfaz sin tipos de SoPlex. */
#ifndef __IBEX_CAJA_SPX_H__
#define __IBEX_CAJA_SPX_H__
namespace ibex {
class CajaSPx {
public:
    CajaSPx();
    ~CajaSPx();
    void cargar(int n, int m, const double* lb, const double* ub, const double* lhs,
                const double* rhs, const double* A, double tol, bool escalar);
    void cota(int i, double lb, double ub);
    /* 0 optimo, 1 infactible con rayo, 2 infactible sin rayo, 3 sin optimo */
    int resolver(int k, bool maximizar, bool frio, double* lam, double* x, int* iters);
    long n_fact() const;
private:
    void* p;
    int n_, m_;
};
}
#endif
