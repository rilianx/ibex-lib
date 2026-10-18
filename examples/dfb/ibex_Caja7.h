/* Caja negra con SoPlex 7 (2026-10-17). SoPlex 7 se compila con el namespace
 * renombrado (`-Dsoplex=soplex7`) para convivir con el SoPlex 4.0.2 de Ibex;
 * esta interfaz no expone tipos de SoPlex. */
#ifndef __IBEX_CAJA7_H__
#define __IBEX_CAJA7_H__
namespace ibex {
class Caja7 {
public:
    Caja7();
    ~Caja7();
    /* A en filas (m x n, densa). Las filas con los dos lados infinitos se cargan vacias. */
    void cargar(int n, int m, const double* lb, const double* ub, const double* lhs,
                const double* rhs, const double* A, double tol, int iterlimit);
    void cota(int i, double lb, double ub);
    /* 0 optimo, 1 infactible con rayo, 2 infactible sin rayo, 3 sin optimo */
    int resolver(int k, bool maximizar, bool frio, double* lam, double* x, int* iters);
private:
    void* sx;
};
}
#endif
