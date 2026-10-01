/* Equivalencia entre el tableau denso de DFB y la base factorizada.
 *
 * Reproduce a mano la eliminacion del camino denso (make_column_identity_f)
 * sobre `[Abar | I]` y, con la MISMA secuencia de pivotes, consulta la base
 * factorizada. Si la afirmacion «el tableau transformado es B^-1 [Abar | I]»
 * es correcta, todas las filas deben coincidir.
 *
 * Es la prueba que hay que pasar antes de cablear DFBBasis dentro de CtcDFB:
 * un error aca seria invisible mas adelante. */
#include "ibex_DFBBasis.h"
#include "ibex_Matrix.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace ibex;

namespace {

const double TOL = 1e-9;

/* Eliminacion del camino denso, copiada de make_column_identity_f. */
bool dense_identity(Matrix& Af, int col, int row) {
    const int m = Af.nb_rows(), nc = Af.nb_cols();
    const double piv = Af[row][col];
    if (std::fabs(piv) < 1e-10) return false;
    const double inv = 1.0 / piv;
    for (int c = 0; c < nc; ++c) Af[row][c] *= inv;
    Af[row][col] = 1.0;
    for (int i = 0; i < m; ++i) {
        if (i == row) continue;
        const double f = Af[i][col];
        if (f == 0.0) continue;
        for (int c = 0; c < nc; ++c) Af[i][c] -= f * Af[row][c];
        Af[i][col] = 0.0;
    }
    return true;
}

double frand() { return (double)rand() / RAND_MAX * 2.0 - 1.0; }

} /* anonymous */

int main(int argc, char** argv) {
    const int reps = (argc > 1) ? atoi(argv[1]) : 40;
    srand(12345);

    int fallos = 0, casos = 0;
    double peor = 0.0;

    for (int rep = 0; rep < reps; ++rep) {
        const int m  = 3 + rand() % 12;
        const int nv = 2 + rand() % 8;
        const int na = nv + m;
        const double dens = 0.3 + 0.7 * ((double)rand() / RAND_MAX);

        Matrix Abar(m, na);
        for (int j = 0; j < m; ++j)
            for (int c = 0; c < na; ++c)
                Abar[j][c] = ((double)rand() / RAND_MAX < dens) ? frand() : 0.0;
        /* la parte de holguras de la linealizacion es -I, como en refA */
        for (int j = 0; j < m; ++j)
            for (int c = nv; c < na; ++c) Abar[j][c] = (c - nv == j) ? -1.0 : 0.0;

        const int k = rand() % nv;

        /* --- camino denso --- */
        Matrix Af(m, na + m);
        for (int j = 0; j < m; ++j) {
            for (int c = 0; c < na; ++c) Af[j][c] = Abar[j][c];
            for (int c = 0; c < m; ++c)  Af[j][na + c] = (j == c) ? 1.0 : 0.0;
        }
        /* interchange: la fila objetivo es la primera con |Abar[j][k]| usable */
        int obj = -1;
        for (int j = 0; j < m; ++j)
            if (std::fabs(Abar[j][k]) >= 1e-10) { obj = j; break; }
        if (obj < 0) continue;
        if (obj != 0) for (int c = 0; c < na + m; ++c) std::swap(Af[0][c], Af[obj][c]);
        if (!dense_identity(Af, k, 0)) continue;

        /* --- base factorizada --- */
        DFBBasis B;
        B.load_matrix(Abar);
        if (!B.init_basis(k, 1e-10)) { printf("rep %d: init_basis fallo\n", rep); ++fallos; continue; }
        if (B.objective_row() != obj) {
            printf("rep %d: obj_row %d != %d\n", rep, B.objective_row(), obj); ++fallos; continue;
        }

        /* La fila del denso que corresponde a la fila r del factorizado:
         * el denso permuto 0 <-> obj, el factorizado no. */
        std::vector<int> dense_of(m);
        for (int r = 0; r < m; ++r) dense_of[r] = r;
        if (obj != 0) { dense_of[0] = obj; dense_of[obj] = 0; }

        std::vector<double> lam(m), g(na);

        /* --- misma secuencia de pivotes en los dos --- */
        const int npiv = 1 + rand() % 5;
        for (int p = 0; p <= npiv; ++p) {

            /* comparar todas las filas */
            for (int r = 0; r < m; ++r) {
                B.btran_row(r, &lam[0]);
                B.row_of_A(&lam[0], &g[0]);
                const int dr = dense_of[r];
                for (int c = 0; c < na; ++c) {
                    const double d = std::fabs(g[c] - Af[dr][c]);
                    if (d > peor) peor = d;
                    if (d > TOL) {
                        printf("rep %d piv %d: fila %d col %d  factor=%g denso=%g\n",
                               rep, p, r, c, g[c], Af[dr][c]);
                        ++fallos;
                    }
                    ++casos;
                }
                for (int c = 0; c < m; ++c) {   /* parte lambda */
                    const double d = std::fabs(lam[c] - Af[dr][na + c]);
                    if (d > peor) peor = d;
                    if (d > TOL) {
                        printf("rep %d piv %d: fila %d lambda %d  factor=%g denso=%g\n",
                               rep, p, r, c, lam[c], Af[dr][na + c]);
                        ++fallos;
                    }
                    ++casos;
                }
            }
            if (p == npiv) break;

            /* elegir un pivote valido: fila j != obj_row, columna i no basica
             * con |Af[j][i]| razonable */
            int j = -1, col = -1;
            for (int tries = 0; tries < 60 && col < 0; ++tries) {
                const int jj = rand() % m;
                if (jj == B.objective_row()) continue;
                const int cc = rand() % (na + m);
                if (B.is_basic(cc)) continue;
                if (std::fabs(Af[dense_of[jj]][cc]) < 1e-3) continue;
                j = jj; col = cc;
            }
            if (col < 0) break;

            /* denso: fila 0 += alpha*fila j (alpha del ratio test simplificado),
             * luego identidad de la columna entrante en la fila j */
            const int dj = dense_of[j];
            const double alpha = -(Af[0][col] / Af[dj][col]);
            for (int c = 0; c < na + m; ++c) Af[0][c] += alpha * Af[dj][c];
            Af[0][col] = 0.0;
            if (!dense_identity(Af, col, dj)) break;
            Af[0][k] = 1.0;

            /* factorizado: solo cambia la columna basica de la fila j */
            if (!B.change_basis(j, col)) { printf("rep %d: change_basis fallo\n", rep); ++fallos; break; }
        }
    }

    printf("\ncasos comparados: %d   fallos: %d   peor diferencia: %.3g\n",
           casos, fallos, peor);
    printf("factorizaciones=%ld updates=%ld btran=%ld ftran=%ld\n",
           DFBBasis::n_factorizations, DFBBasis::n_updates,
           DFBBasis::n_btran, DFBBasis::n_ftran);
    return fallos ? 1 : 0;
}
