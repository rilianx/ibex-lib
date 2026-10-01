/* Equivalencia bit a bit entre el tableau denso de DFB y el disperso.
 *
 * Las dos representaciones deben dar el MISMO resultado, no uno parecido: el
 * orden de las operaciones por entrada es el mismo (columnas crecientes) y los
 * valores son los mismos, asi que la unica diferencia posible seria un error de
 * implementacion. Por eso la tolerancia es 0.
 *
 * Se reproduce la secuencia completa: init (make_identity de la columna k en la
 * fila 0, con interchange), y luego pivotes fila0 += alpha*fila_j seguidos de
 * make_identity de la columna entrante. */
#include "ibex_DFBSparse.h"
#include "ibex_Matrix.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace ibex;

namespace {

const double TOL = 1e-10;   /* tolerancia de pivote, igual que DFB_F_TOL */

bool dense_identity(Matrix& Af, int col, int row, bool interchange) {
    const int m = Af.nb_rows(), nc = Af.nb_cols();
    if (interchange && std::fabs(Af[row][col]) < TOL) {
        int found = -1;
        for (int i = 0; i < m; ++i)
            if (i != row && std::fabs(Af[i][col]) >= TOL) { found = i; break; }
        if (found < 0) return false;
        for (int c = 0; c < nc; ++c) std::swap(Af[row][c], Af[found][c]);
    }
    const double piv = Af[row][col];
    if (std::fabs(piv) < TOL) return false;
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
    const int reps = (argc > 1) ? atoi(argv[1]) : 50;
    srand(9876);

    long casos = 0; int fallos = 0;
    double peor = 0.0;
    long nnz_denso = 0, nnz_disperso = 0;

    for (int rep = 0; rep < reps; ++rep) {
        const int m  = 4 + rand() % 20;
        const int nv = 2 + rand() % 10;
        const int na = nv + m;
        const int nc = na + m;
        const double dens = 0.05 + 0.5 * ((double)rand() / RAND_MAX);

        Matrix Abar(m, na);
        for (int j = 0; j < m; ++j) {
            for (int c = 0; c < na; ++c)
                Abar[j][c] = ((double)rand() / RAND_MAX < dens) ? frand() : 0.0;
            for (int c = nv; c < na; ++c) Abar[j][c] = (c - nv == j) ? -1.0 : 0.0;
        }
        const int k = rand() % nv;

        Matrix Af(m, nc);
        for (int j = 0; j < m; ++j) {
            for (int c = 0; c < na; ++c) Af[j][c] = Abar[j][c];
            for (int c = 0; c < m; ++c)  Af[j][na + c] = (j == c) ? 1.0 : 0.0;
        }

        DFBSparseTableau S;
        S.build(Abar);

        const bool okd = dense_identity(Af, k, 0, true);
        const bool oks = S.make_identity(k, 0, true, TOL);
        if (okd != oks) { printf("rep %d: init denso=%d disperso=%d\n", rep, okd, oks); ++fallos; continue; }
        if (!okd) continue;

        const int npiv = 1 + rand() % 6;
        for (int p = 0; p <= npiv; ++p) {

            for (int j = 0; j < m; ++j) {
                for (int c = 0; c < nc; ++c) {
                    const double d = std::fabs(S.entry(j, c) - Af[j][c]);
                    if (d > peor) peor = d;
                    if (d != 0.0) {
                        if (fallos < 12)
                            printf("rep %d piv %d: (%d,%d) disperso=%.17g denso=%.17g\n",
                                   rep, p, j, c, S.entry(j, c), Af[j][c]);
                        ++fallos;
                    }
                    ++casos;
                    if (Af[j][c] != 0.0) ++nnz_denso;
                }
            }
            nnz_disperso += S.total_nnz();
            if (p == npiv) break;

            /* pivote valido: fila j >= 1, columna con entrada usable */
            int j = -1, col = -1;
            for (int tries = 0; tries < 80 && col < 0; ++tries) {
                const int jj = 1 + rand() % (m - 1);
                const int cc = rand() % nc;
                if (std::fabs(Af[jj][cc]) < 1e-3) continue;
                if (std::fabs(Af[0][cc]) < 1e-12) continue;
                j = jj; col = cc;
            }
            if (col < 0) break;

            const double alpha = -(Af[0][col] / Af[j][col]);

            for (int c = 0; c < nc; ++c) Af[0][c] += alpha * Af[j][c];
            Af[0][col] = 0.0;
            S.add_to_row0(alpha, j);
            if (col < na) S.g[col] = 0.0; else S.lam[col - na] = 0.0;

            const bool a = dense_identity(Af, col, j, false);
            const bool b = S.make_identity(col, j, false, TOL);
            if (a != b) { printf("rep %d: pivote denso=%d disperso=%d\n", rep, a, b); ++fallos; break; }
            if (!a) break;
            Af[0][k] = 1.0; S.g[k] = 1.0;
        }
    }

    printf("\ncasos comparados: %ld   fallos: %d   peor diferencia: %.3g\n",
           casos, fallos, peor);
    printf("entradas no nulas del denso: %ld   almacenadas por el disperso: %ld  (x%.2f)\n",
           nnz_denso, nnz_disperso, nnz_denso ? (double)nnz_disperso / nnz_denso : 0.0);
    return fallos ? 1 : 0;
}
