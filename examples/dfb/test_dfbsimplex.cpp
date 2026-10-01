/* Valida el simplex dual de DFB contra SoPlex.
 *
 * Sobre relajaciones aleatorias con la misma estructura que las de DFB (la
 * parte b es -I), resuelve `min z_k` con DFBSimplex y con LPSolver, y compara.
 * Es la prueba que decide si el simplex alcanza el optimo del LP, que es todo
 * el objetivo del ejercicio.
 *
 * Tambien verifica la dualidad fuerte con los multiplicadores devueltos:
 * f(y) = -sum_{i!=k} ub(gamma_i*[z_i]) con gamma = c - d debe coincidir con el
 * optimo. Si coincide, `lambda()` sirve para certificar la cota en intervalos.
 */
#include "ibex_DFBSimplex.h"
#include "ibex_LPSolver.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace ibex;

namespace {

double frand() { return (double)rand() / RAND_MAX * 2.0 - 1.0; }

/* Cota dual a partir de gamma: -sum_{i!=k} ub(gamma_i * [z_i]), dividida por
 * gamma_k. Es exactamente lo que hace la certificacion de DFB. */
double dual_bound(const std::vector<double>& gamma, const IntervalVector& z,
                  int k) {
    const double gk = gamma[k];
    if (std::fabs(gk) < 1e-12) return -1e300;
    double S_ub = 0.0;
    for (size_t i = 0; i < gamma.size(); ++i) {
        if ((int)i == k) continue;
        const double g = gamma[i];
        if (g == 0.0) continue;
        const double a = g * z[i].lb(), b = g * z[i].ub();
        S_ub += (a > b) ? a : b;
    }
    /* gk*z_k = -sum_{i!=k} g_i z_i  =>  z_k >= -S_ub/gk  si gk > 0 */
    return (gk > 0) ? (-S_ub / gk) : (-S_ub / gk);
}

} /* anonymous */

int main(int argc, char** argv) {
    const int reps = (argc > 1) ? atoi(argv[1]) : 200;
    /* Con SX_ESCALA=1 se usan escalas extremas (coeficientes hasta 1e36, cajas
     * hasta 1e8), que es el regimen de Brown-* y donde el simplex FALLA: 217
     * discrepancias de 400 en frio y 107 de 641 en tibio. Sin escalado previo
     * de filas y columnas no es robusto ahi. Queda como caso de prueba
     * reproducible; ver MEDICIONES_PODA.md §7.12. */
    const bool escala_extrema = (getenv("SX_ESCALA") != NULL);
    /* SX_ENCOGE=1: la caja se aprieta entre cotas, como en la propagacion. Es
     * la dimension que faltaba: el camino tibio sobre caja FIJA pasa
     * 1066/1066, y sobre caja cambiante falla al maximizar (ver §8.4). */
    const bool encoge = (getenv("SX_ENCOGE") != NULL);
    printf("escalas: %s\n", escala_extrema ?
           "calibradas al banco (filas 1..1e6, cajas 1..1e6)" : "moderadas");
    srand(4242);

    int casos = 0, ok = 0, dif = 0, infeas_ok = 0, iterlim = 0, singular = 0;
    int dual_ok = 0, dual_mal = 0, dual_mal2 = 0, inv_ok = 0, con_cota = 0;
    double peor = 0.0;

    for (int rep = 0; rep < reps; ++rep) {
        const int nx = 2 + rand() % 6;
        const int m  = 2 + rand() % 6;
        const int na = nx + m;
        const double dens = 0.4 + 0.6 * ((double)rand() / RAND_MAX);

        /* Mismas escalas extremas que en el bloque tibio: sin esto el test no
         * cubre el regimen de Brown-* (coeficientes ~1e72, cajas de radio 1e9),
         * que es donde aparecen los problemas numericos reales. */
        /* Regimen FIEL a las instancias reales: medido con `spread`, la
         * dispersion de los terminos a*z DENTRO de una fila es ~8-11 en
         * Brown-* y hasta 1e8 en Katsura/Virasoro. Lo que varia mucho es la
         * escala ENTRE filas y el ancho de las cajas. Generar coeficientes
         * dispares dentro de una misma fila —como hacia la version anterior de
         * este test— inventa una patologia que no ocurre. */
        const double escz_w = escala_extrema ? std::pow(10.0, ((rep / 4) % 3) * 3.0) : 1.0;
        Matrix Abar(m, na);
        for (int j = 0; j < m; ++j) {
            for (int c = 0; c < na; ++c) Abar[j][c] = 0.0;
            const double rj = escala_extrema ?
                std::pow(10.0, (double)((j * 7 + rep) % 5) * 1.5) : 1.0;
            for (int c = 0; c < nx; ++c)
                if ((double)rand() / RAND_MAX < dens) Abar[j][c] = frand() * rj;
            Abar[j][nx + j] = -1.0;
        }

        IntervalVector z(na);
        for (int c = 0; c < nx; ++c) {
            const double a = frand() * 5.0 * escz_w;
            const double w = (0.5 + 3.0 * ((double)rand() / RAND_MAX)) * escz_w;
            z[c] = Interval(a, a + w);
        }
        /* cotas de las b: la evaluacion de la fila, encogida al azar para que a
         * veces el LP sea infactible, igual que en el uso real */
        for (int j = 0; j < m; ++j) {
            Interval e(0.0);
            for (int c = 0; c < nx; ++c) e += Interval(Abar[j][c]) * z[c];
            const double shrink = 0.3 + 0.9 * ((double)rand() / RAND_MAX);
            const double mid = e.mid(), rad = e.rad() * shrink;
            z[nx + j] = Interval(mid - rad, mid + rad);
        }

        const int k = rand() % nx;
        const bool maxi = (rand() % 2) == 1;   /* la mitad de los contractores
                                                  de DFB piden la cota superior */

        DFBSimplex sx;
        sx.load(Abar, nx);
        DFBSimplex::Status st = sx.solve(k, maxi, z, 500);

        /* referencia: SoPlex */
        LPSolver lp(na, LPSolver::Mode::NotCertified, 1e-9, 10.0, 100000);
        for (int j = 0; j < m; ++j) {
            Vector row(na, 0.0);
            for (int c = 0; c < na; ++c) row[c] = Abar[j][c];
            lp.add_constraint(0.0, row, 0.0);
        }
        lp.set_bounds(z);
        Vector cost(na, 0.0); cost[k] = maxi ? -1.0 : 1.0;
        lp.set_cost(cost);
        LPSolver::Status lst = lp.minimize();

        ++casos;
        if (lst == LPSolver::Status::Infeasible ||
            lst == LPSolver::Status::InfeasibleProved) {
            if (st == DFBSimplex::INFEASIBLE) ++infeas_ok;
            else if (st == DFBSimplex::ITER_LIMIT) ++iterlim;
            else if (st == DFBSimplex::SINGULAR) ++singular;
            else {
                printf("rep %d: SoPlex infactible, simplex dice OPTIMAL bnd=%g\n",
                       rep, sx.bound());
                ++dif;
            }
            continue;
        }
        if (lst != LPSolver::Status::Optimal && lst != LPSolver::Status::OptimalProved)
            continue;   /* unbounded/timeout: no es un caso de prueba */

        const double ref = lp.minimum().lb();
        if (st == DFBSimplex::ITER_LIMIT) { ++iterlim; continue; }
        if (st == DFBSimplex::SINGULAR)   { ++singular; continue; }
        if (st == DFBSimplex::INFEASIBLE) {
            printf("rep %d: simplex dice INFEASIBLE, SoPlex da %g\n", rep, ref);
            ++dif; continue;
        }

        const double zk_ref = maxi ? -ref : ref;
        const double err = std::fabs(sx.bound() - zk_ref) /
                           (1.0 + std::fabs(ref));
        if (err > peor) peor = err;
        if (err < 1e-7) ++ok;
        else {
            if (dif < 8)
                printf("rep %d: nx=%d m=%d k=%d  simplex=%.12g  SoPlex=%.12g  "
                       "err=%.3g  (%d pivotes)\n", rep, nx, m, k,
                       sx.bound(), zk_ref, err, sx.iterations());
            ++dif;
        }

        /* Invariante que hace certificable la cota: gamma = y^T*Abar = c - d.
         * Si no se cumple, lambda() no sirve para certificar. */
        std::vector<double> gamma(na, 0.0);
        const double ck = maxi ? -1.0 : 1.0;
        for (int c = 0; c < na; ++c) gamma[c] = ((c == k) ? ck : 0.0) - sx.reduced_costs()[c];
        double err_inv = 0.0;
        for (int c = 0; c < na; ++c) {
            double acc = 0.0;
            for (int j = 0; j < m; ++j) acc += sx.lambda()[j] * Abar[j][c];
            const double e = std::fabs(acc - gamma[c]);
            if (e > err_inv) err_inv = e;
        }
        if (err_inv > 1e-9) {
            if (dual_mal < 5)
                printf("   rep %d: gamma != y^T*Abar, error %.3g\n", rep, err_inv);
            ++dual_mal;
        } else ++inv_ok;

        /* Dualidad fuerte, solo donde hay cota que certificar (gamma_k != 0):
         * si gamma_k = 0 el optimo es la propia cota de z_k y no hay
         * contraccion, que es el caso en que la certificacion de DFB devuelve
         * false. */
        if (std::fabs(gamma[k]) > 1e-9) {
            ++con_cota;
            const double db = dual_bound(gamma, z, k);
            if (std::fabs(db - zk_ref) / (1.0 + std::fabs(zk_ref)) < 1e-7) ++dual_ok;
            else {
                if (dual_mal2 < 5)
                    printf("   rep %d: cota dual %.12g != optimo %.12g (gamma_k=%g, max=%d)\n",
                           rep, db, zk_ref, gamma[k], (int)maxi);
                ++dual_mal2;
            }
        }
    }

    /* ---- camino TIBIO: varias cotas sobre la misma base ----
     * Es el camino que no estaba probado, y el que hace util el warm start. */
    printf("\n=== warm start: secuencias de cotas sobre una sola base ===\n");
    int w_casos = 0, w_ok = 0, w_dif = 0; double w_peor = 0.0;
    srand(777);
    for (int rep = 0; rep < reps / 4; ++rep) {
        const int nx = 3 + rand() % 5;
        const int m  = 2 + rand() % 5;
        const int na = nx + m;
        const double dens = 0.5 + 0.5 * ((double)rand() / RAND_MAX);
        /* Escalas extremas a proposito: las linealizaciones de Brown-* tienen
         * coeficientes de ~1e72 sobre cajas de radio 1e9, y es ahi donde el
         * camino tibio falla en el banco real (una cota superior devuelve la
         * cota INFERIOR de la variable declarando optimalidad). El test con
         * escalas moderadas no lo reproduce. */
        /* Regimen FIEL a las instancias reales: medido con `spread`, la
         * dispersion de los terminos a*z DENTRO de una fila es ~8-11 en
         * Brown-* y hasta 1e8 en Katsura/Virasoro. Lo que varia mucho es la
         * escala ENTRE filas y el ancho de las cajas. Generar coeficientes
         * dispares dentro de una misma fila —como hacia la version anterior de
         * este test— inventa una patologia que no ocurre. */
        const double escz_w = escala_extrema ? std::pow(10.0, ((rep / 4) % 3) * 3.0) : 1.0;
        Matrix Abar(m, na);
        for (int j = 0; j < m; ++j) {
            for (int c = 0; c < na; ++c) Abar[j][c] = 0.0;
            const double rj = escala_extrema ?
                std::pow(10.0, (double)((j * 7 + rep) % 5) * 1.5) : 1.0;
            for (int c = 0; c < nx; ++c)
                if ((double)rand() / RAND_MAX < dens) Abar[j][c] = frand() * rj;
            Abar[j][nx + j] = -1.0;
        }
        IntervalVector z(na);
        const double escz = escala_extrema ? std::pow(10.0, ((rep / 4) % 3) * 3.0) : 1.0;
        for (int c = 0; c < nx; ++c) {
            const double a = frand() * 5.0 * escz;
            const double w = (0.5 + 3.0 * ((double)rand() / RAND_MAX)) * escz;
            z[c] = Interval(a, a + w);
        }
        /* Cotas de las b: en el uso real son `lhs_rhs & (fila . caja)`, o sea
         * a menudo el enclosure COMPLETO de la fila, no una version encogida.
         * Con SX_BANCHAS=1 se usa el enclosure entero, que es el caso real. */
        const bool b_anchas = (getenv("SX_BANCHAS") != NULL);
        for (int j = 0; j < m; ++j) {
            Interval e(0.0);
            for (int c = 0; c < nx; ++c) e += Interval(Abar[j][c]) * z[c];
            const double f = b_anchas ? 1.0 : 0.8;
            z[nx + j] = Interval(e.mid() - e.rad() * f, e.mid() + e.rad() * f);
        }

        DFBSimplex sx;
        sx.load(Abar, nx);
        /* Las 2*nx cotas en secuencia, reusando la base. La caja se APRIETA
         * un poco entre cotas, que es lo que pasa en la propagacion real y es
         * la dimension que el test no cubria: antes `z` era fijo. */
        for (int kk = 0; kk < nx; ++kk)
            for (int up = 0; up < 2; ++up) {
                if (encoge && (kk + up) > 0) {
                    const int v = (kk * 2 + up) % nx;
                    const double r = z[v].rad() * 0.05;
                    if (r > 0) z[v] = Interval(z[v].lb() + r, z[v].ub() - r);
                }
                DFBSimplex::Status st = sx.solve(kk, up == 1, z, 500, true);
                LPSolver lp(na, LPSolver::Mode::NotCertified, 1e-9, 10.0, 100000);
                for (int j = 0; j < m; ++j) {
                    Vector row(na, 0.0);
                    for (int c = 0; c < na; ++c) row[c] = Abar[j][c];
                    lp.add_constraint(0.0, row, 0.0);
                }
                lp.set_bounds(z);
                Vector cost(na, 0.0); cost[kk] = (up == 1) ? -1.0 : 1.0;
                lp.set_cost(cost);
                LPSolver::Status lst = lp.minimize();
                if (lst != LPSolver::Status::Optimal &&
                    lst != LPSolver::Status::OptimalProved) continue;
                if (st != DFBSimplex::OPTIMAL) continue;
                const double zk = (up == 1) ? -lp.minimum().lb() : lp.minimum().lb();
                const double e = std::fabs(sx.bound() - zk) / (1.0 + std::fabs(zk));
                ++w_casos;
                if (e > w_peor) w_peor = e;
                if (e < 1e-7) ++w_ok;
                else {
                    if (w_dif < 8)
                        printf("  tibio rep %d k=%d up=%d: simplex=%.12g SoPlex=%.12g err=%.3g\n",
                               rep, kk, up, sx.bound(), zk, e);
                    ++w_dif;
                }
            }
    }
    printf("  cotas comparadas: %d   coinciden: %d   DISCREPANCIAS: %d  (peor %.3g)\n",
           w_casos, w_ok, w_dif, w_peor);
    printf("  tibias=%ld frias=%ld pivotes=%ld\n",
           DFBSimplex::n_warm, DFBSimplex::n_cold, DFBSimplex::n_pivots);

    printf("\ncasos: %d\n", casos);
    printf("  optimo coincide con SoPlex : %d\n", ok);
    printf("  DISCREPANCIAS              : %d   (peor error relativo %.3g)\n", dif, peor);
    printf("  infactibilidad detectada   : %d\n", infeas_ok);
    printf("  limite de iteraciones      : %d\n", iterlim);
    printf("  pivote singular            : %d\n", singular);
    printf("  invariante gamma = y^T*Abar : %d ok, %d mal\n", inv_ok, dual_mal);
    printf("  casos con cota certificable  : %d de %d (gamma_k != 0)\n", con_cota, ok);
    printf("  dualidad fuerte en esos      : %d ok, %d mal\n", dual_ok, dual_mal2);
    printf("  estados: opt=%ld infact=%ld iterlim=%ld\n",
           DFBSimplex::n_optimal, DFBSimplex::n_infeasible, DFBSimplex::n_iterlimit);
    printf("  pivotes totales=%ld  solves=%ld\n", DFBSimplex::n_pivots, DFBSimplex::n_solves);
    return (dif || w_dif) ? 1 : 0;
}
