#include "ibex_LSmearDFB.h"
#include "ibex_LSmear.h"
#include <cstdio>
#include <string>
#include <cmath>
#include <cstdlib>

namespace ibex {

long LSmearDFB::n_llamadas = 0, LSmearDFB::n_sin_sens = 0, LSmearDFB::n_fuera = 0,
     LSmearDFB::n_dfb = 0, LSmearDFB::n_iguales = 0;

LSmearDFB::LSmearDFB(ExtendedSystem& sys, const Vector& prec, LargestFirst& lf, Bsc* comparar) :
    SmearSumRelative(sys, prec, lf), esys(sys), gv(sys.goal_var()), comparar(comparar) {}

LSmearDFB::LSmearDFB(System& sys, double prec) :
    SmearSumRelative(sys, prec), esys(sys), gv(-1), comparar(NULL) {}

BisectionPoint LSmearDFB::choose_var(const Cell& cell) {
    ++n_llamadas;
    const IntervalVector& box = cell.box;
    const BxpDFBSens* p = (const BxpDFBSens*) cell.prop[BxpDFBSens::id()];
    int lvar = -1;
    if (!p || !p->ok) { ++n_sin_sens; }
    else {
        bool dentro = (p->caja.size() == box.size());
        for (int j = 0; dentro && j < box.size(); ++j)
            if (!p->caja[j].is_superset(box[j])) dentro = false;
        if (!dentro) { ++n_fuera; }
        else {
            /* `DFBH_SENS=jac`: los multiplicadores de DFB por restriccion, con
             * signo, proyectados sobre la jacobiana del medio como en LSmear. */
            static const bool jac = getenv("DFBH_SENS") && std::string(getenv("DFBH_SENS")) == "jac";
            IntervalMatrix J(1,1);
            const bool usa_jac = jac && (int)p->lam_ctr.size() == esys.f_ctrs.image_dim();
            if (usa_jac) {
                J.resize(esys.f_ctrs.image_dim(), esys.nb_var);
                IntervalVector bm(box.mid()); bm.inflate(1e-8); bm &= box;
                esys.f_ctrs.jacobian(bm, J);
            }
            double mejor = 0.0; int k = 0;
            for (int j = 0; j < esys.nb_var; ++j) {
                if (too_small(box, j)) continue;
                if (!(box[j].mag() < 1 || box[j].diam()/box[j].mag() >= prec(j))) continue;
                double sc;
                if (j == gv) sc = p->contrajo ? 1.0 : 0.0;
                else if (usa_jac) {
                    double acc = 0.0;
                    for (int c = 0; c < esys.f_ctrs.image_dim(); ++c) acc += p->lam_ctr[c] * J[c][j].mid();
                    sc = std::fabs(acc);
                } else sc = p->sens[j];
                sc *= box[j].diam();
                if (getenv("IBEX_LSMEAR_DUMP")) {
                    static long nd = 0;
                    if (nd < 8*esys.nb_var) { ++nd; fprintf(stderr, "[dfb-dump] x%d score=%.3g diam=%.3g\n", j, sc, box[j].diam()); }
                }
                if (sc > 1e-10) { ++k; if (sc > mejor) { mejor = sc; lvar = j; } }
            }
            if (k == 1 && lvar == gv) lvar = -1;
            if (lvar >= 0) ++n_dfb;
        }
    }
    if (comparar) {
        const int otro = comparar->choose_var(cell).var;
        if (lvar >= 0 && otro == lvar) ++n_iguales;
        if ((n_llamadas % 2000) == 0)
            fprintf(stderr, "[lsmeardfb] llamadas=%ld | sin sensibilidades=%ld fuera de la caja=%ld | "
                    "eligio DFB=%ld | coincide con LSmear=%ld (%.1f%% de las que eligio DFB)\n",
                    n_llamadas, n_sin_sens, n_fuera, n_dfb, n_iguales, n_dfb ? 100.0*n_iguales/n_dfb : 0.0);
    }
    if (lvar >= 0 && box[lvar].is_bisectable()) return BisectionPoint(lvar, lf->ratio, true);
    return SmearSumRelative::choose_var(cell);
}

/* Sin celda no hay sensibilidades: respaldo. */
int LSmearDFB::var_to_bisect(IntervalMatrix& J, const IntervalVector& box) const {
    return SmearSumRelative::var_to_bisect(J, box);
}

} /* namespace ibex */
