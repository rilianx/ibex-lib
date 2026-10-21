/* Reproductor de pasadas (2026-10-18). Carga un problema, construye el mismo
 * contractor lineal que `ibexopt --lr=dfbhull` (sistema extendido, linealizador
 * afin, HC4 interno, variable objetivo) o que `ibexsolve_dfb --lr=dfbhull`
 * (sistema original, linealizador afin), y lo llama sobre cada caja capturada
 * con `DFBH_CAPTURA`, en el mismo orden. Las variantes se eligen por variables
 * de entorno, como en las corridas; con `DFBH_ABCAJA=1` cada pasada se resuelve
 * tambien con la caja negra sobre los mismos LPs.
 * uso: ibexreplay problema.bch opt|sol cajas.txt */
#include "ibex.h"
#include "ibex_CtcDFBHull.h"
#include "ibex_LinearizerAffine2.h"
#include "ibex_DFBSimplex.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>
using namespace ibex;

static double cpu() { struct timespec ts; clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "uso: ibexreplay problema.bch opt|sol cajas.txt\n"); return 1; }
    System sys(argv[1]);
    const bool opt = std::string(argv[2]) == "opt";
    if (opt) {
        /* las mismas omisiones que fija `ibexopt` en su main */
        if (!getenv("DFB_SX_PIVOTEO")) {
            DFBSimplex::pivoteo_dual = false;
            DFBSimplex::dual_puro = false;
            if (!getenv("DFB_SX_REVPRIMAL")) DFBSimplex::rev_primal = true;
            if (!getenv("DFB_SX_MOTOR_M")) DFBSimplex::motor_umbral = 15;
        }
        if (!getenv("DFBH_ORDEN")) setenv("DFBH_ORDEN", "cerca", 0);
    } else if (!getenv("DFBH_ORDEN")) setenv("DFBH_ORDEN", "costo_o_cerca", 0);
    {
        const char* tol = "1e-9";
        if (!getenv("DFB_SX_RAYO"))           setenv("DFB_SX_RAYO", "fila", 0);
        if (!getenv("DFB_SX_PISO"))           setenv("DFB_SX_PISO", "error", 0);
        if (!getenv("DFB_SX_FLIPINC"))        setenv("DFB_SX_FLIPINC", "1", 0);
        if (!getenv("DFBH_FILAS_PARALELAS"))  setenv("DFBH_FILAS_PARALELAS", "paralelas", 0);
        if (!getenv("DFB_SX_FEASTOL"))        setenv("DFB_SX_FEASTOL", tol, 0);
        if (!getenv("DFB_SX_DUALTOL"))        setenv("DFB_SX_DUALTOL", tol, 0);
        if (!getenv("DFB_SX_TOLANCHO"))       setenv("DFB_SX_TOLANCHO", "0", 0);
        if (!getenv("DFBH_CAJA_TOPE"))        setenv("DFBH_CAJA_TOPE", "igual", 0);
    }
    System* s = &sys;
    ExtendedSystem* ext = NULL;
    if (opt) { ext = new ExtendedSystem(sys, 1e-8); s = ext; }
    LinearizerAffine2* lr = new LinearizerAffine2(*s);
    CtcHC4* hc4 = opt ? new CtcHC4(s->ctrs, 0.01, true) : NULL;
    CtcDFBHull* ctc = opt ? new CtcDFBHull(*lr, ext->goal_var(), LPSolver::default_tolerance, hc4)
                          : new CtcDFBHull(*lr);
    FILE* f = fopen(argv[3], "r");
    if (!f) { fprintf(stderr, "no se pudo abrir %s\n", argv[3]); return 1; }
    long n_cajas = 0, n_vacias = 0;
    double t = 0.0, contr = 0.0;
    int n;
    while (fscanf(f, "%d", &n) == 1) {
        IntervalVector box(n);
        bool ok = true;
        for (int i = 0; i < n && ok; ++i) {
            char a[64], b[64];
            if (fscanf(f, "%63s %63s", a, b) != 2) { ok = false; break; }
            box[i] = Interval(strtod(a, 0), strtod(b, 0));
        }
        if (!ok || n != s->nb_var) break;
        IntervalVector entrada(box);
        const double t0 = cpu();
        ctc->contract(box);
        t += cpu() - t0;
        ++n_cajas;
        if (box.is_empty()) { ++n_vacias; contr += 1.0; continue; }
        double per0 = 0, per1 = 0;
        for (int i = 0; i < n; ++i) if (entrada[i].diam() < 1e18) { per0 += entrada[i].diam(); per1 += box[i].diam(); }
        if (per0 > 0) contr += 1.0 - per1 / per0;
    }
    fclose(f);
    printf("[replay] cajas=%ld vacias=%ld contraccion_media=%.5f cpu_contractor=%.4f s us_por_llamada=%.1f LPs=%ld pivotes=%ld\n",
           n_cajas, n_vacias, n_cajas ? contr / n_cajas : 0.0, t, n_cajas ? 1e6 * t / n_cajas : 0.0,
           CtcDFBHull::n_lps, CtcDFBHull::n_pivots);
    delete ctc;
    return 0;
}
