#include <cfenv>
#include <map>
#include <sstream>
//============================================================================
//                                  I B E X
// File        : ibex_CtcDFBHull.cpp
//============================================================================

#include "affine/ibex_LinearizerAffine2.h"
#include <limits>
#include <climits>
#include "ibex_CtcDFBHull.h"
#include "ibex_Caja7.h"
#include "ibex_CajaHighs.h"
#include "ibex_CajaGurobi.h"
#include "ibex_CajaSPx.h"
#include "ibex_LPSolver.h"
#include "ibex_LinearizerCompo.h"

#include <cmath>
#include <cstdlib>
#include <string>
#include <ctime>

/* `getenv` recorre todo el entorno; las sondas que se consultan dentro del
 * bucle de LPs lo pagaban en cada LP (perf: ~1 % del total). Esto lo lee una
 * sola vez por punto de llamada. */
#define DFB_ENV(n) ([]() { static const bool v_ = (getenv(n) != NULL); return v_; }())

namespace {
/* Estado de `DFBH_ADAPT`: umbral, promedio movil de cambios de lado por LP
 * tibio (-1 = sin datos todavia) y nodos resueltos con cada estrategia. */
double adapt_umbral() {
    static const double u = getenv("DFBH_ADAPT") ? atof(getenv("DFBH_ADAPT")) : 0.0;
    return u;
}
double& adapt_ema() { static double e = -1.0; return e; }
/* DIAGNOSTICO DEL PRIMER LP DE LA PASADA (`DFBH_PRIMERDIAG=1`, 2026-10-22):
 * pivotes y contraccion del LP del paso 0 contra los demas; histograma de
 * pivotes del primero. Se imprime al salir. */
bool primerdiag_on() { static const bool v = getenv("DFBH_PRIMERDIAG") != NULL; return v; }
struct PrimerDiag {
    long n[2] = {0,0}, piv[2] = {0,0}, corte[2] = {0,0}, nc[2] = {0,0};
    double ctr[2] = {0,0};
    long hist[7] = {0,0,0,0,0,0,0};   /* 0, 1-5, 6-10, 11-20, 21-40, 41-60, >60 */
    long hpiv[7] = {0,0,0,0,0,0,0};
    ~PrimerDiag() {
        if (!n[0] && !n[1]) return;
        fprintf(stderr, "[primerdiag] primero: LPs=%ld piv/LP=%.1f cortados=%ld contraccion media=%.3f | resto: LPs=%ld piv/LP=%.1f contraccion media=%.3f | primero = %.1f%% de los pivotes\n",
                n[0], (double)piv[0]/std::max(1L,n[0]), corte[0], ctr[0]/std::max(1L,nc[0]),
                n[1], (double)piv[1]/std::max(1L,n[1]), ctr[1]/std::max(1L,nc[1]),
                100.0*piv[0]/std::max(1L,piv[0]+piv[1]));
        const char* et[7] = {"0","1-5","6-10","11-20","21-40","41-60",">60"};
        fprintf(stderr, "[primerdiag] hist primero (LPs%%/pivotes%%):");
        for (int i = 0; i < 7; ++i) fprintf(stderr, " %s:%.0f/%.0f", et[i], 100.0*hist[i]/std::max(1L,n[0]), 100.0*hpiv[i]/std::max(1L,piv[0]));
        fprintf(stderr, "\n");
    }
};
PrimerDiag& primerdiag() { static PrimerDiag d; return d; }
void primerdiag_piv(bool primero, int it, bool cortado) {
    PrimerDiag& d = primerdiag(); const int i = primero ? 0 : 1;
    ++d.n[i]; d.piv[i] += it; if (cortado) ++d.corte[i];
    if (primero) {
        const int b = it == 0 ? 0 : it <= 5 ? 1 : it <= 10 ? 2 : it <= 20 ? 3 : it <= 40 ? 4 : it <= 60 ? 5 : 6;
        ++d.hist[b]; d.hpiv[b] += it;
    }
}
void primerdiag_ctr(bool primero, double c) { PrimerDiag& d = primerdiag(); const int i = primero ? 0 : 1; ++d.nc[i]; d.ctr[i] += c; }
/* Estado de `DFBH_WARM=medir` (2026-10-21): costo por LP (promedio movil) de
 * cada arranque tibio, 0 = dual (reubica) y 1 = primal con devex; nodos con
 * cada uno; estrategia del nodo en curso. */
double g_ws_ema[2] = {-1.0, -1.0};
long g_ws_nodos[2] = {0, 0};
long g_ws_total = 0;
int g_ws_estr = -1;
bool warm_medir() { static const bool v = getenv("DFBH_WARM") && std::string(getenv("DFBH_WARM")) == "medir"; return v; }
long& adapt_nodos(bool dual) { static long n[2] = {0, 0}; return n[dual ? 1 : 0]; }
}

namespace ibex {

long CtcDFBHull::n_lps = 0;
long CtcDFBHull::n_pivots = 0;
long CtcDFBHull::n_congela = 0;
long CtcDFBHull::n_sinpar = 0;
long CtcDFBHull::n_acum_gana = 0;
long CtcDFBHull::n_m_tot = 0, CtcDFBHull::n_m_llam = 0;
long CtcDFBHull::n_orig = 0, CtcDFBHull::n_orig_mejor = 0, CtcDFBHull::n_orig_peor = 0;
long CtcDFBHull::n_orig_igual = 0, CtcDFBHull::n_orig_vacio = 0;
double CtcDFBHull::acc_orig = 0.0;
long CtcDFBHull::n_fila_tot = 0, CtcDFBHull::n_fila_vieja = 0;
long CtcDFBHull::n_apr = 0, CtcDFBHull::n_apr_nada = 0;
double CtcDFBHull::acc_apr = 0.0;
long CtcDFBHull::n_chk = 0, CtcDFBHull::n_chk_ident = 0;
long CtcDFBHull::n_chk_cond = 0, CtcDFBHull::n_chk_viola = 0;
double CtcDFBHull::chk_peor = 0.0, CtcDFBHull::chk_peor_ident = 0.0;
long CtcDFBHull::n_barata = 0;
long CtcDFBHull::n_ctr_evaluadas = 0;
long CtcDFBHull::n_diag = 0, CtcDFBHull::n_diag_vieja = 0, CtcDFBHull::n_diag_nueva = 0, CtcDFBHull::n_diag_igual = 0;
double CtcDFBHull::acc_diag = 0.0, CtcDFBHull::acc_corr_v = 0.0, CtcDFBHull::acc_corr_n = 0.0;
long CtcDFBHull::n_ctr_total = 0;
/* Acumuladores de la sonda `DFBH_ABCONG`, a nivel de archivo para que el
 * destructor pueda emitir el resumen aunque la instancia tenga pocos nodos. */
/* Reparto del tiempo del contractor (`DFBH_PERFIL=1`): linealizacion, pasada
 * del simplex y resto. Lo que decide si congelar `A` puede pagar es cuanto pesa
 * la linealizacion, porque congelar NO la ahorra: la rama congelada tambien
 * evalua todas las restricciones en afin para poder corregir las cotas. */
static double pf_lin = 0.0, pf_pas = 0.0, pf_tot = 0.0;
static long   pf_n = 0;
static void pf_resumen() {
    if (!pf_n || pf_tot <= 0.0) return;
    fprintf(stderr, "[perfil] llamadas=%ld  total=%.3fs | linealizacion %.1f%%  "
            "pasada del simplex %.1f%%  resto %.1f%% | us por llamada: lin %.1f  pasada %.1f\n",
            pf_n, pf_tot, 100.0*pf_lin/pf_tot, 100.0*pf_pas/pf_tot,
            100.0*(pf_tot-pf_lin-pf_pas)/pf_tot,
            1e6*pf_lin/pf_n, 1e6*pf_pas/pf_n);
}
static long   ab_n = 0, ab_ganaF = 0, ab_ganaC = 0, ab_ig = 0, ab_vacF = 0, ab_vacC = 0;
/* Sonda `DFBH_AB`: contraccion de `CtcPolytopeHull` contra DFB sobre la MISMA
 * caja, con la misma relajacion, ya que los dos comparten el linearizador. */
static long   ph_n = 0, ph_vac_prod = 0, ph_vac_dfb = 0, ph_gana_prod = 0, ph_gana_dfb = 0;
static double ph_acc_prod = 0.0, ph_acc_dfb = 0.0;
static void ph_resumen() {
    if (!ph_n) return;
    fprintf(stderr, "[ab] nodos=%ld | contraccion media produccion=%.6f DFB=%.6f razon=%.4f | "
            "gana produccion=%ld DFB=%ld iguales=%ld | vacios prod=%ld DFB=%ld razon=%.4f\n",
            ph_n, ph_acc_prod/ph_n, ph_acc_dfb/ph_n,
            ph_acc_prod > 0 ? ph_acc_dfb/ph_acc_prod : 0.0,
            ph_gana_prod, ph_gana_dfb, ph_n-ph_gana_prod-ph_gana_dfb,
            ph_vac_prod, ph_vac_dfb, ph_vac_prod ? (double)ph_vac_dfb/ph_vac_prod : 0.0);
}
static long   ab_LF = 0, ab_LC = 0, ab_PF = 0, ab_PC = 0;
/* Sonda `DFBH_ABCAJA`: contadores por motor (0 = DFB, 1 = caja negra):
 * LPs, optimos con cota certificada, optimos sin certificado (gamma_k
 * contiene al cero), vacios probados por el certificado, LPs no optimos,
 * cotas saltadas por Achterberg. */
static long abc_cnt[2][6] = {{0,0,0,0,0,0},{0,0,0,0,0,0}};
static long abc_n = 0, abc_vac[2] = {0,0}, abc_gana[2] = {0,0};
static double abc_contr[2] = {0.0, 0.0};
static double abc_t[2] = {0.0, 0.0};   /* segundos en la pasada, por motor */
static long abc_piv[2] = {0, 0};        /* pivotes (iteraciones de SoPlex en la caja) */
static void abc_resumen() {
    if (!abc_n) return;
    for (int e = 0; e < 2; ++e)
        fprintf(stderr, "[abcaja] %s: LPs=%ld | optimo certificado=%ld  optimo sin certificado=%ld  vacio probado=%ld  no optimo=%ld | saltadas=%ld\n",
                e ? "caja " : "DFB  ", abc_cnt[e][0], abc_cnt[e][1], abc_cnt[e][2], abc_cnt[e][3], abc_cnt[e][4], abc_cnt[e][5]);
    fprintf(stderr, "[abcaja] nodos=%ld | contraccion media DFB=%.5f caja=%.5f | nodos vaciados DFB=%ld caja=%ld | gana DFB=%ld caja=%ld\n",
            abc_n, abc_contr[0]/abc_n, abc_contr[1]/abc_n, abc_vac[0], abc_vac[1], abc_gana[0], abc_gana[1]);
    for (int e = 0; e < 2; ++e)
        fprintf(stderr, "[abcaja] costo %s: pasada %.3f s  pivotes=%ld  us/LP=%.2f  piv/LP=%.2f  us/piv=%.3f\n",
                e ? "caja " : "DFB  ", abc_t[e], abc_piv[e], 1e6*abc_t[e]/std::max(1L, abc_cnt[e][0]),
                (double)abc_piv[e]/std::max(1L, abc_cnt[e][0]), 1e6*abc_t[e]/std::max(1L, abc_piv[e]));
}

static double ab_accF = 0.0, ab_accC = 0.0;

/* Sonda sombra: pivotes y tiempo de SC y de la caja negra en el mismo LP. */
static long sm_n = 0, sm_menos = 0, sm_igual = 0, sm_mas = 0, sm_piv[2] = {0,0}, sm_opt[2] = {0,0};
static double sm_t[2] = {0.0, 0.0};
static long sm_hist[2][8] = {{0}};   /* pivotes por LP: 0,1,2,3,4-5,6-9,10-19,20+ */
static long sm_cero_bb[8] = {0};   /* pivotes de SC en los LPs que la caja resuelve sin pivotes */
static int sm_bin(int p) { return p <= 3 ? p : p <= 5 ? 4 : p <= 9 ? 5 : p <= 19 ? 6 : 7; }
static void sombra_anotar(int it_sc, int it_bb, double t_sc, double t_bb, bool o_sc, bool o_bb) {
    ++sm_n; sm_piv[0] += it_sc; sm_piv[1] += it_bb; sm_t[0] += t_sc; sm_t[1] += t_bb;
    sm_opt[0] += o_sc; sm_opt[1] += o_bb;
    if (it_bb < it_sc) ++sm_menos; else if (it_bb == it_sc) ++sm_igual; else ++sm_mas;
    ++sm_hist[0][sm_bin(it_sc)]; ++sm_hist[1][sm_bin(it_bb)];
    if (it_bb == 0) ++sm_cero_bb[sm_bin(it_sc)];
}
static long sz_lps = 0, sz_piv = 0, sz_prim = 0, sz_cero = 0;
/* infactibilidades de partida de SC por clase de LP: [0] SC>10 y caja<=1, [1] el resto */
static long si_n[2] = {0,0}, si_p[2] = {0,0}, si_d[2] = {0,0}, si_piv[2] = {0,0};
static long sb_n[2] = {0,0}, sb_ksc[2] = {0,0}, sb_xsc[2] = {0,0}, sb_kbb[2] = {0,0}, sb_xbb[2] = {0,0}, sb_nbb[2] = {0,0};
static void sombra_base(int it_sc, int it_bb, int ksc, int xsc, int kbb, int xbb) {
    const int c = (it_sc > 10 && it_bb <= 1) ? 0 : 1;
    ++sb_n[c]; sb_ksc[c] += ksc; sb_xsc[c] += xsc;
    if (kbb >= 0) { ++sb_nbb[c]; sb_kbb[c] += kbb; sb_xbb[c] += xbb; }
}
static long sf_n = 0, sf_b = 0, sf_c = 0;
static void sombra_fin(int b, int c) { ++sf_n; sf_b += b; sf_c += c; }
static void sombra_inicio(int it_sc, int it_bb, int pinf, int dinf) {
    const int c = (it_sc > 10 && it_bb <= 1) ? 0 : 1;
    ++si_n[c]; si_p[c] += pinf; si_d[c] += dinf; si_piv[c] += it_sc;
}
static void sombra_cero(int it_sc, long prim, long cero) { ++sz_lps; sz_piv += it_sc; sz_prim += prim; sz_cero += cero; }
static void sombra_resumen() {
    if (DFBSimplex::n_inst_ok + DFBSimplex::n_inst_fallo)
        fprintf(stderr, "[sombra] bases de SoPlex instaladas en SC: %ld, fallidas %ld\n", DFBSimplex::n_inst_ok, DFBSimplex::n_inst_fallo);
    if (sf_n) fprintf(stderr, "[sombra] al terminar cada LP, SC deja %.2f estructurales basicas, %.2f de ellas en una cota (%.0f%%)\n",
                      (double)sf_b/sf_n, (double)sf_c/sf_n, sf_b ? 100.0*sf_c/sf_b : 0.0);
    for (int c = 0; c < 2; ++c) if (sb_n[c] && sb_nbb[c])
        fprintf(stderr, "[sombra] %s: k basica al empezar SC %.1f%% caja %.1f%% | estructurales basicas SC %.2f caja %.2f\n",
                c == 0 ? "SC>10 y caja<=1" : "resto          ", 100.0*sb_ksc[c]/sb_n[c], 100.0*sb_kbb[c]/sb_nbb[c], (double)sb_xsc[c]/sb_n[c], (double)sb_xbb[c]/sb_nbb[c]);
    for (int c = 0; c < 2; ++c) if (si_n[c])
        fprintf(stderr, "[sombra] %s: LPs %ld | al empezar SC tiene %.2f basicas fuera de cota y %.2f costos reducidos de signo malo | piv/LP SC %.2f\n",
                c == 0 ? "SC>10 y caja<=1" : "resto          ", si_n[c], (double)si_p[c]/si_n[c], (double)si_d[c]/si_n[c], (double)si_piv[c]/si_n[c]);
    if (sz_lps) fprintf(stderr, "[sombra] LPs con caja=0 y SC>0: %ld, pivotes SC %ld (primales %ld, de paso cero %ld)\n", sz_lps, sz_piv, sz_prim, sz_cero);
    if (!sm_n) return;
    fprintf(stderr, "[sombra] LPs=%ld | piv/LP SC %.2f caja %.2f | us/LP SC %.2f caja %.2f | caja pivotea menos %.1f%% igual %.1f%% mas %.1f%% | optimos SC %ld caja %ld\n",
            sm_n, (double)sm_piv[0]/sm_n, (double)sm_piv[1]/sm_n, 1e6*sm_t[0]/sm_n, 1e6*sm_t[1]/sm_n,
            100.0*sm_menos/sm_n, 100.0*sm_igual/sm_n, 100.0*sm_mas/sm_n, sm_opt[0], sm_opt[1]);
    const char* et[8] = {"0","1","2","3","4-5","6-9","10-19","20+"};
    for (int e = 0; e < 2; ++e) {
        fprintf(stderr, "[sombra] %s pivotes por LP:", e ? "caja" : "SC  ");
        for (int b = 0; b < 8; ++b) fprintf(stderr, " %s:%.1f%%", et[b], 100.0*sm_hist[e][b]/sm_n);
        fprintf(stderr, "\n");
    }
    long nz = 0; for (int b = 0; b < 8; ++b) nz += sm_cero_bb[b];
    if (nz) {
        fprintf(stderr, "[sombra] en los %ld LPs sin pivotes de la caja, SC hace:", nz);
        for (int b = 0; b < 8; ++b) fprintf(stderr, " %s:%.1f%%", et[b], 100.0*sm_cero_bb[b]/nz);
        fprintf(stderr, "\n");
    }
}
long CtcDFBHull::n_nodes = 0;

/* Reloj de la caja negra (`DFBH_CAJA_RELOJ=1`, 2026-10-16): reparte el tiempo de
 * `resolver_caja` en cargar cotas y objetivo, `minimize` (dentro, el `solve`
 * de SoPlex y la copia de la solucion, si ademas `IBEX_PH_STATS=1`) y leer
 * multiplicadores. `DFBH_CAJA_COTAS=cambio` carga solo las cotas que cambiaron
 * desde el LP anterior de la misma relajacion, en vez de las n cada vez. */
extern double ph_t_solve, ph_t_extrae;
static double cr_t[5] = {0.0, 0.0, 0.0, 0.0, 0.0};   /* cotas+obj, minimize, lectura, solve, copia */
static long cr_n = 0, cr_cotas = 0;
static void cr_resumen() {
    if (!cr_n) return;
    fprintf(stderr, "[cajareloj] LPs=%ld  us/LP: cotas+obj %.2f  minimize %.2f (solve %.2f, copia %.2f)  lectura %.2f  | cotas cargadas por LP %.2f\n",
            cr_n, 1e6*cr_t[0]/cr_n, 1e6*cr_t[1]/cr_n, 1e6*cr_t[3]/cr_n, 1e6*cr_t[4]/cr_n, 1e6*cr_t[2]/cr_n, (double)cr_cotas/cr_n);
}
long CtcDFBHull::n_goal_lp = 0;
long CtcDFBHull::n_sens_extra = 0;
long CtcDFBHull::n_filas_fusionadas = 0;
long CtcDFBHull::n_filas_paralelas = 0;
long CtcDFBHull::n_filas_quitadas = 0;
long CtcDFBHull::n_inf_falso = 0;
long CtcDFBHull::n_bounds = 0;
long CtcDFBHull::n_empty = 0;
long CtcDFBHull::n_cortes = 0;
long CtcDFBHull::n_gtest = 0;
long CtcDFBHull::n_gvacio = 0;
long CtcDFBHull::n_gfuera = 0;
static double t_gtest = 0.0;

void CtcDFBHull::reset_counters() {
    n_lps = n_pivots = n_nodes = n_bounds = n_empty = n_cortes = 0;
    n_congela = n_sinpar = n_acum_gana = 0;
    n_barata = n_ctr_evaluadas = n_ctr_total = 0;
    n_m_tot = n_m_llam = 0;
    n_diag = n_diag_vieja = n_diag_nueva = n_diag_igual = 0;
    acc_diag = acc_corr_v = acc_corr_n = 0.0;
    n_gtest = n_gvacio = n_gfuera = 0;
}

namespace {

double ahora_s() {
    struct timespec ts;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
    return ts.tv_sec + 1e-9*ts.tv_nsec;
}

/* ===================== SONDA DEL DUELO POR COTA =====================
 *
 * `ex7_3_4` gasta x2.58 celdas con un costo por nodo igual de bueno: la perdida
 * es PURA poda. Esta sonda mide de donde sale, resolviendo el MISMO LP —misma
 * relajacion `Af`, misma caja `z` del momento— con SoPlex certificado, y
 * anotando cota por cota lo que cada uno logra.
 *
 * Lo que se registra por cota, una linea por LP:
 *   nodo cota k lado  anchoz  st_dfb piv  st_spx it
 *   lo_z hi_z            dominio antes
 *   nb_dfb               cota certificada por gamma = lambda^T Abar
 *   nb_spx               cota certificada por Neumaier-Shcherbina sobre SoPlex
 *   bf_dfb               optimo FLOTANTE de DFB (sin certificar)
 *   red_dfb red_spx      reduccion relativa del diametro de z_k
 *   vac_dfb vac_spx      cada uno probo vacio?
 *
 * Con eso se separan las cuatro causas posibles: DFB no llega al optimo,
 * el salteo lo da por hecho, el lambda es malo aunque el flotante sea igual,
 * o la certificacion en intervalos pierde mas que Neumaier-Shcherbina. */
struct Duelo {
    FILE* f;
    long cada;
    Duelo() : f(NULL), cada(0) {
        const char* c = getenv("DFBH_DUELO");
        if (!c) return;
        cada = atol(c);
        if (cada <= 0) { cada = 0; return; }
        const char* o = getenv("DFBH_DUELO_OUT");
        f = o ? fopen(o, "w") : stderr;
        if (!f) f = stderr;
        fprintf(f, "nodo,k,lado,st_dfb,piv_dfb,st_spx,it_spx,"
                   "lo_z,hi_z,nb_lo,nb_hi,spx_lo,spx_hi,bf_dfb,spx_f,gk_lo,gk_hi,"
                   "red_dfb,red_spx,vac_dfb,vac_spx,anchoz,t_spx_us,m\n");
    }
};

Duelo& duelo() { static Duelo d; return d; }

/* Reduccion relativa del diametro de z_k al aplicar la cota `nb`. 1 = vacio. */
double reduccion(const Interval& zk, const Interval& nb) {
    if (nb.is_empty() || !zk.intersects(nb)) return 1.0;
    const Interval t = zk & nb;
    const double d0 = zk.diam();
    if (!(d0 > 0.0) || d0 == POS_INFINITY)
        return (t.diam() < d0) ? 1.0 : 0.0;
    return 1.0 - t.diam()/d0;
}

/* Cota de z_k desde gamma, en INTERVALOS. `gamma.z = 0` vale para todo z
 * factible, asi que si gamma_k no contiene al cero
 *     z_k in (-sum_{i!=k} gamma_i*[z_i]) / gamma_k
 * El cociente de intervalos cubre los dos signos de gamma_k. */
Interval cota_desde_gamma(const IntervalVector& gamma, const IntervalVector& z,
                          int k) {
    if (gamma[k].contains(0.0)) return Interval::ALL_REALS;
    /* `z` puede ser mas largo que `gamma` cuando se trabaja con el sistema
     * TRUNCADO (ver m_uso): el producto se hace sobre las `na` componentes de
     * `gamma`, que son las que el sistema en uso tiene. */
    const int na = gamma.size();
    const Interval gk = gamma[k];
    Interval acc(0.0);
    for (int c = 0; c < na; ++c) {
        if (c == k) continue;
        acc += gamma[c] * z[c];
    }
    return (-acc) / gk;
}

/* Dos sensibilidades posibles para el bisector (`DFBH_SENS=redcost|abs`):
 *  redcost: |gamma_j / gamma_goal|, el costo reducido normalizado. Es cero en
 *           toda variable BASICA, asi que ordena solo las no basicas: medido,
 *           hace explotar el arbol.
 *  abs:     sum_i |lambda_i| |A_ij| / |gamma_goal|, la smear ponderada por los
 *           multiplicadores, que es lo que LSmear calcula en los hechos. */
static bool sens_abs() { static const bool v = getenv("DFBH_SENS") && std::string(getenv("DFBH_SENS")) == "abs"; return v; }

/* Multiplicadores por RESTRICCION del sistema (suma con signo de los de sus
 * filas), normalizados por |gamma_goal|; para proyectarlos sobre la jacobiana
 * que usa LSmear (`DFBH_SENS=jac`). Requiere el linearizador afin. */
static void llenar_lam_ctr(std::vector<double>& out, const std::vector<double>& lam, Linearizer& lr,
                           const IntervalVector& gamma, int gv) {
    out.clear();
    LinearizerAffine2* la = dynamic_cast<LinearizerAffine2*>(&lr);
    if (!la || la->fila_ctr.size() != lam.size()) return;
    const double gk = std::fabs(gamma[gv].mid());
    out.assign(la->sistema().nb_ctr, 0.0);
    for (size_t r = 0; r < lam.size(); ++r) out[la->fila_ctr[r]] += lam[r] / gk;
}

static void llenar_sens(std::vector<double>& out, const std::vector<double>& lam, const IntervalMatrix& A,
                        const IntervalVector& gamma, int nb_var, int gv) {
    const double gk = std::fabs(gamma[gv].mid());
    out.assign(nb_var, 0.0);
    if (!sens_abs()) {
        for (int j = 0; j < nb_var; ++j) if (j != gv) out[j] = std::fabs(gamma[j].mid()) / gk;
        return;
    }
    const int m = A.nb_rows();
    for (int i = 0; i < m; ++i) {
        if (lam[i] == 0.0) continue;
        const double li = std::fabs(lam[i]);
        for (int j = 0; j < nb_var; ++j) if (j != gv) out[j] += li * A[i][j].mag();
    }
    for (int j = 0; j < nb_var; ++j) out[j] /= gk;
}


/* SONDA DEL TESTIGO (`DFBH_TESTIGO=archivo`): un punto que se sabe solucion.
 * Si DFB lo deja fuera de la caja, algo no es riguroso. */
static const std::vector<double>& testigo() {
    static std::vector<double> t; static bool cargado = false;
    if (!cargado) {
        cargado = true;
        const char* fn = getenv("DFBH_TESTIGO");
        if (fn) { FILE* fp = fopen(fn, "r"); double v; if (fp) { while (fscanf(fp, "%lf", &v) == 1) t.push_back(v); fclose(fp); } }
    }
    return t;
}
static bool testigo_en(const IntervalVector& z, int nb_var) {
    const std::vector<double>& t = testigo();
    if ((int)t.size() != nb_var) return false;
    for (int i = 0; i < nb_var; ++i) if (!z[i].contains(t[i])) return false;
    return true;
}

/* Contadores de la sonda del corte del objetivo (`DFBH_GOALPROBE=1`). Se
 * registran tambien en la rama que vacia la caja: antes se contaba despues de
 * ella, asi que la sonda nunca veia un nodo muerto por la cota del objetivo. */
struct GoalRep {
    static long n, muere, cruza, pv_tot, pv_cruce, pv_muere, ahorro, paso_tot, restantes, pv_total_run;
    static void anota(bool vacia, int iters, int obs, int paso, int restan) {
        ++n; pv_tot += iters; paso_tot += paso;
        if (vacia) {
            ++muere; pv_muere += iters;
            if (obs >= 0) { ++cruza; pv_cruce += obs; ahorro += iters - obs; restantes += restan; }
        }
    }
    static void volcar(long piv_run) {
        if (!n) return;
        fprintf(stderr, "[goalprobe] cotas del objetivo=%ld | matan el nodo %.1f%% | de esas, el flotante cruza antes "
                "del final en %.1f%%, en el pivote %.2f de %.2f | pivotes ahorrables=%ld (%.2f%% de todos los del "
                "contractor) | cotas de la pasada que no harian falta %.1f por nodo muerto\n",
                n, 100.0*muere/n, muere ? 100.0*cruza/muere : 0.0,
                cruza ? (double)pv_cruce/cruza : 0.0, muere ? (double)pv_muere/muere : 0.0,
                ahorro, piv_run ? 100.0*ahorro/piv_run : 0.0, cruza ? (double)restantes/cruza : 0.0);
    }
};
long GoalRep::n = 0, GoalRep::muere = 0, GoalRep::cruza = 0, GoalRep::pv_tot = 0, GoalRep::pv_cruce = 0,
     GoalRep::pv_muere = 0, GoalRep::ahorro = 0, GoalRep::paso_tot = 0, GoalRep::restantes = 0, GoalRep::pv_total_run = 0;

/* Contadores de la sonda del lagrangiano (`DFBH_LAGR=1`), ver `una_pasada`. */
struct LagrRep {
    static long n, n_ctrl_mal, n_bp_mejor, n_fix, n_c_mejor, n_c_mejor_a, n_mata, n_mata_bp, n_mata_c_solo, n_gap;
    static long n_fallo_b, n_fallo_bp, n_fallo_grad, n_fallo_c;
    static double max_ctrl, s_bp, s_cbp, s_ca, s_fix;
    static void volcar() {
        if (!n) return;
        fprintf(stderr, "[lagr] LPs del objetivo=%ld (fallos: b=%ld b'=%ld grad=%ld c=%ld) | control |b-a| rel max %.1e, >1e-7 en %ld | "
                "b' (rho en la caja contraida) mejora a en %ld (%.1f%%), media acotada (b'-a)/brecha %.3f | "
                "monotonia: fija alguna variable en %ld (%.1f%%), %.2f por LP; c mejora b' en %ld (%.1f%%), "
                "media acotada (c-b')/brecha %.3f; c mejora a en %ld, media acotada (c-a)/brecha %.3f | "
                "con brecha finita %ld: matan el nodo y a no -> b' %ld, c %ld, de los cuales solo c (b' no) %ld\n",
                n, n_fallo_b, n_fallo_bp, n_fallo_grad, n_fallo_c, max_ctrl, n_ctrl_mal,
                n_bp_mejor, 100.0*n_bp_mejor/n, n_gap ? s_bp/n_gap : 0.0,
                n_fix, 100.0*n_fix/n, s_fix/n, n_c_mejor, 100.0*n_c_mejor/n, n_gap ? s_cbp/n_gap : 0.0,
                n_c_mejor_a, n_gap ? s_ca/n_gap : 0.0,
                n_gap, n_mata_bp, n_mata, n_mata_c_solo);
    }
};
long LagrRep::n = 0, LagrRep::n_ctrl_mal = 0, LagrRep::n_bp_mejor = 0, LagrRep::n_fix = 0,
     LagrRep::n_c_mejor = 0, LagrRep::n_c_mejor_a = 0, LagrRep::n_mata = 0, LagrRep::n_mata_bp = 0,
     LagrRep::n_mata_c_solo = 0, LagrRep::n_gap = 0,
     LagrRep::n_fallo_b = 0, LagrRep::n_fallo_bp = 0, LagrRep::n_fallo_grad = 0, LagrRep::n_fallo_c = 0;
double LagrRep::max_ctrl = 0.0, LagrRep::s_bp = 0.0, LagrRep::s_cbp = 0.0, LagrRep::s_ca = 0.0, LagrRep::s_fix = 0.0;

/** \brief Engancha/desengancha la traza de `lambda` por pivote. Solo sonda. */
void sk_traza(DFBSimplex& sk, std::vector<std::vector<double> >* t) {
    sk.traza_y = t;
}

} /* anonymous */

CtcDFBHull* CtcDFBHull::ultimo = NULL;

CtcDFBHull::CtcDFBHull(Linearizer& lr, int goal_var, double eps, Ctc* sonda_ctc) :
    Ctc(lr.nb_var()), lr(lr), nb_var(lr.nb_var()), goal_var(goal_var),
    sonda_ctc(sonda_ctc), caja_lin(1), caja_sens(1), sens_ok(false), sens_contrajo(false), calcular_sens(false), sens_pase(false), jac_cortes_ok(false), m_afin(0), lin_valida(false), g_apagado(false), sin_recarga(false),
    mylinearsolver(lr.nb_var(), LPSolver::Mode::NotCertified, eps),
    A(1,1), z(1), Af(1,1), gamma(1), hull_primal(1), hull_vacio(true), caja_acum(1), ph_ab(NULL) {
    m_uso = 0;
    hc4_actuo = false;
    {   const char* pr = getenv("DFBH_PRIMERA");
        regla_primera = (pr && std::string(pr)=="ancha")   ? 1 :
                        (pr && std::string(pr)=="angosta") ? 2 :
                        (pr && std::string(pr)=="azar")    ? 3 : 0;
    }
    traza_cotas = false;
    nodo_nuevo = true;
    traza_resumen = (getenv("DFBH_SINCOTA") != NULL);
    if (DFB_ENV("DFBH_GVAR"))
        fprintf(stderr, "[dfb] nb_var=%d  goal_var=%d\n", nb_var, goal_var);
    ultimo = this;
    if (DFB_ENV("DFBH_AB")) ph_ab = new CtcPolytopeHull(lr);
}

CtcDFBHull::~CtcDFBHull() {
    ibex::perf_volcar();
    if (DFB_ENV("DFBH_ENACID"))
        fprintf(stderr, "[gtest] llamadas=%ld  rechazadas=%ld (%.1f%%)  "
                "vacios=%ld (%.2f%% de las aplicables)%s\n",
                n_gtest, n_gfuera, n_gtest ? 100.0*n_gfuera/n_gtest : 0.0,
                n_gvacio, (n_gtest-n_gfuera) ? 100.0*n_gvacio/(n_gtest-n_gfuera) : 0.0,
                g_apagado ? "  [AUTOAPAGADO]" : "");
    /* Que estrategia corrio, para que cualquier medicion lo diga sin que haya
     * que reconstruirlo del entorno. */
    fprintf(stderr, "[dfb] pivoteo=%s%s\n",
            !DFBSimplex::pivoteo_dual ? "primal (dual en la cota 0, primal de la 1 en adelante)"
            : DFBSimplex::dual_puro ? "dual (puro: reubica todo, k incluida)"
            : "mixta (reubica todo menos k)",
            DFBSimplex::n_entrak ? "  +entrada explicita de k" : "");
    fprintf(stderr, "[dfb] nodos=%ld  vacios probados por DFB=%ld (%.1f%% de los nodos)  "
            "LPs=%ld  pivotes=%ld\n", n_nodes, n_empty,
            n_nodes ? 100.0*n_empty/n_nodes : 0.0, n_lps, n_pivots);
    fprintf(stderr, "[dfb] piv/LP=%.2f  degen=%.1f%%  degen-enmascarados=%.1f%%  "
            "LPs perturbados=%ld  pivotes de limpieza=%ld\n",
            n_lps ? (double)n_pivots/n_lps : 0.0,
            n_pivots ? 100.0*DFBSimplex::n_piv_degen/n_pivots : 0.0,
            n_pivots ? 100.0*DFBSimplex::n_piv_degen_pert/n_pivots : 0.0,
            DFBSimplex::n_pert_lps, DFBSimplex::n_pert_limpieza);
    if (n_orig)
        fprintf(stderr,
          "[corig] cotas %ld | refrescada mas angosta %.2f%%  igual %.2f%%  mas ancha %.2f%%"
          " | ancho refrescada/original x%.4f | vacia solo la refrescada %ld\n",
          n_orig, 100.0*n_orig_mejor/n_orig, 100.0*n_orig_igual/n_orig,
          100.0*n_orig_peor/n_orig,
          (n_orig-n_orig_vacio) ? std::exp(acc_orig/(n_orig-n_orig_vacio)) : 0.0,
          n_orig_vacio);
    if (n_fila_tot)
        fprintf(stderr,
          "[refc] filas-llamada %ld | sin refresco de c: %ld (%.2f%%)"
          " | ancho corr/triv medio %.4f | no aprieta nada: %.2f%%\n",
          n_fila_tot, n_fila_vieja, 100.0*n_fila_vieja/n_fila_tot,
          n_apr ? acc_apr/n_apr : 0.0, n_apr ? 100.0*n_apr_nada/n_apr : 0.0);
    if (n_chk)
        fprintf(stderr,
          "[chkc] muestras %ld | identidad A.x-a'.x in dif: %ld violaciones (peor %.3g)"
          " | puntos que cumplen rango': %ld, de ellos fuera de corr: %ld (peor %.3g)\n",
          n_chk, n_chk_ident, chk_peor_ident, n_chk_cond, n_chk_viola, chk_peor);
    if (n_m_tot)
        fprintf(stderr, "[dfb] filas de la relajacion: %.2f por linealizacion (%ld llamadas)\n",
                (double)n_m_tot/n_m_llam, n_m_llam);
        if (n_filas_paralelas || n_filas_quitadas)
            fprintf(stderr, "[dfb] filas paralelas: fundidas=%ld quitadas=%ld (%.2f y %.2f por linealizacion)\n", n_filas_paralelas, n_filas_quitadas, (double)n_filas_paralelas/n_m_llam, (double)n_filas_quitadas/n_m_llam);
    if (DFB_ENV("DFBH_AFALLOC")) {
        extern long af_n_alloc, af_n_double;
        fprintf(stderr, "[afalloc] asignaciones de la afin=%ld  doubles=%ld  (%.1f por linealizacion)\n",
                af_n_alloc, af_n_double, n_nodes ? (double)af_n_alloc/n_nodes : 0.0);
    }
    if (DFB_ENV("DFBH_LINPERF")) lin_perf_volcar();
    pf_resumen();
    ph_resumen();
    abc_resumen();
    cr_resumen();
    sombra_resumen();
    if (ab_n)
        fprintf(stderr, "[abcong] nodos=%ld | contraccion media fresca=%.6f congelada=%.6f "
                "razon=%.4f | gana fresca=%ld congelada=%ld iguales=%ld | vacios F=%ld C=%ld | "
                "LPs/nodo F=%.2f C=%.2f | pivotes/nodo F=%.2f C=%.2f\n",
                ab_n, ab_accF/ab_n, ab_accC/ab_n,
                ab_accF > 0 ? ab_accC/ab_accF : 0.0,
                ab_ganaF, ab_ganaC, ab_ig, ab_vacF, ab_vacC,
                (double)ab_LF/ab_n, (double)ab_LC/ab_n,
                (double)ab_PF/ab_n, (double)ab_PC/ab_n);
    if (n_congela)
        fprintf(stderr, "[dfb] matriz congelada: congelaciones=%ld de %ld llamadas (%.1f%%)  filas sin par=%ld\n",
                n_congela, n_nodes, n_nodes ? 100.0*n_congela/n_nodes : 0.0, n_sinpar);
    if (n_congela)
        fprintf(stderr, "[dfb] cotas de b que la acumulacion aprieta=%ld  vueltas sin linealizar=%ld\n", n_acum_gana, n_barata);
    if (n_diag)
        fprintf(stderr, "[congdiag] cotas comparadas=%ld | gana la VIEJA %.1f%%  la NUEVA %.1f%%  iguales %.1f%% | ancho vieja/nueva x%.4f | termino de correccion relativo: vieja %.4f  nueva %.4f\n",
                n_diag, 100.0*n_diag_vieja/n_diag, 100.0*n_diag_nueva/n_diag,
                100.0*n_diag_igual/n_diag, std::exp(acc_diag/n_diag),
                acc_corr_v/n_diag, acc_corr_n/n_diag);
    if (n_ctr_total)
        fprintf(stderr, "[dfb] linealizacion selectiva: %ld de %ld restricciones evaluadas (%.1f%%)\n",
                n_ctr_evaluadas, n_ctr_total, 100.0*n_ctr_evaluadas/n_ctr_total);
    if (DFBSimplex::n_base_heredada || DFBSimplex::n_base_rechazada)
        fprintf(stderr, "[dfb] bases heredadas entre nodos=%ld  rechazadas por tableau=%ld\n",
                DFBSimplex::n_base_heredada, DFBSimplex::n_base_rechazada);
    if (DFBSimplex::n_hereda)
        fprintf(stderr, "[dfb] cargas que heredaron la ubicacion de las no basicas=%ld\n",
                DFBSimplex::n_hereda);
    if (DFBSimplex::n_anytime)
        fprintf(stderr, "[dfb] cortes anytime por brecha=%ld (%.1f%% de los LPs)\n",
                DFBSimplex::n_anytime, n_lps ? 100.0*DFBSimplex::n_anytime/n_lps : 0.0);
    fprintf(stderr, "[dfb] vueltas del bucle=%ld  (%.2f por resolucion, contra %.2f pivotes)\n",
            DFBSimplex::n_vueltas,
            DFBSimplex::n_solves ? (double)DFBSimplex::n_vueltas/DFBSimplex::n_solves : 0.0,
            DFBSimplex::n_solves ? (double)n_pivots/DFBSimplex::n_solves : 0.0);
    fprintf(stderr, "[dfb] compute_basics=%ld  (%.2f por pivote)\n",
            DFBSimplex::n_basics, n_pivots ? (double)DFBSimplex::n_basics/n_pivots : 0.0);
    fprintf(stderr, "[dfb] estados: solves=%ld  optimal=%.1f%%  infeasible=%.1f%%  iterlimit=%.1f%%  "
            "otros(singular/unbounded)=%.1f%%  | warm=%ld cold=%ld fact_fallo=%ld\n",
            DFBSimplex::n_solves,
            DFBSimplex::n_solves ? 100.0*DFBSimplex::n_optimal/DFBSimplex::n_solves : 0.0,
            DFBSimplex::n_solves ? 100.0*DFBSimplex::n_infeasible/DFBSimplex::n_solves : 0.0,
            DFBSimplex::n_solves ? 100.0*DFBSimplex::n_iterlimit/DFBSimplex::n_solves : 0.0,
            DFBSimplex::n_solves ? 100.0*(DFBSimplex::n_solves-DFBSimplex::n_optimal-DFBSimplex::n_infeasible-DFBSimplex::n_iterlimit)/DFBSimplex::n_solves : 0.0,
            DFBSimplex::n_warm, DFBSimplex::n_cold, DFBSimplex::n_fact_fallo);
    LagrRep::volcar();
    GoalRep::volcar(n_pivots);
    if (DFB_ENV("DFBH_INFSTATS"))
        fprintf(stderr, "[infstats] filas donde el infinito de SoPlex (1e100) habria cortado=%ld\n", n_inf_falso);
    fprintf(stderr, "[dfb] LPs de la cota inferior del objetivo con certificado=%ld (%.2f por nodo); LPs extra para sensibilidades=%ld\n",
            n_goal_lp, n_nodes ? (double)n_goal_lp/n_nodes : 0.0, n_sens_extra);
    if (DFBSimplex::n_nodos_rev)
        fprintf(stderr, "[dfb] motor por nodo: revisado=%ld tableau=%ld\n", DFBSimplex::n_nodos_rev, DFBSimplex::n_nodos_tab);
        if (n_hib_caja + n_hib_sc > 0)
            fprintf(stderr, "[dfb] hibrido: pasadas con SoPlex directo=%ld con SC=%ld\n", n_hib_caja, n_hib_sc);
        if (warm_medir())
            fprintf(stderr, "[dfb] warm medir: nodos dual=%ld primal=%ld | costo por LP (prom. movil) dual %.1f us, primal %.1f us\n",
                    g_ws_nodos[0], g_ws_nodos[1], 1e6*g_ws_ema[0], 1e6*g_ws_ema[1]);
        if (adapt_umbral() > 0.0)
            fprintf(stderr, "[dfb] adaptativo (umbral %.2f): nodos dual=%ld primal=%ld, cambios de lado por LP (prom. movil)=%.2f\n",
                    adapt_umbral(), adapt_nodos(true), adapt_nodos(false), adapt_ema());
    if (DFB_ENV("DFBH_LUSTATS"))
        fprintf(stderr, "[lustats] factorizaciones=%ld updates=%ld btran=%ld ftran=%ld\n",
                DFBBasis::n_factorizations, DFBBasis::n_updates, DFBBasis::n_btran, DFBBasis::n_ftran);
    if (n_lagr_aplic)
        fprintf(stderr, "[dfb] cota lagrangiana: aplicada en %ld cotas, mejoro la del LP en %ld (%.1f%%), "
                "vacios que solo ella prueba=%ld, variables fijadas por cota %.2f\n",
                n_lagr_aplic, n_lagr_mejora, 100.0*n_lagr_mejora/n_lagr_aplic, n_lagr_vacio,
                (double)n_lagr_fijas/n_lagr_aplic);
    if (DFBSimplex::sonda_deriva && DFBSimplex::dv_refrescos)
        fprintf(stderr, "[deriva] refrescos=%ld, con algun veredicto distinto=%ld (%.2f%%) | d: %ld signos distintos de %ld no basicas, "
                "max|dif| %.2e | x_B: %ld factibilidades distintas de %ld basicas, max dif rel %.2e\n",
                DFBSimplex::dv_refrescos, DFBSimplex::dv_refr_dist,
                100.0*DFBSimplex::dv_refr_dist/DFBSimplex::dv_refrescos,
                DFBSimplex::dv_d_signo, DFBSimplex::dv_d_nobas, DFBSimplex::dv_d_max,
                DFBSimplex::dv_xb_fact, DFBSimplex::dv_xb_bas, DFBSimplex::dv_xb_max);
    if (DFBSimplex::lu_tableau)
        fprintf(stderr, "[dfb] LU tableau: refrescos por tableau=%ld, guardia fallo en d=%ld x_B=%ld, lambda final=%ld\n",
                DFBSimplex::n_guardia_ok, DFBSimplex::n_guardia_fallo_d,
                DFBSimplex::n_guardia_fallo_xb, DFBSimplex::n_guardia_fallo_fin);
    if (DFBSimplex::lu_diferida)
        fprintf(stderr, "[dfb] LU diferida: refactorizaciones=%ld (%.2f por resolucion)\n",
                DFBSimplex::n_refact_dif,
                DFBSimplex::n_solves ? (double)DFBSimplex::n_refact_dif/DFBSimplex::n_solves : 0.0);
    if (DFBSimplex::n_verif)
        fprintf(stderr, "[dfb] verif: optimos verificados=%ld  primal-infactibles=%.1f%%  dual-infactibles=%.1f%%\n",
                DFBSimplex::n_verif, 100.0*DFBSimplex::n_verif_pinf/DFBSimplex::n_verif,
                100.0*DFBSimplex::n_verif_dinf/DFBSimplex::n_verif);
    delete ph_ab;
}

int CtcDFBHull::linearize(const IntervalVector& box, ContractContext& context) {
    static const int cong = getenv("DFBH_CONG") ? atoi(getenv("DFBH_CONG")) : 0;
    if (cong > 0 && dynamic_cast<LinearizerAffine2*>(&lr) != NULL)
        return linearize_congelada(box, cong);
    return linearize_fresca(box, context);
}

int CtcDFBHull::linearize_fresca(const IntervalVector& box, ContractContext& context) {
    reusa_A = false;
    mylinearsolver.clear_constraints();
    if (box.is_unbounded()) return 0;

    const int m0 = lr.linearize(box, mylinearsolver, context.prop);
    if (m0 <= 0) return m0;

    const Matrix rows0 = mylinearsolver.rows();
    const IntervalVector lhs_rhs0 = mylinearsolver.lhs_rhs();

    /* FUSION DE FILAS GEMELAS. `linearize` parte cada ecuacion en dos filas
     * unilateras con la MISMA normal (una por lado). En la forma homogenea
     * basta una fila con rango [lhs, rhs]: el poliedro es el mismo y `m` baja
     * a la mitad en los sistemas cuadrados (Bratu-0150: 300 -> 150), lo que
     * achica la base del motor revisado. Se fusionan solo filas CONSECUTIVAS
     * identicas, que es lo que emite el linearizador. La caja negra recibe la
     * misma relajacion fusionada (`DFBH_FUSIONA=0` lo apaga). */
    static const bool fusiona = !(getenv("DFBH_FUSIONA") && std::string(getenv("DFBH_FUSIONA")) == "0");
    int m = m0;
    std::vector<int> fuente;   /* fila del LPSolver que da los coeficientes de cada fila fusionada */
    IntervalVector lados(m0);
    for (int i = 0; i < m0; ++i) lados[i] = lhs_rhs0[nb_var+i];
    if (fusiona) {
        int k = 0;
        for (int i = 0; i < m0; ++i) {
            if (k > 0) {
                const int p = fuente[k-1];
                bool igual = true;
                for (int j = 0; j < nb_var && igual; ++j) igual = (rows0[nb_var+i][j] == rows0[nb_var+p][j]);
                const Interval& a = lados[k-1]; const Interval& b = lhs_rhs0[nb_var+i];
                const double lo = std::max(a.lb(), b.lb()), hi = std::min(a.ub(), b.ub());
                /* lados cruzados: se dejan las dos filas, que el LP pruebe el vacio */
                if (igual && lo <= hi) {
                    lados[k-1] = Interval(lo, hi);
                    ++n_filas_fusionadas;
                    continue;
                }
            }
            fuente.push_back(i); lados[k] = lhs_rhs0[nb_var+i]; ++k;
        }
        m = k;
    } else {
        for (int i = 0; i < m0; ++i) fuente.push_back(i);
    }

    /* FILAS PARALELAS (`DFBH_FILAS_PARALELAS=paralelas|quita`, 2026-10-20). La relajacion
     * `compo` suma a la fila afin de cada restriccion las de XTaylor, casi
     * paralelas a ella, y el LP se vuelve muy degenerado (Katsurabis-14 con
     * compo: 53 % de pivotes degenerados y 64 piv/LP contra 10 con el afin).
     * `paralelas`: dos filas PROPORCIONALES dentro de `DFBH_FILAS_TOL`
     * (omision 1e-9, relativa) se funden en una: se escalan los lados de la
     * segunda y se intersecan con los de la primera; el poliedro es el mismo.
     * `quita`: ademas, de dos filas con coseno > 1 - `DFBH_FILAS_COS` (omision
     * 1e-6) se deja la mas apretada sobre la caja (ancho relativo de la holgura
     * sobre el rango de la fila); esto si cambia la relajacion. Se compara solo
     * entre filas con el mismo soporte. Vale para los dos motores. */
    /* (Se llamaba `DFBH_FILAS`, que ya era la COSECHA de cotas de otras filas:
     * encenderla encendia tambien esa cosecha, cara y sin relacion.) */
    static const char* filas_env = getenv("DFBH_FILAS_PARALELAS");
    static const bool filas_par = filas_env && (std::string(filas_env) == "paralelas" || std::string(filas_env) == "quita");
    static const bool filas_quita = filas_env && std::string(filas_env) == "quita";
    static const double filas_tol = getenv("DFBH_FILAS_TOL") ? atof(getenv("DFBH_FILAS_TOL")) : 1e-9;
    static const double filas_cos = getenv("DFBH_FILAS_COS") ? atof(getenv("DFBH_FILAS_COS")) : 1e-6;
    if (filas_par && m > 1) {
        std::vector<char> viva(m, 1);
        std::map<std::vector<int>, std::vector<int> > grupos;
        std::vector<double> norma(m, 0.0), apriete(m, 0.0);
        std::vector<std::vector<int> > sop(m);
        for (int i = 0; i < m; ++i) {
            const Vector& r = rows0[nb_var+fuente[i]];
            Interval rango(0.0);
            for (int j = 0; j < nb_var; ++j)
                if (r[j] != 0.0) { sop[i].push_back(j); norma[i] += r[j]*r[j]; rango += Interval(r[j]) * box[j]; }
            norma[i] = std::sqrt(norma[i]);
            const Interval h = lados[i] & rango;
            apriete[i] = (h.is_empty() || rango.diam() <= 0.0) ? 0.0 : h.diam() / rango.diam();   /* 0 = la mas apretada */
        }
        for (int i = 0; i < m; ++i) {
            if (sop[i].empty()) continue;
            std::vector<int>& g = grupos[sop[i]];
            for (size_t t = 0; t < g.size() && viva[i]; ++t) {
                const int q = g[t];
                if (!viva[q]) continue;
                const Vector& a = rows0[nb_var+fuente[q]];
                const Vector& b = rows0[nb_var+fuente[i]];
                const int j0 = sop[i][0];
                const double c = b[j0] / a[j0];
                bool prop = true; double dot = 0.0;
                for (size_t u = 0; u < sop[i].size(); ++u) {
                    const int j = sop[i][u];
                    dot += a[j]*b[j];
                    if (std::fabs(b[j] - c*a[j]) > filas_tol * std::max(std::fabs(b[j]), std::fabs(c*a[j]))) prop = false;
                }
                if (prop) {
                    const Interval inter = lados[q] & (lados[i] / c);
                    if (!inter.is_empty()) { lados[q] = inter; viva[i] = 0; ++n_filas_paralelas; }
                } else if (filas_quita && std::fabs(dot) / (norma[q]*norma[i]) > 1.0 - filas_cos) {
                    if (apriete[i] < apriete[q]) { viva[q] = 0; g[t] = i; } else viva[i] = 0;
                    ++n_filas_quitadas;
                }
            }
            if (viva[i]) g.push_back(i);
        }
        int k = 0;
        for (int i = 0; i < m; ++i)
            if (viva[i]) { fuente[k] = fuente[i]; lados[k] = lados[i]; ++k; }
        m = k;
    }
    /* con `DFBH_HIBRIDO` el motor se decide despues de linealizar: se recarga siempre */
    if (m < m0 && (motor_caja || DFB_ENV("DFBH_ABCAJA") || DFB_ENV("DFBH_SOMBRA") || DFB_ENV("DFBH_HIBRIDO"))) {
        /* La caja negra resuelve sobre `mylinearsolver`: se recarga con las
         * filas fusionadas para que los multiplicadores sigan siendo fila a
         * fila los de `A`. */
        mylinearsolver.clear_constraints();
        for (int i = 0; i < m; ++i) {
            Vector row(nb_var);
            for (int j = 0; j < nb_var; ++j) row[j] = rows0[nb_var+fuente[i]][j];
            mylinearsolver.add_constraint(lados[i].lb(), row, lados[i].ub());
        }
    }

    caja_lin = box;          /* la caja para la que valen esta relajacion y sus gamma */
    lin_valida = false;      /* recien es valida cuando haya gamma calculados */
    n_m_tot += m; ++n_m_llam;
    m_afin = m;

    /* Cortes del pozo de la celda que valen para esta caja (`DFBH_CORTES`). */
    std::vector<const CorteDFB*> vc;
    if (modo_cortes > 0) {
        const BxpDFBCortes* pc = (const BxpDFBCortes*) context.prop[BxpDFBCortes::id()];
        if (pc)
            for (size_t t = 0; t < pc->pozo.size(); ++t)
                if (pc->pozo[t].caja.size() == box.size() && pc->pozo[t].caja.is_superset(box)
                    && (int)pc->pozo[t].a.size() == nb_var)
                    vc.push_back(&pc->pozo[t]);
        if (!vc.empty()) { n_cortes_filas += vc.size(); ++n_cortes_llam; }
    }
    const int mt = m + (int)vc.size();

    A.resize(mt, nb_var + mt);
    A.clear();
    z.resize(nb_var + mt);
    for (int i = 0; i < nb_var; ++i) z[i] = box[i];

    /* b_i = fila_i . x, asi que su rango valido es la evaluacion por intervalos
     * de la fila sobre la caja, intersectada con la cota que reporta el solver.
     * Las dos son validas para toda solucion, de modo que su interseccion
     * tambien; si fuera vacia se cae al enclosure, que siempre es valido, para
     * no propagar vacuidad desde la linealizacion. */
    for (int i = 0; i < m; ++i) {
        Interval enclosure(0.0);
        for (int j = 0; j < nb_var; ++j)
            enclosure += Interval(rows0[nb_var+fuente[i]][j]) * box[j];
        /* SoPlex representa el lado infinito de una fila como +-1e100
         * (`soplex::infinity`). Hay que leerlo como infinito: tomado como cota
         * finita, con coeficientes del orden de 1e113 (`Brown-15`, cajas de
         * 1e8 con productos de muchas variables) corta soluciones reales, y DFB
         * certificaba contracciones y vacios falsos (2026-10-05). */
        Interval ladob = lados[i];
        if ((ladob.lb() <= -1e100 && enclosure.lb() < -1e100) || (ladob.ub() >= 1e100 && enclosure.ub() > 1e100)) ++n_inf_falso;
        if (ladob.lb() <= -1e100) ladob = Interval(NEG_INFINITY, ladob.ub());
        if (ladob.ub() >=  1e100) ladob = Interval(ladob.lb(), POS_INFINITY);
        const Interval bi = ladob & enclosure;
        /* `DFBH_HOLGURA=lados` (2026-10-19, diagnostico): la holgura usa solo el
         * lado del solver, como la caja negra, sin intersectarlo con la evaluacion
         * por intervalos de la fila. Valido igual (el lado vale para toda
         * solucion); el certificado puede quedar mas debil si un lado es infinito. */
        static const bool holgura_lados = getenv("DFBH_HOLGURA") && std::string(getenv("DFBH_HOLGURA")) == "lados";
        z[nb_var+i] = holgura_lados ? ladob : (bi.is_empty() ? enclosure : bi);
    }
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < nb_var; ++j) A[i][j] = rows0[nb_var+fuente[i]][j];
        A[i][nb_var+i] = Interval(-1.0);   /* la parte b de Abar es -I */
    }
    for (int r = 0; r < (int)vc.size(); ++r) {
        const int i = m + r;
        Interval enclosure(0.0);
        for (int j = 0; j < nb_var; ++j) {
            A[i][j] = Interval(vc[r]->a[j]);
            enclosure += Interval(vc[r]->a[j]) * box[j];
        }
        A[i][nb_var+i] = Interval(-1.0);
        const Interval bi = vc[r]->rango & enclosure;
        z[nb_var+i] = bi.is_empty() ? enclosure : bi;
    }

    return mt;
}

/* MATRIZ CONGELADA. Ver `cong_rows` en la cabecera.
 *
 * Se queda con las filas de la linealizacion anterior y recalcula solo el rango
 * de cada `b_i`, usando una linealizacion FRESCA de la misma restriccion:
 *
 *     b_i = A_i.x = (A_i - a'_i).x + a'_i.x   con  a'_i.x en rango'_i
 *     ==> b_i  en  rango'_i + (A_i - a'_i).X
 *
 * La fila diferencia se evalua por intervalos sobre la caja, que para una forma
 * lineal es exacto, asi que la cancelacion de los simbolos de ruido de la afin
 * queda contemplada. Se intersecta ademas con `A_i.X`, siempre valido.
 *
 * Como `A` no cambia, el simplex no recarga: tableau, base y factorizacion
 * sobreviven, y entre llamadas solo se mueven las cotas de las `b`. */
/* SONDA DEL ACOTAMIENTO AFIN DE `c` (`DFBH_CHKC=1`).
 *
 * La cota congelada de `b_i` descansa en dos afirmaciones. La primera,
 * `a'_i.x in rango'_i` para todo `x` factible, es la del linearizador afin y ya
 * quedo verificada aparte: congelando en cada llamada el camino reproduce bit a
 * bit al fresco. La segunda es la nueva, y es la que se prueba aqui:
 *
 *   identidad   `A_i.x - a'_i.x in dif`  para TODO `x` de la caja, sin condicion;
 *   implicacion `a'_i.x in rango'_i  ==>  A_i.x in corr`.
 *
 * Se muestrean puntos de la caja (centro, esquinas al azar, interiores) y se
 * evalua todo en intervalos, de modo que una violacion es un defecto real y no
 * redondeo. La implicacion se prueba contra `corr`, antes de intersectar con
 * `b_acum`: `b_acum` viene de cajas anteriores y un punto muestreado no tiene
 * por que respetarlo. */
/* Rango capturado al congelar la fila `i`. Solo lo usa la sonda `DFBH_CORIG`
 * cuando la acumulacion se corta y hay que volver al punto de partida. */
Interval CtcDFBHull::gn_cong(int i, const std::vector<int>&, const std::vector<Interval>&) {
    return ((int)b_orig_cong.size() == (int)cong_rows.size()) ? b_orig_cong[i]
                                                              : Interval::ALL_REALS;
}

void CtcDFBHull::chequear_c(int i, const Vector& Ai, const Vector& ap,
                            const Interval& rango_p, const IntervalVector& box,
                            const Interval& dif, const Interval& corr) {
    static const int chk = getenv("DFBH_CHKC") ? atoi(getenv("DFBH_CHKC")) : 0;
    if (chk <= 0) return;
    if ((int)Ai.size() != nb_var || (int)ap.size() != nb_var) return;
    for (int t = 0; t < chk; ++t) {
        Interval pA(0.0), pa(0.0);
        std::vector<double> xs(nb_var);
        for (int j = 0; j < nb_var; ++j) {
            double xj;
            if (t == 0) xj = box[j].mid();
            else if (box[j].is_bisectable() && (rand() & 1)) {
                const double u = (double)rand()/(double)RAND_MAX;
                xj = box[j].lb() + u*box[j].diam();
                if (!std::isfinite(xj)) xj = box[j].mid();
            } else xj = (rand() & 1) ? box[j].lb() : box[j].ub();
            if (!std::isfinite(xj)) xj = box[j].mid();
            xs[j] = xj;
            const Interval x(xj);
            pA += Interval(Ai[j]) * x;
            pa += Interval(ap[j])  * x;
        }
        if (pA.is_empty() || pa.is_empty()) continue;
        ++n_chk;
        /* `A.x - a'.x` NO se puede evaluar como `pA - pa`: los dos comparten
         * `x`, y al restarlos como intervalos independientes la cancelacion
         * infla el resultado hasta hacerlo inutil. Se acumula termino a
         * termino, en precision extendida, que es lo que `dif` acota. */
        long double ld = 0.0L;
        for (int j = 0; j < nb_var; ++j)
            ld += ((long double)Ai[j] - (long double)ap[j]) * (long double)xs[j];
        const double dv = (double)ld;
        if (std::isfinite(dv) && !dif.contains(dv)) {
            ++n_chk_ident;
            const double e = std::max(dif.lb()-dv, dv-dif.ub());
            if (e > chk_peor_ident) chk_peor_ident = e;
        }
        if (rango_p.is_superset(pa)) {
            ++n_chk_cond;
            if (!corr.is_superset(pA)) {
                ++n_chk_viola;
                const double e = std::max(corr.lb()-pA.lb(), pA.ub()-corr.ub());
                if (e > chk_peor) chk_peor = e;
            }
        }
    }
    (void)i;
}

int CtcDFBHull::linearize_congelada(const IntervalVector& box, int cong) {
    if (box.is_unbounded()) { reusa_A = false; return 0; }
    LinearizerAffine2* la = static_cast<LinearizerAffine2*>(&lr);

    /* VUELTA BARATA (`DFBH_CONGBARATO=1`).
     *
     * Con `A` congelada, una vuelta del punto fijo puede saltearse la
     * linealizacion entera. La informacion de las restricciones ya esta en la
     * cota ACUMULADA de cada `b_i`, valida mientras la caja solo se encoja; lo
     * unico nuevo que aporta esta vuelta es que la caja es mas chica, y eso lo
     * captura `A_i.X`. Queda `b_i = b_acum ∩ A_i.X`, que cuesta `O(m*n)`
     * operaciones de intervalo y ninguna evaluacion afin.
     *
     * Importa porque la linealizacion afin pesa entre el 10 % y el 70 % del
     * contractor segun la instancia, y congelar `A` no la ahorraba. */
    static const bool barato = (getenv("DFBH_CONGBARATO") != NULL);
    if (barato && !nodo_nuevo && !cong_ids.empty()
        && (int)b_acum.size() == (int)cong_rows.size()
        && caja_acum.size() == nb_var) {
        bool nested = true;
        for (int j = 0; nested && j < nb_var; ++j)
            if (!caja_acum[j].is_superset(box[j])) nested = false;
        if (nested) {
            const int m2 = (int)cong_rows.size();
            reusa_A = true;
            caja_lin = box; lin_valida = false;
            z.resize(nb_var + m2);
            for (int i = 0; i < nb_var; ++i) z[i] = box[i];
            for (int i = 0; i < m2; ++i) {
                const Vector& Ai = cong_rows[i];
                Interval triv(0.0);
                for (int j = 0; j < nb_var; ++j) triv += Interval(Ai[j]) * box[j];
                const Interval bi = triv & b_acum[i];
                if (bi.is_empty()) { reusa_A = false; return -1; }
                z[nb_var+i] = bi;
                b_acum[i] = bi;
            }
            caja_acum = box;
            ++n_barata;
            return m2;
        }
    }

    /* LINEALIZACION SELECTIVA (`DFBH_CONGVAR=r`).
     *
     * Solo se re-evaluan las restricciones alguna de cuyas variables encogio
     * mas de `r` en diametro relativo desde su ultima evaluacion. Las demas
     * conservan `a'_i` y `rango'_i`: siguen siendo validos porque la caja solo
     * se encoge, y la correccion `(A_i - a'_i).X` hasta mejora al evaluarse
     * sobre la caja nueva. Es la logica de una cola de propagacion aplicada a
     * la linealizacion, y apunta al 10-70 % del tiempo del contractor que se va
     * en la afin. */
    static const double congvar = getenv("DFBH_CONGVAR")
                                ? atof(getenv("DFBH_CONGVAR")) : -1.0;
    std::vector<char> mascara;
    const bool selectivo = (congvar >= 0.0) && !nodo_nuevo && !cong_ids.empty()
                        && ult_rows.size() == cong_rows.size();
    if (selectivo) {
        mascara.assign(la->nb_ctr(), 0);
        for (size_t i = 0; i < cong_ids.size(); ++i) {
            const int ctr = cong_ids[i]/2;
            if (ctr < 0 || ctr >= (int)mascara.size() || mascara[ctr]) continue;
            if ((int)ult_diam[i].size() != nb_var) { mascara[ctr] = 1; continue; }
            /* CONTENCION: sin esto la cota guardada puede no valer en esta caja. */
            if (i >= ult_caja.size() || ult_caja[i].size() != nb_var) { mascara[ctr] = 1; continue; }
            bool dentro = true;
            for (int j = 0; dentro && j < nb_var; ++j)
                if (!ult_caja[i][j].is_superset(box[j])) dentro = false;
            if (!dentro) { mascara[ctr] = 1; continue; }
            const std::vector<int>& uv = la->used_vars(ctr);
            for (size_t t = 0; t < uv.size(); ++t) {
                const int j = uv[t];
                if (j < 0 || j >= nb_var) continue;
                const double d0 = ult_diam[i][j], d1 = box[j].diam();
                if (!(d0 > 0.0) || d0 >= 1e18) { mascara[ctr] = 1; break; }
                if (d1 < d0*(1.0 - congvar)) { mascara[ctr] = 1; break; }
            }
        }
        long ev = 0; for (size_t t = 0; t < mascara.size(); ++t) if (mascara[t]) ++ev;
        n_ctr_evaluadas += ev; n_ctr_total += (long)mascara.size();
    }

    std::vector<Vector> rn; std::vector<Interval> gn; std::vector<int> idn;
    const int mm = la->linearize_id(box, rn, gn, idn, selectivo ? &mascara : NULL);
    if (mm == -1) { reusa_A = false; return -1; }

    /* Se congela al empezar un nodo. `cong` ya no es una antiguedad sino
     * solo el interruptor: el alcance del congelamiento es el nodo. */
    /* `DFBH_CONGSIEMPRE=1`: re-congela en CADA llamada. Es un control, no una
     * variante: congelando siempre, la relajacion tiene que ser identica a la
     * del camino normal, porque `A_i = a'_i` y el termino de correccion es
     * exactamente cero. Si los resultados difieren, el error esta en este
     * camino y no en la antiguedad de las filas. */
    static const bool cong_siempre = (getenv("DFBH_CONGSIEMPRE") != NULL);
    const bool congelar = (cong_siempre || nodo_nuevo || cong_ids.empty()
                        || (int)cong_rows.size() == 0);
    if (congelar) {
        if (mm <= 0) { reusa_A = false; cong_ids.clear(); return 0; }
        cong_rows = rn; cong_ids = idn;
        ult_rows = rn; ult_rango = gn;
        ult_diam.assign(rn.size(), std::vector<double>(nb_var, 0.0));
        ult_caja.assign(rn.size(), IntervalVector(nb_var));
        for (size_t i = 0; i < rn.size(); ++i) {
            for (int j = 0; j < nb_var; ++j) ult_diam[i][j] = box[j].diam();
            for (int j = 0; j < nb_var; ++j) ult_caja[i][j] = box[j];
        }
        reusa_A = false;
        b_acum.clear();          /* la acumulacion arranca de cero */
        /* `c` ORIGINAL: el rango capturado en ESTE momento, que es de donde
         * parte la variante sin refresco afin. Se guarda aparte porque
         * `ult_rango` se sobreescribe en cada vuelta con el rango fresco. */
        b_orig = gn; b_orig_cong = gn;
        ++n_congela;
    } else {
        reusa_A = true;
    }
    (void)cong;

    const int m = (int)cong_rows.size();
    if (m <= 0) { reusa_A = false; return 0; }
    n_m_tot += m; ++n_m_llam;
    /* `DFBH_CHKC=k` enciende las sondas del acotamiento afin de `c`: la
     * verificacion por muestreo (k puntos por fila) y las estadisticas de
     * cobertura y apriete. */
    static const bool sondas_c = (getenv("DFBH_CHKC") != NULL);
    static const bool sonda_orig = (getenv("DFBH_CORIG") != NULL);

    caja_lin = box;
    lin_valida = false;
    z.resize(nb_var + m);
    for (int i = 0; i < nb_var; ++i) z[i] = box[i];

    /* La cota acumulada sigue valiendo si la caja de ahora esta CONTENIDA en
     * aquella en que se acumulo: toda solucion de esta lo era de aquella. */
    bool acumular = ((int)b_acum.size() == m && caja_acum.size() == nb_var);
    for (int j = 0; acumular && j < nb_var; ++j)
        if (!caja_acum[j].is_superset(box[j])) acumular = false;

    for (int i = 0; i < m; ++i) {
        const Vector& Ai = cong_rows[i];
        Interval triv(0.0);
        for (int j = 0; j < nb_var; ++j) triv += Interval(Ai[j]) * box[j];
        Interval bi = triv;
        size_t q = idn.size();
        for (size_t t = 0; t < idn.size(); ++t)
            if (idn[t] == cong_ids[i]) { q = t; break; }
        if (sondas_c) {
            ++n_fila_tot;
            if (q >= idn.size()) ++n_fila_vieja;   /* su `c` no se pudo refrescar */
        }
        if (q < idn.size() && i < ult_rows.size()) {       /* se re-evaluo */
            /* DIAGNOSTICO (`DFBH_CONGDIAG=1`): con la fila congelada `A_i`, la
             * cota de `b_i` se puede armar con la linealizacion NUEVA o con la
             * que se conservaba. Las dos son validas. La nueva tiene mejor
             * `rango'`, la vieja tiene `a'` mas cerca de `A_i` —se congelaron
             * en la misma caja— asi que su termino de correccion
             * `(A_i - a'_i).X` es mas chico. Cual gana no es obvio y decide si
             * re-evaluar menos es un mecanismo o solo suerte de trayectoria. */
            static const bool congdiag = (getenv("DFBH_CONGDIAG") != NULL);
            if (congdiag && ult_rows[i].size() == nb_var) {
                Interval dv(0.0), dn(0.0);
                for (int j = 0; j < nb_var; ++j) {
                    dv += (Interval(Ai[j]) - Interval(ult_rows[i][j])) * box[j];
                    dn += (Interval(Ai[j]) - Interval(rn[q][j]))       * box[j];
                }
                const Interval cv = (ult_rango[i] + dv) & triv;
                const Interval cn = (gn[q]        + dn) & triv;
                if (!cv.is_empty() && !cn.is_empty() && triv.diam() > 0.0) {
                    const double a = cv.diam(), b = cn.diam();
                    ++n_diag;
                    if (a < b*(1.0-1e-12))      ++n_diag_vieja;
                    else if (b < a*(1.0-1e-12)) ++n_diag_nueva;
                    else                        ++n_diag_igual;
                    if (a > 0 && b > 0) acc_diag += std::log(a/b);
                    /* y cuanto pesa cada termino */
                    acc_corr_v += dv.diam()/triv.diam();
                    acc_corr_n += dn.diam()/triv.diam();
                }
            }
            ult_rows[i] = rn[q]; ult_rango[i] = gn[q];
            if ((int)ult_diam[i].size() == nb_var)
                for (int j = 0; j < nb_var; ++j) ult_diam[i][j] = box[j].diam();
            if (i < ult_caja.size() && ult_caja[i].size() == nb_var)
                for (int j = 0; j < nb_var; ++j) ult_caja[i][j] = box[j];
        }
        if (i < ult_rows.size() && ult_rows[i].size() == nb_var) {
            /* `A_i.x = a'_i.x + (A_i - a'_i).x`, asi que de `a'_i.x in rango'_i`
             * se sigue `A_i.x in rango'_i + (A_i - a'_i).X`. La resta va en
             * INTERVALOS: `Interval(Ai[j] - ult_rows[i][j])` la hacia en punto
             * flotante y envolvia el resultado ya redondeado, perdiendo hasta
             * un ulp de la diferencia. */
            Interval dif(0.0);
            for (int j = 0; j < nb_var; ++j)
                dif += (Interval(Ai[j]) - Interval(ult_rows[i][j])) * box[j];
            const Interval corr = (ult_rango[i] + dif) & triv;
            if (corr.is_empty()) { reusa_A = false; return -1; }
            /* Sondas: apagadas por omision y sin costo en el camino caliente. */
            if (sondas_c) {
                chequear_c(i, Ai, ult_rows[i], ult_rango[i], box, dif, corr);
                if (triv.diam() > 0.0 && std::isfinite(triv.diam())
                    && std::isfinite(corr.diam())) {
                    ++n_apr; acc_apr += corr.diam()/triv.diam();
                    if (corr.diam() >= triv.diam()*(1.0-1e-12)) ++n_apr_nada;
                }
            }
            bi = corr;
        } else ++n_sinpar;
        /* SONDA `DFBH_CORIG=1`: cuanto aporta refrescar `c` por afin, contra
         * dejar el `c` ORIGINAL del momento de congelar. La cota original se
         * propaga sola —`b_orig` arranca en el rango capturado al congelar y
         * cada vuelta se intersecta con `A_i.X`, que es lo unico nuevo que
         * aporta una caja mas chica— y se lleva en paralelo, sin mirar nunca la
         * refrescada, para que las dos trayectorias de acumulacion sean
         * independientes. Se compara en la MISMA caja, asi que el resultado no
         * tiene el ruido de arbol. */
        if (sonda_orig && (int)b_orig.size() == m) {
            /* La acumulacion de la cota original obedece la MISMA condicion de
             * anidamiento que la refrescada: sin ella se guardaria una cota de
             * una caja que ya no contiene a esta, y la comparacion seria
             * contra algo invalido. */
            const Interval bo = acumular ? (b_orig[i] & triv) : (gn_cong(i, idn, gn) & triv);
            b_orig[i] = bo;
            const Interval br = acumular ? (bi & b_acum[i]) : bi;
            if (!bo.is_empty() && !br.is_empty()
                && std::isfinite(bo.diam()) && std::isfinite(br.diam())) {
                ++n_orig;
                if (bo.diam() > 0.0) acc_orig += std::log(std::max(br.diam(),1e-300)
                                                        / bo.diam());
                if (br.diam() < bo.diam()*(1.0-1e-12))      ++n_orig_mejor;
                else if (bo.diam() < br.diam()*(1.0-1e-12)) ++n_orig_peor;
                else                                        ++n_orig_igual;
            } else if (br.is_empty() && !bo.is_empty()) { ++n_orig; ++n_orig_vacio; }
        }
        if (acumular) {
            const Interval t2 = bi & b_acum[i];
            if (t2.is_empty()) { reusa_A = false; return -1; }
            if (t2.diam() < bi.diam()) ++n_acum_gana;
            bi = t2;
        }
        z[nb_var+i] = bi;
    }
    b_acum.assign(m, Interval::ALL_REALS);
    for (int i = 0; i < m; ++i) b_acum[i] = z[nb_var+i];
    caja_acum = box;

    if (!reusa_A) {
        A.resize(m, nb_var + m);
        A.clear();
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < nb_var; ++j) A[i][j] = cong_rows[i][j];
            A[i][nb_var+i] = Interval(-1.0);
        }
    }
    return m;
}

/* Igual que `una_pasada` pero sin recargar: el simplex ya tiene la
 * relajacion y bases utiles de la llamada anterior del punto fijo. */
long CtcDFBHull::n_lagr_aplic = 0, CtcDFBHull::n_lagr_mejora = 0,
     CtcDFBHull::n_lagr_vacio = 0, CtcDFBHull::n_lagr_fijas = 0;

Interval CtcDFBHull::cota_lagrangiana(const IntervalVector& gamma, int k, int lado) {
    static const int rondas = getenv("DFBH_LAGRR") ? atoi(getenv("DFBH_LAGRR")) : 2;
    LinearizerAffine2* la = dynamic_cast<LinearizerAffine2*>(&lr);
    const int m = A.nb_rows();
    if (!la || (int)la->fila_ctr.size() != m || gamma[k].contains(0.0)) return Interval::ALL_REALS;
    const System& sys = la->sistema();
    const Interval gk = gamma[k];
    const bool inferior = (lado == 0);

    /* filas activas: lado de la restriccion y signo del multiplicador */
    std::vector<char> act(m, 0);
    std::vector<double> aj(m * nb_var);
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < nb_var; ++i) aj[j*nb_var + i] = A[j][i].mid();
        const CmpOp op = sys.ctrs[la->fila_ctr[j]].op;
        const Interval cb = -gamma[nb_var + j] / gk;
        if (op == EQ) { act[j] = 1; continue; }
        Interval enc(0.0);
        for (int i = 0; i < nb_var; ++i) enc += Interval(aj[j*nb_var + i]) * caja_lin[i];
        const Interval& bj = z[nb_var + j];
        const bool leq = (op == LEQ || op == LT);
        /* cota inferior: hace falta cb*g >= 0; superior: cb*g <= 0 */
        const bool signo_ok = inferior ? (leq ? cb.ub() <= 0.0 : cb.lb() >= 0.0)
                                       : (leq ? cb.lb() >= 0.0 : cb.ub() <= 0.0);
        const bool lado_ctr = leq ? (bj.ub() < enc.ub()) : (bj.lb() > enc.lb());
        if (signo_ok && lado_ctr) act[j] = 1;
    }

    IntervalVector Bc = z.subvector(0, nb_var - 1);
    int nfix = 0;
    /* gradiente de N' por intervalos, y caras */
    std::vector<IntervalVector> gradc;   /* por restriccion, cache dentro de la ronda */
    for (int r = 0; r < rondas; ++r) {
        IntervalVector gN(nb_var);
        for (int i = 0; i < nb_var; ++i) gN[i] = (i == k) ? Interval(0.0) : -gamma[i];
        gradc.assign(sys.nb_ctr, IntervalVector(1));
        std::vector<char> hecho(sys.nb_ctr, 0);
        for (int j = 0; j < m; ++j) {
            const Interval gb = gamma[nb_var + j];
            if (!act[j]) {
                for (int i = 0; i < nb_var; ++i) gN[i] -= gb * Interval(aj[j*nb_var + i]);
                continue;
            }
            const int c = la->fila_ctr[j];
            if (!hecho[c]) { gradc[c] = sys.ctrs[c].f.gradient(Bc); hecho[c] = 1; }
            if (gradc[c].is_empty()) return Interval::ALL_REALS;
            for (int i = 0; i < nb_var; ++i)
                gN[i] += gb * (gradc[c][i] - Interval(aj[j*nb_var + i]));
        }
        bool cambio = false;
        for (int i = 0; i < nb_var; ++i) {
            if (i == k) continue;
            const double esp = 1e-12 * (1.0 + Bc[i].mag());
            if (Bc[i].diam() <= 2.0*esp || Bc[i].is_unbounded()) continue;
            /* signo de dPhi/dx_i = signo(dN'_i) * signo(gk) */
            const Interval d = gN[i] * gk;
            int dir = 0;                       /* +1: Phi crece con x_i */
            if (d.lb() > 0.0) dir = 1; else if (d.ub() < 0.0) dir = -1;
            if (dir == 0) continue;
            /* minimizo Phi (cota inferior) o lo maximizo (superior) */
            const bool abajo = inferior ? (dir > 0) : (dir < 0);
            Bc[i] = abajo ? Interval(Bc[i].lb(), Bc[i].lb() + esp) : Interval(Bc[i].ub() - esp, Bc[i].ub());
            cambio = true; ++nfix;
        }
        if (!cambio) break;
    }
    n_lagr_fijas += nfix;

    /* N' sobre la caja reducida, con rho_j en afin */
    Interval N(0.0);
    for (int i = 0; i < nb_var; ++i) if (i != k) N -= gamma[i] * Bc[i];
    Affine2Variables var(Bc);
    std::vector<char> hecho(sys.nb_ctr, 0);
    std::vector<Affine2> afc(sys.nb_ctr);
    for (int j = 0; j < m; ++j) {
        const Interval gb = gamma[nb_var + j];
        if (!act[j]) {
            Interval lin(0.0);
            for (int i = 0; i < nb_var; ++i) lin += Interval(aj[j*nb_var + i]) * Bc[i];
            N -= gb * lin;
            continue;
        }
        const int c = la->fila_ctr[j];
        if (!hecho[c]) { afc[c] = la->forma_afin(c, Bc); hecho[c] = 1; }
        if (afc[c].is_empty()) return Interval::ALL_REALS;
        Affine2 rho = afc[c];
        for (int i = 0; i < nb_var; ++i)
            if (aj[j*nb_var + i] != 0.0) rho = rho - var[i] * aj[j*nb_var + i];
        N += gb * rho.itv();
    }
    const Interval phi = N / gk;
    ++n_lagr_aplic;
    return inferior ? Interval(phi.lb(), POS_INFINITY) : Interval(NEG_INFINITY, phi.ub());
}

void CtcDFBHull::una_pasada_reusada(IntervalVector& box) {
    sin_recarga = true;
    una_pasada(box);
    sin_recarga = false;
}

namespace {
/* Mide el costo por LP de cada pasada y lo carga a la estrategia de arranque
 * tibio del nodo (`DFBH_WARM=medir`). */
struct GuardiaWarm {
    bool on; double t0; long lp0;
    GuardiaWarm(bool o) : on(o), t0(o ? ahora_s() : 0.0), lp0(o ? DFBSimplex::n_solves : 0) {}
    ~GuardiaWarm() {
        if (!on || g_ws_estr < 0) return;
        const long lps = DFBSimplex::n_solves - lp0;
        if (lps <= 0) return;
        const double c = (ahora_s() - t0) / lps;
        double& e = g_ws_ema[g_ws_estr];
        e = (e < 0.0) ? c : 0.9 * e + 0.1 * c;
    }
};
}
void CtcDFBHull::una_pasada(IntervalVector& box) {
    if (warm_medir() && !motor_caja) {
        static bool puesto = false;
        if (!puesto) {   /* devex primal en los dos motores, antes del primer paso primal */
            if (!getenv("DFB_SX_DEVEXP"))  setenv("DFB_SX_DEVEXP", "1", 0);
            if (!getenv("DFB_SX_DEVEXPT")) setenv("DFB_SX_DEVEXPT", "1", 0);
            puesto = true;
        }
    }
    GuardiaWarm _gw(warm_medir() && !motor_caja);
    const int m = (m_uso > 0 && m_uso < A.nb_rows()) ? m_uso : A.nb_rows();
    const int na = nb_var + m;
    const int nc = 2*nb_var;
    const bool hay_testigo = testigo_en(z, nb_var);
    if (hay_testigo) {
        const std::vector<double>& t = testigo();
        int malas = 0;
        for (int j = 0; j < m; ++j) {
            Interval v(0.0);
            for (int c = 0; c < nb_var; ++c) v += A[j][c] * Interval(t[c]);
            if (!v.intersects(z[nb_var+j])) {
                if (malas++ < 5)
                    fprintf(stderr, "[testigo] FILA %d VIOLADA: a.x* = [%.17g,%.17g]  rango b = [%.17g,%.17g]\n",
                            j, v.lb(), v.ub(), z[nb_var+j].lb(), z[nb_var+j].ub());
            }
        }
        static long n_chk = 0;
        if (malas || (++n_chk % 50) == 1)
            fprintf(stderr, "[testigo] pasada con el testigo adentro: m=%d filas violadas=%d\n", m, malas);
    }
    /* Tope de pivotes POR LP. Distinto del de §20, que era por nodo y por eso
     * dejaba sin nada a las ultimas cotas. En el arbol el simplex necesita ~30
     * pivotes por LP (contra ~3 en la caja raiz), asi que aca un tope si muerde.
     * Cortar es valido: la cota se certifica desde el lambda que haya. */
    static const int lpcap = getenv("DFBH_LPCAP") ? atoi(getenv("DFBH_LPCAP")) : 0;
    /* `DFBH_LPCAP_M=k` (2026-10-21): tope de k*m pivotes por LP, proporcional al
     * tamano de la relajacion (cualquier base certifica: cortar da una cota
     * valida, mas debil). */
    static const double lpcap_m = getenv("DFBH_LPCAP_M") ? atof(getenv("DFBH_LPCAP_M")) : 0.0;
    const int tope = lpcap_m > 0.0 ? std::max(1, (int)(lpcap_m * m)) : (lpcap > 0 ? lpcap : 20*m + 200);
    /* Base COMPARTIDA: un solo simplex para las 2n cotas, en vez de 2n
     * ejemplares sembrados. En la caja raiz la siembra ganaba x0.75; con LPs
     * mucho mas duras el balance puede cambiar. */
    const double PREC = 1e-8;   /* misma precision relativa que usa Achterberg */
    /* Perillas de diagnostico: el duelo sin propagacion dijo que el orden de la
     * mas lejana y el salteo eran lo mejor, pero en el arbol el resultado se da
     * vuelta, asi que hay que poder aislarlos. */
    static const char* ord_env = getenv("DFBH_ORDEN");
    /* ORDEN POR INFACTIBILIDADES DUALES (`DFBH_ORDEN=costo`, orden 3).
     *
     * Medido con el camino primal, el orden «cerca» divide por dos los
     * pivotes por LP contra «lejos» pero resuelve mas LPs por nodo y contrae
     * menos: el vertice lejano compra salteos, el cercano compra LPs baratos.
     * La distancia primal |x*_j - borde| es solo un sustituto del costo; el
     * costo real de arrancar la cota j desde la base actual es cuantas no
     * basicas tienen el signo del costo reducido violado para +-e_j, que se
     * lee de la fila de x_j en el tableau. Con el hull decidiendo el SALTEO y
     * este conteo decidiendo el ORDEN, se busca conservar los salteos de
     * «lejos» con los pivotes de «cerca». Desempate: mayor ganancia relativa
     * (con hull) o mayor distancia (sin hull).
     *
     * ES EL DEFAULT. Medido sobre 137 instancias y 4 semillas contra
     * Achterberg + «lejos»: pivotes por LP x0.78, cpu x0.903 (> 0.1 s),
     * x0.924 (> 1 s), x0.960 (> 5 s), 128 pares mas rapidos contra 21, celdas
     * x0.97, 0 optimos incompatibles. Solo rinde MAS que combinado con el
     * hull (x0.918 / x0.943 / x0.965): el hull compra celdas pagando LPs por
     * nodo, y los dos ahorros no se suman. `DFBH_ORDEN=lejos` restaura el
     * default anterior. */
    /* El valor se llama `costo` —no `dual`— para no confundirlo con
     * `DFB_SX_PIVOTEO`, que elige la ESTRATEGIA DE PIVOTEO. Son dos ejes
     * independientes: esto decide QUE cota va despues, aquello COMO se resuelve.
     * `dual` se acepta como sinonimo por compatibilidad. */
    /* `DFBH_ORDEN=costo_o_cerca` (orden 4, 2026-10-20): el conteo donde hay
     * tableau y la MAS CERCANA con el motor revisado. El orden 3 cae a la mas
     * lejana con el revisado (abajo), que en relajaciones grandes es el peor
     * orden para el warm start dual: con compo, Bratu-0120 pasa de 61 a 11
     * piv/LP y Katsurabis-14 de 64 a 48 con la mas cercana. */
    static const int   orden = (ord_env && std::string(ord_env)=="cerca") ? 0 :
                               (ord_env && std::string(ord_env)=="lejos") ? 1 :
                               (ord_env && std::string(ord_env)=="nat")   ? 2 :
                               (ord_env && std::string(ord_env)=="costo_o_cerca") ? 4 : 3;
    /* La caja negra no tiene tableau: el orden 3 lee el de DFB, que en la sonda
     * `DFBH_ABCAJA` es el de la pasada anterior. Por omision usa el de la mas
     * lejana; `DFBH_ORDEN_CAJA=cerca|lejos|nat` le da un orden propio sin
     * tocar el de DFB (los brazos del banco fijan `DFBH_ORDEN` y no cambian). */
    static const char* ordc_env = getenv("DFBH_ORDEN_CAJA");
    static const int   orden_caja = !ordc_env ? ((orden == 3 || orden == 4) ? 1 : orden)
                               : std::string(ordc_env)=="cerca" ? 0
                               : std::string(ordc_env)=="nat"   ? 2 : 1;
    const int ord = motor_caja ? orden_caja : (orden == 4 ? (sx.es_revisado() ? 0 : 3) : orden);
    static const bool  no_salteo = (getenv("DFBH_NOSKIP") != NULL);

    /* --- CORTE TEMPRANO (anytime) ---
     *
     * Con el orden de la mas lejana la ganancia esta muy adelantada (§13): al
     * 50 % del presupuesto ya hay 15.63 de los 34.17 puntos, y el ultimo cuarto
     * compra 4.7. Entonces no conviene un tope fijo de pivotes sino un criterio
     * de RENDIMIENTO DECRECIENTE sobre el barrido ya ordenado:
     *
     *   - una cota «rinde» si contrae z_k mas de `tol` relativo;
     *   - se corta la pasada tras `pac` cotas consecutivas que no rindieron.
     *
     * La paciencia hace falta porque una cota sola puede dar cero y la
     * siguiente mucho. Cortar es SIEMPRE valido: la caja que hay es correcta, y
     * el CtcFixPoint de afuera vuelve a llamar si conviene.
     *
     * `DFBH_MAXPIV` es la alternativa cruda —tope de pivotes por nodo— para
     * comparar contra el criterio de rendimiento.
     *
     * Riesgo a vigilar: cortar antes tambien pierde oportunidades de PROBAR
     * VACIO, que es lo que de verdad poda el arbol (§19). */
    static const double tol = getenv("DFBH_TOL") ? atof(getenv("DFBH_TOL")) : 0.0;
    static const int    pac = getenv("DFBH_PAC") ? atoi(getenv("DFBH_PAC")) : 3;
    static const int    maxpiv = getenv("DFBH_MAXPIV") ? atoi(getenv("DFBH_MAXPIV")) : 0;
    /* Criterios atados a lo que la busqueda USA, en vez de a cuanto trabajo se
     * hace. Lo unico que poda un nodo es la cota inferior de la variable
     * objetivo contra el loup; las demas solo mejoran la caja para bisecar.
     *   DFBH_GOAL=1     : resolver esa cota PRIMERO, antes del orden habitual.
     *   DFBH_SOLOGOAL=1 : ademas, parar ahi. Es el extremo de ese eje. */
    static const bool goal_first = (getenv("DFBH_GOAL") != NULL)
                                || (getenv("DFBH_SOLOGOAL") != NULL);
    static const bool solo_goal  = (getenv("DFBH_SOLOGOAL") != NULL);
    const bool hay_goal = goal_first && goal_var >= 0 && goal_var < nb_var;

    int flojas = 0;
    long piv_nodo = 0;
    double acum_contr = 1.0;   /* producto de d1/d0 desde la ultima propagacion */
    /* Variables apretadas desde la ultima propagacion. Entre dos llamadas a HC4
     * pasan varias cotas, asi que son varias: el conjunto se acumula y se
     * reinicia al propagar. */
    BitSet tocadas(nb_var);
    tocadas.clear();
    /* Contraccion acumulada POR VARIABLE desde la ultima propagacion. HC4
     * propaga DESDE una variable, asi que lo que decide si vale la pena no es
     * un total repartido entre muchas sino que ALGUNA se haya movido lo
     * suficiente. */
    std::vector<double> acum_var(nb_var, 1.0);
    bool hubo_optimo = false;   /* alguna cota de esta pasada probo que P != vacio */
    int est_previos[8] = {0,0,0,0,0,0,0,0};   /* estados de las cotas ya resueltas */

    Af.resize(m, na);
    if (!motor_caja) for (int i = 0; i < m; ++i)
        for (int j = 0; j < na; ++j) {
            /* `mid()` de gaol guarda y restaura el modo de redondeo aun en un
             * intervalo degenerado; la mayoria lo son (ceros, la parte -I). */
            const double l = A[i][j].lb(), u = A[i][j].ub();
            Af[i][j] = (l == u) ? l : A[i][j].mid();
        }

    /* CUANTO CAMBIA LA RELAJACION ENTRE LLAMADAS (`DFBH_ACHG=1`).
     *
     * Decide si tiene sentido heredar la base de un nodo al siguiente. Si los
     * coeficientes de `Abar` cambian mucho, la base heredada es una conjetura
     * sobre OTRO poliedro y no sirve por mas bien implementada que este. Con
     * `LinearizerXTaylor` cambiaban en cada nodo porque expande en una esquina
     * aleatoria; con la afin, que es el default de hoy, hay que remedirlo. */
    if (DFB_ENV("DFBH_ACHG")) {
        static Matrix Aprev(1,1); static bool hay = false;
        static long n = 0, iguales = 0, chico = 0, grande = 0, dim = 0;
        static double acc = 0.0;
        if (hay && Aprev.nb_rows() == m && Aprev.nb_cols() == na) {
            long e = 0, ig = 0, ch = 0, gr = 0; double sum = 0.0;
            for (int i = 0; i < m; ++i)
                for (int j = 0; j < na; ++j) {
                    const double a = Aprev[i][j], b = Af[i][j];
                    if (a == 0.0 && b == 0.0) continue;
                    ++e;
                    const double den = std::max(std::fabs(a), std::fabs(b));
                    const double r = (den > 0.0) ? std::fabs(a-b)/den : 0.0;
                    sum += r;
                    if (r <= 1e-12) ++ig; else if (r <= 1e-2) ++ch; else ++gr;
                }
            if (e > 0) {
                ++n; iguales += 100*ig/e; chico += 100*ch/e; grande += 100*gr/e;
                acc += sum/e;
            }
        } else if (hay) ++dim;
        Aprev = Af; hay = true;
        if (n > 0 && (n % 300) == 0)
            fprintf(stderr, "[achg] comparaciones=%ld (%ld con dimension distinta) | "
                    "entradas no nulas: identicas %.1f%%  cambio <=1%% %.1f%%  "
                    "cambio >1%% %.1f%% | cambio relativo medio %.4f\n",
                    n, dim, (double)iguales/n, (double)chico/n, (double)grande/n, acc/n);
    }

    /* Reparto del costo por nodo, para saber cuanto se va en la SIEMBRA
     * (cargar la relajacion en los 2n ejemplares) contra el pivoteo. */
    static const bool cronometro = (getenv("DFBH_CRONO") != NULL);
    static double t_load = 0.0, t_resto = 0.0, t_solve = 0.0, t_cert = 0.0;
    static long nn2 = 0, n_lp2 = 0;
    static long cert_tot = 0, cert_nulo = 0, cert_flojo = 0, cert_bueno = 0, cert_vacio = 0;
    static long cert_brecha = 0, cert_hay_margen = 0;
    double t_ini = 0.0, t_carga = 0.0;

    /* Con `sin_recarga` el simplex ya tiene la relajacion y una base util de la
     * llamada anterior del punto fijo: no se recarga ni se siembra. */
    if (!sin_recarga) {
        t_ini = cronometro ? ahora_s() : 0.0;
        /* Los `2n` ejemplares cargan la MISMA `Af`: se construye una vez y los
         * demas copian lo ya armado. `load` se llevaba del 16 % al 45 % del
         * tiempo del simplex (§43). */
        /* ARRANQUE TIBIO ADAPTATIVO (`DFBH_ADAPT=u`). El dual reubica lados
         * gratis y el primal paga una iteracion por cada no basica que cambia
         * de cota entre cotas consecutivas (§ESTADO 2026-10-14). Se lleva un
         * promedio movil de esos cambios por LP tibio y, al cargar cada nodo,
         * se elige dual (umbral de motor 30) si superan `u`, primal (15) si no. */
        if (adapt_umbral() > 0.0 && !motor_caja) {
            static bool def_leido = false, def_rev_primal = false;
            if (!def_leido) { def_rev_primal = DFBSimplex::rev_primal; def_leido = true; }
            const bool dual = adapt_ema() < 0.0 || adapt_ema() > adapt_umbral();
            DFBSimplex::pivoteo_dual = dual;
            DFBSimplex::dual_puro = dual;
            DFBSimplex::rev_primal = dual ? def_rev_primal : true;
            DFBSimplex::motor_umbral = dual ? 30 : 15;
            ++adapt_nodos(dual);
        }
        if (warm_medir() && !motor_caja) {
            static bool def_leido = false, def_rev_primal = false;
            if (!def_leido) { def_rev_primal = DFBSimplex::rev_primal; def_leido = true; }
            static const long W = getenv("DFBH_WARM_EXPLORA") ? atol(getenv("DFBH_WARM_EXPLORA")) : 5;
            static const long R = getenv("DFBH_WARM_CADA") ? atol(getenv("DFBH_WARM_CADA")) : 25;
            const long nn = ++g_ws_total;
            int est;
            if (nn <= 2*W || g_ws_ema[0] < 0.0 || g_ws_ema[1] < 0.0) est = (int)(nn % 2);
            else {
                const int mejor = (g_ws_ema[0] <= g_ws_ema[1]) ? 0 : 1;
                est = (nn % R == 0) ? 1 - mejor : mejor;
            }
            g_ws_estr = est;
            DFBSimplex::pivoteo_dual = (est == 0);
            DFBSimplex::dual_puro = (est == 0);
            DFBSimplex::rev_primal = (est == 0) ? def_rev_primal : true;
            ++g_ws_nodos[est];
        }
        if (!motor_caja) sx.load(Af, nb_var);
        t_carga = cronometro ? ahora_s() : 0.0;

    } else {
        t_ini = t_carga = cronometro ? ahora_s() : 0.0;
    }

    hecho.assign(nc, 0);
    /* OJO: `gs`/`gs_ok` son miembros y persisten entre nodos. Sin reiniciarlos,
     * la sonda usaba `gamma` de un nodo anterior —otra relajacion, a veces otra
     * dimension— y daba 100 % de trozos vacios, ademas de romper por indice. */
    gs.assign(nc, IntervalVector(na));
    gs_ok.assign(nc, 0);
    lin_valida = true;
    xstar.assign(na, 0.0);
    jac_cortes_ok = false;
    {   /* El hull se reinicia SIEMPRE. Es una caja de puntos factibles y sirve
         * de techo exacto solo si todos siguen dentro de la caja actual; entre
         * llamadas la caja se encogio, asi que conservarlo seria un techo
         * invalido. */
        hull_primal.resize(nb_var);
        hull_vacio = true;
        /* y los certificados de fila, que valen solo para ESTA relajacion */
        gfila.assign(nb_var, IntervalVector(na));
        gfila_ok.assign(nb_var, 0);
    }
    gamma.resize(na);
    sens_pase = false;
    n_bounds += nc;

    /* BARRIDO CIRCULAR con punto fijo INTERNO.
     *
     * En vez de envolver el contractor en un `CtcFixPoint` que re-entra y
     * RE-LINEALIZA, el ciclo vive adentro: se siguen recorriendo las cotas y se
     * termina cuando pasa una vuelta entera de `2n` sin que el HC4 intercalado
     * haya contraido. Cada contraccion de HC4 corre la meta hacia adelante.
     *
     * Es el criterio de una cola de propagacion aplicado a las cotas: sabe
     * DONDE paso algo en vez de reiniciar todo. Y no re-linealiza —la
     * relajacion sigue valida porque la caja solo se encoge— ni pierde la base,
     * asi que las vueltas extra cuestan mucho menos que la primera.
     *
     * `DFBH_CIRC=k`: tope de `k` vueltas (0 = apagado, comportamiento de
     * siempre: una sola pasada de `2n`). */
    static const int circ = getenv("DFBH_CIRC") ? atoi(getenv("DFBH_CIRC")) : 0;
    const int tope_pasos = circ > 0 ? circ*nc : nc;
    int ultimo_hc4 = 0;
    hc4_actuo = false;
    for (int paso = 0; paso < tope_pasos && (paso - ultimo_hc4) < nc; ++paso) {
        int idx = -1;
        /* Conteo que la heuristica de orden PREDIJO para la cota elegida, para
         * contrastarlo despues con los pivotes que realmente costo
         * (`DFBH_PREDICE=1`). La pregunta es si el mismo predictor sirve en los
         * dos caminos de pivoteo: cuenta columnas con el signo violado, que en
         * el camino primal hay que meter a la base de a un pivote, pero que en
         * el dual se arreglan REUBICANDO, sin pivotear. */
        int cnt_elegido = -1;
        if (paso == 0) {
            /* QUE COTA VA PRIMERO.
             *
             * En el paso 0 la heuristica por costo no puede opinar: la base
             * viene en frio, todas las `x` estan no basicas en su cota
             * inferior, asi que el conteo vale 0 para toda cota `min` y 1 para
             * toda `max`. Y el salteo de Achterberg tampoco, porque `xstar`
             * esta en ceros y marcaria como hechas las cotas que valen 0. Por
             * eso esta cota se elige aparte, y hasta ahora se elegia `idx = 0`
             * —la variable 0, minimo— que es una decision arbitraria.
             *
             * Importa mas de lo que su peso en la pasada sugiere: es la unica
             * que se resuelve en frio, o sea la mas cara (37.7 pivotes contra
             * 4.9 en `ex5_3_2`), y de ella salen la base tibia y el primer
             * `x*` con que trabajan todas las demas.
             *
             * `DFBH_PRIMERA=0|ancha|angosta|azar`. `DFBH_GOAL=1` tiene
             * prioridad y elige la variable objetivo. */
            const int primera = regla_primera;
            if (hay_goal) idx = 2*goal_var;
            else if (primera == 0) idx = 0;
            else if (primera == 3) {
                static unsigned sem = 987654321u;
                sem = sem*1664525u + 1013904223u;
                idx = (int)((sem >> 8) % (unsigned)nc);
            } else {
                int mejor_j = 0; double mejor_d = -1.0;
                for (int j = 0; j < nb_var; ++j) {
                    const double dj = z[j].diam();
                    if (!(dj > 0.0) || dj >= POS_INFINITY) continue;
                    if (mejor_d < 0.0
                        || (primera == 1 ? (dj > mejor_d) : (dj < mejor_d))) {
                        mejor_d = dj; mejor_j = j;
                    }
                }
                idx = 2*mejor_j;
            }
        } else if (solo_goal && hay_goal) {
            ++n_cortes; break;
        } else {
            /* SELECCION POR HULL DE PUNTOS PRIMALES.
             *
             * `gan` = cuanto puede contraer esta cota COMO MAXIMO, dado que
             * todos los `x*` vistos son factibles: para (j,min) es
             * `H.lb_j - X.lb_j`, para (j,max) es `X.ub_j - H.ub_j`. Se saltea
             * si la ganancia relativa cae bajo `tau` y se elige la de `gan`
             * maximo; la pasada termina cuando ninguna pendiente supera `tau`
             * —corte por exito, no por haber recorrido las `2n`—.
             *
             * `tau` nunca baja de 1e-8 para no chocar con el residuo de los
             * testigos primales. */
            static const bool por_hull = (getenv("DFBH_HULL") != NULL);
            static const double tau = std::max(1e-8,
                getenv("DFBH_HULLTAU") ? atof(getenv("DFBH_HULLTAU")) : 1e-6);
            if (por_hull && !hull_vacio) {
                double mejor_gan = -1.0;
                int    mejor_cnt = INT_MAX;
                for (int j = 0; j < nb_var; ++j)
                    for (int lado = 0; lado < 2; ++lado) {
                        const int q = 2*j + lado;
                        if (hecho[q]) continue;
                        const double gan = lado ? (z[j].ub() - hull_primal[j].ub())
                                                : (hull_primal[j].lb() - z[j].lb());
                        const double dd = z[j].diam();
                        const double rel = (dd > 0.0 && dd < POS_INFINITY)
                                         ? gan/dd : gan;
                        if (!(rel > tau)) { hecho[q] = 1; continue; }
                        /* Se ordena por la ganancia RELATIVA, la misma
                         * magnitud con que se saltea. Con la absoluta, las
                         * variables de dominio grande ganaban siempre aunque su
                         * ganancia relativa fuera despreciable, y en instancias
                         * con escalas dispares eso desordenaba la pasada
                         * —`alkylbis` se iba de 70 a 700 celdas—. */
                        if (ord == 3 && !sx.es_revisado()) {
                            const int cnt = sx.dual_infeasibilities(j, lado == 1);
                            if (idx < 0 || cnt < mejor_cnt
                                || (cnt == mejor_cnt && rel > mejor_gan)) {
                                mejor_cnt = cnt; mejor_gan = rel; idx = q;
                                cnt_elegido = cnt;
                            }
                        } else if (rel > mejor_gan) { mejor_gan = rel; idx = q; }
                    }
                if (idx < 0) { ++n_cortes; break; }   /* nadie supera tau */
                hecho[idx] = 1;
                goto elegida;
            }

            /* Salteo de Achterberg + orden de la mas lejana. */
            double mejor = POS_INFINITY;
            int    mejor_cnt = INT_MAX;
            for (int j = 0; j < nb_var; ++j)
                for (int lado = 0; lado < 2; ++lado) {
                    const int q = 2*j + lado;
                    if (hecho[q]) continue;
                    const double bd = lado ? z[j].ub() : z[j].lb();
                    const double dj = std::fabs(xstar[j] - bd);
                    /* `DFBH_SALTEO=ancho`: tolerancia relativa al ANCHO de la
                     * caja. La de Achterberg es relativa a |cota| y, en cajas
                     * angostas lejos del cero, salteaba cotas que todavia se
                     * podian contraer (sonda `DFBH_ABHULL`). */
                    static const bool salteo_ancho = getenv("DFBH_SALTEO") && std::string(getenv("DFBH_SALTEO")) == "ancho";
                    const double wj = z[j].diam();
                    const double rel = (salteo_ancho && wj > 0.0 && wj < POS_INFINITY) ? dj / wj
                                     : (std::fabs(bd) < 1.0) ? dj
                                                             : std::fabs(dj/bd);
                    if (!no_salteo && rel < PREC) { hecho[q] = 1; ++abc_cnt[motor_caja][5]; continue; }
                    /* Con el motor revisado el tableau no esta al dia: orden de
                     * la mas lejana (`orden == 1`). */
                    if (ord == 3 && !sx.es_revisado()) {
                        const int cnt = sx.dual_infeasibilities(j, lado == 1);
                        if (idx < 0 || cnt < mejor_cnt
                            || (cnt == mejor_cnt && -dj < mejor)) {
                            mejor_cnt = cnt; mejor = -dj; idx = q;
                            cnt_elegido = cnt;
                        }
                        continue;
                    }
                    const double key = (ord == 0) ? dj
                                     : (ord == 1 || ord == 3) ? -dj
                                                    : (double)q;
                    if (key < mejor) { mejor = key; idx = q; }
                }
        }
        if (idx < 0) break;
        hecho[idx] = 1;
        elegida: ;
        const int k = idx/2, lado = idx%2;
        /* Predictor alternativo, sobre la cota ya elegida: cuesta un pivote, se
         * mide sin cambiar la seleccion (`DFBH_PREDICE2=1`). */
        static const bool sonda_pinf = (getenv("DFBH_PREDICE2") != NULL);
        int cnt_pinf = -1;
        if (sonda_pinf)
            cnt_pinf = sx.primal_infeasibilities(k, lado == 1, z);

        if (maxpiv > 0 && piv_nodo >= maxpiv) { ++n_cortes; break; }

        /* PRESUPUESTO GRADUADO.
         *
         * Todo lo probado antes fue uniforme: todo exacto (§25), todo truncado
         * (§12, §20) o solo la cota del objetivo (§20). Esto es graduado, y lo
         * motiva el §35: el nodo muere en la cota 0.8a de 18, 2a de 46, 8a de 78
         * —el 5-11 % del recorrido—, asi que ESE tramo tiene que ser exacto,
         * porque de ahi sale la prueba de vacio. Las demas solo mejoran la caja
         * para bisecar, y ahi podria alcanzar con ser parcial.
         *
         * `DFBH_EXACTAS=k`  primeras k cotas al optimo (k<0: fraccion de 2n)
         * `DFBH_RESTO=b`    presupuesto de pivotes para el resto */
        static const int n_exactas = getenv("DFBH_EXACTAS")
                                   ? atoi(getenv("DFBH_EXACTAS")) : -1;
        static const int b_resto   = getenv("DFBH_RESTO")
                                   ? atoi(getenv("DFBH_RESTO")) : 0;
        int tope_paso = tope;
        if (b_resto > 0) {
            /* Negativo = porcentaje de `2n`, para que escale con la instancia. */
            int lim = (n_exactas >= 0) ? n_exactas : (nc*(-n_exactas))/100;
            if (lim < 2) lim = 2;
            if (paso >= lim) tope_paso = b_resto;
        }
        /* TOPE SOLO PARA EL PRIMER LP DE LA PASADA (`DFBH_CAP0=c`, 2026-10-22).
         * Es el mas largo: arranca de la base que dejo el nodo padre (o en frio)
         * y no de la de la cota anterior. Cortarlo da una cota valida (cualquier
         * base certifica) y deja una base intermedia para la cota siguiente. */
        static const int cap0 = getenv("DFBH_CAP0") ? atoi(getenv("DFBH_CAP0")) : 0;
        if (cap0 > 0 && paso == 0) tope_paso = std::min(tope_paso, cap0);

        /* SONDA DE SALIDA TEMPRANA POR EXITO.
         *
         * Hoy el vacio se detecta cuando el LP TERMINA. Pero `lambda` existe en
         * cada pivote, asi que `0 not-in gamma.z` se puede evaluar durante el
         * bucle. Si el nodo muere en el pivote 3, hoy se pagan los 30 restantes
         * MAS todas las cotas siguientes del nodo: lo que se ahorra no es un LP,
         * es el resto del nodo.
         *
         * Antes de implementar nada, la pregunta: ¿en que pivote se vuelve
         * detectable? Si es en el ultimo, la idea muere. */
        static const bool sonda_v = (getenv("DFBH_VACIO_PROBE") != NULL);
        std::vector<std::vector<double> > trz;
        if (sonda_v) sk_traza(sx, &trz);

        DFBSimplex& sk = sx;
        ++n_lps;
        /* Volcado de LPs DUROS. `DFBH_DUMP=prefijo` guarda la relajacion y la
         * caja de las resoluciones que superan `DFBH_DUMP_MIN` pivotes, para
         * poder resolverlas despues con SoPlex y con DFB y comparar. Hay que
         * copiar el estado ANTES de resolver. */
        static const char* dump_pref = getenv("DFBH_DUMP");
        static const int dump_min = getenv("DFBH_DUMP_MIN")
                                  ? atoi(getenv("DFBH_DUMP_MIN")) : 40;
        static int n_dump = 0;
        IntervalVector z_antes(1);
        if (dump_pref && n_dump < 12) z_antes = z;

        /* SONDA DE GENERACION PEREZOSA DE FILAS.
         *
         * Solo el 16-21 % de las filas esta activo en el optimo. La idea:
         * resolver primero con un SUBCONJUNTO —las de afin— y agregar las
         * demas despues; las que no cortan no cuestan un pivote, y la cota
         * intermedia ya es valida, asi que se puede parar en cualquier ronda.
         *
         * Una fila se «desactiva» dandole a su componente `b` la evaluacion por
         * INTERVALOS de la fila sobre la caja: es finita —no puede romper la
         * certificacion como haria un infinito— y no corta por construccion,
         * porque contiene todo valor alcanzable de `A_j.x`.
         *
         * RESULTADO: la fase 2 cuesta 1-8 pivotes, o sea que agregar las filas
         * a un optimo ya alcanzado ES barato, como predecia el conteo de filas
         * activas (DFBH_ACTIVAS: solo el 16-21 % esta activo). Pero el total
         * sale 34-56 % PEOR que resolver compo de una vez, por dos defectos de
         * ESTA implementacion, no de la idea:
         *
         *  (a) El cambio se hace POR COTA, asi que el LP alterna entre dos
         *      poliedros en cada una de las `2n`. La cota k arranca del optimo
         *      de la k-1 con el poliedro completo y se le cambia debajo: el
         *      arranque tibio secuencial, que es de donde sale la ventaja de
         *      DFB, queda destruido. La fase 1 sola cuesta tanto como compo
         *      entero. Habria que hacer dos PASADAS completas sobre las 2n.
         *
         *  (b) Ensanchar la cota `b` vuelve la fila no cortante pero NO la saca
         *      de la matriz: el tableau sigue con `m` filas y el costo por
         *      pivote sigue siendo O(m). Para cobrar el beneficio las filas
         *      tienen que estar AUSENTES, o sea cargar una matriz mas chica en
         *      la fase 1 y sembrar la base hacia la completa en la fase 2.
         *
         * Estimado con lo medido, la version correcta daria ~10.7 pivotes
         * equivalentes contra 18.85, pero pide cargar dos matrices por nodo. */
        static const bool lazy = (getenv("DFBH_LAZY") != NULL);
        int piv_f1 = 0;
        if (lazy) {
            const LinearizerCompo* lc = dynamic_cast<const LinearizerCompo*>(&lr);
            const int corte = lc ? lc->nb_l1 : 0;
            if (corte > 0 && corte < m) {
                std::vector<Interval> guardadas(m - corte);
                for (int j = corte; j < m; ++j) {
                    guardadas[j-corte] = z[nb_var+j];
                    Interval enc(0.0);
                    for (int c = 0; c < nb_var; ++c) enc += A[j][c] * z[c];
                    z[nb_var+j] = enc;               /* ancha y valida */
                }
                sk.solve(k, lado == 1, z, tope_paso, true);
                piv_f1 = sk.iterations();
                for (int j = corte; j < m; ++j) z[nb_var+j] = guardadas[j-corte];
            }
        }

        /* TECHO DE LA GENERACION PEREZOSA, con oraculo.
         *
         * Se resuelve normal, se miran que filas quedaron ACTIVAS (lambda != 0)
         * y se vuelve a resolver con SOLO esas —las demas ensanchadas—. El
         * segundo solve da el mismo optimo, porque el conjunto activo lo
         * determina, y dice cuantos pivotes costaria si uno supiera de antemano
         * que filas hacen falta. Es el techo de cualquier esquema perezoso:
         * ningun metodo real puede hacerlo mejor que el oraculo.
         *
         * OJO: esto NO mide el ahorro por pivote, porque ensanchar no saca la
         * fila de la matriz. Mide solo el ahorro en CANTIDAD de pivotes. El
         * ahorro por pivote seria adicional y proporcional a la reduccion de
         * filas, que se reporta al lado. */
        static const bool techo = (getenv("DFBH_TECHO") != NULL);
        const double ts0 = cronometro ? ahora_s() : 0.0;
        /* COTA LIBRE: el certificado que ya existe con CERO pivotes.
         *
         * Si `k` esta en la base, `y` es la fila de `B^-1` y hay cota
         * certificada antes de pivotear. `solve` con presupuesto 0 la entrega
         * sin tocar la base —ITER_LIMIT no la invalida—, o sea que es la
         * propiedad anytime usada literalmente. La pregunta que decide si la
         * idea vale: de la contraccion que consigue el LP resuelto, cuanta ya
         * estaba disponible gratis. `DFBH_LIBRE=1`. */
        static const bool sonda_libre = (getenv("DFBH_LIBRE") != NULL);
        Interval nb_libre(Interval::ALL_REALS);
        bool hay_libre = false;
        if (sonda_libre && sk.es_basica(k)) {
            sk.solve(k, lado == 1, z, 0, true);
            IntervalVector g0(na);
            for (int c = 0; c < na; ++c) g0[c] = Interval(0.0);
            const std::vector<double>& lam0 = sk.lambda();
            for (int j = 0; j < m; ++j) {
                if (lam0[j] == 0.0) continue;
                const Interval lj(lam0[j]);
                for (int c = 0; c < na; ++c)
                    if (A[j][c] != Interval(0.0)) g0[c] += lj * A[j][c];
            }
            nb_libre = cota_desde_gamma(g0, z, k);
            hay_libre = true;
        }
        /* Techo del hull para el corte anytime (`DFBH_ANYTAU=tau`). El hull
         * acota EXACTAMENTE la contraccion alcanzable de esta cota, asi que si
         * la cota flotante del dual ya llego a ese techo menos `tau*diam`, lo
         * que queda por ganar esta por debajo de la tolerancia. Es el mismo
         * criterio con que el hull saltea cotas enteras, aplicado por pivote. */
        static const double anytau = getenv("DFBH_ANYTAU")
                                   ? atof(getenv("DFBH_ANYTAU")) : 0.0;
        sk.tope_anytime = std::numeric_limits<double>::quiet_NaN();
        if (anytau > 0.0 && !hull_vacio && k < nb_var) {
            const double dd = z[k].diam();
            const double mar = (dd > 0.0 && dd < POS_INFINITY) ? anytau*dd : 0.0;
            const double techo = lado ? hull_primal[k].ub() : hull_primal[k].lb();
            if (techo == techo && techo > NEG_INFINITY && techo < POS_INFINITY)
                sk.tope_anytime = lado ? (techo + mar) : (techo - mar);
        }
        /* SONDA DEL TECHO DE LA IDEA «cortar cuando el objetivo cruza el loup»
         * (`DFBH_GOALPROBE=1`). El optimizador guarda el loup en la cota
         * superior de la variable objetivo, asi que el nodo muere por
         * acotamiento cuando la cota inferior certificada de `z_goal` la supera.
         * Se observa en que pivote lo cruza el valor flotante, SIN cortar, para
         * saber cuanto se ahorraria antes de implementar nada. */
        static const bool goalprobe = (getenv("DFBH_GOALPROBE") != NULL);
        const bool es_goal = goalprobe && k == goal_var && lado == 0;
        sk.obs_tope = es_goal ? z[k].ub() : std::numeric_limits<double>::quiet_NaN();
        static const bool lp_largo2 = (getenv("DFBH_LPLARGO") != NULL);
        static long lps_pasada = 0; ++lps_pasada;
        if (lp_largo2 && (lps_pasada % 1000) == 0)
            fprintf(stderr, "[pasada] LPs=%ld paso=%d de %d (nc=%d) k=%d lado=%d nodos=%ld\n", lps_pasada, paso, tope_pasos, nc, k, lado, n_nodes);
        static const bool sombra = DFB_ENV("DFBH_SOMBRA");
        if (sombra) DFBSimplex::medir_inicio = true;
        /* `DFBH_SOMBRA_BASE=1` (2026-10-19): SC arranca cada LP desde la base final de
         * SoPlex del LP anterior de la sombra (estructurales y holguras basicas, cota de
         * cada no basica), para separar el efecto de la cadena de bases del algoritmo. */
        static const bool base_de_soplex = DFB_ENV("DFBH_SOMBRA_BASE");
        if (sombra && base_de_soplex && !motor_caja && sx_directa && sx_directa->hasBasis()
            && sx_directa->numColsReal() == nb_var && sx_directa->numRowsReal() == A.nb_rows()
            && sk.nb_x() == nb_var) {
            const int mr = sx_directa->numRowsReal();
            std::vector<soplex::SPxSolver::VarStatus> br(mr), bc(nb_var);
            sx_directa->getBasis(br.data(), bc.data());
            std::vector<int> bas; std::vector<char> atup(nb_var + mr, 0);
            for (int i = 0; i < nb_var; ++i) {
                if (bc[i] == soplex::SPxSolver::BASIC) bas.push_back(i);
                else atup[i] = (bc[i] == soplex::SPxSolver::ON_UPPER);
            }
            for (int j = 0; j < mr; ++j) {
                if (br[j] == soplex::SPxSolver::BASIC) bas.push_back(nb_var + j);
                else atup[nb_var + j] = (br[j] == soplex::SPxSolver::ON_UPPER);
            }
            static long f_tam = 0, f_otro = 0, n_vis = 0;
            if ((int)bas.size() == mr) { if (!sk.instalar_base(bas, atup)) ++f_otro; }
            else { ++DFBSimplex::n_inst_fallo; ++f_tam; }
            if (++n_vis % 2000 == 0) fprintf(stderr, "[sombrabase] vistos %ld | tamano de base distinto %ld (ultima %d de %d) | otros fallos %ld\n", n_vis, f_tam, (int)bas.size(), mr, f_otro);
        }
        const double sm0 = (sombra && !motor_caja) ? ahora_s() : 0.0;
        const long pc0 = DFBSimplex::n_primal_cero, pt0 = DFBSimplex::n_primal_total;
        DFBSimplex::Status st = motor_caja ? resolver_caja(k, lado == 1, z)
                                           : sk.solve(k, lado == 1, z, tope_paso, true);
        if (primerdiag_on() && !motor_caja) primerdiag_piv(paso == 0, sk.iterations(), st == DFBSimplex::ITER_LIMIT);
        if (sombra && !motor_caja) {
            /* SOMBRA (`DFBH_SOMBRA=1`, 2026-10-19): el MISMO LP (cota, sentido, caja)
             * se resuelve tambien con la caja negra, en la misma cadena de LPs, para
             * comparar pivotes y tiempo LP a LP. No toca el estado de SC. */
            const double sm1 = ahora_s();
            const int it_sc = sk.iterations();
            const std::vector<double> lam_g = lam_caja, x_g = x_caja;
            int bb_kbas = -1, bb_nxbas = -1;
            if (sx_directa && sx_directa->hasBasis() && sx_directa->numColsReal() == nb_var) {
                std::vector<soplex::SPxSolver::VarStatus> br(sx_directa->numRowsReal()), bc(nb_var);
                sx_directa->getBasis(br.data(), bc.data());
                bb_kbas = bc[k] == soplex::SPxSolver::BASIC;
                bb_nxbas = 0; for (int i = 0; i < nb_var; ++i) bb_nxbas += bc[i] == soplex::SPxSolver::BASIC;
            }
            const DFBSimplex::Status st_bb = resolver_caja(k, lado == 1, z);
            const double sm2 = ahora_s();
            lam_caja = lam_g; x_caja = x_g;
            sombra_anotar(it_sc, it_caja, sm1 - sm0, sm2 - sm1, st == DFBSimplex::OPTIMAL, st_bb == DFBSimplex::OPTIMAL);
            if (it_caja == 0 && it_sc > 0) sombra_cero(it_sc, DFBSimplex::n_primal_total - pt0, DFBSimplex::n_primal_cero - pc0);
            sombra_inicio(it_sc, it_caja, DFBSimplex::inicio_pinf, DFBSimplex::inicio_dinf);
            sombra_base(it_sc, it_caja, DFBSimplex::inicio_kbasica, DFBSimplex::inicio_nxbas, bb_kbas, bb_nxbas);
            sombra_fin(DFBSimplex::fin_nxbas, DFBSimplex::fin_nxcota);
        }
        /* `DFBH_FRIO_SI_NO_CERTIFICA` (encendido por omision; `=0` lo apaga):
         * un INFEASIBLE cuyo rayo no certifica se vuelve a resolver EN FRIO.
         * En CountercurrentReactors2-26 (semilla 2, orden
         * propio) la busqueda se estanca en una zona donde el simplex en
         * caliente declara vacios espurios (rayos de 3 filas, el cero en medio
         * de la suma) que el certificado rechaza, y el nodo no se poda; la
         * misma caja resuelta desde cero se prueba vacia. A 60 s es neutro
         * (x1.00, celdas identicas en el 96 % de las corridas). Vale igual
         * para la caja negra (se borra la base de SoPlex), que ademas reintenta
         * cuando declara infactible sin entregar rayo. */
        static const bool reintenta = []() {
            const char* e = getenv("DFBH_FRIO_SI_NO_CERTIFICA");
            return !(e && e[0] == '0'); }();
        if (reintenta && (st == DFBSimplex::INFEASIBLE || (motor_caja && caja_inf_sin_rayo))) {
            bool certifica = false;
            if (st == DFBSimplex::INFEASIBLE) {
                const std::vector<double>& l0 = motor_caja ? lam_caja : sk.lambda();
                IntervalVector g0(na, Interval(0.0));
                for (int j = 0; j < m; ++j) {
                    if (l0[j] == 0.0) continue;
                    const Interval lj(l0[j]);
                    for (int c = 0; c < na; ++c)
                        if (A[j][c] != Interval(0.0)) g0[c] += lj * A[j][c];
                }
                Interval s0(0.0);
                for (int c = 0; c < na; ++c) s0 += g0[c] * z[c];
                certifica = !s0.contains(0.0);
            }
            if (!certifica) {
                st = motor_caja ? resolver_caja(k, lado == 1, z, true)
                                : sk.solve(k, lado == 1, z, tope_paso, false);
                static long nr = 0, nopt = 0, ninf = 0;
                ++nr; if (st == DFBSimplex::OPTIMAL) ++nopt; else if (st == DFBSimplex::INFEASIBLE) ++ninf;
                if (nr % 1000 == 0)
                    fprintf(stderr, "[frio_si_no_certifica] %ld en frio: OPTIMAL %ld, INFEASIBLE %ld, otros %ld\n",
                            nr, nopt, ninf, nr - nopt - ninf);
            }
        }
        ++abc_cnt[motor_caja][0];
        if (techo) {
            const std::vector<double>& lam0 = sk.lambda();
            std::vector<char> act(m, 0);
            int nact = 0;
            for (int j = 0; j < m; ++j) if (lam0[j] != 0.0) { act[j] = 1; ++nact; }
            if (nact > 0 && nact < m) {
                const int piv_full = sk.iterations();
                std::vector<Interval> guard(m);
                for (int j = 0; j < m; ++j) {
                    guard[j] = z[nb_var+j];
                    if (act[j]) continue;
                    Interval enc(0.0);
                    for (int c = 0; c < nb_var; ++c) enc += A[j][c] * z[c];
                    z[nb_var+j] = enc;
                }
                sk.solve(k, lado == 1, z, tope_paso, false);   /* en FRIO, sin herencia */
                const int piv_or = sk.iterations();
                for (int j = 0; j < m; ++j) z[nb_var+j] = guard[j];
                static long n = 0, pf = 0, po = 0, sm = 0, sa = 0;
                ++n; pf += piv_full; po += piv_or; sm += m; sa += nact;
                if ((n % 500) == 0)
                    fprintf(stderr, "[techo] LPs=%ld | filas %.1f -> activas %.1f (%.0f%%) | "
                            "pivotes: normal %.2f, oraculo en frio %.2f | "
                            "techo del costo por LP ~x%.2f\n",
                            n, (double)sm/n, (double)sa/n, 100.0*sa/sm,
                            (double)pf/n, (double)po/n,
                            ((double)po/n)/((double)pf/n) * ((double)sa/sm));
                sk.solve(k, lado == 1, z, tope_paso, false);   /* restaurar estado */
            }
        }
        if (lazy && piv_f1 >= 0) {
            static long n = 0, a1 = 0, a2 = 0;
            ++n; a1 += piv_f1; a2 += sk.iterations();
            if ((n % 500) == 0)
                fprintf(stderr, "[lazy] LPs=%ld | fase1 (solo afin) %.2f piv | "
                        "fase2 (agregando xtaylor) %.2f piv | total %.2f\n",
                        n, (double)a1/n, (double)a2/n, (double)(a1+a2)/n);
        }
        if (dump_pref && n_dump < 12 && sk.iterations() >= dump_min) {
            char nom[512];
            snprintf(nom, sizeof(nom), "%s_%02d.lp", dump_pref, n_dump++);
            FILE* f = fopen(nom, "w");
            if (f) {
                fprintf(f, "# m na nx k lado pivotes_dfb\n%d %d %d %d %d %d\n",
                        m, na, nb_var, k, lado, sk.iterations());
                fprintf(f, "# Abar (m x na), en doubles\n");
                for (int i = 0; i < m; ++i) {
                    for (int c = 0; c < na; ++c) fprintf(f, "%.17g ", Af[i][c]);
                    fprintf(f, "\n");
                }
                fprintf(f, "# caja z (na intervalos)\n");
                for (int c = 0; c < na; ++c)
                    fprintf(f, "%.17g %.17g\n", z_antes[c].lb(), z_antes[c].ub());
                fclose(f);
            }
        }
        if (cronometro) { t_solve += ahora_s() - ts0; ++n_lp2; }
        /* Desenganchar YA: si el bucle hace `continue` mas abajo, el puntero
         * quedaria apuntando a un vector destruido y la proxima iteracion
         * escribiria en memoria liberada. */
        if (sonda_v) sk_traza(sk, NULL);
        n_pivots += motor_caja ? it_caja : sk.iterations();
        piv_nodo += motor_caja ? it_caja : sk.iterations();

        /* LOS DOS PREDICTORES CONTRA LA REALIDAD.
         *
         * `dual_infeasibilities` cuenta columnas a meter en la base: modelo de
         * costo del PRIMAL. `primal_infeasibilities` cuenta filas que quedan
         * fuera de rango tras reubicar: modelo de costo del DUAL. Se registran
         * los dos contra los pivotes que la cota costo de verdad, para saber
         * cual predice mejor en cada camino. */
        if (DFB_ENV("DFBH_PREDICE") && cnt_elegido >= 0) {
            static long nrec[2][9], piv[2][9], tot = 0;
            const int br = DFBSimplex::pivoteo_dual ? 1 : 0;
            const int b = cnt_elegido > 8 ? 8 : cnt_elegido;
            ++nrec[br][b]; piv[br][b] += sk.iterations(); ++tot;
            if ((tot % 5000) == 0) {
                for (int q2 = 0; q2 < 2; ++q2) {
                    if (!nrec[q2][0] && !nrec[q2][1]) continue;
                    fprintf(stderr, "[predice] %s conteo DUAL (columnas)  ",
                            q2 ? "pivoteo dual  |" : "pivoteo primal|");
                    for (int b2 = 0; b2 <= 8; ++b2)
                        if (nrec[q2][b2])
                            fprintf(stderr, "  %d%s:%.2f(n=%ld)", b2, b2==8?"+":"",
                                    (double)piv[q2][b2]/nrec[q2][b2], nrec[q2][b2]);
                    fprintf(stderr, "\n");
                }
            }
        }
        if (DFB_ENV("DFBH_PREDICE2") && cnt_pinf >= 0) {
            static long nrec[2][9], piv[2][9], tot = 0;
            const int br = DFBSimplex::pivoteo_dual ? 1 : 0;
            const int b = cnt_pinf > 8 ? 8 : cnt_pinf;
            ++nrec[br][b]; piv[br][b] += sk.iterations(); ++tot;
            if ((tot % 5000) == 0)
                for (int q2 = 0; q2 < 2; ++q2) {
                    if (!nrec[q2][0] && !nrec[q2][1]) continue;
                    fprintf(stderr, "[predice] %s conteo PRIMAL (filas)    ",
                            q2 ? "pivoteo dual  |" : "pivoteo primal|");
                    for (int b2 = 0; b2 <= 8; ++b2)
                        if (nrec[q2][b2])
                            fprintf(stderr, "  %d%s:%.2f(n=%ld)", b2, b2==8?"+":"",
                                    (double)piv[q2][b2]/nrec[q2][b2], nrec[q2][b2]);
                    fprintf(stderr, "\n");
                }
        }

        /* PIVOTES POR POSICION EN LA PASADA.
         *
         * Con el simplex secuencial, la cota del `paso` 0 arranca de `B = -I`
         * —en frio, porque `load` reinicia la base— y las siguientes heredan la
         * base de la anterior. La hipotesis es que el grueso del costo esta en
         * esa primera. Si es asi, lo que hay que atacar es el arranque del
         * NODO, no el de cada cota. */
        if (DFB_ENV("DFBH_PORPASO")) {
            /* Por motor (0 = DFB, 1 = caja negra), para compararlos con
             * `DFBH_ABCAJA` sobre las mismas cajas. */
            static std::vector<long> piv_m[2], cnt_m[2];
            static long nn_m[2] = {0, 0};
            const int e = motor_caja ? 1 : 0;
            std::vector<long>& piv = piv_m[e]; std::vector<long>& cnt = cnt_m[e];
            if ((int)piv.size() < nc) { piv.resize(nc, 0); cnt.resize(nc, 0); }
            piv[paso] += motor_caja ? it_caja : sk.iterations(); ++cnt[paso];
            long& nn = nn_m[e];
            if (++nn % 200 == 0) {
                long tot = 0; for (int q = 0; q < nc; ++q) tot += piv[q];
                fprintf(stderr, "[porpaso %s] nc=%d  pivotes por posicion:", e ? "caja" : "DFB ", nc);
                for (int q = 0; q < nc && q < 10; ++q)
                    fprintf(stderr, " %d:%.2f", q, cnt[q] ? (double)piv[q]/cnt[q] : 0.0);
                double resto = 0.0; long cr = 0;
                for (int q = 10; q < nc; ++q) { resto += piv[q]; cr += cnt[q]; }
                if (cr) fprintf(stderr, "  10+:%.2f", resto/cr);
                fprintf(stderr, " | la cota 0 es el %.0f%% de los pivotes del nodo\n",
                        tot ? 100.0*piv[0]/tot : 0.0);
            }
        }

        /* BASE COMPARTIDA ENTRE COTAS CONSECUTIVAS.
         *
         * El arranque tibio dual conserva la base y reubica lados; el primal
         * conserva el vertice. Si el optimo de la cota `paso` tiene casi la
         * misma base que el de `paso-1`, el dual llega sin pivotear. Se cuenta
         * cuantas columnas basicas del optimo actual no lo eran en el anterior:
         * es una cota inferior de los pivotes que necesita cualquier simplex. */
        if (adapt_umbral() > 0.0 && !motor_caja) {
            static std::vector<char> prev_a;
            static bool prev_a_ok = false;
            static long n_a = 0;
            std::vector<char> cur(na, 0);
            for (int c = 0; c < na; ++c)
                cur[c] = sk.es_basica(c) ? 1 : (sk.en_superior(c) ? 3 : 2);
            if (paso > 0 && prev_a_ok && (int)prev_a.size() == na) {
                long lado = 0;
                for (int c = 0; c < na; ++c)
                    if (cur[c] != 1 && prev_a[c] != 1 && prev_a[c] != cur[c]) ++lado;
                double& e = adapt_ema();
                e = (e < 0.0) ? (double)lado : 0.98 * e + 0.02 * (double)lado;
                if ((++n_a % 20000) == 0 && getenv("DFBH_ADAPT_TRAZA"))
                    fprintf(stderr, "[adapt] LPs=%ld  cambios de lado (prom. movil)=%.2f  "
                            "nodos dual=%ld primal=%ld\n", n_a, e,
                            adapt_nodos(true), adapt_nodos(false));
            }
            prev_a.swap(cur); prev_a_ok = true;
        }

        if (DFB_ENV("DFBH_BASECOMP") && !motor_caja) {
            static std::vector<char> prev;
            static bool prev_ok = false;
            static long nlp = 0, sum_dif = 0, sum_m = 0, sum_piv = 0, iguales = 0;
            static long sum_lado = 0, sum_nb = 0;
            /* 1 = basica, 2 = no basica en la inferior, 3 = en la superior */
            std::vector<char> cur(na, 0);
            for (int c = 0; c < na; ++c)
                cur[c] = sk.es_basica(c) ? 1 : (sk.en_superior(c) ? 3 : 2);
            if (paso > 0 && prev_ok && (int)prev.size() == na) {
                long dif = 0, nb = 0, lado = 0, nnb = 0;
                for (int c = 0; c < na; ++c) {
                    if (cur[c] == 1) { ++nb; if (prev[c] != 1) ++dif; }
                    else if (prev[c] != 1) { ++nnb; if (prev[c] != cur[c]) ++lado; }
                }
                ++nlp; sum_dif += dif; sum_m += nb; sum_piv += sk.iterations();
                sum_lado += lado; sum_nb += nnb;
                if (dif == 0) ++iguales;
                /* Por tipo de columna: variable (cota de caja) u holgura de fila;
                 * finita en las dos cotas o no; y lo angosta que es la franja de
                 * cada fila: ancho de su holgura sobre el rango de la fila en la caja. */
                static long fl_var = 0, fl_hol = 0, nb_var_c = 0, nb_hol = 0, nb_fin = 0, fl_inf = 0;
                static double sum_rel = 0.0, sum_rel_fl = 0.0; static long n_rel = 0, n_rel_fl = 0;
                for (int c = 0; c < na; ++c) {
                    if (cur[c] == 1 || prev[c] == 1) continue;
                    const bool flip = (prev[c] != cur[c]);
                    const bool fin = z[c].lb() > -1e300 && z[c].ub() < 1e300;
                    if (fin) ++nb_fin; else if (flip) ++fl_inf;
                    if (c < nb_var) { ++nb_var_c; if (flip) ++fl_var; }
                    else {
                        ++nb_hol; if (flip) ++fl_hol;
                        const int j = c - nb_var;
                        if (fin && j < m) {
                            Interval enc(0.0);
                            for (int q = 0; q < nb_var; ++q) enc += A[j][q] * z[q];
                            const double r = enc.diam() > 0 ? std::min(1.0, z[c].diam()/enc.diam()) : 1.0;
                            sum_rel += r; ++n_rel; if (flip) { sum_rel_fl += r; ++n_rel_fl; }
                        }
                    }
                }
                if ((nlp % (atol(getenv("DFBH_BASECOMP")) > 1 ? atol(getenv("DFBH_BASECOMP")) : 2000)) == 0)
                    fprintf(stderr, "[basecomp2] no basicas: variables %.1f (cambian %.2f), holguras %.1f (cambian %.2f) | "
                            "finitas en las dos cotas %.0f%%, cambios en columnas con una cota infinita %ld | "
                            "franja (ancho holgura / rango de la fila) media %.3g, en las que cambian %.3g\n",
                            (double)nb_var_c/nlp, (double)fl_var/nlp, (double)nb_hol/nlp, (double)fl_hol/nlp,
                            100.0*nb_fin/std::max(1L, nb_var_c+nb_hol), fl_inf,
                            n_rel ? sum_rel/n_rel : 0.0, n_rel_fl ? sum_rel_fl/n_rel_fl : 0.0);
                static const long cada = atol(getenv("DFBH_BASECOMP")) > 1 ? atol(getenv("DFBH_BASECOMP")) : 2000;
                if ((nlp % cada) == 0)
                    fprintf(stderr, "[basecomp] LPs tibios=%ld  m=%.1f  columnas basicas nuevas=%.2f "
                            "(%.1f%% de la base)  base identica en el %.1f%%  pivotes=%.2f  | no basicas que siguen: %.1f, "
                            "cambian de lado: %.2f (%.1f%%)\n",
                            nlp, (double)sum_m/nlp, (double)sum_dif/nlp,
                            100.0*sum_dif/sum_m, 100.0*iguales/nlp, (double)sum_piv/nlp,
                            (double)sum_nb/nlp, (double)sum_lado/nlp,
                            sum_nb ? 100.0*sum_lado/sum_nb : 0.0);
            }
            prev.swap(cur); prev_ok = true;
        }

        /* FILAS ACTIVAS POR ORIGEN.
         *
         * Una fila esta activa en el optimo si su multiplicador es no nulo. Si
         * las filas de xtaylor activas son POCAS, conviene resolver primero con
         * las de afin y agregar las otras despues con pasos duales: las que no
         * cortan no cuestan un solo pivote. Si son la mitad, no hay negocio. */
        if (DFB_ENV("DFBH_ACTIVAS")) {
            const LinearizerCompo* lc = dynamic_cast<const LinearizerCompo*>(&lr);
            if (lc && lc->nb_l1 > 0 && lc->nb_l1 < m) {
                const std::vector<double>& lam = sk.lambda();
                long a1 = 0, a2 = 0;
                for (int j = 0; j < m; ++j)
                    if (lam[j] != 0.0) { if (j < lc->nb_l1) ++a1; else ++a2; }
                static long n = 0, s1 = 0, s2 = 0, t1 = 0, t2 = 0;
                ++n; s1 += a1; s2 += a2; t1 += lc->nb_l1; t2 += m - lc->nb_l1;
                if ((n % 500) == 0)
                    fprintf(stderr, "[activas] LPs=%ld | afin: %.1f de %.1f filas activas "
                            "(%.0f%%) | xtaylor: %.1f de %.1f (%.0f%%) | "
                            "las de xtaylor son el %.0f%% de las activas\n",
                            n, (double)s1/n, (double)t1/n, 100.0*s1/t1,
                            (double)s2/n, (double)t2/n, 100.0*s2/t2,
                            100.0*s2/(s1+s2));
            }
        }

        /* certificacion: gamma = lambda^T Abar, en intervalos. Se hace SIEMPRE,
         * no solo cuando el simplex llega al optimo: `lambda` sale de la
         * factorizacion pase lo que pase, y sirve igual para el test de vacio. */
        const double tc0 = cronometro ? ahora_s() : 0.0;
        for (int c = 0; c < na; ++c) gamma[c] = Interval(0.0);
        const std::vector<double>& lam = motor_caja ? lam_caja : sk.lambda();
        for (int j = 0; j < m; ++j) {
            if (lam[j] == 0.0) continue;
            const Interval lj(lam[j]);
            for (int c = 0; c < na; ++c)
                if (A[j][c] != Interval(0.0)) gamma[c] += lj * A[j][c];
        }

        /* PRUEBA DE VACIO, y es la que faltaba.
         *
         * `gamma = lambda^T Abar` y `Abar z = 0` para toda solucion, asi que
         * `0 in gamma.z` es obligatorio. Si esa suma NO contiene al cero, la
         * caja no tiene ninguna solucion de la relajacion: es un certificado de
         * Farkas verificado en intervalos, del mismo tipo que el test de
         * infactibilidad de Neumaier-Shcherbina que usa CtcPolytopeHull.
         *
         * Sin esto DFB tiraba TODAS las pruebas de vacio: medido en el arbol,
         * PolytopeHull vaciaba 163 de cada 200 nodos y DFB ninguno, con lo que
         * el arbol crecia x1.6 en celdas aunque las cotas fueran igual de
         * buenas. Es la diferencia entre la medicion aislada y la del arbol. */
        if (cronometro) t_cert += ahora_s() - tc0;

        /* SONDA DE DERIVA (`DFB_SX_DERIVA=1`): el mismo certificado con el
         * `lambda` del tableau, sobre la misma caja. Mide cuanto perderia la
         * cota si se sacara la LU del camino. No cambia la trayectoria. */
        if (DFBSimplex::sonda_deriva && (int)sk.y_tab.size() == m) {
            static IntervalVector g2(1);
            if (g2.size() != na) g2.resize(na);
            for (int c = 0; c < na; ++c) g2[c] = Interval(0.0);
            double ymax = 0.0, dmax = 0.0;
            for (int j = 0; j < m; ++j) {
                const double a = lam[j], b = sk.y_tab[j];
                if (std::fabs(a) > ymax) ymax = std::fabs(a);
                if (std::fabs(a-b) > dmax) dmax = std::fabs(a-b);
                if (b == 0.0) continue;
                const Interval lj(b);
                for (int c = 0; c < na; ++c)
                    if (A[j][c] != Interval(0.0)) g2[c] += lj * A[j][c];
            }
            static long n = 0, n_nz = 0, peor6 = 0, peor3 = 0, mejor6 = 0, vac_lu = 0, vac_tab = 0,
                        sin_cert_lu = 0, sin_cert_tab = 0;
            static double s_rel = 0.0, max_rel = 0.0, s_perd = 0.0, s_gan = 0.0;
            ++n;
            if (ymax > 0.0) {
                ++n_nz;
                const double rel = dmax / ymax;
                s_rel += rel; if (rel > max_rel) max_rel = rel;
            }
            Interval f1(0.0), f2(0.0);
            for (int c = 0; c < na; ++c) { f1 += gamma[c] * z[c]; f2 += g2[c] * z[c]; }
            const bool v1 = !f1.contains(0.0), v2 = !f2.contains(0.0);
            if (v1) ++vac_lu;
            if (v2) ++vac_tab;
            if (!v1 && !v2) {
                const Interval c1 = cota_desde_gamma(gamma, z, k), c2 = cota_desde_gamma(g2, z, k);
                const double w = z[k].diam();
                if (c1 == Interval::ALL_REALS) ++sin_cert_lu;
                if (c2 == Interval::ALL_REALS) ++sin_cert_tab;
                if (w > 0.0 && w < POS_INFINITY) {
                    /* ganancia de cada certificado sobre la cota que se contrae, en fraccion del ancho */
                    double g1 = lado ? (z[k].ub() - c1.ub()) / w : (c1.lb() - z[k].lb()) / w;
                    double g2v = lado ? (z[k].ub() - c2.ub()) / w : (c2.lb() - z[k].lb()) / w;
                    if (!(g1 > 0.0)) g1 = 0.0;
                    if (!(g2v > 0.0)) g2v = 0.0;
                    if (g1 > 1.0) g1 = 1.0;
                    if (g2v > 1.0) g2v = 1.0;
                    if (g2v < g1 - 1e-6) ++peor6;
                    if (g2v < g1 - 1e-3) ++peor3;
                    if (g2v > g1 + 1e-6) ++mejor6;
                    if (g2v < g1) s_perd += g1 - g2v; else s_gan += g2v - g1;
                }
            }
            if ((n % 2000) == 0)
                fprintf(stderr, "[deriva] LPs=%ld | lambda: dif rel media %.2e max %.2e (%ld con y!=0) | "
                        "vacio LU=%ld tab=%ld | sin cert LU=%ld tab=%ld | tab peor por >1e-6 del ancho: %ld (%.3f%%), "
                        ">1e-3: %ld, mejor >1e-6: %ld | perdida media %.2e ganancia media %.2e (del ancho)\n",
                        n, n_nz ? s_rel/n_nz : 0.0, max_rel, n_nz, vac_lu, vac_tab, sin_cert_lu, sin_cert_tab,
                        peor6, 100.0*peor6/n, peor3, mejor6, s_perd/n, s_gan/n);
        }

        /* SONDA DEL LAGRANGIANO (`DFBH_LAGR=1`). Solo lee.
         *
         * El certificado del LP del objetivo dice  y = -(sum_{i!=y} g_i z_i +
         * sum_j g_{b_j} b_j)/g_y  con `b_j = a_j.x`. El LP acota cada `b_j` por
         * su rango. El lagrangiano sustituye en las filas ACTIVAS —las que el
         * LP acoto por el lado de la restriccion— `b_j = g_c(x) - rho_j(x)`,
         * con `rho_j = g_c - a_j.x` el resto no lineal, y descarta `cb_j g_c`
         * porque su signo lo hace >= 0 en todo factible. Queda
         *
         *     Phi(x) = sum_{i!=y} cx_i x_i - sum_{act} cb_j rho_j(x) + sum_{inact} cb_j a_j.x
         *
         * y `min_caja Phi` es cota valida de `y`. Evaluada en afin sobre la
         * misma caja, `rho_j` da el mismo intervalo que uso el LP: (b) tiene
         * que coincidir con (a), es el control. Lo nuevo es (c): si `dPhi/dx_i`
         * tiene signo constante, el minimo esta en una cara, y la forma afin de
         * `rho_j` sobre la cara tiene menos error. (b') aisla el otro efecto,
         * `rho_j` sobre la caja ya contraida en vez de la de linealizacion. */
        static const bool sonda_lagr = (getenv("DFBH_LAGR") != NULL);
        if (sonda_lagr && k == goal_var && lado == 0 && st == DFBSimplex::OPTIMAL
            && !gamma[k].contains(0.0)) {
            LinearizerAffine2* la = dynamic_cast<LinearizerAffine2*>(&lr);
            if (la && (int)la->fila_ctr.size() == m) {
                const System& sys = la->sistema();
                const double gy = gamma[k].mid();
                std::vector<double> cx(nb_var, 0.0), cb(m, 0.0);
                for (int i = 0; i < nb_var; ++i) if (i != k) cx[i] = -gamma[i].mid()/gy;
                for (int j = 0; j < m; ++j) cb[j] = -gamma[nb_var+j].mid()/gy;
                const IntervalVector B = z.subvector(0, nb_var-1);

                Interval Sa(0.0);
                for (int i = 0; i < nb_var; ++i) if (i != k) Sa += Interval(cx[i]) * z[i];
                for (int j = 0; j < m; ++j) Sa += Interval(cb[j]) * z[nb_var+j];
                const double lb_a = Sa.lb();

                std::vector<char> act(m, 0);
                std::vector<Vector> aj(m, Vector(nb_var));
                for (int j = 0; j < m; ++j) {
                    for (int i = 0; i < nb_var; ++i) aj[j][i] = A[j][i].mid();
                    const CmpOp op = sys.ctrs[la->fila_ctr[j]].op;
                    Interval enc(0.0);
                    for (int i = 0; i < nb_var; ++i) enc += Interval(aj[j][i]) * caja_lin[i];
                    const Interval& bj = z[nb_var+j];
                    if (op == EQ) act[j] = 1;
                    else if ((op == LEQ || op == LT) && cb[j] < 0.0 && bj.ub() < enc.ub()) act[j] = 1;
                    else if ((op == GEQ || op == GT) && cb[j] > 0.0 && bj.lb() > enc.lb()) act[j] = 1;
                }
                bool ok = true;
                /* cota inferior de Phi: parte lineal sobre X, rho_j sobre R */
                struct Phi {
                    static double lb(const IntervalVector& X, const IntervalVector& R, int k, int nb_var, int m,
                                     const std::vector<double>& cx, const std::vector<double>& cb,
                                     const std::vector<char>& act, const std::vector<Vector>& aj,
                                     LinearizerAffine2* la, bool& ok) {
                        Interval S(0.0);
                        for (int i = 0; i < nb_var; ++i) if (i != k) S += Interval(cx[i]) * X[i];
                        Affine2Variables var(R);
                        for (int j = 0; j < m; ++j) {
                            if (!act[j]) {
                                Interval lin(0.0);
                                for (int i = 0; i < nb_var; ++i) lin += Interval(aj[j][i]) * X[i];
                                S += Interval(cb[j]) * lin;
                                continue;
                            }
                            Affine2 rho = la->forma_afin(la->fila_ctr[j], R);
                            if (rho.is_empty()) { ok = false; return NEG_INFINITY; }
                            for (int i = 0; i < nb_var; ++i)
                                if (aj[j][i] != 0.0) rho = rho - var[i] * aj[j][i];
                            S += Interval(-cb[j]) * rho.itv();
                        }
                        return S.lb();
                    }
                };
                const double lb_b  = Phi::lb(B, caja_lin, k, nb_var, m, cx, cb, act, aj, la, ok);
                if (!ok) ++LagrRep::n_fallo_b;
                const double lb_bp = ok ? Phi::lb(B, B,   k, nb_var, m, cx, cb, act, aj, la, ok) : NEG_INFINITY;
                if (!ok && lb_b > NEG_INFINITY) ++LagrRep::n_fallo_bp;

                /* (c): fijar las variables en que Phi es monotona, hasta dos rondas */
                IntervalVector Bc = B; int nfix = 0;
                for (int ronda = 0; ronda < 2 && ok; ++ronda) {
                    IntervalVector grad(nb_var);
                    for (int i = 0; i < nb_var; ++i) grad[i] = Interval(cx[i]);
                    for (int j = 0; j < m; ++j) {
                        if (!act[j]) { for (int i = 0; i < nb_var; ++i) grad[i] += Interval(cb[j]*aj[j][i]); continue; }
                        IntervalVector gj = sys.ctrs[la->fila_ctr[j]].f.gradient(Bc);
                        if (gj.is_empty()) { ok = false; ++LagrRep::n_fallo_grad; break; }
                        for (int i = 0; i < nb_var; ++i) grad[i] += Interval(-cb[j]) * (gj[i] - Interval(aj[j][i]));
                    }
                    bool cambio = false;
                    for (int i = 0; i < nb_var && ok; ++i) {
                        if (i == k || Bc[i].diam() <= 2e-12 * (1.0 + Bc[i].mag())) continue;
                        /* La cara se toma con un espesor relativo de 1e-12 y no
                         * degenerada: sigue conteniendo al minimo, y evita el
                         * caso especial de radio 0 en las variables afines. */
                        const double esp = 1e-12 * (1.0 + Bc[i].mag());
                        if (grad[i].lb() > 0.0)      { Bc[i] = Interval(Bc[i].lb(), Bc[i].lb() + esp) & Bc[i]; cambio = true; ++nfix; }
                        else if (grad[i].ub() < 0.0) { Bc[i] = Interval(Bc[i].ub() - esp, Bc[i].ub()) & Bc[i]; cambio = true; ++nfix; }
                        if (cambio && Bc[i].diam() <= 0.0) Bc[i] = Interval(Bc[i].lb() - esp, Bc[i].ub() + esp);
                    }
                    if (!cambio) break;
                }
                const bool ok_grad = ok;
                const double lb_c = ok ? Phi::lb(Bc, Bc, k, nb_var, m, cx, cb, act, aj, la, ok) : NEG_INFINITY;
                if (ok_grad && !ok) ++LagrRep::n_fallo_c;

                if (ok) {
                    ++LagrRep::n;
                    const double esc = 1.0 + std::fabs(lb_a), tol = 1e-9*esc;
                    const double ctrl = std::fabs(lb_b - lb_a)/esc;
                    if (ctrl > LagrRep::max_ctrl) LagrRep::max_ctrl = ctrl;
                    if (ctrl > 1e-7) ++LagrRep::n_ctrl_mal;
                    const double ub_y = z[k].ub();
                    const double gap = (ub_y < POS_INFINITY && ub_y > lb_a) ? ub_y - lb_a : -1.0;
                    if (lb_bp > lb_a + tol) ++LagrRep::n_bp_mejor;
                    if (nfix > 0) ++LagrRep::n_fix;
                    LagrRep::s_fix += nfix;
                    if (lb_c > lb_bp + tol) ++LagrRep::n_c_mejor;
                    if (lb_c > lb_a + tol) ++LagrRep::n_c_mejor_a;
                    if (gap > 0.0) {
                        ++LagrRep::n_gap;
                        LagrRep::s_bp  += std::min(1.0, std::max(0.0, lb_bp - lb_a)/gap);
                        LagrRep::s_cbp += std::min(1.0, std::max(0.0, lb_c - lb_bp)/gap);
                        LagrRep::s_ca  += std::min(1.0, std::max(0.0, lb_c - lb_a)/gap);
                        if (lb_bp >= ub_y) ++LagrRep::n_mata_bp;
                        if (lb_c >= ub_y) { ++LagrRep::n_mata; if (lb_bp < ub_y) ++LagrRep::n_mata_c_solo; }
                    }
                    if ((LagrRep::n % 2000) == 0) LagrRep::volcar();
                }
            }
        }

        /* SONDA DEL DUELO: el mismo LP, resuelto tambien con SoPlex. Se engancha
         * ANTES del corte por vacio para poder anotar los dos veredictos. */
        if (duelo().cada > 0 && (n_nodes % duelo().cada) == 0) {
            const Interval zk_antes = z[k];
            const Interval nb_dfb = cota_desde_gamma(gamma, z, k);
            Interval far(0.0);
            for (int c = 0; c < na; ++c) far += gamma[c] * z[c];
            const int vac_dfb = far.contains(0.0) ? 0 : 1;

            LPSolver lp(nb_var, LPSolver::Mode::Certified, 1e-9, 10.0, 100000);
            for (int j = 0; j < m; ++j) {
                Vector row(nb_var, 0.0);
                for (int c = 0; c < nb_var; ++c) row[c] = Af[j][c];
                lp.add_constraint(z[nb_var+j].lb(), row, z[nb_var+j].ub());
            }
            lp.set_bounds(z.subvector(0, nb_var-1));
            lp.set_cost(k, lado ? -1.0 : 1.0);
            /* SOLO el minimize(), no el armado: lo que interesa es el costo
             * por ITERACION de SoPlex, para contrastarlo con el costo por
             * pivote de DFB, que crece proporcional a `m`. */
            const double t_spx0 = ahora_s();
            const LPSolver::Status st_spx = lp.minimize();
            const double t_spx_us = 1e6*(ahora_s() - t_spx0);
            const int it_spx = lp.mysoplex->numIterations();
            Interval nb_spx = Interval::ALL_REALS;
            int vac_spx = 0;
            if (st_spx == LPSolver::Status::OptimalProved) {
                const Interval o = lp.minimum();
                nb_spx = lado ? Interval(NEG_INFINITY, (-o).ub())
                              : Interval(o.lb(), POS_INFINITY);
            } else if (st_spx == LPSolver::Status::Infeasible) {
                vac_spx = 1;
            }
            /* El optimo FLOTANTE de SoPlex, SIN el ensanche de
             * Neumaier-Shcherbina: es el unico termino comparable contra el
             * flotante de DFB para decidir si el simplex resolvio o no. */
            double spx_f = 0.0;
            if (st_spx == LPSolver::Status::OptimalProved
             || st_spx == LPSolver::Status::Optimal) {
                const double v = lp.mysoplex->objValueReal();
                spx_f = lado ? -v : v;
            }

            FILE* fo = duelo().f;
            fprintf(fo, "%ld,%d,%d,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,"
                        "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.6g,%.6g,%d,%d,%.6g,%.3f,%d\n",
                    n_nodes, k, lado, (int)st, sk.iterations(),
                    (int)st_spx, it_spx,
                    zk_antes.lb(), zk_antes.ub(),
                    nb_dfb.is_empty() ? 0.0 : nb_dfb.lb(),
                    nb_dfb.is_empty() ? 0.0 : nb_dfb.ub(),
                    nb_spx.lb(), nb_spx.ub(),
                    (st == DFBSimplex::OPTIMAL) ? sk.bound() : 0.0,
                    spx_f, gamma[k].lb(), gamma[k].ub(),
                    vac_dfb ? 1.0 : reduccion(zk_antes, nb_dfb),
                    vac_spx ? 1.0 : reduccion(zk_antes, nb_spx),
                    vac_dfb, vac_spx, z.max_diam(), t_spx_us, m);
        }

        if (sonda_v) {
            /* primer pivote en que el certificado ya prueba vacio */
            int primero = -1;
            IntervalVector gg(na);
            for (size_t t = 0; t < trz.size(); ++t) {
                for (int c = 0; c < na; ++c) gg[c] = Interval(0.0);
                for (int r = 0; r < m; ++r) {
                    if (trz[t][r] == 0.0) continue;
                    const Interval lr(trz[t][r]);
                    for (int c = 0; c < na; ++c)
                        if (A[r][c] != Interval(0.0)) gg[c] += lr * A[r][c];
                }
                Interval sm(0.0);
                for (int c = 0; c < na; ++c) sm += gg[c] * z[c];
                if (!sm.contains(0.0)) { primero = (int)t; break; }
            }
            static long n_lp_v = 0, n_lp_vac = 0, acc_prim = 0, acc_tot = 0;
            static long n_nodo_v = 0, acc_pos = 0, acc_2n = 0;
            ++n_lp_v;
            if (primero >= 0) {
                ++n_lp_vac;
                acc_prim += primero;
                acc_tot  += (long)trz.size() - 1;
                ++n_nodo_v; acc_pos += paso; acc_2n += nc;
            }
            if ((n_lp_v % 2000) == 0 && n_lp_vac)
                fprintf(stderr, "[vacio] LPs=%ld  con vacio detectable=%ld (%.1f%%) | "
                        "primer pivote %.2f de %.2f (%.0f%%) | "
                        "la cota que lo detecta es la %.1f de %.1f del nodo (%.0f%%)\n",
                        n_lp_v, n_lp_vac, 100.0*n_lp_vac/n_lp_v,
                        (double)acc_prim/n_lp_vac, (double)acc_tot/n_lp_vac,
                        acc_tot ? 100.0*acc_prim/acc_tot : 0.0,
                        (double)acc_pos/n_nodo_v, (double)acc_2n/n_nodo_v,
                        acc_2n ? 100.0*acc_pos/acc_2n : 0.0);
        }
        {
            /* APRETAR LA PARTE `b` ANTES DE CERTIFICAR.
             *
             * El simplex termina `INFEASIBLE` y entrega el rayo de Farkas ya en
             * la primera cota, pero al levantar ese `lambda` a la matriz de
             * intervalos `gamma = lambda^T Abar` se ensancha y `0 in gamma.z`
             * vuelve a cumplirse: la prueba se pierde por aritmetica, no por
             * falta de informacion. Medido, solo el 52-54 % de los vacios se
             * certifica en la cota 0; el resto necesita ~3 lambdas mas.
             *
             * `b_j = A_j.x`, asi que la componente `b` admite reevaluarse por
             * intervalos sobre la `x` actual e intersecarse con la cota de la
             * linealizacion. Es una restriccion valida —no un recorte
             * heuristico— y angosta los terminos de la suma, que es justo lo
             * que decide si el certificado sobrevive.
             *
             * `DFBH_NOAPRIETA=1` lo apaga. */
            static const bool aprieta = (getenv("DFBH_APRIETAB") != NULL);  /* banco: neutro */
            if (aprieta) {
                for (int j = 0; j < m; ++j) {
                    Interval enc(0.0);
                    for (int c = 0; c < nb_var; ++c)
                        if (A[j][c] != Interval(0.0)) enc += A[j][c] * z[c];
                    const Interval t = z[nb_var+j] & enc;
                    if (t.is_empty()) { ++n_empty; box.set_empty(); return; }
                    z[nb_var+j] = t;
                }
            }

            Interval suma(0.0);
            for (int c = 0; c < na; ++c) suma += gamma[c] * z[c];
            /* `DFBH_RAYODBG=k` (2026-10-20): cada k rayos que NO certifican, los
             * terminos de la suma por ancho, la violacion y el lado de la basica. */
            {
                static const long rdbg = getenv("DFBH_RAYODBG") ? atol(getenv("DFBH_RAYODBG")) : 0;
                if (rdbg > 0 && st == DFBSimplex::INFEASIBLE && !motor_caja) {
                    static long nf = 0, nok = 0;
                    if (suma.contains(0.0)) ++nf; else ++nok;
                    if (suma.contains(0.0) && ((nf - 1) % rdbg) == 0) {
                        std::vector<std::pair<double,int> > t;
                        double wtot = 0.0;
                        for (int c = 0; c < na; ++c) { const Interval tc = gamma[c] * z[c]; const double w = tc.diam(); wtot += w; if (w > 0) t.push_back(std::make_pair(w, c)); }
                        std::sort(t.rbegin(), t.rend());
                        int nzl = 0; double lmax = 0.0;
                        const std::vector<double>& lm = sk.lambda();
                        for (int j = 0; j < m; ++j) if (lm[j] != 0.0) { ++nzl; lmax = std::max(lmax, std::fabs(lm[j])); }
                        fprintf(stderr, "[rayodbg] falla #%ld (certifican %ld) | m=%d nnz(lambda)=%d max|lambda|=%.3g | suma=[%.6g, %.6g] ancho total %.3g | mayores terminos:",
                                nf, nok, m, nzl, lmax, suma.lb(), suma.ub(), wtot);
                        for (size_t q = 0; q < t.size() && q < 5; ++q) {
                            const int c = t[q].second;
                            fprintf(stderr, " %s%d gamma=[%.4g,%.4g] z=[%.6g,%.6g] ancho %.3g;", c < nb_var ? "x" : "s", c < nb_var ? c : c - nb_var,
                                    gamma[c].lb(), gamma[c].ub(), z[c].lb(), z[c].ub(), t[q].first);
                        }
                        fprintf(stderr, "\n");
                    }
                }
            }
            if (motor_caja && st == DFBSimplex::INFEASIBLE && getenv("DFBH_CAJA_DBG")) {
                static int dbg2 = 0;
                if (dbg2++ < 3) fprintf(stderr, "[caja] suma=[%.6g,%.6g]\n", suma.lb(), suma.ub());
            }

            /* CUANTO FALTA para certificar, en los casos donde el simplex YA
             * probo infactibilidad pero el certificado en intervalos no
             * sobrevive. `suma` contiene el cero; el margen que falta es
             * min(-lo, hi), y lo que dice si el caso esta al alcance es ese
             * margen RELATIVO al ancho de la suma. Cerca de 0 => la suma apenas
             * cruza el cero y una mejora chica del ancho lo recupera. Cerca de
             * 0.5 => el cero esta en el medio y no hay nada que hacer. */
            /* DENSIDAD DE `lambda`: el ancho de `gamma_i = sum_j lambda_j
             * A[j][i]` se acumula con la cantidad de terminos, asi que un rayo
             * de Farkas mas disperso deberia certificar mas seguido. Se compara
             * nnz(lambda) entre los que certifican y los que no. */
            if (DFB_ENV("DFBH_MARGEN") && st == DFBSimplex::INFEASIBLE) {
                const std::vector<double>& lm = sk.lambda();
                int nz = 0;
                for (int j = 0; j < m; ++j) if (lm[j] != 0.0) ++nz;
                static long nok = 0, nfa = 0, zok = 0, zfa = 0;
                if (suma.contains(0.0)) { ++nfa; zfa += nz; } else { ++nok; zok += nz; }
                if ((nok+nfa) % 500 == 0)
                    fprintf(stderr, "[lambda] m=%d | CERTIFICAN %ld con nnz=%.2f | "
                            "FALLAN %ld con nnz=%.2f\n", m,
                            nok, nok ? (double)zok/nok : 0.0,
                            nfa, nfa ? (double)zfa/nfa : 0.0);
            }
            /* BUSQUEDA DE OTRO CERTIFICADO, sin pivotear.
             *
             * Cuando el rayo del simplex no sobrevive a la verificacion en
             * intervalos, se prueban las `m` filas de `B^-1` como
             * multiplicadores alternativos. Cada una cuesta un BTRAN y ningun
             * pivote, y el test `0 not-in gamma.z` es valido para CUALQUIER
             * `lambda`, asi que no hace falta que sea un rayo del simplex.
             *
             * Medido: recupera el 100 % de los fallos, probando 2.6 de 8 filas
             * y 6.6 de 12.
             *
             * Y NO SIRVE, por eso esta apagado. Las pruebas no se perdian, se
             * POSTERGABAN: la caja moria igual unas 3 cotas despues, mismo
             * nodo y mismo arbol. Lo unico que ahorra son ~3 LPs dentro de los
             * nodos que mueren, que son una fraccion chica del tiempo: medido
             * en 7 instancias da celdas identicas y tiempos dentro del ruido.
             *
             * Vale como dato sobre el metodo: el certificado en intervalos es
             * valido para CUALQUIER lambda, y cuando el rayo del simplex no
             * sobrevive siempre hay otra fila de B^-1 que si. */
            if (DFB_ENV("DFBH_OTRORAYO") && st == DFBSimplex::INFEASIBLE
                && suma.contains(0.0)) {
                std::vector<double> lam2;
                IntervalVector g2(na);
                int probadas = 0, exito = -1;
                for (int r = 0; r < m && exito < 0; ++r) {
                    if (!sk.rayo_de_fila(r, lam2)) break;
                    ++probadas;
                    for (int cc = 0; cc < na; ++cc) g2[cc] = Interval(0.0);
                    for (int jj = 0; jj < m; ++jj) {
                        if (lam2[jj] == 0.0) continue;
                        const Interval lj(lam2[jj]);
                        for (int cc = 0; cc < na; ++cc)
                            if (A[jj][cc] != Interval(0.0)) g2[cc] += lj * A[jj][cc];
                    }
                    Interval s2(0.0);
                    for (int cc = 0; cc < na; ++cc) s2 += g2[cc] * z[cc];
                    if (!s2.contains(0.0)) exito = r;
                }
                if (DFB_ENV("DFBH_OTRORAYO")) {
                    static long n = 0, rec = 0, acc = 0;
                    ++n;
                    if (exito >= 0) { ++rec; acc += probadas; }
                    if (n % 500 == 0)
                        fprintf(stderr, "[otrorayo] fallos=%ld | recuperados: %ld (%.0f%%) | "
                                "filas probadas: %.1f (de m=%d)\n",
                                n, rec, 100.0*rec/n, rec ? (double)acc/rec : 0.0, m);
                }
                if (exito >= 0) { ++n_empty; box.set_empty(); return; }
            }

            if (DFB_ENV("DFBH_MARGEN") && st == DFBSimplex::INFEASIBLE
                && suma.contains(0.0)) {
                const double lo = suma.lb(), hi = suma.ub();
                const double anc = hi - lo;
                if (anc > 0.0 && anc < 1e300) {
                    const double margen = std::min(-lo, hi) / anc;
                    static long n = 0, b1 = 0, b5 = 0, b20 = 0, bres = 0;
                    ++n;
                    /* `DFBH_DUMPFALLO=N`: vuelca la caja (x) del N-esimo fallo, para
                     * resolverla aparte con otro orden. */
                    static const long dump_n = getenv("DFBH_DUMPFALLO") ? atol(getenv("DFBH_DUMPFALLO")) : 0;
                    if (dump_n > 0 && n + 1 == dump_n) {
                        fprintf(stderr, "[dumpfallo]");
                        for (int c = 0; c < nb_var; ++c) fprintf(stderr, " %.17g %.17g", z[c].lb(), z[c].ub());
                        fprintf(stderr, "\n");
                    }
                    if (margen < 0.01) ++b1;
                    else if (margen < 0.05) ++b5;
                    else if (margen < 0.20) ++b20;
                    else ++bres;
                    if (n % 500 == 0)
                        fprintf(stderr, "[margen] fallos de certificacion=%ld | "
                                "margen relativo <1%%: %.0f%%  <5%%: %.0f%%  "
                                "<20%%: %.0f%%  resto: %.0f%%\n",
                                n, 100.0*b1/n, 100.0*b5/n, 100.0*b20/n, 100.0*bres/n);
                }
            }
            if ((int)st >= 0 && (int)st < 8) ++est_previos[(int)st];
            static const bool sin_vacio = (getenv("DFBH_NOVACIO") != NULL);
            if (!suma.contains(0.0) && !sin_vacio) {
                if (testigo_en(z, nb_var)) {
                    const std::vector<double>& t = testigo();
                    Interval gz(0.0);
                    for (int c = 0; c < nb_var; ++c) gz += gamma[c] * Interval(t[c]);
                    Interval gb(0.0);
                    for (int j = 0; j < m; ++j) gb += gamma[nb_var+j] * z[nb_var+j];
                    fprintf(stderr, "[testigo] VACIO POR FARKAS con el testigo adentro: paso=%d k=%d lado=%d estado=%d "
                            "suma=[%.6g,%.6g] gamma_x.x*=[%.6g,%.6g] gamma_b.b=[%.6g,%.6g]\n",
                            paso, k, lado, (int)st, suma.lb(), suma.ub(), gz.lb(), gz.ub(), gb.lb(), gb.ub());
                }
                ++n_empty; ++abc_cnt[motor_caja][3];
                if (traza_cotas)
                    fprintf(stderr, "    paso=%-3d x%-3d %s  estado=%d piv=%-3d | VACIO por Farkas, suma=[%.6g,%.6g]\n",
                            paso, k, lado ? "max" : "min", (int)st, sk.iterations(), suma.lb(), suma.ub());
                /* EN QUE POSICION DE LA PASADA aparece el vacio.
                 *
                 * Si el poliedro sobre la caja no es vacio, NINGUN `lambda`
                 * puede probar vacio: `gamma.z* = lambda^T Abar z* = 0` para
                 * cualquier solucion. Asi que con la caja FIJA el vacio tendria
                 * que aparecer en la cota 0 o nunca. Aparece despues solo
                 * porque la pasada aprieta `z[k] &= nb` sobre la marcha: lo que
                 * poda es la contraccion ACUMULADA, no la cota de turno. */
                if (DFB_ENV("DFBH_POSVACIO")) {
                    /* Si alguna cota previa de esta pasada llego a OPTIMAL, el
                     * poliedro era NO VACIO; y como `z[k] &= nb` recorta solo
                     * fuera de el, no puede haberse vaciado despues. Un vacio
                     * en esa situacion seria incoherente. */
                    static long v0 = 0, vresto = 0, incoherentes = 0;
                    if (paso == 0) ++v0; else { ++vresto; if (hubo_optimo) ++incoherentes; }
                    /* estado con que termino ESTA resolucion, para ver si el
                     * simplex esta llegando al optimo o quedandose corto */
                    static long ac[8] = {0,0,0,0,0,0,0,0};
                    if (paso > 0)
                        for (int q = 0; q < 8; ++q) ac[q] += est_previos[q];
                    if ((v0 + vresto) % 200 == 0)
                        fprintf(stderr, "[posvacio] vacios=%ld | en la cota 0: %ld (%.0f%%) | "
                                "despues: %ld (%.0f%%) | de esos, CON un optimo previo: "
                                "%ld (%.0f%% de los tardios)\n",
                                v0+vresto, v0, 100.0*v0/(v0+vresto),
                                vresto, 100.0*vresto/(v0+vresto),
                                incoherentes, vresto ? 100.0*incoherentes/vresto : 0.0),
                        fprintf(stderr, "           estados de las cotas PREVIAS en los vacios "
                                "tardios: OPT=%ld INFEAS=%ld UNBOUND=%ld ITERLIM=%ld SING=%ld\n",
                                ac[0], ac[1], ac[2], ac[3], ac[4]);
                }
                box.set_empty();
                return;
            }
        }

        if (st != DFBSimplex::OPTIMAL) {
            ++abc_cnt[motor_caja][4];
            if (traza_cotas)
                fprintf(stderr, "    paso=%-3d x%-3d %s  estado=%d piv=%-3d | sin aplicar (no OPTIMAL)\n",
                        paso, k, lado ? "max" : "min", (int)st, sk.iterations());
            continue;
        }
        hubo_optimo = true;
        if (motor_caja) { for (int i = 0; i < nb_var && i < (int)x_caja.size(); ++i) xstar[i] = x_caja[i]; }
        else sk.primal_solution(z, xstar);
        if (modo_cortes > 0) generar_cortes(xstar);
        /* H |= x*, con el punto de ESTE ejemplar, no solo del ultimo. */
        if (hull_vacio) {
            for (int i = 0; i < nb_var; ++i) hull_primal[i] = Interval(xstar[i]);
            hull_vacio = false;
        } else {
            for (int i = 0; i < nb_var; ++i) hull_primal[i] |= Interval(xstar[i]);
        }

        gs[idx] = gamma; gs_ok[idx] = 1;

        Interval nb = cota_desde_gamma(gamma, z, k);
        if (lado ? (nb.ub() >= POS_INFINITY) : (nb.lb() <= NEG_INFINITY)) ++abc_cnt[motor_caja][2];
        else ++abc_cnt[motor_caja][1];

        /* COTA LAGRANGIANA CON MONOTONIA (`DFBH_LAGRC=1|2`), ver la cabecera. */
        static const int lagrc = getenv("DFBH_LAGRC") ? atoi(getenv("DFBH_LAGRC")) : 0;
        if (lagrc > 0 && !nb.is_empty() && (lagrc >= 2 || (k == goal_var && lado == 0))) {
            const Interval lag = cota_lagrangiana(gamma, k, lado);
            if (lag != Interval::ALL_REALS) {
                if (!lag.intersects(nb) || !lag.intersects(z[k])) {
                    /* dos cotas validas sin punto comun: no hay factible en la caja.
                     * Se cuenta como vacio PROPIO solo si la del LP no lo probaba. */
                    if (z[k].intersects(nb)) ++n_lagr_vacio;
                    ++n_empty;
                    box.set_empty();
                    return;
                }
                const Interval antes = nb;
                nb &= lag;
                if (nb != antes) ++n_lagr_mejora;
            }
        }

        if (nb.is_empty()) { if (++flojas >= pac && tol > 0) { ++n_cortes; break; } continue; }
        if (!z[k].intersects(nb)) {
            /* `gamma.z = 0` vale para toda solucion de la relajacion, asi que
             * esto prueba que la caja no contiene ninguna: es riguroso. */
            ++n_empty;
            if (es_goal) {
                int q2 = 0; for (int t = 0; t < nc; ++t) if (!hecho[t]) ++q2;
                GoalRep::anota(true, sk.iterations(), sk.obs_pivote, paso, q2);
            }
            if (traza_cotas)
                fprintf(stderr, "    paso=%-3d x%-3d %s  OPTIMAL piv=%-3d | VACIO por cota: z=[%.12g,%.12g] cota=[%.12g,%.12g]\n",
                        paso, k, lado ? "max" : "min", sk.iterations(), z[k].lb(), z[k].ub(), nb.lb(), nb.ub());
            box.set_empty();
            return;
        }
        /* COTAS SIN CERTIFICADO. `cota_desde_gamma` devuelve el lado infinito
         * cuando `gamma_k` contiene al cero, y eso pasa exactamente cuando `k`
         * quedo NO BASICA en el optimo: entonces `c_B = 0`, `y = 0` y no hay
         * nada que certificar. Se cuenta por separado en cada estrategia
         * porque es el sintoma que distingue al camino dual (`DFBH_ABDUAL`). */
        static long n_sin_cota[2] = {0, 0}, n_cota[2] = {0, 0};
        if (traza_resumen) {
            const int br = DFBSimplex::pivoteo_dual ? 1 : 0;
            const bool sin = lado ? (nb.ub() >= POS_INFINITY) : (nb.lb() <= NEG_INFINITY);
            ++n_cota[br];
            if (sin) ++n_sin_cota[br];
            if (((n_cota[0] + n_cota[1]) % 2000) == 0)
                fprintf(stderr, "[sincota] PRIMAL: %ld de %ld cotas sin certificado (%.1f%%) | "
                        "DUAL: %ld de %ld (%.1f%%)\n",
                        n_sin_cota[0], n_cota[0], n_cota[0] ? 100.0*n_sin_cota[0]/n_cota[0] : 0.0,
                        n_sin_cota[1], n_cota[1], n_cota[1] ? 100.0*n_sin_cota[1]/n_cota[1] : 0.0);
        }
        if (es_goal) GoalRep::anota(false, sk.iterations(), sk.obs_pivote, paso, 0);
        /* Sensibilidades para `LSmearDFB`: del LP del objetivo, cota inferior. */
        if (k == goal_var && lado == 0 && !gamma[k].contains(0.0)) {
            ++n_goal_lp;
            llenar_sens(sens_goal, sk.lambda(), A, gamma, nb_var, k);
            llenar_lam_ctr(lam_ctr_goal, sk.lambda(), lr, gamma, k);
            caja_sens = caja_lin;
            sens_contrajo = nb.lb() > z[k].lb();
            sens_ok = true; sens_pase = true;
        }
        /* Sin objetivo (solver): no hay un LP privilegiado, asi que se SUMAN las
         * sensibilidades `abs` de todos los LPs certificados de la pasada, cada
         * uno normalizado por su |gamma_k|. En el optimizador la de un solo
         * vertice concentraba el peso en 2-4 filas activas y el arbol
         * explotaba; sumar las `2n` cotas da el peso ancho que sirve para
         * bisecar. */
        if (calcular_sens && goal_var < 0 && !gamma[k].contains(0.0)) {
            const std::vector<double>& lam = sk.lambda();
            const double gk = std::fabs(gamma[k].mid());
            if ((int)sens_goal.size() != nb_var || !sens_pase) sens_goal.assign(nb_var, 0.0);
            for (int i = 0; i < (int)lam.size() && i < A.nb_rows(); ++i) {
                if (lam[i] == 0.0) continue;
                const double li = std::fabs(lam[i]) / gk;
                for (int j = 0; j < nb_var; ++j) if (j != k) sens_goal[j] += li * A[i][j].mag();
            }
            lam_ctr_goal.clear();
            caja_sens = caja_lin;
            sens_contrajo = false;
            sens_ok = true; sens_pase = true;
        }
        const double d0 = z[k].diam();
        const Interval zk_pre = z[k];
        if (hay_libre) {
            const double dl = (z[k] & nb_libre).is_empty() ? 0.0 : (z[k] & nb_libre).diam();
            const double ds = (z[k] & nb).is_empty()       ? 0.0 : (z[k] & nb).diam();
            const double gl = d0 - dl, gs = d0 - ds;
            static long n = 0, sin_gan = 0, todo = 0, casi = 0, medio = 0, algo = 0, nada = 0;
            static long piv_todo = 0, piv_casi = 0, piv_tot = 0;
            ++n; piv_tot += sk.iterations();
            if (!(gs > 0.0)) ++sin_gan;
            else {
                const double r = gl/gs;
                if (r >= 0.999)     { ++todo;  piv_todo += sk.iterations(); }
                else if (r >= 0.9)  { ++casi;  piv_casi += sk.iterations(); }
                else if (r >= 0.5)    ++medio;
                else if (r > 0.0)     ++algo;
                else                  ++nada;
            }
            if ((n % 500) == 0)
                fprintf(stderr, "[libre] LPs con k basica=%ld | el LP no contrae: %.0f%% | "
                        "de los que contraen, la cota libre da: todo %.0f%%  >=90%% %.0f%%  "
                        ">=50%% %.0f%%  algo %.0f%%  nada %.0f%% | pivotes ahorrables "
                        "(todo+>=90%%) = %.0f%% de los de este tramo\n",
                        n, 100.0*sin_gan/n,
                        100.0*todo/std::max(1L, n-sin_gan), 100.0*casi/std::max(1L, n-sin_gan),
                        100.0*medio/std::max(1L, n-sin_gan), 100.0*algo/std::max(1L, n-sin_gan),
                        100.0*nada/std::max(1L, n-sin_gan),
                        piv_tot ? 100.0*(piv_todo+piv_casi)/piv_tot : 0.0);
        }
        if (testigo_en(z, nb_var) && !nb.contains(testigo()[k]))
            fprintf(stderr, "[testigo] LA COTA EXCLUYE AL TESTIGO: paso=%d k=%d lado=%d x*_k=%.17g cota=[%.17g,%.17g] z_k=[%.17g,%.17g]\n",
                    paso, k, lado, testigo()[k], nb.lb(), nb.ub(), z[k].lb(), z[k].ub());
        z[k] &= nb;
        const double d1 = z[k].diam();
        if (primerdiag_on() && !motor_caja && d0 > 0.0 && d0 < POS_INFINITY) primerdiag_ctr(paso == 0, 1.0 - d1/d0);
        /* DIAGNOSTICO DEL PRIMER PIVOTE (`DFBH_PIVDIAG`, 2026-10-21): tras
         * contraer la cota k, HC4 sobre una copia de la caja; se agrupa por la
         * violacion relativa del primer pivote dual del LP y se cuenta en cuantos
         * casos la propagacion vacia la caja. */
        static const bool pivdiag_on = getenv("DFBH_PIVDIAG") != NULL;
        if (pivdiag_on && hc4_diag && !motor_caja && k < nb_var) {
            static long cnt[8][4] = {{0}};   /* [cubeta][0 LPs que contraen, 1 HC4 vacia, 2 sale holgura, 3 HC4 vacia y sale holgura] */
            static long ctr_cnt[6][2] = {{0}}; /* por contraccion relativa de la cota: [cubeta][LPs, vacia] */
            static long ntot = 0;
            if (d1 < d0) {
                IntervalVector xb(nb_var);
                for (int i = 0; i < nb_var; ++i) xb[i] = z[i];
                bool vacia = false;
                try { hc4_diag->contract(xb); vacia = xb.is_empty(); } catch (...) { vacia = false; }
                const int pc = sk.primer_columna(); const double pv = sk.primer_violacion();
                const int b = pc < 0 ? 0 : pv < 1e-3 ? 1 : pv < 1e-2 ? 2 : pv < 1e-1 ? 3 : pv < 1.0 ? 4 : pv < 10.0 ? 5 : 6;
                const bool holg = pc >= nb_var;
                cnt[b][0]++; if (vacia) cnt[b][1]++; if (holg) cnt[b][2]++; if (vacia && holg) cnt[b][3]++;
                const double rc = d0 > 0.0 ? (d0 - d1) / d0 : 0.0;
                const int bc = rc < 1e-3 ? 0 : rc < 1e-2 ? 1 : rc < 1e-1 ? 2 : rc < 0.5 ? 3 : rc < 0.9 ? 4 : 5;
                ctr_cnt[bc][0]++; if (vacia) ctr_cnt[bc][1]++;
                if ((++ntot % 5000) == 0) {
                    const char* nb_[] = {"sin pivote","<1e-3","<1e-2","<1e-1","<1",  "<10",">=10"};
                    fprintf(stderr, "[pivdiag] %ld LPs que contraen | por violacion del primer pivote (LPs, %% HC4 vacia, %% sale holgura):", ntot);
                    for (int q = 0; q < 7; ++q) if (cnt[q][0]) fprintf(stderr, " %s: %ld %.1f%% %.0f%%;", nb_[q], cnt[q][0], 100.0*cnt[q][1]/cnt[q][0], 100.0*cnt[q][2]/cnt[q][0]);
                    const char* nc_[] = {"<0.1%","<1%","<10%","<50%","<90%",">=90%"};
                    fprintf(stderr, " | por contraccion de la cota (LPs, %% vacia):");
                    for (int q = 0; q < 6; ++q) if (ctr_cnt[q][0]) fprintf(stderr, " %s: %ld %.1f%%;", nc_[q], ctr_cnt[q][0], 100.0*ctr_cnt[q][1]/ctr_cnt[q][0]);
                    fprintf(stderr, "\n");
                }
            }
        }
        if (traza_cotas) {
            const char* en[] = {"OPTIMAL","INFEASIBLE","UNBOUNDED","ITER_LIMIT","SINGULAR","?","?","?"};
            fprintf(stderr, "    paso=%-3d x%-3d %s  %-10s piv=%-3d | z antes=[%.12g,%.12g] "
                    "cota=[%.12g,%.12g] despues=[%.12g,%.12g]  gano=%.3g\n",
                    paso, k, lado ? "max" : "min",
                    ((int)st >= 0 && (int)st < 5) ? en[(int)st] : "?",
                    sk.iterations(), zk_pre.lb(), zk_pre.ub(),
                    nb.is_empty() ? 0.0 : nb.lb(), nb.is_empty() ? 0.0 : nb.ub(),
                    z[k].lb(), z[k].ub(), d0 - d1);
        }
        if (d1 < d0) {
            tocadas.add(k);
            if (d0 > 0.0 && d0 < POS_INFINITY) acum_var[k] *= (d1/d0);
        }

        /* COSECHA DE LAS FILAS DE `B^-1` DE LAS OTRAS BASICAS (`DFBH_FILAS=1`).
         *
         * Al terminar la cota `k`, cada OTRA variable `j` que quedo en la base
         * tiene su propia fila de `B^-1`, que es un `lambda` valido y DISTINTO
         * del que certifica a `k`. Como `0 in gamma.z` vale para cualquier
         * `lambda`, de ahi sale una cota certificada para `z_j` sin pivotear:
         * son las cotas que el pivoteo ya contrajo de paso.
         *
         * No es lo mismo que las colaterales del historico, que usaban el
         * `gamma` del objetivo: ese `gamma_j` vale exactamente 0 en toda basica
         * —su costo reducido es 0— asi que no daba cota para ninguna de ellas,
         * solo fijacion por costo reducido para las NO basicas. Aca hay un
         * BTRAN por basica y un `gamma` propio.
         *
         * Solo se cosechan las cotas que la pasada todavia no dio por hechas,
         * para no pagar el `gamma` de lo que ya no se va a usar. */
        static const bool cosecha_filas = (getenv("DFBH_FILAS") != NULL);
        /* Cada cuantas resoluciones se CONSTRUYEN certificados nuevos. Con 1 se
         * reconstruye siempre, que es la version cara. Con k > 1 se construye
         * cada k y en las demas solo se REEVALUA lo cacheado sobre la caja ya
         * apretada, que cuesta `O(na)` por certificado en vez de
         * `O(nnz(lambda)*na)` y ademas mejora sola: la misma `gamma` sobre una
         * caja mas chica da una cota mas fuerte. */
        static const int filas_cada = getenv("DFBH_FILASCADA")
                                    ? atoi(getenv("DFBH_FILASCADA")) : 1;
        if (cosecha_filas && st == DFBSimplex::OPTIMAL
            && (int)gfila_ok.size() == nb_var) {
            static long n_sol = 0;
            const bool construir = ((n_sol++ % filas_cada) == 0);
            static long n_con = 0, n_re = 0, n_gan = 0, n_int = 0;
            std::vector<double> lam_j;

            if (construir) {
                for (int j = 0; j < nb_var; ++j) {
                    if (j == k) continue;
                    if (hecho[2*j] && hecho[2*j+1]) continue;
                    const int r = sk.fila_de(j);
                    if (r < 1) continue;
                    if (!sk.rayo_de_fila(r - 1, lam_j)) continue;
                    IntervalVector& gj = gfila[j];
                    if (gj.size() != na) gj.resize(na);
                    for (int c = 0; c < na; ++c) gj[c] = Interval(0.0);
                    for (int jj = 0; jj < m; ++jj) {
                        if (lam_j[jj] == 0.0) continue;
                        const Interval lj(lam_j[jj]);
                        for (int c = 0; c < na; ++c)
                            if (A[jj][c] != Interval(0.0)) gj[c] += lj * A[jj][c];
                    }
                    gfila_ok[j] = 1;
                    ++n_con;
                }
            }

            /* Reevaluacion: barata y valida siempre, venga el certificado de
             * esta base o de una anterior. */
            for (int j = 0; j < nb_var; ++j) {
                if (!gfila_ok[j] || j == k) continue;
                if (hecho[2*j] && hecho[2*j+1]) continue;
                if (gfila[j].size() != na) continue;
                const Interval cj = cota_desde_gamma(gfila[j], z, j);
                ++n_re; ++n_int;
                if (cj.is_empty()) continue;
                if (!z[j].intersects(cj)) { ++n_empty; box.set_empty(); return; }
                const double e0 = z[j].diam();
                z[j] &= cj;
                const double e1 = z[j].diam();
                if (e1 < e0) {
                    ++n_gan;
                    tocadas.add(j);
                    if (e0 > 0.0 && e0 < POS_INFINITY) acum_var[j] *= (e1/e0);
                }
            }
            if (DFB_ENV("DFBH_FILASTAT") && (n_int % 20000) == 0)
                fprintf(stderr, "[filas] evaluadas=%ld | construidas=%ld (%.2f por evaluacion) | "
                        "contraen %.1f%%\n", n_re, n_con, n_re ? (double)n_con/n_re : 0.0,
                        100.0*n_gan/std::max(1L, n_re));
        }

        /* COLATERALES COMO ALIMENTO DEL HC4.
         *
         * `gamma.z = 0` vale para toda coordenada, asi que UNA cota da cota
         * para toda `j` con `0 not-in gamma_j`. Solas no sirven —medido: no
         * aportan contraccion, porque el LP llega a lo mismo al resolver la
         * cota `j`, ni ahorran pivotes—. Pero ese argumento no vale cuando hay
         * HC4 intercalado: la colateral no da contraccion nueva, da la MISMA
         * antes, y HC4 la propaga por las restricciones NO lineales, que la
         * pasada no puede replicar llegue cuando llegue.
         *
         * Cada `j` que contraiga se marca en `tocadas` para que entre en el
         * `impact` de la proxima propagacion. `DFBH_COLAT=1`. */
        static const bool colat_hc4 = (getenv("DFBH_COLAT") != NULL);
        if (colat_hc4) {
            for (int j = 0; j < nb_var; ++j) {
                if (j == k || gamma[j].contains(0.0)) continue;
                const Interval cj = cota_desde_gamma(gamma, z, j);
                if (cj.is_empty()) continue;
                if (!z[j].intersects(cj)) { ++n_empty; box.set_empty(); return; }
                const double e0 = z[j].diam();
                z[j] &= cj;
                const double e1 = z[j].diam();
                if (e1 < e0) {
                    tocadas.add(j);
                    if (e0 > 0.0 && e0 < POS_INFINITY) acum_var[j] *= (e1/e0);
                }
            }
        }

        /* HC4 INTERCALADO, sobre el sistema ORIGINAL.
         *
         * Es lo unico que aporta informacion que la relajacion lineal no puede
         * derivar: propaga por las restricciones NO lineales. Las colaterales y
         * el apriete de `b` ya estaban implicados por el LP y por eso no
         * rindieron; esto no.
         *
         * Se inyecta DESPUES de aplicar la cota, no antes: asi propaga la
         * contraccion recien obtenida y la cota siguiente corre sobre la caja
         * ya apretada. Con el orden inverso cada pasada va un paso atrasada y
         * la contraccion de la ultima cota no se propaga nunca.
         *
         * Se inyecta ENTRE cotas para que las restantes corran sobre una caja
         * mas chica —y sobre todo para que `gamma.z` se evalue mas angosto, que
         * es lo que decide si el certificado sobrevive—. La cantidad que deberia
         * moverse es la fraccion de nodos que DFB mata, que tiene apalancamiento
         * x2.5 a x10 sobre el arbol.
         *
         * `DFBH_HC4MED=k`: correr HC4 cada `k` cotas. */
        static const int hc4_cada = getenv("DFBH_HC4MED")
                                  ? atoi(getenv("DFBH_HC4MED")) : 0;
        /* DISPARO POR CONTRACCION en vez de por contador.
         *
         * Un contador fijo gasta HC4 en las cotas que no contrajeron nada, donde
         * no hay informacion nueva que propagar. La magnitud ya esta calculada
         * —`d0` y `d1` deciden si la cota «rindio»— asi que condicionarlo es un
         * `if`. `DFBH_HC4RED=r`: correr HC4 cuando la cota redujo el diametro
         * mas de `r` relativo. */
        static const double hc4_red = getenv("DFBH_HC4RED")
                                    ? atof(getenv("DFBH_HC4RED")) : -1.0;
        /* ACUMULADO, que es la forma correcta del criterio anterior.
         *
         * Una cota sola rara vez contrae un 10 %, pero cinco del 2 % si suman
         * algo que vale propagar; con umbral POR COTA ninguna dispara y la
         * informacion se acumula sin usarse. Se lleva el producto de `d1/d0`
         * desde la ultima llamada —que es la contraccion total real, venga de
         * una cota o de seis— y se resetea al propagar.
         *
         * `DFBH_HC4ACUM=r`. */
        static const double hc4_acum = getenv("DFBH_HC4ACUM")
                                     ? atof(getenv("DFBH_HC4ACUM")) : -1.0;
        if (hc4_acum >= 0.0 && d0 > 0.0 && d0 < POS_INFINITY && d1 <= d0)
            acum_contr *= (d1/d0);

        /* `DFBH_HC4VAR=r`: disparar si ALGUNA variable acumulo mas de `r`. */
        static const double hc4_var = getenv("DFBH_HC4VAR")
                                    ? atof(getenv("DFBH_HC4VAR")) : -1.0;
        bool disparar = false;
        if (hc4_var >= 0.0) {
            for (int i = 0; i < nb_var && !disparar; ++i)
                if ((1.0 - acum_var[i]) > hc4_var) disparar = true;
        }
        else if (hc4_acum >= 0.0)
            disparar = ((1.0 - acum_contr) > hc4_acum);
        else if (hc4_red >= 0.0)
            disparar = (d0 > 0.0 && d0 < POS_INFINITY && (1.0 - d1/d0) > hc4_red);
        else if (hc4_cada > 0)
            disparar = ((paso % hc4_cada) == (hc4_cada-1));
        Ctc* ctc_mid = sonda_ctc ? sonda_ctc : hc4_diag;   /* en el solver de sistemas, el HC4 del enganche */
        if (disparar && ctc_mid) {
            IntervalVector bx(nb_var);
            for (int i = 0; i < nb_var; ++i) bx[i] = z[i];
            /* PROPAGACION INCREMENTAL. `hc4` se construye con
             * `incremental=true`, asi que `CtcPropag` siembra la agenda SOLO
             * con las restricciones que contienen las variables marcadas en
             * `impact`. Llamarlo completo recorre todas para revisar unas
             * pocas. `DFBH_HC4FULL=1` vuelve al modo completo. */
            /* QUE variables entran al `impact`.
             *
             * Con el criterio por variable se propaga SOLO desde las que
             * cruzaron el umbral, y las demas quedan PENDIENTES para la
             * proxima: el `impact` queda chico y enfocado —que es donde la
             * propagacion incremental rinde— y ninguna contraccion se pierde,
             * simplemente espera a acumular lo suficiente por si sola.
             *
             * Con el contador o el acumulado global no hay criterio por
             * variable, asi que entran todas las tocadas. */
            BitSet prop(nb_var);
            prop.clear();
            if (hc4_var >= 0.0) {
                for (int i = 0; i < nb_var; ++i)
                    if (tocadas[i] && (1.0 - acum_var[i]) > hc4_var) prop.add(i);
            } else {
                prop = tocadas;
            }

            /* Nota: con criterio por variable y `prop` vacio se propagaba
             * completo. Se probo la variante estricta —no propagar— y da
             * resultados IDENTICOS en las 5 instancias de control: con umbral
             * alto el caso no ocurre, siempre hay al menos una que cruzo. */
            static const bool hc4_full = (getenv("DFBH_HC4FULL") != NULL);
            if (hc4_full || prop.empty()) {
                ctc_mid->contract(bx);
                tocadas.clear();
                std::fill(acum_var.begin(), acum_var.end(), 1.0);
            } else {
                ContractContext ctx(bx);
                ctx.impact = prop;
                ctc_mid->contract(bx, ctx);
                /* se limpian SOLO las propagadas */
                for (int i = 0; i < nb_var; ++i)
                    if (prop[i]) { tocadas.remove(i); acum_var[i] = 1.0; }
            }
            if (bx.is_empty()) { ++n_empty; box.set_empty(); return; }
            bool hc4_contrajo = false;
            for (int i = 0; i < nb_var; ++i) {
                if (bx[i].diam() < z[i].diam()) hc4_contrajo = true;
                z[i] = bx[i];
            }
            if (hc4_contrajo) hc4_actuo = true;
            if (hc4_contrajo && circ > 0) {
                /* La contraccion corre la meta, y reabre las cotas salteadas:
                 * el salteo de Achterberg decidio mirando el punto primal de
                 * antes, y esas decisiones quedaron obsoletas. */
                ultimo_hc4 = paso;
                hecho.assign(nc, 0);
                hecho[idx] = 1;
            }
            acum_contr = 1.0;
        }

        const bool rindio = (d0 > 0.0) && (d1 < d0*(1.0 - tol));
        if (tol > 0.0) {
            if (rindio) flojas = 0;
            else if (++flojas >= pac) { ++n_cortes; break; }
        }
    }

    if (cronometro) {
        t_load  += t_carga - t_ini;
        t_resto += ahora_s() - t_carga;
        if ((++nn2 % 50) == 0)
        {
            const double T = t_load + t_resto;
            fprintf(stderr, "[crono] n=%d m=%d na=%d | nodos=%ld LPs=%ld | "
                    "carga %.0f%%  simplex %.0f%%  certif %.0f%%  otro %.0f%% | "
                    "LP/nodo %.1f (de %d)  piv/LP %.2f  degen %.0f%%  us/LP %.1f "
                    "(total %.2fs)\n",
                    nb_var, m, na, nn2, n_lp2,
                    100*t_load/T, 100*t_solve/T, 100*t_cert/T,
                    100*(T-t_load-t_solve-t_cert)/T,
                    (double)n_lp2/nn2, 2*nb_var, (double)n_pivots/n_lp2,
                    n_pivots ? 100.0*DFBSimplex::n_piv_degen/n_pivots : 0.0,
                    1e6*t_solve/n_lp2, T);
            if (cert_tot)
                fprintf(stderr, "        certificados=%ld | nulos %.0f%%  flojos %.0f%%  "
                        "buenos %.0f%% | el flotante ve margen en %.0f%%, y de esos el "
                        "certificado NO lo captura en %.0f%% (= %.0f%% del total)\n", cert_tot,
                        100.0*cert_nulo/cert_tot, 100.0*cert_flojo/cert_tot,
                        100.0*cert_bueno/cert_tot,
                        100.0*cert_hay_margen/cert_tot,
                        cert_hay_margen ? 100.0*cert_brecha/cert_hay_margen : 0.0,
                        100.0*cert_brecha/cert_tot);
        }
    }

    /* SENSIBILIDADES DEL OBJETIVO CUANDO LA PASADA NO LAS DEJO. Si la cota
     * inferior de `y` se salteo, o termino con `y` no basica (`y = 0`, sin
     * informacion), se resuelve `min y` con la cota inferior de `y` SUELTA:
     * es lo que hace `LSmear` (`set_bounds(goal, [-1e10, 1e10])`), y fuerza a
     * `y` a la base, de donde salen multiplicadores no triviales. La cota que
     * da es valida (el poliedro solo se agrando) pero no se usa para
     * contraer: solo se guarda `lambda`. */
    if (calcular_sens && !sens_pase && goal_var >= 0 && goal_var < nb_var && !box.is_empty()) {
        const Interval zg = z[goal_var];
        z[goal_var] = Interval(NEG_INFINITY, zg.ub());
        const DFBSimplex::Status st2 = sx.solve(goal_var, false, z, tope, true);
        ++n_sens_extra; ++n_lps; n_pivots += sx.iterations();
        if (st2 == DFBSimplex::OPTIMAL) {
            const std::vector<double>& lam = sx.lambda();
            for (int c = 0; c < na; ++c) gamma[c] = Interval(0.0);
            for (int j = 0; j < m; ++j) {
                if (lam[j] == 0.0) continue;
                const Interval lj(lam[j]);
                for (int c = 0; c < na; ++c)
                    if (A[j][c] != Interval(0.0)) gamma[c] += lj * A[j][c];
            }
            if (!gamma[goal_var].contains(0.0)) {
                llenar_sens(sens_goal, lam, A, gamma, nb_var, goal_var);
                llenar_lam_ctr(lam_ctr_goal, lam, lr, gamma, goal_var);
                caja_sens = caja_lin;
                z[goal_var] = zg;
                sens_contrajo = cota_desde_gamma(gamma, z, goal_var).lb() > zg.lb();
                sens_ok = true; sens_pase = true;
            }
        }
        z[goal_var] = zg;
    }

    for (int i = 0; i < nb_var; ++i) box[i] = z[i];

}

bool CtcDFBHull::gamma_prueba_vacio(const IntervalVector& box) const {
    /* Cuanto cuesta el test EN SI, para separarlo del efecto que tiene sobre la
     * trayectoria de la busqueda: son dos cosas distintas y conviene no
     * confundirlas al decir «es caro». */
    const double tg0 = ahora_s();
    if (!lin_valida || gs.empty() || m_uso > 0) { t_gtest += ahora_s()-tg0; return false; }
    const int na = A.nb_cols();
    if (na <= nb_var || (int)z.size() != na) { t_gtest += ahora_s()-tg0; return false; }
    if (box.size() != nb_var || caja_lin.size() != nb_var) { t_gtest += ahora_s()-tg0; return false; }

    /* AUTOAPAGADO. El test gana donde acierta mucho (17 % en `alkylbis` y
     * `launch`, que mejoran) y pierde donde acierta poco (0.8 % en `ex5_3_2`,
     * 0.2 % en `house`, que empeoran). `m` no los separa —`ex5_3_2` tiene m=49 y
     * pierde—, pero la tasa de acierto si, y se mide sola: tras un periodo de
     * prueba, si no rinde se apaga para el resto de la corrida. Es la misma
     * logica adaptativa que usa ACID para elegir cuantas variables rebanar. */
    /* Apagado por omision: probado, sale peor. `launch` tiene una tasa temprana
     * baja aunque la global sea 16.8 %, asi que se autoapaga y termina en 832
     * celdas contra las 552 de base; y apagarlo no restaura el comportamiento
     * base, porque las primeras llamadas ya desviaron la trayectoria. */
    static const long calentar = getenv("DFBH_GWARM") ? atol(getenv("DFBH_GWARM")) : 0;
    static const double gmin = getenv("DFBH_GMIN") ? atof(getenv("DFBH_GMIN")) : 0.05;
    if (g_apagado) { t_gtest += ahora_s()-tg0; return false; }

    ++n_gtest;
    if (calentar > 0 && n_gtest == calentar) {
        const long apl = n_gtest - n_gfuera;
        if (apl < 50 || (double)n_gvacio/apl < gmin) g_apagado = true;
    }
    /* Los `gamma` valen para la caja en la que se linealizo. Si la de ahora no
     * esta contenida en aquella, la relajacion no es valida aca y el test
     * NO se aplica: el objeto se comparte en todo el arbol y estos `gamma`
     * pueden venir de otra rama. */
    for (int i = 0; i < nb_var; ++i)
        if (!caja_lin[i].is_superset(box[i])) { ++n_gfuera; t_gtest += ahora_s()-tg0; return false; }

    /* DISTANCIA a la caja de linealizacion, para saber si la herencia PADRE->HIJO
     * serviria. Un padre esta a una biseccion: una variable al doble y el resto
     * igual. Un ancestro lejano esta mucho mas lejos. `DFBH_GDIST=1` reporta la
     * tasa de acierto por cubeta de distancia; si la cubeta «una biseccion»
     * acierta poco, heredar del padre no vale la pena por mas que la contencion
     * pase siempre. */
    static const bool gdist = (getenv("DFBH_GDIST") != NULL);
    int cub = -1;
    if (gdist) {
        double r = 1.0; int ndif = 0;
        for (int i = 0; i < nb_var; ++i) {
            const double a = caja_lin[i].diam(), b = box[i].diam();
            if (!(a > 0.0) || !(b > 0.0) || a >= 1e18 || b >= 1e18) continue;
            if (a > b*(1.0+1e-9)) { ++ndif; if (a/b > r) r = a/b; }
        }
        /* cubetas: 0 = misma caja, 1 = una biseccion (una variable, factor ~2),
         * 2 = una variable pero mas de 2, 3 = pocas variables, 4 = muchas */
        cub = (ndif == 0) ? 0
            : (ndif == 1 && r <= 2.5) ? 1
            : (ndif == 1) ? 2
            : (ndif <= 3) ? 3 : 4;
    }
    struct GRep {
        static void anota(int cub, bool mata) {
            static long n[5] = {0,0,0,0,0}, v[5] = {0,0,0,0,0}, tot = 0;
            if (cub < 0) return;
            ++n[cub]; if (mata) ++v[cub];
            if ((++tot % 200) == 0) {
                static const char* et[5] = {"misma caja","1 var x<=2.5","1 var x>2.5","2-3 vars","4+ vars"};
                fprintf(stderr, "[gdist] aplicables=%ld |", tot);
                for (int q = 0; q < 5; ++q)
                    if (n[q]) fprintf(stderr, "  %s: %.1f%% de %ld", et[q], 100.0*v[q]/n[q], n[q]);
                fprintf(stderr, "\n");
            }
        }
    };

    IntervalVector zs(z);
    for (int i = 0; i < nb_var; ++i) zs[i] = box[i];

    /* La parte `b` tambien se aprieta con el trozo, y esto es gratis en pivotes.
     *
     * `b_j = A_j.x`, asi que sobre una caja mas angosta su rango es mas angosto:
     * se intersecta la cota guardada —deducida para la caja ancha, valida pero
     * floja— con la evaluacion por intervalos de la fila sobre `box`. El test
     * queda estrictamente mas fuerte por `m*nx` operaciones de intervalo y
     * ningun pivote. Usar la cota vieja desperdiciaba justamente la informacion
     * del trozo en `m` de las `na` componentes. */
    const int m = A.nb_rows();
    for (int j = 0; j < m; ++j) {
        Interval enc(0.0);
        for (int c = 0; c < nb_var; ++c) enc += A[j][c] * box[c];
        const Interval t = z[nb_var+j] & enc;
        if (t.is_empty()) { ++n_gvacio; GRep::anota(cub,true); t_gtest += ahora_s()-tg0; return true; }
        zs[nb_var+j] = t;
    }

    const int nc = (int)gs.size();
    for (int q = 0; q < nc; ++q) {
        if (!gs_ok[q]) continue;
        Interval suma(0.0);
        for (int c = 0; c < na; ++c) suma += gs[q][c] * zs[c];
        if (!suma.contains(0.0)) { ++n_gvacio; GRep::anota(cub,true); t_gtest += ahora_s()-tg0; return true; }
    }
    GRep::anota(cub,false);
    t_gtest += ahora_s()-tg0;
    return false;
}


void CtcDFBHull::contract(IntervalVector& box) {
    ContractContext context(box);
    contract(box, context);
}

long BxpDFBSens::id() { static const long i = next_id(); return i; }

/* Tope del pozo: `DFBH_CORTES_TOPE` veces las filas afines (omision 1). */
static int tope_cortes(int m_afin) {
    static const double f = getenv("DFBH_CORTES_TOPE") ? atof(getenv("DFBH_CORTES_TOPE")) : 1.0;
    return std::max(1, (int)(f * m_afin + 0.5));
}
int CtcDFBHull::modo_cortes = getenv("DFBH_CORTES") ? atoi(getenv("DFBH_CORTES")) : 0;
long CtcDFBHull::n_cortes_gen = 0, CtcDFBHull::n_cortes_filas = 0, CtcDFBHull::n_cortes_llam = 0;
long BxpDFBCortes::id() { static const long i = next_id(); return i; }
Bxp* BxpDFBCortes::copy(const IntervalVector& box, const BoxProperties& prop) const {
    /* Las hijas estan contenidas en la caja del padre, asi que los cortes
     * siguen siendo validos; solo se heredan con `DFBH_CORTES=2`. */
    if (CtcDFBHull::modo_cortes == 2) return new BxpDFBCortes(*this);
    if (CtcDFBHull::modo_cortes == 3) {
        /* la hija aplica los del padre: los que el padre aplicaba y los que
         * genero; los suyos propios iran a sus hijas */
        BxpDFBCortes* h = new BxpDFBCortes();
        h->tope = tope;
        h->pozo = pozo;
        h->pozo.insert(h->pozo.end(), pendientes.begin(), pendientes.end());
        /* sin tope el pozo acumularia los cortes de todos los ancestros (todos
         * siguen valiendo): se quedan los mas nuevos */
        if (tope > 0 && h->pozo.size() > tope)
            h->pozo.erase(h->pozo.begin(), h->pozo.begin() + (h->pozo.size() - tope));
        return h;
    }
    return new BxpDFBCortes();
}

long CtcDFBHull::n_hib_caja = 0, CtcDFBHull::n_hib_sc = 0;
Ctc* CtcDFBHull::hc4_diag = NULL;
bool CtcDFBHull::motor_caja = getenv("DFBH_MOTOR") && std::string(getenv("DFBH_MOTOR")) == "caja";

__attribute__((noinline)) void CtcDFBHull::pasada_ab_caja(IntervalVector& b) { una_pasada(b); }
__attribute__((noinline)) void CtcDFBHull::pasada_ab_sc(IntervalVector& b) { una_pasada(b); }

/* CAJA NEGRA DIRECTA (`DFBH_CAJA_DIRECTA=1`, 2026-10-16). El `LPSolver` de Ibex
 * guarda las cotas de las n variables tambien como n filas (para los duales de
 * Neumaier-Shcherbina), asi que SoPlex trabaja con m+n filas en vez de m. Aca
 * la caja negra usa su propio SoPlex, con los mismos parametros que
 * `mylinearsolver` (se copian todos), las mismas filas y lados de la relajacion
 * y las cotas solo como cotas de columna. Las cotas se cargan solo si cambiaron. */
/* CAJA NEGRA CON SOPLEX 7 (`DFBH_CAJA_SOPLEX7=1`, 2026-10-17): como la caja directa
 * (mismas filas y lados, cotas solo de columna, filas libres vacias, cotas
 * cargadas solo si cambian), pero resuelta por SoPlex 7.1.6 (ver ibex_Caja7.h). */
DFBSimplex::Status CtcDFBHull::resolver_caja7(int k, bool maximizar, const IntervalVector& z, bool frio) {
    it_caja = 0;
    caja_inf_sin_rayo = false;
    const int m = A.nb_rows();
    if (!caja7) caja7 = new Caja7();
    if (caja7_ep != n_m_llam || (int)caja_cotas.size() != nb_var) {
        const Matrix R = mylinearsolver.rows();
        const IntervalVector LR = mylinearsolver.lhs_rhs();
        std::vector<double> lb(nb_var), ub(nb_var), lhs(m), rhs(m), a((size_t)m * nb_var);
        for (int i = 0; i < nb_var; ++i) { lb[i] = z[i].lb(); ub[i] = z[i].ub(); }
        for (int j = 0; j < m; ++j) {
            lhs[j] = LR[nb_var + j].lb(); rhs[j] = LR[nb_var + j].ub();
            for (int i = 0; i < nb_var; ++i) a[(size_t)j * nb_var + i] = R[nb_var + j][i];
        }
        soplex::SoPlex& lp = *mylinearsolver.mysoplex;
        caja7->cargar(nb_var, m, &lb[0], &ub[0], &lhs[0], &rhs[0], &a[0],
                      lp.realParam(soplex::SoPlex::FEASTOL), lp.intParam(soplex::SoPlex::ITERLIMIT));
        caja_cotas = IntervalVector(nb_var);
        for (int i = 0; i < nb_var; ++i) caja_cotas[i] = z[i];
        caja7_ep = n_m_llam;
    } else {
        for (int i = 0; i < nb_var; ++i)
            if (!(caja_cotas[i] == z[i])) { caja7->cota(i, z[i].lb(), z[i].ub()); caja_cotas[i] = z[i]; }
    }
    lam_caja.assign(m, 0.0);
    x_caja.assign(nb_var, 0.0);
    const int r = caja7->resolver(k, maximizar, frio, &lam_caja[0], &x_caja[0], &it_caja);
    if (r == 0) return DFBSimplex::OPTIMAL;
    if (r == 1) return DFBSimplex::INFEASIBLE;
    if (r == 2) caja_inf_sin_rayo = true;
    lam_caja.assign(m, 0.0);
    return DFBSimplex::SINGULAR;
}

DFBSimplex::Status CtcDFBHull::resolver_cajahighs(int k, bool maximizar, const IntervalVector& z, bool frio) {
    it_caja = 0;
    caja_inf_sin_rayo = false;
    const int m = A.nb_rows();
    if (!cajah) cajah = new CajaHighs();
    if (cajah_ep != n_m_llam || (int)caja_cotas.size() != nb_var) {
        const Matrix R = mylinearsolver.rows();
        const IntervalVector LR = mylinearsolver.lhs_rhs();
        std::vector<double> lb(nb_var), ub(nb_var), lhs(m), rhs(m), a((size_t)m * nb_var);
        for (int i = 0; i < nb_var; ++i) { lb[i] = z[i].lb(); ub[i] = z[i].ub(); }
        for (int j = 0; j < m; ++j) {
            lhs[j] = LR[nb_var + j].lb(); rhs[j] = LR[nb_var + j].ub();
            for (int i = 0; i < nb_var; ++i) a[(size_t)j * nb_var + i] = R[nb_var + j][i];
        }
        soplex::SoPlex& lp = *mylinearsolver.mysoplex;
        cajah->cargar(nb_var, m, &lb[0], &ub[0], &lhs[0], &rhs[0], &a[0],
                      lp.realParam(soplex::SoPlex::FEASTOL),
                      !(getenv("DFBH_CAJA_SIMP") && std::string(getenv("DFBH_CAJA_SIMP")) == "0"));
        caja_cotas = IntervalVector(nb_var);
        for (int i = 0; i < nb_var; ++i) caja_cotas[i] = z[i];
        cajah_ep = n_m_llam;
    } else {
        for (int i = 0; i < nb_var; ++i)
            if (!(caja_cotas[i] == z[i])) { cajah->cota(i, z[i].lb(), z[i].ub()); caja_cotas[i] = z[i]; }
    }
    lam_caja.assign(m, 0.0);
    x_caja.assign(nb_var, 0.0);
    const int r = cajah->resolver(k, maximizar, frio, &lam_caja[0], &x_caja[0], &it_caja);
    if (r == 0) return DFBSimplex::OPTIMAL;
    if (r == 1) return DFBSimplex::INFEASIBLE;
    if (r == 2) caja_inf_sin_rayo = true;
    lam_caja.assign(m, 0.0);
    return DFBSimplex::SINGULAR;
}

DFBSimplex::Status CtcDFBHull::resolver_cajagurobi(int k, bool maximizar, const IntervalVector& z, bool frio) {
    it_caja = 0;
    caja_inf_sin_rayo = false;
    const int m = A.nb_rows();
    if (!cajag) cajag = new CajaGurobi();
    if (cajag_ep != n_m_llam || (int)caja_cotas.size() != nb_var) {
        const Matrix R = mylinearsolver.rows();
        const IntervalVector LR = mylinearsolver.lhs_rhs();
        std::vector<double> lb(nb_var), ub(nb_var), lhs(m), rhs(m), a((size_t)m * nb_var);
        for (int i = 0; i < nb_var; ++i) { lb[i] = z[i].lb(); ub[i] = z[i].ub(); }
        for (int j = 0; j < m; ++j) {
            lhs[j] = LR[nb_var + j].lb(); rhs[j] = LR[nb_var + j].ub();
            for (int i = 0; i < nb_var; ++i) a[(size_t)j * nb_var + i] = R[nb_var + j][i];
        }
        soplex::SoPlex& lp = *mylinearsolver.mysoplex;
        cajag->cargar(nb_var, m, &lb[0], &ub[0], &lhs[0], &rhs[0], &a[0],
                      lp.realParam(soplex::SoPlex::FEASTOL),
                      !(getenv("DFBH_CAJA_SIMP") && std::string(getenv("DFBH_CAJA_SIMP")) == "0"));
        caja_cotas = IntervalVector(nb_var);
        for (int i = 0; i < nb_var; ++i) caja_cotas[i] = z[i];
        cajag_ep = n_m_llam;
    } else {
        for (int i = 0; i < nb_var; ++i)
            if (!(caja_cotas[i] == z[i])) { cajag->cota(i, z[i].lb(), z[i].ub()); caja_cotas[i] = z[i]; }
    }
    lam_caja.assign(m, 0.0);
    x_caja.assign(nb_var, 0.0);
    const int r = cajag->resolver(k, maximizar, frio, &lam_caja[0], &x_caja[0], &it_caja);
    if (r == 0) return DFBSimplex::OPTIMAL;
    if (r == 1) return DFBSimplex::INFEASIBLE;
    if (r == 2) caja_inf_sin_rayo = true;
    lam_caja.assign(m, 0.0);
    return DFBSimplex::SINGULAR;
}

DFBSimplex::Status CtcDFBHull::resolver_cajaspx(int k, bool maximizar, const IntervalVector& z, bool frio) {
    it_caja = 0;
    caja_inf_sin_rayo = false;
    const int m = A.nb_rows();
    if (!cajas) cajas = new CajaSPx();
    if (cajas_ep != n_m_llam || (int)caja_cotas.size() != nb_var) {
        const Matrix R = mylinearsolver.rows();
        const IntervalVector LR = mylinearsolver.lhs_rhs();
        std::vector<double> lb(nb_var), ub(nb_var), lhs(m), rhs(m), a((size_t)m * nb_var);
        for (int i = 0; i < nb_var; ++i) { lb[i] = z[i].lb(); ub[i] = z[i].ub(); }
        for (int j = 0; j < m; ++j) {
            lhs[j] = LR[nb_var + j].lb(); rhs[j] = LR[nb_var + j].ub();
            for (int i = 0; i < nb_var; ++i) a[(size_t)j * nb_var + i] = R[nb_var + j][i];
        }
        soplex::SoPlex& lp = *mylinearsolver.mysoplex;
        cajas->cargar(nb_var, m, &lb[0], &ub[0], &lhs[0], &rhs[0], &a[0],
                      lp.realParam(soplex::SoPlex::FEASTOL),
                      !(getenv("DFBH_CAJA_SCALER") && std::string(getenv("DFBH_CAJA_SCALER")) == "0"));
        caja_cotas = IntervalVector(nb_var);
        for (int i = 0; i < nb_var; ++i) caja_cotas[i] = z[i];
        cajas_ep = n_m_llam;
    } else {
        for (int i = 0; i < nb_var; ++i)
            if (!(caja_cotas[i] == z[i])) { cajas->cota(i, z[i].lb(), z[i].ub()); caja_cotas[i] = z[i]; }
    }
    lam_caja.assign(m, 0.0);
    x_caja.assign(nb_var, 0.0);
    const int r = cajas->resolver(k, maximizar, frio, &lam_caja[0], &x_caja[0], &it_caja);
    if (r == 0) return DFBSimplex::OPTIMAL;
    if (r == 1) return DFBSimplex::INFEASIBLE;
    if (r == 2) caja_inf_sin_rayo = true;
    lam_caja.assign(m, 0.0);
    return DFBSimplex::SINGULAR;
}

DFBSimplex::Status CtcDFBHull::resolver_caja_directa(int k, bool maximizar, const IntervalVector& z, bool frio) {
    it_caja = 0;
    caja_inf_sin_rayo = false;
    static const bool reloj = DFB_ENV("DFBH_CAJA_RELOJ");
    const double cr0 = reloj ? ahora_s() : 0.0;
    const int m = A.nb_rows();
    soplex::SoPlex& lp = *mylinearsolver.mysoplex;
    /* `DFBH_CAJA_NUEVO=1`: un objeto SoPlex nuevo por LP (diagnostico del estado que SoPlex arrastra) */
    if (sx_directa && DFB_ENV("DFBH_CAJA_NUEVO")) { delete sx_directa; sx_directa = nullptr; sx_directa_ep = -1; }
    if (!sx_directa) sx_directa = new soplex::SoPlex();
    soplex::SoPlex& sx = *sx_directa;
    static const bool reb = DFB_ENV("DFBH_CAJA_DIRECTA_REB");
    if (reb || sx_directa_ep != n_m_llam || sx.numRowsReal() != m || sx.numColsReal() != nb_var) {
        /* los parametros van DESPUES de `clearLPReal`: el clear devuelve el LP
         * al sentido de fabrica (maximizar) sin volver a aplicar OBJSENSE */
        sx.clearLPReal();
        for (int i = 0; i < soplex::SoPlex::INTPARAM_COUNT; ++i)
            sx.setIntParam((soplex::SoPlex::IntParam) i, lp.intParam((soplex::SoPlex::IntParam) i), true);
        for (int i = 0; i < soplex::SoPlex::REALPARAM_COUNT; ++i)
            sx.setRealParam((soplex::SoPlex::RealParam) i, lp.realParam((soplex::SoPlex::RealParam) i), true);
        for (int i = 0; i < soplex::SoPlex::BOOLPARAM_COUNT; ++i)
            sx.setBoolParam((soplex::SoPlex::BoolParam) i, lp.boolParam((soplex::SoPlex::BoolParam) i), true);
        sx.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MINIMIZE, true);
        /* `DFBH_CAJA_TOPE=igual`: el mismo tope de iteraciones que SC (20m+200)
         * en vez del de Ibex (100), para comparar en igualdad de condiciones. */
        static const bool tope_igual = getenv("DFBH_CAJA_TOPE") && std::string(getenv("DFBH_CAJA_TOPE")) == "igual";
        if (tope_igual) sx.setIntParam(soplex::SoPlex::ITERLIMIT, 20*m + 200, true);
        if (DFB_ENV("DFBH_CAJA_VERB")) sx.setIntParam(soplex::SoPlex::VERBOSITY, soplex::SoPlex::VERBOSITY_NORMAL, true);
        soplex::LPColSetReal cols(nb_var);
        soplex::DSVectorReal vacio(0);
        for (int i = 0; i < nb_var; ++i) cols.add(0.0, z[i].lb(), vacio, z[i].ub());
        sx.addColsReal(cols);
        soplex::LPRowSetReal rows(m);
        const Matrix R = mylinearsolver.rows();
        const IntervalVector LR = mylinearsolver.lhs_rhs();
        for (int j = 0; j < m; ++j) {
            soplex::DSVectorReal fila(nb_var);
            for (int i = 0; i < nb_var; ++i) if (R[nb_var + j][i] != 0.0) fila.add(i, R[nb_var + j][i]);
            /* `DFBH_CAJA_FILTRA=t`: la fila se deja libre (lados infinitos) si sus coeficientes
             * no nulos abarcan mas de t ordenes de magnitud; sigue siendo una relajacion valida
             * y el multiplicador de esa fila sale cero. */
            static const double filtra = getenv("DFBH_CAJA_FILTRA") ? atof(getenv("DFBH_CAJA_FILTRA")) : 0.0;
            bool libre = false;
            if (filtra > 0) {
                double amax = 0, amin = std::numeric_limits<double>::infinity();
                for (int i = 0; i < nb_var; ++i) { const double v = std::fabs(R[nb_var + j][i]); if (v > 0) { amax = std::max(amax, v); amin = std::min(amin, v); } }
                libre = amax > 0 && std::log10(amax / amin) > filtra;
            }
            /* `DFBH_CAJA_SINLIBRES=1`: una fila libre (los dos lados infinitos para SoPlex) no
             * restringe nada; se le pasa vacia para que sus coeficientes no entren al LP. */
            static const bool sin_libres = DFB_ENV("DFBH_CAJA_SINLIBRES");
            const double inf_sx = sx.realParam(soplex::SoPlex::INFTY);
            if (sin_libres && LR[nb_var + j].lb() <= -inf_sx && LR[nb_var + j].ub() >= inf_sx) {
                soplex::DSVectorReal vacia(0);
                rows.add(-soplex::infinity, vacia, soplex::infinity);
            }
            else if (libre) rows.add(-soplex::infinity, fila, soplex::infinity);
            else rows.add(LR[nb_var + j].lb(), fila, LR[nb_var + j].ub());
            if (false) {
                static long n_raras = 0;
                const double lo = LR[nb_var + j].lb(), up = LR[nb_var + j].ub();
                bool rara = !(lo <= up) || std::isnan(lo) || std::isnan(up) || LR[nb_var + j].is_empty();
                for (int i = 0; i < nb_var; ++i) if (std::isnan(R[nb_var + j][i]) || std::isinf(R[nb_var + j][i])) rara = true;
                if (rara && ++n_raras <= 5)
                    fprintf(stderr, "[cajafila] fila %d lados [%g, %g] vacia=%d\n", j, lo, up, (int)LR[nb_var + j].is_empty());
            }
        }
        if (DFB_ENV("DFBH_CAJA_DBG")) {
            static long n_cr = 0;
            for (int i = 0; i < nb_var; ++i)
                if ((!(z[i].lb() <= z[i].ub()) || std::isnan(z[i].lb()) || std::isnan(z[i].ub())) && ++n_cr <= 5)
                    fprintf(stderr, "[cajacol] columna %d cotas [%g, %g]\n", i, z[i].lb(), z[i].ub());
        }
        sx.addRowsReal(rows);
        caja_cotas = IntervalVector(nb_var);
        for (int i = 0; i < nb_var; ++i) caja_cotas[i] = z[i];
        sx_directa_ep = n_m_llam;
    } else {
        for (int i = 0; i < nb_var; ++i)
            if (!(caja_cotas[i] == z[i])) { sx.changeBoundsReal(i, z[i].lb(), z[i].ub()); caja_cotas[i] = z[i]; }
    }
    sx.changeObjReal(k, maximizar ? -1.0 : 1.0);
    if (frio || DFB_ENV("DFBH_CAJA_FRIO")) sx.clearBasis();   /* DFBH_CAJA_FRIO=1: la caja directa parte en frio en cada LP */
    /* `DFBH_CAJA_DUMP=dir`: guarda cada LP ANTES de resolverlo (MPS, base de
     * partida y parametros) y conserva los que SoPlex no termina (hasta 40) y
     * uno de cada 2000 que si (control). Despues de una falla el LP interno de
     * SoPlex queda inconsistente y ya no se puede escribir. */
    static const char* dump_dir = getenv("DFBH_CAJA_DUMP");
    std::string cur;
    if (dump_dir) {
        /* formato propio, exacto (hexadecimal): n m k max / cotas de columnas / filas (lhs rhs coeficientes) */
        FILE* t = fopen((std::string(dump_dir) + "/cur.txtlp").c_str(), "w");
        if (t) {
            const Matrix R = mylinearsolver.rows();
            const IntervalVector LR = mylinearsolver.lhs_rhs();
            fprintf(t, "%d %d %d %d\n", nb_var, m, k, (int)maximizar);
            for (int i = 0; i < nb_var; ++i) fprintf(t, "%a %a\n", z[i].lb(), z[i].ub());
            for (int j = 0; j < m; ++j) {
                fprintf(t, "%a %a", LR[nb_var + j].lb(), LR[nb_var + j].ub());
                for (int i = 0; i < nb_var; ++i) fprintf(t, " %a", R[nb_var + j][i]);
                fprintf(t, "\n");
            }
            fclose(t);
        }
    }
    if (dump_dir) {
        cur = std::string(dump_dir) + "/cur";
        static long n_esc_mal = 0, n_esc = 0;
        ++n_esc;
        try {
            sx.writeFileReal((cur + ".mps").c_str());
            sx.saveSettingsFile((cur + ".set").c_str());
            std::remove((cur + ".bas").c_str());
            if (sx.hasBasis()) sx.writeBasisFile((cur + ".bas").c_str());
        } catch (const soplex::SPxException& ex) {
            std::remove((cur + ".mps").c_str());
            if (n_esc_mal < 3) fprintf(stderr, "[cajadump] excepcion al escribir: %s\n", ex.what().c_str());
            ++n_esc_mal;
        } catch (...) {
            std::remove((cur + ".mps").c_str());
            if (++n_esc_mal % 1000 == 1) fprintf(stderr, "[cajadump] LP %ld no se pudo escribir (%ld de %ld): LP interno de SoPlex inconsistente\n", n_esc, n_esc_mal, n_esc);
        }
    }
    const double cr1 = reloj ? ahora_s() : 0.0;
    static const int modo_red = getenv("DFBH_CAJA_REDONDEO") ? atoi(getenv("DFBH_CAJA_REDONDEO")) : -1;
    static int red_vista = -2;
    const int red0 = fegetround();
    if (red_vista == -2 && DFB_ENV("DFBH_CAJA_DBG")) {
        red_vista = red0;
        fprintf(stderr, "[cajad] modo de redondeo al llamar a SoPlex: %s\n",
                red0 == FE_TONEAREST ? "mas cercano" : red0 == FE_UPWARD ? "hacia arriba" : red0 == FE_DOWNWARD ? "hacia abajo" : "hacia cero");
    }
    if (modo_red == 1) fesetround(FE_TONEAREST);   /* `DFBH_CAJA_REDONDEO=1`: SoPlex con redondeo al mas cercano */
    sx.optimize();
    sx.ignoreUnscaledViolations();   /* como `LPSolver::minimize` */
    if (modo_red == 1) fesetround(red0);
    const soplex::SPxSolver::Status st = sx.status();
    if (dump_dir) {
        static long n_mal = 0, n_bien = 0, n_tot = 0; ++n_tot;
        const bool mal = !(st == soplex::SPxSolver::OPTIMAL || st == soplex::SPxSolver::INFEASIBLE);
        if ((mal && n_mal < 40) || (!mal && n_tot % 2000 == 0 && n_bien < 20)) {
            const long id = mal ? n_mal++ : n_bien++;
            char base[512];
            snprintf(base, sizeof base, "%s/%s%03ld", dump_dir, mal ? "mal" : "bien", id);
            const std::string f(base);
            std::rename((std::string(dump_dir) + "/cur.txtlp").c_str(), (f + ".txtlp").c_str());
            std::rename((cur + ".mps").c_str(), (f + ".mps").c_str());
            std::rename((cur + ".set").c_str(), (f + ".set").c_str());
            std::rename((cur + ".bas").c_str(), (f + ".bas").c_str());
            FILE* g = fopen((f + ".txt").c_str(), "w");
            if (g) { fprintf(g, "estado %d iteraciones %d k %d max %d\n", (int)st, sx.numIterations(), k, (int)maximizar); fclose(g); }
        }
    }
fin_dump:
    const double cr2 = reloj ? ahora_s() : 0.0;
    if (reloj) { cr_t[0] += cr1 - cr0; cr_t[1] += cr2 - cr1; cr_t[3] += cr2 - cr1; ++cr_n; }
    it_caja = sx.numIterations();
    lam_caja.assign(m, 0.0);
    DFBSimplex::Status res = DFBSimplex::SINGULAR;
    if (st == soplex::SPxSolver::OPTIMAL) {
        soplex::DVectorReal y(m), x(nb_var);
        if (sx.getDualReal(y) && sx.getPrimalReal(x)) {
            for (int j = 0; j < m; ++j) lam_caja[j] = y[j];
            x_caja.assign(nb_var, 0.0);
            for (int i = 0; i < nb_var; ++i) x_caja[i] = x[i];
            res = DFBSimplex::OPTIMAL;
        }
    } else if (st == soplex::SPxSolver::INFEASIBLE) {
        soplex::DVectorReal f(m);
        if (sx.getDualFarkasReal(f)) {
            for (int j = 0; j < m; ++j) lam_caja[j] = f[j];
            res = DFBSimplex::INFEASIBLE;
        } else caja_inf_sin_rayo = true;
    }
    sx.changeObjReal(k, 0.0);
    if (reloj) cr_t[2] += ahora_s() - cr2;
    if (DFB_ENV("DFBH_CAJA_DBG")) {
        /* estados de SoPlex en la caja directa: optimo, infactible con/sin rayo, no acotado, tope de iteraciones, otro */
        static long c[6] = {0,0,0,0,0,0}, tot = 0;
        int i = 5;
        if (st == soplex::SPxSolver::OPTIMAL) i = 0;
        else if (st == soplex::SPxSolver::INFEASIBLE) i = (res == DFBSimplex::INFEASIBLE) ? 1 : 2;
        else if (st == soplex::SPxSolver::UNBOUNDED) i = 3;
        else if (st == soplex::SPxSolver::ABORT_ITER) i = 4;
        ++c[i]; ++tot;
        static std::map<int,long> hist; ++hist[(int)st];
        if (tot % 2000 == 0) {
            fprintf(stderr, "[cajad] LPs=%ld optimo=%ld inf_con_rayo=%ld inf_sin_rayo=%ld noacotado=%ld tope_iter=%ld otro=%ld | estados SoPlex:",
                    tot, c[0], c[1], c[2], c[3], c[4], c[5]);
            for (std::map<int,long>::const_iterator it = hist.begin(); it != hist.end(); ++it) fprintf(stderr, " %d:%ld", it->first, it->second);
            fprintf(stderr, "\n");
        }
    }
    return res;
}

DFBSimplex::Status CtcDFBHull::resolver_caja(int k, bool maximizar, const IntervalVector& z, bool frio) {
    it_caja = 0;
    caja_inf_sin_rayo = false;
    /* Tolerancias de SoPlex en la caja negra (`DFBH_CAJA_TOL`, omision 1e-6, la
     * de fabrica de SoPlex). Con la de Ibex (1e-9) SoPlex declara infactibles
     * LPs factibles cuando la relajacion afin tiene casi-igualdades de ancho
     * ~1e-14 en cajas grandes (Brown); el rigor lo pone el certificado. */
    static const double tol_caja = getenv("DFBH_CAJA_TOL") ? atof(getenv("DFBH_CAJA_TOL")) : 1e-6;
    static bool tol_puesta = false;
    if (!tol_puesta) {
        mylinearsolver.set_tolerance(tol_caja); tol_puesta = true;
        /* `DFBH_CAJA_SIMP=0` apaga el presolve de SoPlex, `DFBH_CAJA_SCALER=0`
         * el escalado: parametros del solver de caja negra, para elegir su
         * mejor configuracion. */
        if (getenv("DFBH_CAJA_SIMP") && std::string(getenv("DFBH_CAJA_SIMP")) == "0")
            mylinearsolver.mysoplex->setIntParam(soplex::SoPlex::SIMPLIFIER, soplex::SoPlex::SIMPLIFIER_OFF);
        if (getenv("DFBH_CAJA_SCALER") && std::string(getenv("DFBH_CAJA_SCALER")) == "0")
            mylinearsolver.mysoplex->setIntParam(soplex::SoPlex::SCALER, soplex::SoPlex::SCALER_OFF);
        if (getenv("DFBH_CAJA_REP") && std::string(getenv("DFBH_CAJA_REP")) == "col")
            mylinearsolver.mysoplex->setIntParam(soplex::SoPlex::REPRESENTATION, soplex::SoPlex::REPRESENTATION_COLUMN);
        if (getenv("DFBH_CAJA_REP") && std::string(getenv("DFBH_CAJA_REP")) == "row")
            mylinearsolver.mysoplex->setIntParam(soplex::SoPlex::REPRESENTATION, soplex::SoPlex::REPRESENTATION_ROW);
        if (getenv("DFBH_CAJA_ALG") && std::string(getenv("DFBH_CAJA_ALG")) == "primal")
            mylinearsolver.mysoplex->setIntParam(soplex::SoPlex::ALGORITHM, soplex::SoPlex::ALGORITHM_PRIMAL);
    }
    if (DFB_ENV("DFBH_CAJA_SOPLEX7")) return resolver_caja7(k, maximizar, z, frio);
    if (DFB_ENV("DFBH_CAJA_HIGHS")) return resolver_cajahighs(k, maximizar, z, frio);
    if (DFB_ENV("DFBH_CAJA_GUROBI")) return resolver_cajagurobi(k, maximizar, z, frio);
    if (DFB_ENV("DFBH_CAJA_SPX")) return resolver_cajaspx(k, maximizar, z, frio);
    static const bool directa = DFB_ENV("DFBH_CAJA_DIRECTA");
    if (directa) return resolver_caja_directa(k, maximizar, z, frio);

    static const bool reloj = DFB_ENV("DFBH_CAJA_RELOJ");
    static const bool cotas_cambio = getenv("DFBH_CAJA_COTAS") && std::string(getenv("DFBH_CAJA_COTAS")) == "cambio";
    const double cr0 = reloj ? ahora_s() : 0.0;
    if (cotas_cambio && caja_cotas_ep == n_m_llam && (int)caja_cotas.size() == nb_var) {
        for (int i = 0; i < nb_var; ++i)
            if (!(caja_cotas[i] == z[i])) { mylinearsolver.set_bounds(i, z[i]); caja_cotas[i] = z[i]; ++cr_cotas; }
    } else {
        IntervalVector xb(nb_var);
        for (int i = 0; i < nb_var; ++i) xb[i] = z[i];
        mylinearsolver.set_bounds(xb);
        cr_cotas += nb_var;
        if (cotas_cambio) { caja_cotas = xb; caja_cotas_ep = n_m_llam; }
    }
    mylinearsolver.set_cost(k, maximizar ? -1.0 : 1.0);
    /* en frio: se borra la base que SoPlex guarda del LP anterior */
    if (frio) mylinearsolver.mysoplex->clearBasis();
    const double cr1 = reloj ? ahora_s() : 0.0;
    const double ps0 = ph_t_solve, pe0 = ph_t_extrae;
    const LPSolver::Status st = mylinearsolver.minimize();
    const double cr2 = reloj ? ahora_s() : 0.0;
    if (reloj) { cr_t[0] += cr1 - cr0; cr_t[1] += cr2 - cr1; cr_t[3] += ph_t_solve - ps0; cr_t[4] += ph_t_extrae - pe0; ++cr_n; }
    if (reloj && getenv("DFBH_CAJA_STATS") && cr_n % atol(getenv("DFBH_CAJA_STATS")) == 0) {
        std::ostringstream os; mylinearsolver.mysoplex->printSolvingStatistics(os);
        fprintf(stderr, "[cajastats] LP %ld\n%s", cr_n, os.str().c_str());
    }
    it_caja = mylinearsolver.mysoplex->numIterations();
    const int m = A.nb_rows();
    lam_caja.assign(m, 0.0);
    DFBSimplex::Status res = DFBSimplex::SINGULAR;

    /* OJO: los multiplicadores se leen ANTES de devolver el costo a cero.
     * Cambiar el objetivo invalida la solucion de SoPlex, rayo de Farkas
     * incluido (el dual lo copia el LPSolver dentro de `minimize`, el rayo no
     * se pedia desde aca). El LPSolver de Ibex pone las cotas de las variables
     * como las primeras `nb_var` filas (ver `linearize_fresca`): los
     * multiplicadores de la relajacion son las componentes siguientes. */
    if (st == LPSolver::Status::Optimal || st == LPSolver::Status::OptimalProved) {
        const Vector y = mylinearsolver.not_proved_dual_sol();
        if (y.size() >= nb_var + m) {
            for (int j = 0; j < m; ++j) lam_caja[j] = y[nb_var + j];
            const Vector x = mylinearsolver.not_proved_primal_sol();
            x_caja.assign(x.size(), 0.0);
            for (int i = 0; i < x.size(); ++i) x_caja[i] = x[i];
            res = DFBSimplex::OPTIMAL;
        }
    } else if (st == LPSolver::Status::Infeasible || st == LPSolver::Status::InfeasibleProved) {
        soplex::DVectorReal f(mylinearsolver.mysoplex->numRowsReal());
        f.clear();
        if (mylinearsolver.mysoplex->getDualFarkasReal(f) && f.dim() >= nb_var + m) {
            for (int j = 0; j < m; ++j) lam_caja[j] = f[nb_var + j];
            res = DFBSimplex::INFEASIBLE;
        } else caja_inf_sin_rayo = true;
    }
    mylinearsolver.set_cost(k, 0.0);
    if (reloj) cr_t[2] += ahora_s() - cr2;
    if (DFB_ENV("DFBH_CAJA_CMP")) {
        /* Compara con la caja directa en el mismo LP: estados y valor optimo. */
        static long n = 0, dif_st = 0, st_ibex[8] = {0}, st_dir[8] = {0}, dif_obj = 0;
        const std::vector<double> lam0 = lam_caja, x0 = x_caja; const int it0 = it_caja; const bool sr0 = caja_inf_sin_rayo;
        const DFBSimplex::Status r2 = resolver_caja_directa(k, maximizar, z, frio);
        ++n; ++st_ibex[(int)res & 7]; ++st_dir[(int)r2 & 7];
        if (r2 != res) ++dif_st;
        else if (res == DFBSimplex::OPTIMAL && (int)x0.size() > k && (int)x_caja.size() > k
                 && std::fabs(x0[k] - x_caja[k]) > 1e-7 * (1 + std::fabs(x0[k]))) {
            if (++dif_obj <= 4) {
                fprintf(stderr, "[cajacmp-dif] k=%d max=%d z_k=[%g,%g] ibex x_k=%.9g directa x_k=%.9g | m=%d filas sx=%d cols sx=%d filas ibex=%d\n",
                        k, (int)maximizar, z[k].lb(), z[k].ub(), x0[k], x_caja[k], A.nb_rows(), sx_directa->numRowsReal(), sx_directa->numColsReal(), mylinearsolver.nb_rows());
                for (int j = 0; j < sx_directa->numRowsReal() && j < 4; ++j) fprintf(stderr, "   fila %d lados sx [%g,%g] ibex [%g,%g]\n", j, sx_directa->lhsReal(j), sx_directa->rhsReal(j), mylinearsolver.lhs_rhs()[nb_var+j].lb(), mylinearsolver.lhs_rhs()[nb_var+j].ub());
            }
        }
        static const long cada = std::max(1L, atol(getenv("DFBH_CAJA_CMP")));
        if (n % cada == 0) {
            fprintf(stderr, "[cajacmp] LPs=%ld estado distinto=%ld optimo distinto=%ld | ibex:", n, dif_st, dif_obj);
            for (int i = 0; i < 8; ++i) fprintf(stderr, " %ld", st_ibex[i]);
            fprintf(stderr, " | directa:");
            for (int i = 0; i < 8; ++i) fprintf(stderr, " %ld", st_dir[i]);
            fprintf(stderr, "\n");
        }
        lam_caja = lam0; x_caja = x0; it_caja = it0; caja_inf_sin_rayo = sr0;
    }
    if (getenv("DFBH_CAJA_DBG")) {
        static long cnt[8] = {0,0,0,0,0,0,0,0}, tot = 0, rayo = 0;
        ++cnt[(int)st]; ++tot; if (res == DFBSimplex::INFEASIBLE) ++rayo;
        if (tot % 2000 == 0)
            fprintf(stderr, "[caja] estados: Opt=%ld Inf=%ld (con rayo %ld) Unb=%ld MaxIt=%ld Unk=%ld\n",
                    cnt[0]+cnt[1], cnt[2]+cnt[3], rayo, cnt[4], cnt[6], cnt[7]);
    }
    return res;
}

/* CORTES EN X* (`DFBH_CORTES`).
 *
 * Para cada restriccion g_c que el `x*` de un LP viola, una fila XTaylor en
 * la esquina `c` de `caja_lin` mas cercana a `x*`, coordenada a coordenada:
 *
 *     g(x) >= g(c) + sum_i a_i (x_i - c_i),   a_i = J_i.lb si c_i = lb, J_i.ub si no
 *
 * vale en toda la caja porque el signo de `x_i - c_i` esta fijo (J en
 * intervalos sobre `caja_lin`). Para `g <= 0` da `a.x <= -g(c) + a.c`; para
 * `g >= 0`, con la cota superior simetrica, `a.x >= ...`. El lado derecho se
 * evalua en intervalos y se toma el extremo seguro, asi que el corte es
 * riguroso. Se guarda solo si separa a `x*`. */
void CtcDFBHull::generar_cortes(const std::vector<double>& xs) {
    LinearizerAffine2* la = dynamic_cast<LinearizerAffine2*>(&lr);
    if (!la || (int)xs.size() < nb_var || caja_lin.size() != nb_var) return;
    const int tope = tope_cortes(m_afin);
    if ((int)cortes_nuevos.size() >= tope) return;
    const System& sys = la->sistema();
    const int nc = sys.f_ctrs.image_dim();
    IntervalVector xp(nb_var);
    for (int i = 0; i < nb_var; ++i) xp[i] = Interval(xs[i]);
    const IntervalVector gx = sys.f_ctrs.eval_vector(xp);
    if (gx.is_empty()) return;
    for (int c = 0; c < nc && (int)cortes_nuevos.size() < tope; ++c) {
        const CmpOp op = sys.ops[c];
        int lado;
        if ((op == LEQ || op == LT || op == EQ) && gx[c].lb() > 0.0) lado = 0;
        else if ((op == GEQ || op == GT || op == EQ) && gx[c].ub() < 0.0) lado = 1;
        else continue;

        /* Solo el gradiente de esta restriccion, una vez por pasada: el
         * jacobiano completo costaba x3.8 por pivote en `Trigexp1-100`. */
        if (!jac_cortes_ok) { grad_cortes.assign(nc, IntervalVector(1)); grad_cortes_ok.assign(nc, 0); jac_cortes_ok = true; }
        if (!grad_cortes_ok[c]) { grad_cortes[c] = sys.f_ctrs[c].gradient(caja_lin); grad_cortes_ok[c] = 1; }
        const IntervalVector& Jc = grad_cortes[c];
        if (Jc.size() != nb_var) continue;
        IntervalVector esq(nb_var);
        std::vector<double> a(nb_var);
        unsigned long clave = 1469598103934665603UL ^ (unsigned long)(2*c + lado);
        bool ok = true;
        for (int i = 0; i < nb_var && ok; ++i) {
            const Interval& Xi = caja_lin[i];
            const Interval& Jci = Jc[i];
            if (Xi.is_unbounded() || Jci.is_unbounded() || Jci.is_empty()) { ok = false; break; }
            const bool abajo = xs[i] <= Xi.mid();
            const double ci = abajo ? Xi.lb() : Xi.ub();
            esq[i] = Interval(ci);
            a[i] = (abajo == (lado == 0)) ? Jci.lb() : Jci.ub();
            clave = (clave ^ (unsigned long)abajo) * 1099511628211UL;
        }
        if (!ok) continue;
        const Interval gc = sys.f_ctrs[c].eval(esq);
        if (gc.is_empty() || gc.is_unbounded()) continue;
        Interval rhs = -gc;
        for (int i = 0; i < nb_var; ++i) rhs += Interval(a[i]) * esq[i];
        double ax = 0.0;
        for (int i = 0; i < nb_var; ++i) ax += a[i] * xs[i];
        const double viol = (lado == 0) ? ax - rhs.ub() : rhs.lb() - ax;
        if (!(viol > 1e-9 * (1.0 + std::fabs(ax)))) continue;
        bool rep = false;
        for (size_t t = 0; t < cortes_nuevos.size() && !rep; ++t) rep = (cortes_nuevos[t].clave == clave);
        if (rep) continue;
        CorteDFB ct;
        ct.a = a;
        ct.rango = (lado == 0) ? Interval(NEG_INFINITY, rhs.ub()) : Interval(rhs.lb(), POS_INFINITY);
        ct.caja = caja_lin;
        ct.clave = clave;
        cortes_nuevos.push_back(ct);
        ++n_cortes_gen;
    }
}

void CtcDFBHull::add_property(const IntervalVector& init_box, BoxProperties& prop) {
    if (calcular_sens && !prop[BxpDFBSens::id()]) prop.add(new BxpDFBSens());
    if (modo_cortes > 0 && !prop[BxpDFBCortes::id()]) prop.add(new BxpDFBCortes());
}

void CtcDFBHull::contract(IntervalVector& box, ContractContext& context) {
    /* CAPTURA DE CAJAS (`DFBH_CAPTURA=archivo`, 2026-10-18): guarda la caja de
     * entrada de cada llamada (hexadecimal exacto, una por linea) para
     * reproducir las mismas pasadas fuera de la busqueda (`ibexreplay`).
     * `DFBH_CAPTURA_CADA=k` guarda una de cada k; `DFBH_CAPTURA_N` es el tope. */
    {
        static const char* cap = getenv("DFBH_CAPTURA");
        static FILE* fcap = cap ? fopen(cap, "w") : NULL;
        static const long cada = getenv("DFBH_CAPTURA_CADA") ? atol(getenv("DFBH_CAPTURA_CADA")) : 1;
        static const long tope = getenv("DFBH_CAPTURA_N") ? atol(getenv("DFBH_CAPTURA_N")) : 2000;
        static long n_llam = 0, n_cap = 0;
        if (fcap && (n_llam++ % cada) == 0 && n_cap < tope) {
            ++n_cap;
            fprintf(fcap, "%d", box.size());
            for (int i = 0; i < box.size(); ++i) fprintf(fcap, " %a %a", box[i].lb(), box[i].ub());
            fprintf(fcap, "\n");
            fflush(fcap);
        }
    }
    cortes_nuevos.clear();
    contract_nucleo(box, context);
    /* Mezcla de los cortes de esta llamada al pozo de la celda: quedan los
     * viejos que siguen valiendo para la caja y los nuevos, sin repetir, hasta
     * el tope (se descartan los mas viejos). */
    if (modo_cortes > 0 && !box.is_empty()) {
        BxpDFBCortes* pc = (BxpDFBCortes*) context.prop[BxpDFBCortes::id()];
        if (pc) {
            std::vector<CorteDFB> pozo;
            const std::vector<CorteDFB>& base = (modo_cortes == 3) ? pc->pendientes : pc->pozo;
            for (size_t t = 0; t < base.size(); ++t)
                if (base[t].caja.size() == box.size() && base[t].caja.is_superset(box))
                    pozo.push_back(base[t]);
            for (size_t t = 0; t < cortes_nuevos.size(); ++t) {
                bool rep = false;
                for (size_t u = 0; u < pozo.size() && !rep; ++u) rep = (pozo[u].clave == cortes_nuevos[t].clave);
                if (!rep) pozo.push_back(cortes_nuevos[t]);
            }
            const size_t tope = (size_t) tope_cortes(m_afin);
            if (pozo.size() > tope) pozo.erase(pozo.begin(), pozo.begin() + (pozo.size() - tope));
            pc->tope = tope;
            if (modo_cortes == 3) pc->pendientes.swap(pozo);
            else pc->pozo.swap(pozo);
        }
    }
    if (!calcular_sens || box.is_empty()) return;
    BxpDFBSens* p = (BxpDFBSens*) context.prop[BxpDFBSens::id()];
    if (p && sens_pase) {
        p->sens = sens_goal; p->lam_ctr = lam_ctr_goal; p->caja = caja_sens; p->ok = true; p->contrajo = sens_contrajo;
    }
}

void CtcDFBHull::contract_nucleo(IntervalVector& box, ContractContext& context) {
    /* `DFBH_LATIDO=1`: una linea cada 200 llamadas, para ver que crece cuando
     * un nodo no termina. */
    static const bool latido = (getenv("DFBH_LATIDO") != NULL);
    static long n_llam = 0;
    if (latido && (++n_llam % 200) == 0) {
        double per = 0.0; int nfin = 0;
        for (int i = 0; i < box.size(); ++i) if (box[i].diam() < 1e300) { per += box[i].diam(); ++nfin; }
        fprintf(stderr, "[latido] llamadas=%ld nodos=%ld LPs=%ld pivotes=%ld iterlimit=%ld | perimetro=%.6g (%d finitas)\n",
                n_llam, n_nodes, n_lps, n_pivots, DFBSimplex::n_iterlimit, per, nfin);
    }
    if (box.is_unbounded()) return;

    /* REUSO DE LA LINEALIZACION DENTRO DEL PUNTO FIJO.
     *
     * `CtcFixPoint` llama a este contractor varias veces por nodo, y entre esas
     * llamadas la caja solo SE ENCOGE. La relajacion `Abar z = 0` calculada para
     * la caja anterior sigue siendo valida para la nueva, asi que se puede
     * reusar: no hay que re-linealizar, ni recargar los `2n` simplex, y los LPs
     * se re-resuelven con el MISMO objetivo sobre una caja apenas mas chica
     * —el caso donde el §25 mostro que el arranque tibio rinde mas, porque la
     * base previa sigue siendo primal-factible y casi optima.
     *
     * Contrapartida: una linealizacion nueva sobre la caja apretada seria mas
     * FUERTE. Se abarata a cambio de un poliedro mas flojo.
     *
     * `DFBH_RELIN=1` lo activa (por omision se re-linealiza, como antes). */
    static const bool reusar = (getenv("DFBH_RELIN") != NULL);
    /* `DFBH_DUMPNODO=N1,N2,...`: vuelca la caja de entrada de esas llamadas al
     * contractor, para resolverla aparte (estancamientos sin vacios espurios). */
    {
        static std::vector<long> dump_nodos;
        static bool dn_ok = false;
        if (!dn_ok) {
            dn_ok = true;
            if (const char* e = getenv("DFBH_DUMPNODO"))
                for (const char* q = e; *q; ) { dump_nodos.push_back(atol(q)); while (*q && *q != ',') ++q; if (*q) ++q; }
        }
        static long n_llam = 0;
        ++n_llam;
        for (size_t q = 0; q < dump_nodos.size(); ++q)
            if (dump_nodos[q] == n_llam) {
                fprintf(stderr, "[dumpnodo] %ld", n_llam);
                for (int c = 0; c < nb_var; ++c) fprintf(stderr, " %.17g %.17g", box[c].lb(), box[c].ub());
                fprintf(stderr, "\n");
            }
    }
    bool reusado = false;
    if (reusar && lin_valida && caja_lin.size() == nb_var && A.nb_rows() > 0) {
        reusado = true;
        for (int i = 0; i < nb_var && reusado; ++i)
            if (!caja_lin[i].is_superset(box[i])) reusado = false;
    }
    if (reusado) {
        /* solo se refresca la parte x de z; la parte b sigue valiendo */
        for (int i = 0; i < nb_var; ++i) z[i] = box[i];
        ++n_nodes;
        una_pasada_reusada(box);
        context.prop.update(BoxEvent(box, BoxEvent::CONTRACT));
        return;
    }

    static const bool perfil = (getenv("DFBH_PERFIL") != NULL);
    const double pf_t0 = perfil ? ahora_s() : 0.0;
    const int m = linearize(box, context);
    if (perfil) { const double t1 = ahora_s(); pf_lin += t1 - pf_t0; ++pf_n; }
    if (m == -1) {                 /* la relajacion es infactible */
        box.set_empty();
        mylinearsolver.clear_constraints();
        return;
    }
    if (m == 0) return;

    ++n_nodes;

    /* Cuantas veces cambia `m` entre nodos consecutivos: si cambia seguido, la
     * herencia de base ni siquiera es aplicable, porque el conjunto de columnas
     * basicas deja de tener sentido. */
    if (DFB_ENV("DFBH_MCHG")) {
        static int m_prev = -1; static long cambios = 0, total = 0;
        if (m_prev >= 0) { ++total; if (m != m_prev) ++cambios; }
        m_prev = m;
        if (total && (total % 100) == 0)
            fprintf(stderr, "[m] nodos=%ld  m cambio respecto del anterior en %ld (%.1f%%)\n",
                    total, cambios, 100.0*cambios/total);
    }

    /* CONTRACCION A IGUALDAD DE CONDICIONES: FRESCA contra CONGELADA
     * (`DFBH_ABCONG=1`).
     *
     * Las dos ramas parten de la MISMA caja y corren su propio bucle de
     * contraccion hasta punto fijo, o sea que se comparan contracciones
     * completas y no una pasada contra otra. La fresca re-linealiza en cada
     * vuelta; la congelada linealiza una vez y despues solo recalcula las cotas
     * de `b`, sin recargar el simplex. Se reportan contraccion, LPs y pivotes
     * de cada una. La busqueda sigue con la fresca, asi que el arbol es el de
     * siempre. */
    static const bool abcong = (getenv("DFBH_ABCONG") != NULL);
    if (abcong && dynamic_cast<LinearizerAffine2*>(&lr) != NULL) {
        const IntervalVector entrada(box);
        const int TOPE = 8;
        double per0 = 0.0;
        for (int i = 0; i < nb_var; ++i)
            if (entrada[i].diam() < 1e18) per0 += entrada[i].diam();
        if (per0 <= 0.0) { mylinearsolver.clear_constraints(); return; }

        long lp0 = n_lps, pv0 = n_pivots;
        IntervalVector bF(entrada);
        for (int it = 0; it < TOPE && !bF.is_empty(); ++it) {
            double a = 0.0; for (int i = 0; i < nb_var; ++i) if (bF[i].diam()<1e18) a += bF[i].diam();
            const int mm = linearize_fresca(bF, context);
            if (mm == -1) { bF.set_empty(); break; }
            if (mm == 0) break;
            sin_recarga = false;
            una_pasada(bF);
            if (bF.is_empty()) break;
            double b2 = 0.0; for (int i = 0; i < nb_var; ++i) if (bF[i].diam()<1e18) b2 += bF[i].diam();
            if (!(b2 < a*(1.0 - 1e-4))) break;
        }
        const long lpF = n_lps-lp0, pvF = n_pivots-pv0;

        lp0 = n_lps; pv0 = n_pivots;
        IntervalVector bC(entrada);
        nodo_nuevo = true; cong_ids.clear(); b_acum.clear();
        for (int it = 0; it < TOPE && !bC.is_empty(); ++it) {
            double a = 0.0; for (int i = 0; i < nb_var; ++i) if (bC[i].diam()<1e18) a += bC[i].diam();
            const int mm = linearize_congelada(bC, 1000000);
            if (mm == -1) { bC.set_empty(); break; }
            if (mm == 0) break;
            sin_recarga = reusa_A;
            una_pasada(bC);
            sin_recarga = false;
            if (bC.is_empty()) break;
            double b2 = 0.0; for (int i = 0; i < nb_var; ++i) if (bC[i].diam()<1e18) b2 += bC[i].diam();
            if (!(b2 < a*(1.0 - 1e-4))) break;
        }
        const long lpC = n_lps-lp0, pvC = n_pivots-pv0;
        nodo_nuevo = true; cong_ids.clear(); b_acum.clear();

        double perF = 0.0, perC = 0.0;
        const bool vF = bF.is_empty(), vC = bC.is_empty();
        if (!vF) for (int i = 0; i < nb_var; ++i) if (entrada[i].diam()<1e18) perF += bF[i].diam();
        if (!vC) for (int i = 0; i < nb_var; ++i) if (entrada[i].diam()<1e18) perC += bC[i].diam();
        const double rF = vF ? 1.0 : 1.0 - perF/per0, rC = vC ? 1.0 : 1.0 - perC/per0;

        ++ab_n; ab_accF += rF; ab_accC += rC;
        ab_LF += lpF; ab_LC += lpC; ab_PF += pvF; ab_PC += pvC;
        if (vF) ++ab_vacF;
        if (vC) ++ab_vacC;
        if (rF > rC + 1e-12) ++ab_ganaF; else if (rC > rF + 1e-12) ++ab_ganaC; else ++ab_ig;

        box = bF;
        mylinearsolver.clear_constraints();
        if (box.is_empty()) return;
        context.prop.update(BoxEvent(box, BoxEvent::CONTRACT));
        return;
    }

    if (ph_ab) {
        /* Misma caja de entrada para los dos. La esquina de la linealizacion es
         * aleatoria, asi que nodo a nodo no son el mismo poliedro, pero sobre
         * miles de nodos la comparacion de medias es valida. */
        const IntervalVector entrada(box);
        IntervalVector copia(entrada);
        double per0 = 0.0;
        for (int i = 0; i < nb_var; ++i)
            if (entrada[i].diam() < 1e18) per0 += entrada[i].diam();
        const double tph0 = ahora_s();
        const int it0 = ph_ab->n_soplex_iterations, ca0 = ph_ab->n_soplex_calls;
        ph_ab->contract(copia);
        static double t_ph = 0.0; static long it_ph = 0, ca_ph = 0;
        t_ph  += ahora_s() - tph0;
        it_ph += ph_ab->n_soplex_iterations - it0;
        ca_ph += ph_ab->n_soplex_calls - ca0;
        const double tdf0 = ahora_s();
        una_pasada(box);
        static double t_df = 0.0; static long lp0 = 0;
        t_df += ahora_s() - tdf0;
        double pph = 0.0, pdf = 0.0;
        bool vph = copia.is_empty(), vdf = box.is_empty();
        if (!vph) for (int i = 0; i < nb_var; ++i)
            if (entrada[i].diam() < 1e18) pph += copia[i].diam();
        if (!vdf) for (int i = 0; i < nb_var; ++i)
            if (entrada[i].diam() < 1e18) pdf += box[i].diam();
        /* VERIFICACION DE RIGOR (`DFBH_ABHULL=1`): donde DFB deja una cota mas
         * apretada que produccion, se resuelve el LP exacto de esa cota sobre
         * la caja de ENTRADA con un SoPlex nuevo. Ninguna contraccion valida
         * puede pasar el hull exacto del poliedro; si DFB lo pasa (mas alla de
         * la tolerancia), corto puntos factibles de la relajacion. Tambien
         * cuenta cuanto se acerca produccion a ese hull. */
        static const bool abhull = DFB_ENV("DFBH_ABHULL");
        static long hv_n = 0, hv_pasa = 0, hv_pasa_tol = 0, hv_dfb_al_hull = 0, hv_ph_al_hull = 0;
        static double hv_peor = 0.0;
        if (abhull && !vph && !vdf) {
            LPSolver lp(nb_var, LPSolver::Mode::NotCertified);
            IntervalVector ent(entrada);
            if (lr.linearize(ent, lp, context.prop) > 0) {
                lp.set_bounds(ent);
                for (int i = 0; i < nb_var; ++i) {
                    const double w = entrada[i].diam();
                    if (!(w > 0) || w > 1e18) continue;
                    for (int lado = 0; lado < 2; ++lado) {
                        const double bd = lado ? box[i].ub() : box[i].lb();
                        const double bp = lado ? copia[i].ub() : copia[i].lb();
                        const bool dfb_mejor = lado ? (bd < bp - 1e-9*w) : (bd > bp + 1e-9*w);
                        const bool ph_mejor  = lado ? (bp < bd - 1e-9*w) : (bp > bd + 1e-9*w);
                        static long hp_n = 0, hp_dfb_hull = 0, hp_ph_hull = 0;
                        if (ph_mejor) {
                            lp.set_cost(i, lado ? -1.0 : 1.0);
                            const LPSolver::Status st2 = lp.minimize();
                            lp.set_cost(i, 0.0);
                            if (st2 == LPSolver::Status::OptimalProved || st2 == LPSolver::Status::Optimal) {
                                const double h2 = lado ? -lp.minimum().mid() : lp.minimum().mid();
                                ++hp_n;
                                if (std::fabs(bd - h2) <= 1e-6*w) ++hp_dfb_hull;
                                if (std::fabs(bp - h2) <= 1e-6*w) ++hp_ph_hull;
                                if ((hp_n % 500) == 0)
                                    fprintf(stderr, "[abhull-ph] cotas donde produccion gana=%ld | produccion en el hull=%ld  DFB en el hull=%ld\n", hp_n, hp_ph_hull, hp_dfb_hull);
                            }
                            continue;
                        }
                        if (!dfb_mejor) continue;
                        lp.set_cost(i, lado ? -1.0 : 1.0);
                        const LPSolver::Status st = lp.minimize();
                        lp.set_cost(i, 0.0);
                        if (st != LPSolver::Status::OptimalProved && st != LPSolver::Status::Optimal) continue;
                        const double h = lado ? -lp.minimum().mid() : lp.minimum().mid();
                        ++hv_n;
                        const double exceso = (lado ? (h - bd) : (bd - h)) / w;   /* > 0: DFB paso el hull */
                        if (exceso > 1e-7) {
                            ++hv_pasa; if (exceso > hv_peor) hv_peor = exceso;
                            /* Mas alla de la tolerancia del oraculo: SoPlex resuelve con
                             * tolerancia absoluta 1e-9; por debajo de eso su hull no es exacto. */
                            if (exceso*w > 1e-9*std::max(1.0, std::fabs(h))) ++hv_pasa_tol;
                            /* Diagnostico: el x* de SoPlex, cumple las filas y la caja? */
                            const Vector xs = lp.not_proved_primal_sol();
                            const Matrix R = lp.rows(); const Vector L = lp.lhs(), U = lp.rhs();
                            double vf = 0.0, vc = 0.0; int rf = -1;
                            for (int r = 0; r < R.nb_rows(); ++r) {
                                double ax = 0.0, sc = 0.0;
                                for (int c = 0; c < nb_var; ++c) { ax += R[r][c]*xs[c]; sc += std::fabs(R[r][c]*xs[c]); }
                                const double v = std::max(std::max(L[r] - ax, ax - U[r]), 0.0) / std::max(1.0, sc);
                                if (v > vf) { vf = v; rf = r; }
                            }
                            for (int c = 0; c < nb_var; ++c) {
                                const double wc = std::max(entrada[c].diam(), 1e-300);
                                vc = std::max(vc, std::max(entrada[c].lb() - xs[c], xs[c] - entrada[c].ub()) / wc);
                            }
                            fprintf(stderr, "[abhull-pasa] var=%d lado=%d ancho=%.3g dfb=%.17g prod=%.17g hull=%.17g "
                                    "exceso=%.3g | x*: viol. fila max=%.3g (fila %d de %d) viol. caja=%.3g\n",
                                    i, lado, w, bd, bp, h, exceso, vf, rf, R.nb_rows(), vc);
                        }
                        if (std::fabs(bd - h) <= 1e-6*w) ++hv_dfb_al_hull;
                        if (std::fabs(bp - h) <= 1e-6*w) ++hv_ph_al_hull;
                    }
                }
            }
            if (hv_n)
                fprintf(stderr, "[abhull] cotas donde DFB gana=%ld | DFB pasa el hull exacto=%ld (peor %.3g del ancho; "
                        "mas alla de la tolerancia del oraculo=%ld) | "
                        "DFB en el hull=%ld  produccion en el hull=%ld\n", hv_n, hv_pasa, hv_peor, hv_pasa_tol, hv_dfb_al_hull, hv_ph_al_hull);
        }
        /* acumuladores a nivel de archivo: el resumen lo emite el destructor */
        if (per0 > 0) {
            ++ph_n;
            ph_acc_prod += vph ? 1.0 : (1.0 - pph/per0);
            ph_acc_dfb  += vdf ? 1.0 : (1.0 - pdf/per0);
            if (vph) ++ph_vac_prod;
            if (vdf) ++ph_vac_dfb;
            if (!vph && !vdf) { if (pph < pdf*(1-1e-9)) ++ph_gana_prod; else if (pdf < pph*(1-1e-9)) ++ph_gana_dfb; }
            if (DFB_ENV("DFBH_ABVERBOSE") && (ph_n % 200) == 0)
            {
                (void)lp0;
                fprintf(stderr, "[ab] nodos=%ld | PH: %ld LPs %ld iters %.2f it/LP %.1f us/LP %.3fs | "
                        "DFB: %ld LPs %ld piv %.2f piv/LP %.1f us/LP %.3fs | "
                        "vacios PH=%ld DFB=%ld\n", ph_n,
                        ca_ph, it_ph, ca_ph? (double)it_ph/ca_ph:0.0,
                        ca_ph? 1e6*t_ph/ca_ph:0.0, t_ph,
                        n_lps, n_pivots, n_lps? (double)n_pivots/n_lps:0.0,
                        n_lps? 1e6*t_df/n_lps:0.0, t_df, ph_vac_prod, ph_vac_dfb);
            }
        }
        mylinearsolver.clear_constraints();
        context.prop.update(BoxEvent(box, BoxEvent::CONTRACT));
        return;
    }

    /* SONDA A/B PRIMAL vs DUAL SOBRE EL MISMO NODO (`DFBH_ABDUAL=umbral`).
     *
     * El barrido dice que el camino dual explora mas celdas (x1.02-1.06), y eso
     * solo puede venir de contraer menos. Pero comparar dos corridas completas
     * no sirve: divergen en el primer nodo distinto y despues ya no ven los
     * mismos nodos. Aca las dos estrategias resuelven LA MISMA caja, con la
     * misma relajacion, desde el mismo estado —`una_pasada` llama a `load`, que
     * reinicia la base, asi que las dos arrancan en frio— y se comparan las
     * cajas de salida. La busqueda sigue con la del default, de modo que el
     * arbol es el del brazo que se quiera medir.
     *
     * Imprime una linea por nodo cuando la diferencia de reduccion de perimetro
     * supera el umbral, con las variables donde las cajas difieren. */
    static const char* abdual = getenv("DFBH_ABDUAL");
    if (abdual) {
        const double umbral = atof(abdual) > 0.0 ? atof(abdual) : 0.01;
        const IntervalVector entrada(box);
        const IntervalVector z_ini(z);
        /* `DFBH_ABPRIM=R` alterna la REGLA DE LA PRIMERA COTA (0 contra R) en vez
         * de la estrategia de pivoteo, sobre la misma caja. */
        static const int abprim = getenv("DFBH_ABPRIM") ? atoi(getenv("DFBH_ABPRIM")) : -1;
        const bool reub_def = DFBSimplex::pivoteo_dual;
        const bool puro_def = DFBSimplex::dual_puro;
        /* `DFBH_ABPURO=1`: el brazo D es el dual PURO (reubica tambien `k`). */
        static const bool abpuro = (getenv("DFBH_ABPURO") != NULL);
        const int  prim_def = regla_primera;

        double per0 = 0.0;
        for (int i = 0; i < nb_var; ++i)
            if (entrada[i].diam() < 1e18) per0 += entrada[i].diam();

        static const long nodo_traza = getenv("DFBH_ABNODO") ? atol(getenv("DFBH_ABNODO")) : -1;
        static long nodo_ab = 0; ++nodo_ab;
        const bool trazar = (nodo_traza >= 0 && nodo_ab == nodo_traza);
        if (trazar) fprintf(stderr, "[abnodo] nodo=%ld  ===== PRIMAL =====\n", nodo_ab);
        traza_cotas = trazar;
        long lp0 = n_lps, pv0 = n_pivots;
        IntervalVector boxP(entrada);
        if (abprim >= 0) regla_primera = 0; else { DFBSimplex::pivoteo_dual = false; DFBSimplex::dual_puro = false; }
        una_pasada(boxP);
        const long lpP = n_lps - lp0, pvP = n_pivots - pv0;

        if (trazar) fprintf(stderr, "[abnodo] nodo=%ld  ===== DUAL =====\n", nodo_ab);
        lp0 = n_lps; pv0 = n_pivots;
        z = z_ini;
        IntervalVector boxD(entrada);
        if (abprim >= 0) regla_primera = abprim; else { DFBSimplex::pivoteo_dual = true; DFBSimplex::dual_puro = abpuro; }
        una_pasada(boxD);
        const long lpD = n_lps - lp0, pvD = n_pivots - pv0;

        DFBSimplex::pivoteo_dual = reub_def;
        DFBSimplex::dual_puro = puro_def;
        regla_primera = prim_def;
        traza_cotas = false;

        double perP = 0.0, perD = 0.0;
        const bool vP = boxP.is_empty(), vD = boxD.is_empty();
        if (!vP) for (int i = 0; i < nb_var; ++i)
            if (entrada[i].diam() < 1e18) perP += boxP[i].diam();
        if (!vD) for (int i = 0; i < nb_var; ++i)
            if (entrada[i].diam() < 1e18) perD += boxD[i].diam();
        const double redP = vP ? 1.0 : (per0 > 0.0 ? 1.0 - perP/per0 : 0.0);
        const double redD = vD ? 1.0 : (per0 > 0.0 ? 1.0 - perD/per0 : 0.0);

        static long nn = 0, gP = 0, gD = 0, ig = 0, vacP = 0, vacD = 0, ndet = 0;
        static double accP = 0.0, accD = 0.0;
        static double peor = 0.0; static long peor_nodo = -1;
        ++nn; accP += redP; accD += redD;
        if (vP) ++vacP;
        if (vD) ++vacD;
        if      (redP > redD + 1e-12) ++gP;
        else if (redD > redP + 1e-12) ++gD;
        else ++ig;
        /* El caso cualitativo: una prueba vacio y la otra no. Es lo unico que
         * puede cambiar el arbol de golpe; una diferencia de diametro solo
         * mueve el punto de biseccion. */
        static long solo_vP = 0, solo_vD = 0;
        if (vP && !vD) { ++solo_vP; fprintf(stderr, "[abdual] nodo=%ld  VACIO SOLO EN P (redD=%.6f)\n", nn, redD); }
        if (vD && !vP) { ++solo_vD; fprintf(stderr, "[abdual] nodo=%ld  VACIO SOLO EN D (redP=%.6f)\n", nn, redP); }
        if ((nn % 200) == 0)
            fprintf(stderr, "[abdual] vacios exclusivos: solo P=%ld  solo D=%ld\n", solo_vP, solo_vD);
        if (redP - redD > peor) { peor = redP - redD; peor_nodo = nn; }

        if (per0 > 0.0 && (redP - redD > umbral || redD - redP > umbral) && ndet < 40) {
            ++ndet;
            fprintf(stderr, "[abdual] nodo=%ld  redP=%.6f redD=%.6f  dif=%+.6f | "
                    "LPs P=%ld D=%ld  piv P=%ld D=%ld | vacio P=%d D=%d\n",
                    nn, redP, redD, redP-redD, lpP, lpD, pvP, pvD, (int)vP, (int)vD);
            if (!vP && !vD)
                for (int i = 0; i < nb_var; ++i) {
                    const double dP = boxP[i].diam(), dD = boxD[i].diam();
                    if (dP >= 1e18 || dD >= 1e18) continue;
                    if (dD <= dP*(1.0+1e-9) && dP <= dD*(1.0+1e-9)) continue;
                    fprintf(stderr, "          x%-3d  entrada=[%.10g,%.10g]  P=[%.10g,%.10g]  D=[%.10g,%.10g]  %s\n",
                            i, entrada[i].lb(), entrada[i].ub(),
                            boxP[i].lb(), boxP[i].ub(), boxD[i].lb(), boxD[i].ub(),
                            dD > dP ? "D PEOR" : "P peor");
                }
        }
        if ((nn % 200) == 0)
            fprintf(stderr, "[abdual] nodos=%ld | reduccion media P=%.5f D=%.5f | "
                    "gana P=%ld  gana D=%ld  iguales=%ld | vacios P=%ld D=%ld | "
                    "peor caso para D: nodo %ld con %+.6f\n",
                    nn, accP/nn, accD/nn, gP, gD, ig, vacP, vacD, peor_nodo, peor);

        box = (abprim >= 0) ? boxP : (reub_def ? boxD : boxP);
        mylinearsolver.clear_constraints();
        if (box.is_empty()) return;
        context.prop.update(BoxEvent(box, BoxEvent::CONTRACT));
        return;
    }

    /* RELAJACION ADAPTATIVA POR NODO.
     *
     * `compo` poda mas que afin (celdas x0.673) y cuesta mas. Pero la poda
     * extra no hace falta en todos los nodos: donde afin ya contrae bien,
     * las filas de xtaylor son gasto puro.
     *
     * Entonces: una pasada con SOLO las filas de afin —matriz mas chica, no
     * solo cotas ensanchadas, asi que el costo por pivote baja con `m`— y
     * si la caja no se contrajo lo suficiente, se rehace con todas. La
     * propiedad anytime es lo que lo vuelve seguro: la caja de la primera
     * pasada ya es valida, la segunda solo puede mejorarla.
     *
     * `DFBH_ADAPT=u`: se agrega la segunda pasada si la reduccion relativa
     * del perimetro quedo por debajo de `u`. Con u=0 nunca, con u=1 siempre. */
    static const char* adapt = getenv("DFBH_ADAPT");
    const LinearizerCompo* lc_ad = adapt ? dynamic_cast<const LinearizerCompo*>(&lr) : NULL;
    if (lc_ad && lc_ad->nb_l1 > 0 && lc_ad->nb_l1 < A.nb_rows()) {
        const double umbral = atof(adapt);
        double per0 = 0.0;
        for (int i = 0; i < nb_var; ++i)
            if (box[i].diam() != POS_INFINITY) per0 += box[i].diam();
        m_uso = lc_ad->nb_l1;
        una_pasada(box);
        m_uso = 0;
        static long n_nodo = 0, n_seg = 0;
        ++n_nodo;
        if (!box.is_empty()) {
            double per1 = 0.0;
            for (int i = 0; i < nb_var; ++i)
                if (box[i].diam() != POS_INFINITY) per1 += box[i].diam();
            const double red = (per0 > 0.0) ? 1.0 - per1/per0 : 0.0;
            if (red < umbral) {
                ++n_seg;
                for (int i = 0; i < nb_var; ++i) z[i] = box[i];
                una_pasada(box);
            }
        }
        if (DFB_ENV("DFBH_CRONO") && (n_nodo % 500) == 0)
            fprintf(stderr, "[adapt] nodos=%ld  segunda pasada en %.0f%%\n",
                    n_nodo, 100.0*n_seg/n_nodo);
    } else {
        /* CICLO EXTERNO PROPIO, en vez del `CtcFixPoint`.
         *
         * El barrido circular interno (DFBH_CIRC) exprime la relajacion sin
         * re-linealizar: barato, pero se queda con el poliedro flojo. El
         * `CtcFixPoint` externo re-lineariza en CADA iteracion, incluidas las
         * que casi no contraen.
         *
         * Aca se combinan: el ciclo interno hace las vueltas baratas y, cuando
         * converge —o sea cuando ya no queda nada que sacarle a esa
         * relajacion—, se re-lineariza sobre la caja nueva y se repite. Se
         * re-lineariza pocas veces y siempre en el momento en que mas rinde.
         *
         * `DFBH_RELINCIC=r`: repetir mientras el ciclo haya reducido el
         * perimetro mas de `r`, con tope de 8 vueltas. */
        /* El disparador es que el HC4 intercalado HAYA ACTUADO en la vuelta:
         * si contrajo, la caja cambio lo bastante como para que una relajacion
         * nueva valga la pena; si no actuo, el ciclo interno ya convergio sobre
         * una caja que no se movio y re-linealizar seria gasto puro.
         * `DFBH_RELINCIC=k`: tope de `k` re-linealizaciones. */
        /* El disparador es que la vuelta HAYA CONTRAIDO, venga de donde venga
         * —del HC4 intercalado o de las propias cotas—, con el mismo criterio
         * por variable que usa HC4.
         *
         * Atarlo a que HC4 disparara acoplaba dos cosas que no deben estarlo:
         * con umbral de HC4 exigente, `hc4_actuo` casi nunca se activaba y el
         * ciclo externo dejaba de re-linealizar, degenerando al caso peor
         * —`ex7_2_3` se iba de 21 882 a 174 826 celdas—. La re-linealizacion
         * tiene que responder a que la caja cambio, no a quien la cambio.
         *
         * `DFBH_RELINCIC=k`: tope de `k` re-linealizaciones. */
        static const int relincic = getenv("DFBH_RELINCIC")
                                  ? atoi(getenv("DFBH_RELINCIC")) : 0;
        static const double umbral_relin = getenv("DFBH_HC4VAR")
                                         ? atof(getenv("DFBH_HC4VAR")) : 0.3;
        IntervalVector antes(box);
        /* Con la matriz congelada el simplex no recarga: conserva tableau, base
         * y factorizacion, y solo ve cotas de `b` distintas. */
        sin_recarga = reusa_A;
        const IntervalVector antes_fp(box);
        const double pf_t1 = perfil ? ahora_s() : 0.0;
        /* SONDA `DFBH_ABCAJA`: la misma relajacion y la misma caja de entrada
         * con los dos motores. La caja negra corre sobre una copia; la busqueda
         * sigue a DFB. */
        static const bool abcaja = DFB_ENV("DFBH_ABCAJA");
        IntervalVector caja_ab(1);
        const IntervalVector entrada_ab(box);
        if (abcaja && !motor_caja) {
            const IntervalVector zsave(z);
            caja_ab = box;
            motor_caja = true;
            const long p0 = n_pivots; const double t0 = ahora_s();
            pasada_ab_caja(caja_ab);
            abc_t[1] += ahora_s() - t0; abc_piv[1] += n_pivots - p0;
            motor_caja = false;
            z = zsave;
        }
        {
            /* CONTRACTOR HIBRIDO (`DFBH_HIBRIDO=T`, 2026-10-21): si la relajacion
             * del nodo anterior tenia mas de T filas, la pasada la resuelve el
             * simplex de SoPlex manejado directamente (SPxSolver, sin la clase
             * SoPlex ni verificacion), con el mismo orden y el mismo certificado;
             * si no, SC. En relajaciones chicas domina el costo fijo por llamada
             * (gana el tableau de SC); en las grandes, los pivotes (gana SoPlex). */
            static const int hib = getenv("DFBH_HIBRIDO") ? atoi(getenv("DFBH_HIBRIDO")) : 0;
            static bool hib_conf = false;
            if (hib > 0 && !hib_conf) {
                setenv("DFBH_CAJA_SPX", "1", 0); setenv("DFBH_CAJA_DIRECTA", "1", 0);
                setenv("DFBH_CAJA_SINLIBRES", "1", 0); setenv("DFBH_ORDEN_CAJA", "cerca", 0);
                hib_conf = true;
            }
            const bool hib_caja = hib > 0 && !abcaja && !motor_caja && m_afin > hib;
            if (hib_caja) { motor_caja = true; ++n_hib_caja; } else if (hib > 0) ++n_hib_sc;
            const long p0 = n_pivots; const double t0 = ahora_s();
            if (abcaja) pasada_ab_sc(box); else una_pasada(box);
            if (abcaja) { abc_t[0] += ahora_s() - t0; abc_piv[0] += n_pivots - p0; }
            if (hib_caja) motor_caja = false;
        }
        if (abcaja && caja_ab.size() == box.size()) {
            double per0 = 0.0, p[2] = {0.0, 0.0};
            const IntervalVector* r[2] = {&box, &caja_ab};
            for (int i = 0; i < nb_var; ++i) if (entrada_ab[i].diam() < 1e18) per0 += entrada_ab[i].diam();
            if (per0 > 0) {
                ++abc_n;
                for (int e = 0; e < 2; ++e) {
                    if (r[e]->is_empty()) { ++abc_vac[e]; abc_contr[e] += 1.0; p[e] = -1.0; continue; }
                    for (int i = 0; i < nb_var; ++i) if (entrada_ab[i].diam() < 1e18) p[e] += (*r[e])[i].diam();
                    abc_contr[e] += 1.0 - p[e]/per0;
                }
                /* Gana el motor que contrae mas. Si los dos vacian el nodo, empatan
                 * (antes ese caso se contaba como victoria de DFB: p = -1 en los dos). */
                if (!(p[0] < 0 && p[1] < 0)) {
                    if (p[0] < p[1] * (1 - 1e-9)) ++abc_gana[0]; else if (p[1] < p[0] * (1 - 1e-9)) ++abc_gana[1];
                }
            }
        }
        if (perfil) { const double t2 = ahora_s(); pf_pas += t2 - pf_t1; pf_tot += t2 - pf_t0; }
        sin_recarga = false;
        /* FRONTERA DE NODO, con el mismo criterio que usa el punto fijo.
         *
         * `CtcFixPoint` vuelve a llamar mientras
         * `old_box.rel_distance(box) > ratio`, con `ratio = relax_ratio = 0.2`
         * por omision (`DFBH_RATIO` lo cambia). Asi que si nuestra contraccion
         * queda por debajo del umbral, el punto fijo se detiene y la proxima
         * llamada es de otro nodo: hay que re-congelar.
         *
         * Es aproximado por abajo: el punto fijo envuelve toda la composicion,
         * no solo este contractor, asi que el conjunto puede contraer mas que
         * nosotros y seguir iterando. En ese caso se re-congela una vez de mas,
         * que es el lado seguro. */
        static const double ratio_fp = getenv("DFBH_RATIO")
                                     ? atof(getenv("DFBH_RATIO")) : 0.2;
        if (box.is_empty()) nodo_nuevo = true;
        else nodo_nuevo = !(antes_fp.rel_distance(box) > ratio_fp);
        for (int vuelta = 1; vuelta <= relincic && !box.is_empty(); ++vuelta) {
            bool contrajo = false;
            for (int i = 0; i < nb_var && !contrajo; ++i) {
                const double a = antes[i].diam(), b = box[i].diam();
                if (a > 0.0 && a < POS_INFINITY && (1.0 - b/a) > umbral_relin)
                    contrajo = true;
            }
            if (!contrajo) break;
            antes = box;
            mylinearsolver.clear_constraints();
            const int m2 = linearize(box, context);
            if (m2 == -1) { box.set_empty(); break; }
            if (m2 == 0) break;
            una_pasada(box);
        }
    }
    mylinearsolver.clear_constraints();

    context.prop.update(BoxEvent(box, BoxEvent::CONTRACT));
}

} /* namespace ibex */
