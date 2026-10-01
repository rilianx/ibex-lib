#include <limits>
#include <string>
#include "ibex_DFBSimplex.h"
#include <algorithm>
#include <ctime>

#include <algorithm>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cstdio>

namespace ibex {

/* --- Perfilado interno del simplex (DFB_SX_PERF=1) ---
 * El simplex es el 57-85 % del tiempo del contractor (§21) pero nunca se miro
 * ADENTRO. La revision del §42 dejo una sola direccion abierta —bajar el costo
 * sin tocar la poda— y para eso hay que saber donde se va ese tiempo. */
namespace perf {
    double t_pivot = 0, t_dual = 0, t_primal = 0, t_refresh = 0,
           t_basics = 0, t_ydes = 0, t_setobj = 0, t_load = 0, t_total = 0,
           t_escala = 0, t_reparar = 0, t_cola = 0, t_frio = 0, t_bucle = 0;
    long   n_solve = 0;
    inline double ahora() {
        struct timespec ts; clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
        return ts.tv_sec + 1e-9*ts.tv_nsec;
    }
    bool on() { static const bool v = (getenv("DFB_SX_PERF") != NULL); return v; }
    /* Tiempo EXCLUSIVO: las fases estan anidadas —`dual_step` llama a `pivot`,
     * `primal_step` a `compute_basics`— asi que cada Cron descuenta lo que
     * consumieron sus hijos y se lo suma al padre. Sin esto los porcentajes
     * pasan de 100 y no se puede saber donde esta el costo. */
    double hijo = 0.0;
    struct Cron {
        double* acc; double t0; double hijo_prev;
        Cron(double* a) : acc(a), t0(on() ? ahora() : 0.0), hijo_prev(hijo) {
            if (on()) hijo = 0.0;
        }
        ~Cron() {
            if (!on()) return;
            const double dt = ahora() - t0;
            *acc += dt - hijo;
            hijo = hijo_prev + dt;
        }
    };
    void volcar() {
        if (!on()) return;
        /* Con tiempo exclusivo, `t_total` es solo lo propio de `solve`; el
         * denominador es la SUMA de todas las fases. `load` va aparte porque se
         * llama desde fuera de `solve`. */
        const double S = t_total + t_pivot + t_dual + t_primal + t_refresh
                       + t_basics + t_ydes + t_setobj + t_escala + t_reparar + t_cola
                       + t_frio + t_bucle;
        if (S <= 0) return;
        fprintf(stderr, "[perf] solves=%ld  solve=%.3fs (load aparte %.3fs) | pivot %.0f%%  "
                "ratio-test dual %.0f%%  primal %.0f%%  refresh %.0f%%  basics %.0f%%  "
                "y_fact %.0f%%  setobj %.0f%% | ESCALAR LA CAJA %.0f%%  reparar cotas %.0f%%  "
                "cola %.0f%%  arranque en frio %.0f%%  bucle %.0f%%  resto %.0f%% | load = %.0f%% de todo\n",
                n_solve, S, t_load,
                100*t_pivot/S, 100*t_dual/S, 100*t_primal/S,
                100*t_refresh/S, 100*t_basics/S, 100*t_ydes/S,
                100*t_setobj/S, 100*t_escala/S, 100*t_reparar/S, 100*t_cola/S,
                100*t_frio/S, 100*t_bucle/S, 100*t_total/S, 100*t_load/(S+t_load));
    }
}

long DFBSimplex::n_solves     = 0;
long DFBSimplex::n_pivots     = 0;
long DFBSimplex::n_optimal    = 0;
long DFBSimplex::n_infeasible = 0;
long DFBSimplex::n_iterlimit  = 0;
long DFBSimplex::n_warm       = 0;
long DFBSimplex::n_cold       = 0;
long DFBSimplex::n_fact_fallo = 0;
long DFBSimplex::n_refact_dif = 0;
/* `pivot` no actualiza la LU salvo con `DFB_SX_LU=actualizada`: vale para la
 * diferida (cpu por celda x0.94) y para el tableau, que es la omision. */
static bool leer_lu_diferida() {
    const char* v = getenv("DFB_SX_LU");
    return v == NULL || std::string(v) != "actualizada";
}
bool DFBSimplex::lu_diferida = leer_lu_diferida();
/* Omision `tableau` desde 2026-09-30 (cpu x0.89 contra la diferida en el
 * banco, ver ESTADO.md §9). `DFB_SX_LU=diferida` o `actualizada` la apagan. */
static bool leer_lu_tableau() {
    const char* v = getenv("DFB_SX_LU");
    return v == NULL || (std::string(v) != "diferida" && std::string(v) != "actualizada");
}
bool DFBSimplex::lu_tableau = leer_lu_tableau();
long DFBSimplex::n_guardia_ok = 0, DFBSimplex::n_guardia_fallo_d = 0,
     DFBSimplex::n_guardia_fallo_xb = 0, DFBSimplex::n_guardia_fallo_fin = 0;
static const double GUARDIA_TOL = getenv("DFB_SX_GUARDIA") ? atof(getenv("DFB_SX_GUARDIA")) : 1e-9;
bool DFBSimplex::sonda_deriva = (getenv("DFB_SX_DERIVA") != NULL);
long DFBSimplex::dv_refrescos = 0, DFBSimplex::dv_d_signo = 0, DFBSimplex::dv_d_nobas = 0,
     DFBSimplex::dv_xb_fact = 0, DFBSimplex::dv_xb_bas = 0, DFBSimplex::dv_refr_dist = 0;
double DFBSimplex::dv_d_max = 0.0, DFBSimplex::dv_xb_max = 0.0;
long DFBSimplex::n_base_heredada = 0;
/* APAGADO por omision. Midio neutro en el test sintetico y **daña** en las
 * instancias reales: Virasoro-icse cae de 100 % a 66.67 % de contraccion. Es la
 * misma leccion que con el regimen de escalas del test: lo que no se mide sobre
 * el banco real no se sabe. Se enciende con DFB_SX_SCALE=1. */
bool DFBSimplex::escalar      = (getenv("DFB_SX_SCALE") != NULL);
long DFBSimplex::n_k_basica     = 0;
long DFBSimplex::n_piv_degen    = 0;
long DFBSimplex::n_pert_lps       = 0;
long DFBSimplex::n_pert_limpieza  = 0;
long DFBSimplex::n_piv_degen_pert = 0;
namespace {
/* `DFB_SX_PIVOTEO=primal|dual|mixta`, ver la cabecera. `mixta` es la que
 * hasta el 2026-09-30 se llamaba `dual`: reubica todo menos `k`, y en los LPs
 * que empiezan con `k` no basica termina pivoteando en primal. */
/* Omision `dual` (puro) desde 2026-09-30, ver ESTADO.md §9 opcion 3.
 * `DFB_SX_PIVOTEO=primal` restaura la estrategia anterior. */
bool leer_pivoteo_dual() {
    const char* v = getenv("DFB_SX_PIVOTEO");
    return v == NULL || std::string(v) == "dual" || std::string(v) == "mixta";
}
bool leer_dual_puro() {
    const char* v = getenv("DFB_SX_PIVOTEO");
    return v == NULL || std::string(v) == "dual";
}
}
bool DFBSimplex::pivoteo_dual = leer_pivoteo_dual();
bool DFBSimplex::dual_puro = leer_dual_puro();
long DFBSimplex::n_entrak = 0;
long DFBSimplex::n_verif      = 0;
long DFBSimplex::n_verif_pinf = 0;
long DFBSimplex::n_verif_dinf = 0;
long DFBSimplex::n_flips        = 0;
long DFBSimplex::n_k_no_basica  = 0;

void DFBSimplex::reset_counters() {
    n_solves = n_pivots = n_optimal = n_infeasible = n_iterlimit = 0;
    n_warm = n_cold = n_fact_fallo = 0;
    n_k_basica = n_k_no_basica = 0;
    n_piv_degen = n_flips = 0;
    n_pert_lps = n_pert_limpieza = n_piv_degen_pert = 0;
    n_verif = n_verif_pinf = n_verif_dinf = 0;
    n_anytime = n_hereda = n_base_rechazada = n_vueltas = 0;
    n_base_heredada = 0;
}

namespace {
const double PIV_TOL  = 1e-9;   /* pivote minimo aceptable */
/* Violacion primal que se considera nula. Escalada por la MAGNITUD de las
 * cotas de la basica. En nodos profundos la caja es diminuta frente a esa
 * magnitud, y entonces la infactibilidad que introduce una reubicacion queda
 * por debajo de la tolerancia: el paso dual no ve nada que reparar, el primal
 * tampoco (la reubicacion ya dejo los signos bien) y el bucle declara OPTIMAL
 * en un punto que no lo es. `DFB_SX_FEASTOL` para medir ese eje. */
const double FEAS_TOL = getenv("DFB_SX_FEASTOL") ? atof(getenv("DFB_SX_FEASTOL")) : 1e-9;
const double DUAL_TOL = 1e-9;
/* TOLERANCIA A ESCALA DE LA CAJA (`DFB_SX_TOLANCHO=tau`, apagada con 0).
 * Medido en `bearing` (2026-09-30): con cajas de ~5e-8 de ancho y cotas de
 * magnitud ~7, `FEAS_TOL*(1+|cota|)` es ~13 % del ancho, y el dual puro
 * declaraba OPTIMAL con `k` no basica —certificado `y = 0`— un poliedro que el
 * primal probaba vacio en intervalos. Se acota la tolerancia por `tau` veces
 * el ancho de la basica, con un piso de redondeo. */
/* Omision 1e-3 desde 2026-09-30; `DFB_SX_TOLANCHO=0` la apaga. */
const double TOL_ANCHO = getenv("DFB_SX_TOLANCHO") ? atof(getenv("DFB_SX_TOLANCHO")) : 1e-3;
inline double tol_primal(double lo, double hi, double esc) {
    double tol = FEAS_TOL * esc;
    if (TOL_ANCHO > 0.0) {
        const double w = hi - lo;
        if (w >= 0.0 && w < POS_INFINITY) {
            const double ta = TOL_ANCHO * w;
            if (ta < tol) tol = ta;
        }
        const double piso = 1e-15 * esc;
        if (tol < piso) tol = piso;
    }
    return tol;
}
long g_ps_sin_entrante = 0, g_ps_pivote_fallido = 0, g_ps_unbounded = 0;
long g_ds_pivote_fallido = 0;   /* costo reducido que se considera nulo */
const double HARRIS_BAND = 1e-3; /* banda relativa del desempate de Harris */
}

bool DFBSimplex::cand_menor(const DFBSimplex::Cand& a, const DFBSimplex::Cand& b) {
    return a.ratio < b.ratio;
}

long DFBSimplex::n_anytime = 0;
long DFBSimplex::n_vueltas = 0;
long DFBSimplex::n_hereda = 0;
long DFBSimplex::n_base_rechazada = 0;
DFBSimplex::DFBSimplex() : traza_y(NULL), tope_anytime(std::numeric_limits<double>::quiet_NaN()),
                           obs_tope(std::numeric_limits<double>::quiet_NaN()), obs_pivote(-1), xb_sucio(true), m(0), nx(0), na(0), nc(0), zs(1),
                           bnd(0.0), iters(0),
                           base_lista(false), fact_ok(false), fact_sucia(false),
                           perturbado(false), pert_eps(0.0), pert_estado(12345u) {
}

/* Equilibrado en norma infinito, alternando columnas y filas, con las escalas
 * redondeadas a potencias de 2 para que multiplicar por ellas sea exacto. Las
 * columnas de b quedan atadas a la fila (C[nx+j] = 1/R[j]) para que la parte b
 * siga siendo -I. */
void DFBSimplex::calcular_escalas(const Matrix& Abar) {
    rscale.assign(m, 1.0);
    cscale.assign(na, 1.0);
    if (!escalar) return;

    for (int pasada = 0; pasada < 3; ++pasada) {
        /* columnas de x */
        for (int c = 0; c < nx; ++c) {
            double mx = 0.0;
            for (int j = 0; j < m; ++j) {
                const double v = std::fabs(rscale[j] * Abar[j][c] * cscale[c]);
                if (v > mx) mx = v;
            }
            if (mx > 0.0) {
                const int e = (int)std::floor(std::log2(mx) + 0.5);
                cscale[c] *= std::ldexp(1.0, -e);
            }
        }
        /* filas; la columna b de la fila aporta rscale[j]*(-1)*cscale[nx+j] = -1 */
        for (int j = 0; j < m; ++j) {
            double mx = 0.0;
            for (int c = 0; c < nx; ++c) {
                const double v = std::fabs(rscale[j] * Abar[j][c] * cscale[c]);
                if (v > mx) mx = v;
            }
            if (mx > 0.0) {
                const int e = (int)std::floor(std::log2(mx) + 0.5);
                rscale[j] *= std::ldexp(1.0, -e);
            }
        }
    }
    for (int j = 0; j < m; ++j) cscale[nx + j] = 1.0 / rscale[j];
}

void DFBSimplex::load(const Matrix& Abar, int nx_) {
    perf::Cron _c(&perf::t_load);
    /* Si la linealizacion conserva las dimensiones, la BASE se conserva como
     * conjetura para el nodo nuevo. Los coeficientes cambian —hay que
     * refactorizar— pero el conjunto de columnas basicas sigue siendo un buen
     * punto de partida, y entre nodos la caja cambio por una biseccion: una
     * variable partida al medio y las demas iguales, que es justo el cambio que
     * un simplex dual repara en pocos pivotes.
     *
     * Sin esto, `load` reseteaba la base en CADA nodo y el warm start entre
     * nodos no ocurria nunca; medido, daba lo mismo que la base compartida. */
    /* HERENCIA DE LA BASE ENTRE NODOS (`DFB_SX_BASENODOS=1`).
     *
     * Estuvo apagada porque `load` restauraba la base y dejaba el TABLEAU
     * vacio: el refresco desde la factorizacion funcionaba, pero el test de
     * razon lee el tableau y no encontraba candidatas. Ahora se reconstruye
     * `T = B^-1 [Abar | I]` con la matriz NUEVA, que es lo que la base heredada
     * significa. Cuesta `m` BTRAN mas `m*nnz(Abar)`, del orden de `m` pivotes.
     *
     * Se exige que la escala este apagada, porque el tableau vive en el espacio
     * escalado y la factorizacion en el sin escalar; con `DFB_SX_SCALE=1` no
     * coinciden. */
    static const bool conservar = (getenv("DFB_SX_BASENODOS") != NULL) && !escalar;
    const bool mismas_dim = conservar && base_lista &&
                            m == Abar.nb_rows() && na == Abar.nb_cols() &&
                            nx == nx_;

    m  = Abar.nb_rows();
    nx = nx_;
    na = Abar.nb_cols();
    nc = na + m;

    calcular_escalas(Abar);

    arow.assign(m, std::vector<std::pair<int,double> >());
    for (int j = 0; j < m; ++j)
        for (int c = 0; c < na; ++c) {
            const double v = rscale[j] * Abar[j][c] * cscale[c];
            if (v != 0.0) arow[j].push_back(std::make_pair(c, v));
        }
    zs.resize(na);

    /* OJO: la base a heredar hay que copiarla ANTES de los `assign` de abajo.
     * Antes se leia `basic` despues de resetearlo a -1, asi que `set_basis`
     * recibia columnas invalidas, fallaba, y la herencia entre nodos no ocurria
     * nunca: `DFB_SX_NODOS=1` daba resultados identicos a no usarlo. */
    std::vector<int>  basic_prev;
    std::vector<int>  rowof_prev;
    std::vector<char> atup_prev;
    if (mismas_dim) { basic_prev = basic; rowof_prev = rowof; atup_prev = atupper; }

    /* HERENCIA DE LAS COTAS DE APOYO ENTRE NODOS (`DFB_SX_HEREDA=1`).
     *
     * Es la version barata de heredar estado: no se hereda la BASE —eso
     * obligaria a reconstruir el tableau, `na` FTRAN— sino solo en cual de sus
     * dos cotas se apoya cada no basica. La base sigue siendo `B = -I`, o sea
     * que el tableau se arma gratis como siempre; lo unico que cambia es el
     * VERTICE del que se parte. Entre nodos la caja cambio por una biseccion,
     * asi que el vertice optimo del nodo anterior es una conjetura razonable.
     *
     * Requiere ademas que la reubicacion del arranque en frio sea MINIMA (ver
     * `set_objective`), porque la agresiva manda todo a la cota inferior y
     * borraria lo heredado. */
    static const bool hereda = (getenv("DFB_SX_HEREDA") != NULL);
    std::vector<char> atup_hered;
    const bool hay_hered = hereda && (int)atupper.size() == nc
                        && m == Abar.nb_rows() && na == Abar.nb_cols() && nx == nx_;
    if (hay_hered) atup_hered = atupper;

    d.assign(na, 0.0);
    ylam.assign(m, 0.0);
    y.assign(m, 0.0);
    tidx.assign(m + 1, std::vector<int>());
    tval.assign(m + 1, std::vector<double>());
    basic.assign(m + 1, -1);
    rowof.assign(nc, -1);
    atupper.assign(nc, 0);
    xN.assign(nc, 0.0);
    xB.assign(m + 1, 0.0);
    colq.assign(m + 1, 0.0);
    lam_buf.assign(m, 0.0);
    vbuf.assign(m, 0.0);
    xhat.assign(na, 0.0);
    cols_buf.assign(m, 0);
    fact.load_matrix(Abar);
    fact_ok = false;
    fact_sucia = false;
    base_lista = false;
    bidx.reserve(nc);
    bval.reserve(nc);
    base_lista = false;

    if (hay_hered && (int)atup_hered.size() == nc) {
        atupper = atup_hered;
        ++n_hereda;
    }

    if (mismas_dim && (int)basic_prev.size() == m+1) {
        /* Se restaura la base anterior y se refactoriza con los coeficientes
         * NUEVOS: la relajacion se recalculo (poliedro fuerte) pero se arranca
         * del conjunto de columnas basicas que ya era optimo para una caja muy
         * parecida. Es lo contrario del §39, que reusaba la relajacion vieja. */
        basic = basic_prev; rowof = rowof_prev; atupper = atup_prev;
        for (int r = 1; r <= m; ++r) cols_buf[r-1] = basic[r];
        fact_ok = fact.set_basis(cols_buf);
        base_lista = fact_ok;
        if (fact_ok) {
            ++n_base_heredada;
            /* RECONSTRUCCION DEL TABLEAU para la base heredada, con la `Abar`
             * nueva. Fila `i` de `T` = (fila `i` de `B^-1`) * [Abar | I]. Se
             * acumula denso sobre `nc` y se comprime. */
            std::vector<double> dens(nc, 0.0);
            for (int r = 1; r <= m && fact_ok; ++r) {
                fact.btran_row(r-1, &lam_buf[0]);
                std::fill(dens.begin(), dens.end(), 0.0);
                for (int j = 0; j < m; ++j) {
                    const double w = lam_buf[j];
                    if (w == 0.0) continue;
                    const std::vector<std::pair<int,double> >& a = arow[j];
                    for (size_t t = 0; t < a.size(); ++t)
                        dens[a[t].first] += w * a[t].second;
                    dens[na + j] += w;          /* la parte identidad */
                }
                tidx[r].clear(); tval[r].clear();
                for (int c = 0; c < nc; ++c)
                    if (dens[c] != 0.0) { tidx[r].push_back(c); tval[r].push_back(dens[c]); }
                /* La columna basica de la fila tiene que valer exactamente 1.
                 * Si no, la factorizacion no corresponde a esta base y es mas
                 * seguro arrancar en frio que pivotear sobre un tableau falso. */
                const int b = basic[r];
                double vb = 0.0;
                for (size_t t = 0; t < tidx[r].size(); ++t)
                    if (tidx[r][t] == b) { vb = tval[r][t]; break; }
                if (!(std::fabs(vb - 1.0) < 1e-7)) fact_ok = false;
            }
            if (!fact_ok) ++n_base_rechazada;
        }
        if (!fact_ok) {
            base_lista = false;
            basic.assign(m+1,-1); rowof.assign(nc,-1); atupper.assign(nc,0);
            for (int r = 1; r <= m; ++r) { tidx[r].clear(); tval[r].clear(); }
        }
    }
}

double DFBSimplex::entry(int r, int col) const {
    if (r == 0) return (col < na) ? d[col] : ylam[col - na];
    const std::vector<int>& ix = tidx[r];
    const std::vector<int>::const_iterator it =
        std::lower_bound(ix.begin(), ix.end(), col);
    if (it == ix.end() || *it != col) return 0.0;
    return tval[r][it - ix.begin()];
}

/* fila_r -= f * fila_rp, las dos dispersas y ordenadas.
 * OJO con `0.0 - f*b` en la rama donde la fila destino no tiene entrada: con
 * redondeo dirigido `-x` y `0-x` no son la misma operacion (ver la trampa 10
 * del plan). Aca no hay una referencia densa que reproducir, pero se mantiene
 * la misma forma por consistencia. */
void DFBSimplex::axpy(int r, int rp, double f, int drop_col) {
    const std::vector<int>&    ai = tidx[r];  const std::vector<double>& av = tval[r];
    const std::vector<int>&    bi = tidx[rp]; const std::vector<double>& bv = tval[rp];
    bidx.clear(); bval.clear();
    size_t p = 0, q = 0;
    while (p < ai.size() || q < bi.size()) {
        int c; double v;
        if (q == bi.size() || (p < ai.size() && ai[p] < bi[q])) {
            c = ai[p]; v = av[p]; ++p;
        } else if (p == ai.size() || bi[q] < ai[p]) {
            c = bi[q]; v = 0.0 - f * bv[q]; ++q;
        } else {
            c = ai[p]; v = av[p] - f * bv[q]; ++p; ++q;
        }
        if (c == drop_col) continue;
        if (v == 0.0) continue;
        bidx.push_back(c); bval.push_back(v);
    }
    tidx[r].swap(bidx); tval[r].swap(bval);
}

static bool leer_dse() { return getenv("DFB_SX_DSE") != NULL && std::string(getenv("DFB_SX_DSE")) != "0"; }
static const bool USA_DSE = leer_dse();

/* Norma al cuadrado de la parte lambda de la fila `r` del tableau: es
 * `||e_r^T B^-1||^2`, el peso de dual steepest edge. */
double DFBSimplex::norma_lambda(int r) const {
    double s2 = 0.0;
    const std::vector<int>& ix = tidx[r]; const std::vector<double>& vl = tval[r];
    for (size_t t = 0; t < ix.size(); ++t) if (ix[t] >= na) s2 += vl[t]*vl[t];
    return s2 > 0.0 ? s2 : 1e-300;
}

/* Hace basica la columna `col` en la fila `r`: normaliza y elimina la columna
 * de las demas filas, incluida la fila objetivo (que es lo que actualiza los
 * costos reducidos y `y`). */
bool DFBSimplex::pivot(int r, int col) {
    perf::Cron _c(&perf::t_pivot);
    /* Todo pivote deja `x_B` obsoleto. Solo el paso dual sabe actualizarlo
     * incrementalmente y lo limpia; cualquier otro camino —`primal_step`, por
     * ejemplo— obliga a recalcular. */
    xb_sucio = true;
    const double p = entry(r, col);
    if (std::fabs(p) < PIV_TOL) return false;
    const double inv = 1.0 / p;

    {   std::vector<double>& vl = tval[r];
        for (size_t t = 0; t < vl.size(); ++t) vl[t] *= inv;
        const std::vector<int>& ix = tidx[r];
        const std::vector<int>::const_iterator it =
            std::lower_bound(ix.begin(), ix.end(), col);
        if (it != ix.end() && *it == col) vl[it - ix.begin()] = 1.0;
    }

    /* fila objetivo */
    const double f0 = (col < na) ? d[col] : ylam[col - na];
    if (f0 != 0.0) {
        const std::vector<int>&    ix = tidx[r];
        const std::vector<double>& vl = tval[r];
        for (size_t t = 0; t < ix.size(); ++t) {
            const int c = ix[t];
            if (c < na) d[c]       -= f0 * vl[t];
            else        ylam[c-na] -= f0 * vl[t];
        }
        if (col < na) d[col] = 0.0; else ylam[col - na] = 0.0;
    }

    /* filas de restriccion. De paso se actualizan los pesos devex: la columna
     * entrante transformada es justo lo que hace falta y ya se recorre. */
    /* Omision encendido desde 2026-09-30 con el dual puro (ESTADO.md §9,
     * opcion 3); `DFB_SX_DEVEX=0` lo apaga. */
    static const bool devex = !(getenv("DFB_SX_DEVEX") && std::string(getenv("DFB_SX_DEVEX")) == "0");
    const double wr = (devex && (int)dw.size() == m+1) ? dw[r] : 0.0;
    const double pr2 = p*p;
    /* `alpha_piv` guarda la columna entrante TAL COMO ESTABA antes de tocar las
     * filas: es justo lo que necesita la actualizacion incremental de `x_B`. */
    if ((int)alpha_piv.size() != m+1) alpha_piv.assign(m+1, 0.0);
    std::fill(alpha_piv.begin(), alpha_piv.end(), 0.0);
    alpha_piv[r] = p;
    for (int i = 1; i <= m; ++i) {
        if (i == r) continue;
        const double f = entry(i, col);
        alpha_piv[i] = f;
        if (devex && wr > 0.0 && f != 0.0) {
            const double cand = (f*f/pr2) * wr;
            if (cand > dw[i]) dw[i] = cand;
        }
        if (f == 0.0) continue;
        axpy(i, r, f, col);
        if (USA_DSE && (int)dse.size() == m+1) dse[i] = norma_lambda(i);
    }
    if (USA_DSE && (int)dse.size() == m+1) dse[r] = norma_lambda(r);

    rowof[basic[r]] = -1;
    if (devex && (int)dw.size() == m+1) {
        const double nw = wr / pr2;
        dw[r] = nw > 1.0 ? nw : 1.0;
    }
    basic[r] = col;
    rowof[col] = r;
    /* La factorizacion sigue el mismo cambio de base, con update tipo
     * Forrest-Tomlin y refactorizacion periodica (ver DFBBasis). */
    if (lu_diferida) fact_sucia = true;
    else if (fact_ok && !fact.change_basis(r - 1, col)) fact_ok = false;
    ++n_pivots;
    return true;
}

/* Valores de las basicas con las no basicas en la cota que indica atupper. */
/* Una no basica tiene que apoyarse en una cota FINITA. Si la que le toca por
 * signo es infinita se usa la otra; si las dos lo son, la variable es libre y se
 * la deja en 0. Sin esto un +-inf se propaga y convierte todo en NaN.
 *
 * Esto estaba DENTRO de `compute_basics`, que asi no era una lectura pura: al
 * saltarla —porque `x_B` ya venia al dia del pivote anterior— la reparacion no
 * ocurria y la cota flotante de una `k` no basica quedaba en +-inf. Daba 8
 * discrepancias de 528 en el test unitario con los MISMOS 810 pivotes, que es la
 * firma de que el estado difiere sin que difiera la trayectoria. */
void DFBSimplex::reparar_cotas_infinitas(const IntervalVector& z) {
    for (int c = 0; c < na; ++c) {
        if (rowof[c] >= 0) continue;
        const double v = atupper[c] ? z[c].ub() : z[c].lb();
        if (v != POS_INFINITY && v != NEG_INFINITY) continue;
        const double o = atupper[c] ? z[c].lb() : z[c].ub();
        if (o != POS_INFINITY && o != NEG_INFINITY) atupper[c] = atupper[c] ? 0 : 1;
    }
}

long DFBSimplex::n_basics = 0;
void DFBSimplex::compute_basics(const IntervalVector& z) {
    perf::Cron _c(&perf::t_basics);
    ++n_basics;
    xb_sucio = false;
    for (int c = 0; c < nc; ++c) {
        if (rowof[c] >= 0) { xN[c] = 0.0; continue; }
        if (c >= na) { xN[c] = 0.0; continue; }  /* columnas lambda: no son variables */
        /* Una no basica tiene que apoyarse en una cota FINITA. Si la que le
         * toca por signo es infinita se usa la otra; si las dos lo son, la
         * variable es libre y se la deja en 0. Sin esto un +-inf se propaga y
         * convierte todo en NaN. */
        double v = atupper[c] ? z[c].ub() : z[c].lb();
        if (v == POS_INFINITY || v == NEG_INFINITY) {
            const double o = atupper[c] ? z[c].lb() : z[c].ub();
            if (o == POS_INFINITY || o == NEG_INFINITY) v = 0.0;
            else { v = o; atupper[c] = atupper[c] ? 0 : 1; }
        }
        xN[c] = v;
    }
    for (int r = 1; r <= m; ++r) {
        double acc = 0.0;
        const std::vector<int>&    ix = tidx[r];
        const std::vector<double>& vl = tval[r];
        for (size_t t = 0; t < ix.size(); ++t) {
            const int c = ix[t];
            if (c >= na) break;
            if (rowof[c] >= 0) continue;
            acc += vl[t] * xN[c];
        }
        xB[r] = -acc;
    }
}

/* Fila objetivo para el costo c = +-e_k sobre la base actual: se arranca de
 * [c | 0] y se eliminan las columnas basicas restando f*fila_r, que es lo que
 * deja d = c - Abar^T y con y = B^-T c_B. Despues cada no basica se apoya en la
 * cota que dicta el signo de su costo reducido, con lo que el punto queda
 * dual-factible sin pivotear. */
void DFBSimplex::set_objective(int k, bool maximize, bool reubicar) {
    perf::Cron _c(&perf::t_setobj);
    xb_sucio = true;
    std::fill(d.begin(), d.end(), 0.0);
    std::fill(ylam.begin(), ylam.end(), 0.0);
    d[k] = maximize ? -1.0 : 1.0;

    for (int r = 1; r <= m; ++r) {
        const int b = basic[r];
        const double f = (b < na) ? d[b] : ylam[b - na];
        if (f == 0.0) continue;
        const std::vector<int>&    ix = tidx[r];
        const std::vector<double>& vl = tval[r];
        for (size_t t = 0; t < ix.size(); ++t) {
            const int c = ix[t];
            if (c < na) d[c]       -= f * vl[t];
            else        ylam[c-na] -= f * vl[t];
        }
        if (b < na) d[b] = 0.0; else ylam[b - na] = 0.0;
    }

    /* Reubicacion del arranque en frio. La forma agresiva manda a la cota
     * INFERIOR todo lo que no viole el signo, incluidas las degeneradas
     * (`d = 0`), que con `c_B = 0` y `c = +-e_k` son todas menos `k`. La minima
     * deja quietas las degeneradas, que es la forma que el camino tibio ya usa
     * y la unica con la que sobrevive una ubicacion heredada. */
    static const bool minima = (getenv("DFB_SX_HEREDA") != NULL);
    if (reubicar)
        for (int c = 0; c < na; ++c) {
            if (rowof[c] >= 0) continue;
            if (minima) {
                if (d[c] < -DUAL_TOL)      atupper[c] = 1;
                else if (d[c] > DUAL_TOL)  atupper[c] = 0;
            } else {
                atupper[c] = (d[c] < -DUAL_TOL) ? 1 : 0;
            }
        }
}

/* Columna q del tableau: se recorre fila por fila. Con el tableau guardado por
 * filas cuesta una busqueda binaria por fila, o sea O(m log nnz), que es
 * despreciable contra el pricing. */
void DFBSimplex::tab_column(int q, std::vector<double>& out) const {
    for (int r = 1; r <= m; ++r) out[r] = entry(r, q);
}


/* Paso dual: sale la basica que mas viola sus cotas (violacion RELATIVA, ver
 * abajo), entra la del menor cociente |d_j/a_rj| con desempate de Harris. */
bool DFBSimplex::dual_step(const IntervalVector& z, bool& infeasible) {
    /* Se recalcula solo si algo lo ensucio; si no, `x_B` viene al dia del
     * pivote anterior. */
    /* Actualizacion incremental de `x_B`: APAGADA por omision. Es correcta
     * —528/528 en el test unitario— pero la acumulacion redondea distinto que el
     * recalculo, y esa deriva mueve las decisiones del arbol: cpu x1.018 y
     * celdas distintas en 25 de 128 instancias. El ahorro era chico porque los
     * bound flips y los pivotes primales ensucian `x_B` seguido. Ver §44.
     * `DFB_SX_XBINC=1` la activa. */
    const bool inc = false;
    (void)inc;
    /* `inc` esta apagado, asi que `xb_sucio` es true siempre que se llega aca:
     * la condicion es equivalente a la anterior y deja que el bucle de `solve`
     * comparta el recalculo cuando hace la comprobacion anytime. */
    if (xb_sucio) compute_basics(z);
    return dual_step_comun(z, infeasible);
}

bool DFBSimplex::dual_step_comun(const IntervalVector& z, bool& infeasible) {
    perf::Cron _c(&perf::t_dual);
    infeasible = false;
    /* DFB_SX_DEVEX=1: eleccion de la fila que sale por devex en vez de Dantzig.
     * Ver el comentario de `dw` en la cabecera. Estuvo forzado a false; se
     * reabre porque el camino dual con el arreglo de `k` paga ~2.4 veces mas
     * pivotes que el primal para el mismo conteo predicho, y Dantzig contra
     * devex es la causa conocida de ese tipo de brecha. */
    /* Omision encendido desde 2026-09-30 con el dual puro (ESTADO.md §9,
     * opcion 3); `DFB_SX_DEVEX=0` lo apaga. */
    static const bool devex_env = !(getenv("DFB_SX_DEVEX") && std::string(getenv("DFB_SX_DEVEX")) == "0");
    /* Con dual steepest edge los pesos exactos reemplazan a los de devex. */
    const bool devex = devex_env || USA_DSE;
    if (USA_DSE && (int)dse.size() != m+1) {
        dse.assign(m+1, 1.0);
        for (int r = 1; r <= m; ++r) dse[r] = norma_lambda(r);
    }
    if (devex && (int)dw.size() != m+1) dw.assign(m+1, 1.0);
    const std::vector<double>& peso = USA_DSE ? dse : dw;
    /* `DFB_SX_DSEABS=1`: con DSE, violacion ABSOLUTA como en el libro, en vez
     * de la relativa a la magnitud de la cota. */
    static const bool dse_abs = USA_DSE && getenv("DFB_SX_DSEABS") != NULL;

    /* Con tolerancia absoluta de 1e-9 sobre variables de magnitud 1e9 —las
     * cajas de Brown-* tienen radio 1e9— el error de redondeo relativo de 1e-16
     * da una violacion absoluta de 1e-7, o sea diez veces la tolerancia: el
     * simplex pivotea sobre ruido hasta quedarse sin candidatas y declara
     * infactible un LP factible. Pasaba en 35 de 71 resoluciones de Brown-20. */
    int r_out = -1; double worst = 0.0; bool below = false;
    for (int r = 1; r <= m; ++r) {
        const int c = basic[r];
        const double lo = z[c].lb(), hi = z[c].ub();
        const double esc = 1.0 + ((std::fabs(lo) > std::fabs(hi)) ?
                                  std::fabs(lo) : std::fabs(hi));
        const double tol = tol_primal(lo, hi, esc);
        if (xB[r] < lo - tol) {
            const double v0 = dse_abs ? (lo - xB[r]) : (lo - xB[r]) / esc;
            const double v = devex ? v0*v0/((int)peso.size()>r ? peso[r] : 1.0) : v0;
            if (v > worst) { worst = v; r_out = r; below = true; }
        } else if (xB[r] > hi + tol) {
            const double v0 = dse_abs ? (xB[r] - hi) : (xB[r] - hi) / esc;
            const double v = devex ? v0*v0/((int)peso.size()>r ? peso[r] : 1.0) : v0;
            if (v > worst) { worst = v; r_out = r; below = false; }
        }
    }
    if (r_out < 0) return false;      /* factible primal */

    /* x_B(r) = -sum_N a_rj x_j, asi que d x_B(r)/d x_j = -a_rj. Si la basica
     * esta por DEBAJO hay que subirla: sirve una no basica en su cota inferior
     * con a_rj < 0 (puede crecer) o en su superior con a_rj > 0 (puede bajar).
     * Al reves si esta por encima. */
    /* --- RATIO TEST DUAL DE PASO LARGO, CON CAMBIO DE COTA ---
     *
     * `DFB_SX_BFRT=1`. Medido dentro del arbol, el **78-88 % de los pivotes son
     * degenerados**: la columna entrante tiene costo reducido ~0, asi que el
     * paso dual no mejora el objetivo y hace falta 2.5-4.5 veces mas pivotes que
     * SoPlex con un costo por pivote comparable (§21).
     *
     * El remedio estandar (Fourer / Maros / Koberstein) es no detenerse en el
     * primer punto de quiebre. Se recorren los quiebres en orden creciente: si
     * la candidata de turno esta ACOTADA por los dos lados, en vez de hacerla
     * basica se la **cambia de cota**, lo que consume `|a_rj|*(u_j-l_j)` de la
     * infactibilidad de la fila que sale; mientras quede infactibilidad del
     * mismo signo, el objetivo dual sigue mejorando y conviene seguir. Se pivotea
     * recien en el quiebre donde la pendiente cambiaria de signo, y se dan asi
     * **un** paso largo en vez de muchos de longitud cero.
     *
     * Nuestro caso es el ideal para esto: en la relajacion de DFB **todas** las
     * variables estan acotadas, incluidas las `b`.
     */
    /* El paso largo gana en relajaciones GRANDES y pierde feo en chicas: medido
     * por franja de `m`, la ventaja en celdas va de -0.150 (m<12) a +0.357
     * (m>=50), y las catastrofes estan todas en m<=12 (`ex9_2_5`, m=10, de 6 a
     * 10 668 celdas). `DFB_SX_BFRT_M=k` lo activa solo cuando `m >= k`;
     * `DFB_SX_BFRT=1` lo activa siempre. */
    /* Se reabre por el mismo motivo que devex: la medicion que lo condeno
     * (`m` chico, celdas x1.35) es anterior al arreglo del pivote fallido, al
     * orden por costo y al arreglo de la reubicacion de `k`. `DFB_SX_BFRT=1`
     * siempre, `DFB_SX_BFRT_M=k` solo si `m >= k`. */
    static const int bfrt_m = getenv("DFB_SX_BFRT_M") ? atoi(getenv("DFB_SX_BFRT_M")) : 0;
    /* Omision encendido desde 2026-09-30 (con el arreglo del ultimo quiebre);
     * `DFB_SX_BFRT=0` lo apaga. */
    static const bool bfrt_on = !(getenv("DFB_SX_BFRT") && std::string(getenv("DFB_SX_BFRT")) == "0");
    const bool bfrt = (bfrt_m > 0) ? (m >= bfrt_m) : bfrt_on;

    int j_in = -1; double best = 0.0, best_piv = 0.0;
    /* Omision: banda relativa. `DFB_SX_HARRISABS=1` activa el Harris de dos
     * pasadas con tolerancia absoluta, que medido en 8 instancias no gana con
     * claridad (ESTADO.md §9, opcion 3). */
    static const bool harris_rel = (getenv("DFB_SX_HARRISABS") == NULL);
    if (bfrt) {
        const std::vector<int>&    ix = tidx[r_out];
        const std::vector<double>& vl = tval[r_out];
        double amax = 0.0;
        for (size_t t = 0; t < ix.size(); ++t) {
            if (ix[t] >= na) break;
            const double av = std::fabs(vl[t]);
            if (av > amax) amax = av;
        }
        const double piv_min = (amax > 0.0) ? amax * 1e-11 : PIV_TOL;

        bfrt_cand.clear();
        for (size_t t = 0; t < ix.size(); ++t) {
            const int c = ix[t];
            if (c >= na) break;
            if (rowof[c] >= 0) continue;
            const double a = vl[t];
            const double av = std::fabs(a);
            if (av < piv_min || av < PIV_TOL) continue;
            const bool sube = atupper[c] ? (a > 0.0) : (a < 0.0);
            if (below != sube) continue;
            Cand cd;
            cd.c = c; cd.a = a;
            cd.ratio = std::fabs(d[c] / a);
            cd.rango = z[c].diam();
            bfrt_cand.push_back(cd);
        }
        static const bool dbg = (getenv("DFB_SX_BFRTDBG") != NULL);
        static long dbg_vacias = 0, dbg_noalcanza = 0, dbg_n = 0;
        if (bfrt_cand.empty()) {
            if (dbg && (++dbg_vacias % 1000) == 1) fprintf(stderr, "[bfrtdbg] sin candidatas (%ld)\n", dbg_vacias);
            infeasible = true; return false;
        }
        std::sort(bfrt_cand.begin(), bfrt_cand.end(), cand_menor);

        /* infactibilidad que hay que consumir, siempre positiva */
        const int cb = basic[r_out];
        double delta = below ? (z[cb].lb() - xB[r_out])
                             : (xB[r_out] - z[cb].ub());
        bfrt_flip.clear();
        size_t parada = bfrt_cand.size();
        for (size_t t = 0; t < bfrt_cand.size(); ++t) {
            const double rango = bfrt_cand[t].rango;
            if (!(rango < 1e18) || rango <= 0.0) { parada = t; break; }
            const double consumo = std::fabs(bfrt_cand[t].a) * rango;
            if (consumo < delta) {
                bfrt_flip.push_back(bfrt_cand[t].c);
                delta -= consumo;
                continue;
            }
            parada = t; break;
        }
        /* NI CAMBIANDO DE COTA TODAS LAS CANDIDATAS ALCANZA. En aritmetica
         * exacta seria infactibilidad primal, pero la prueba no vale con
         * tolerancias: las columnas con |a| bajo el umbral de pivote se
         * excluyeron de las candidatas y aun asi mueven `x_B(r)` (una `b` con
         * rango grande), y las violaciones de este caso son del orden de la
         * tolerancia. Declarar INFEASIBLE aca hacia que el 100 % de los LPs de
         * `ex7_2_6` y `ex9_2_5` terminaran sin cota (§9, 2026-09-30). Se hace
         * lo que la regla normal y las implementaciones robustas: pivotear en
         * el ULTIMO quiebre, cambiando de cota los anteriores, y dejar que la
         * infactibilidad —si la hay— la pruebe el certificado de Farkas.
         * `DFB_SX_BFRTINF=1` restaura el comportamiento anterior. */
        static const bool bfrt_inf_viejo = (getenv("DFB_SX_BFRTINF") != NULL);
        if (parada == bfrt_cand.size() && !bfrt_inf_viejo) {
            parada = bfrt_cand.size() - 1;
            if (!bfrt_flip.empty() && bfrt_flip.back() == bfrt_cand[parada].c) bfrt_flip.pop_back();
        }
        if (parada == bfrt_cand.size()) {
            /* ni cambiando de cota todas las candidatas alcanza: el LP es
             * primal-infactible. */
            if (dbg && (++dbg_noalcanza % 500) == 1 && dbg_n++ < 8) {
                const int cb2 = basic[r_out];
                double tot = 0.0;
                for (size_t t = 0; t < bfrt_cand.size(); ++t) tot += std::fabs(bfrt_cand[t].a) * bfrt_cand[t].rango;
                fprintf(stderr, "[bfrtdbg] no alcanza (%ld): fila %d basica %d en [%g,%g] xB=%g delta0=%g | %zu candidatas, consumo total %g\n",
                        dbg_noalcanza, r_out, cb2, z[cb2].lb(), z[cb2].ub(), xB[r_out],
                        below ? (z[cb2].lb() - xB[r_out]) : (xB[r_out] - z[cb2].ub()),
                        bfrt_cand.size(), tot);
                for (size_t t = 0; t < bfrt_cand.size() && t < 6; ++t)
                    fprintf(stderr, "   cand c=%d a=%g d=%g rango=%g atupper=%d z=[%g,%g] xN=%g\n",
                            bfrt_cand[t].c, bfrt_cand[t].a, d[bfrt_cand[t].c], bfrt_cand[t].rango,
                            (int)atupper[bfrt_cand[t].c], z[bfrt_cand[t].c].lb(), z[bfrt_cand[t].c].ub(), xN[bfrt_cand[t].c]);
            }
            infeasible = true; return false;
        }

        /* Desempate de Harris en el quiebre elegido. */
        best = bfrt_cand[parada].ratio;
        j_in = bfrt_cand[parada].c;
        best_piv = std::fabs(bfrt_cand[parada].a);
        if (harris_rel) {
            for (size_t t = parada+1; t < bfrt_cand.size(); ++t) {
                if (bfrt_cand[t].ratio > best*(1.0 + HARRIS_BAND)) break;
                const double av = std::fabs(bfrt_cand[t].a);
                if (av > best_piv) { best_piv = av; j_in = bfrt_cand[t].c; }
            }
        } else {
            /* Harris de dos pasadas con tolerancia ABSOLUTA sobre `d`: el paso
             * no puede pasar de `theta_max`, asi que ninguna de las que quedan
             * atras viola su signo en mas de `DUAL_TOL`. */
            double theta_max = POS_INFINITY;
            for (size_t t = parada; t < bfrt_cand.size(); ++t) {
                if (bfrt_cand[t].ratio > theta_max) break;
                const double tm = (std::fabs(d[bfrt_cand[t].c]) + DUAL_TOL) / std::fabs(bfrt_cand[t].a);
                if (tm < theta_max) theta_max = tm;
            }
            for (size_t t = parada; t < bfrt_cand.size(); ++t) {
                if (bfrt_cand[t].ratio > theta_max) break;
                const double av = std::fabs(bfrt_cand[t].a);
                if (av > best_piv || t == parada) { best_piv = av; j_in = bfrt_cand[t].c; best = bfrt_cand[t].ratio; }
            }
        }

        if (!bfrt_flip.empty()) {
            for (size_t t = 0; t < bfrt_flip.size(); ++t) {
                const int c = bfrt_flip[t];
                atupper[c] = atupper[c] ? 0 : 1;
            }
            n_flips += (long)bfrt_flip.size();
            xb_sucio = true;
            compute_basics(z);          /* las cotas cambiadas mueven x_B */
        }
    } else {
        const std::vector<int>&    ix = tidx[r_out];
        const std::vector<double>& vl = tval[r_out];

        /* Tolerancia de pivote RELATIVA a la fila: con coeficientes de ~1e72 un
         * umbral absoluto acepta ruido como pivote. */
        double amax = 0.0;
        for (size_t t = 0; t < ix.size(); ++t) {
            if (ix[t] >= na) break;
            const double av = std::fabs(vl[t]);
            if (av > amax) amax = av;
        }
        const double piv_min = (amax > 0.0) ? amax * 1e-11 : PIV_TOL;

        if (harris_rel) {
        for (size_t t = 0; t < ix.size(); ++t) {
            const int c = ix[t];
            if (c >= na) break;
            if (rowof[c] >= 0) continue;
            const double a = vl[t];
            const double av = std::fabs(a);
            if (av < piv_min || av < PIV_TOL) continue;
            const bool sube = atupper[c] ? (a > 0.0) : (a < 0.0);
            if (below != sube) continue;
            const double ratio = std::fabs(d[c] / a);
            /* Desempate de Harris: entre las candidatas dentro de una banda
             * relativa del minimo, la de mayor |a|. Sin esto cicla en
             * instancias degeneradas. */
            if (j_in < 0 || ratio < best * (1.0 - HARRIS_BAND)) {
                best = ratio; best_piv = av; j_in = c;
            } else if (ratio <= best * (1.0 + HARRIS_BAND) && av > best_piv) {
                best_piv = av; j_in = c;
                if (ratio < best) best = ratio;
            }
        }
        } else {
            /* HARRIS DE DOS PASADAS con tolerancia ABSOLUTA sobre `d`.
             * Primera: `theta_max = min (|d_j| + DUAL_TOL)/|a_j|`. Segunda:
             * entre las de cociente <= theta_max, la de mayor |a|. La banda
             * RELATIVA anterior dejaba pasar el paso hasta un 0.1 % mas alla
             * del quiebre minimo, y las columnas que quedaban atras terminaban
             * con el signo violado en hasta 1e-3*|d|: era el ~1 % de pasos
             * primales del dual puro (medido, 2026-09-30).
             * `DFB_SX_HARRISABS=1` la activa; no es la omision. */
            double theta_max = POS_INFINITY;
            for (size_t t = 0; t < ix.size(); ++t) {
                const int c = ix[t];
                if (c >= na) break;
                if (rowof[c] >= 0) continue;
                const double a = vl[t];
                const double av = std::fabs(a);
                if (av < piv_min || av < PIV_TOL) continue;
                const bool sube = atupper[c] ? (a > 0.0) : (a < 0.0);
                if (below != sube) continue;
                const double tm = (std::fabs(d[c]) + DUAL_TOL) / av;
                if (tm < theta_max) theta_max = tm;
            }
            for (size_t t = 0; t < ix.size(); ++t) {
                const int c = ix[t];
                if (c >= na) break;
                if (rowof[c] >= 0) continue;
                const double a = vl[t];
                const double av = std::fabs(a);
                if (av < piv_min || av < PIV_TOL) continue;
                const bool sube = atupper[c] ? (a > 0.0) : (a < 0.0);
                if (below != sube) continue;
                const double ratio = std::fabs(d[c] / a);
                if (ratio > theta_max) continue;
                if (j_in < 0 || av > best_piv) { best_piv = av; j_in = c; best = ratio; }
            }
        }
    }
    if (j_in < 0) { infeasible = true; return false; }

    if (std::fabs(d[j_in]) < DUAL_TOL) ++n_piv_degen;
    if (perturbado && std::fabs(d[j_in]) <= 3.0*pert_eps) ++n_piv_degen_pert;
    const int c_out = basic[r_out];

    /* Valor al que va a parar la basica que sale, y paso de la entrante. */
    const double v_out = below ? z[c_out].lb() : z[c_out].ub();
    const double xb_r  = xB[r_out];
    /* El valor de la entrante hay que leerlo de `xN`, NO recalcularlo de
     * `atupper`: `compute_basics` tiene un tratamiento especial para cotas
     * infinitas —usa la otra cota, o 0— y deja en `xN` el valor realmente
     * usado. Recalcularlo daba otro numero y la actualizacion incremental
     * quedaba desfasada (8 discrepancias de 528 en el test unitario). */
    const double x_in  = xN[j_in];

    if (!pivot(r_out, j_in)) {
        pivote_fallido = true;
        infeasible = false; return false;
    }
    atupper[c_out] = below ? 0 : 1;

    /* Actualizacion incremental de `x_B`, `O(m)` en vez de `O(m*nnz)`. */
    const bool inc = false;
    if (inc && (int)alpha_piv.size() == m+1) {
        const double ar = alpha_piv[r_out];
        if (ar != 0.0 && v_out == v_out && x_in == x_in &&
            v_out != POS_INFINITY && v_out != NEG_INFINITY &&
            x_in  != POS_INFINITY && x_in  != NEG_INFINITY) {
            const double delta = (xb_r - v_out) / ar;
            for (int i = 1; i <= m; ++i)
                if (i != r_out && alpha_piv[i] != 0.0)
                    xB[i] -= alpha_piv[i] * delta;
            xB[r_out] = x_in + delta;
            xN[c_out] = v_out;
            xN[j_in]  = 0.0;
            reparar_cotas_infinitas(z);
            xb_sucio = false;
        } else {
            xb_sucio = true;   /* cotas infinitas o pivote raro: recalcular */
        }
    }
    return true;
}

/* Paso primal: entra la no basica cuyo costo reducido viola el signo de su
 * cota, sale la basica que primero toca una cota. Si lo que primero se agota es
 * el rango de la entrante, es un *bound flip*: no hay pivote, solo cambia de
 * cota. */
bool DFBSimplex::primal_step(const IntervalVector& z, bool& unbounded) {
    perf::Cron _c(&perf::t_primal);
    unbounded = false;

    /* entrante: regla de Dantzig sobre los costos reducidos que violan */
    int q = -1; double bestd = 0.0;
    for (int c = 0; c < na; ++c) {
        if (rowof[c] >= 0) continue;
        const double dc = d[c];
        const double esc = 1.0 + std::fabs(dc);
        (void)esc;
        if (!atupper[c] && dc < -DUAL_TOL) {
            if (-dc > bestd) { bestd = -dc; q = c; }
        } else if (atupper[c] && dc > DUAL_TOL) {
            if (dc > bestd) { bestd = dc; q = c; }
        }
    }
    if (q < 0) return false;   /* factible dual: optimo */

    tab_column(q, colq);
    /* `x_B` SOLO SI ESTA SUCIO.
     *
     * El bucle compuesto llama `dual_step` antes que `primal_step`, y aquel ya
     * recalcula `x_B` al entrar. Si devuelve false es porque no habia
     * infactibilidad primal que reparar, o sea que no toco nada: `x_B` sigue al
     * dia. Recalcularlo aca era hacer el mismo `O(m*nnz)` dos veces por
     * iteracion —medido: 2.06 y 2.25 llamadas por pivote en `ex5_3_2` y
     * `house`, cuando corresponde una—.
     *
     * No tiene nada que ver con la actualizacion incremental del S44, que esta
     * apagada porque acumula redondeo distinto y mueve el arbol. Aca no cambia
     * ninguna cuenta: es la MISMA funcion sobre el MISMO estado, sencillamente
     * no se ejecuta dos veces. Los resultados son identicos bit a bit. */
    static const bool sin_skip = (getenv("DFB_SX_NOXBSKIP") != NULL);
    if (sin_skip || xb_sucio) compute_basics(z);

    /* signo del movimiento: en la cota inferior crece, en la superior baja */
    const double sgn = atupper[q] ? -1.0 : 1.0;

    /* ratio test primal: x_B(r) cambia en -a_rq * delta
     *
     * FILTRO DE ESTABILIDAD. Sin el, se elegia `r_out` por razon minima
     * aceptando cualquier `a != 0`, incluso 1e-18, y despues `pivot` lo
     * rechazaba por `|p| < PIV_TOL`. Ese fallo era indistinguible de «soy
     * optimo» y el bucle lo leia como OPTIMAL. El paso DUAL ya descartaba esos
     * pivotes (Harris); este camino —el que quedo por omision al no reubicar—
     * nunca recibio el mismo filtro. */
    /* ACTIVO POR OMISION. Medido sobre 139 instancias y 4 semillas: en las
     * instancias que pasan de 1 s, celdas x0.935 y cpu x0.805 contra
     * produccion, donde sin el filtro eran x1.139 y x0.890. `DFB_SX_NOPIVFIX=1`
     * restaura el comportamiento anterior para reproducir mediciones viejas. */
    double amax = 0.0;
    for (int r = 1; r <= m; ++r) {
        const double av = std::fabs(colq[r]);
        if (av > amax) amax = av;
    }
    const double piv_min = (amax > 0.0) ? amax*1e-11 : PIV_TOL;
    double lim = POS_INFINITY; int r_out = -1; bool sale_arriba = false;
    double piv_out = 0.0;
    for (int r = 1; r <= m; ++r) {
        const double a = colq[r];
        if (a == 0.0) continue;
        {
            const double av = std::fabs(a);
            if (av < piv_min || av < PIV_TOL) continue;
        }
        const int c = basic[r];
        const double lo = z[c].lb(), hi = z[c].ub();
        const double der = -a * sgn;            /* d x_B(r) / d|delta| */
        double cota;
        bool arriba;
        if (der > 0.0) { cota = (hi - xB[r]) / der; arriba = true; }
        else           { cota = (lo - xB[r]) / der; arriba = false; }
        if (cota < 0.0) cota = 0.0;
        const double av = std::fabs(a);
        if (r_out < 0 || cota < lim*(1.0 - 1e-9)) {
            lim = cota; r_out = r; sale_arriba = arriba; piv_out = av;
        } else if (cota <= lim*(1.0 + 1e-9) && av > piv_out) {
            r_out = r; sale_arriba = arriba; piv_out = av;
            if (cota < lim) lim = cota;
        }
    }

    /* bound flip: la entrante puede recorrer a lo sumo su propio rango */
    const double rango = z[q].diam();
    if (rango < lim) {
        atupper[q] = atupper[q] ? 0 : 1;
        /* Un *bound flip* mueve `x_B` sin pivotear, asi que hay que marcarlo
         * sucio a mano: `pivot()` es quien normalmente lo hace y aca no se
         * llama. Sin esto, la actualizacion incremental de `x_B` (§44) usaba un
         * valor obsoleto y la resolucion devolvia la respuesta de la anterior
         * —los pares (min, max) de una variable salian intercambiados—. Solo se
         * notaba con `reubicar = false`, que es el default del §25, porque ese
         * camino re-optimiza con el simplex PRIMAL y es el unico donde los
         * bound flips ocurren. */
        xb_sucio = true;
        return true;
    }
    if (r_out < 0) { unbounded = true; return false; }

    const int c_out = basic[r_out];
    if (!pivot(r_out, q)) {
        /* No es optimalidad: es un pivote que no se pudo aplicar. Informarlo
         * como tal, para que el bucle no lo confunda con la parada limpia. */
        pivote_fallido = true;
        return false;
    }
    atupper[c_out] = sale_arriba ? 1 : 0;
    return true;
}

/* Con la LU diferida, la refactoriza desde la base del tableau si algun pivote
 * la dejo atras. Una base que no factoriza deja `fact_ok = false`, y los
 * lectores caen al respaldo del tableau igual que con un update fallido. */
bool DFBSimplex::asegurar_fact() {
    if (fact_sucia) {
        fact_sucia = false;
        for (int r = 1; r <= m; ++r) cols_buf[r-1] = basic[r];
        fact_ok = fact.set_basis(cols_buf);
        ++n_refact_dif;
    }
    return fact_ok;
}

bool DFBSimplex::rayo_de_fila(int r, std::vector<double>& out) {
    asegurar_fact();
    if (!fact_ok || m <= 0 || r < 0 || r >= m) return false;
    fact.btran_row(r, &lam_buf[0]);
    out.resize(m);
    for (int j = 0; j < m; ++j) out[j] = lam_buf[j] * rscale[j];
    return true;
}

bool DFBSimplex::y_desde_factorizacion(int k, bool maximize) {
    perf::Cron _c(&perf::t_ydes);
    if (m <= 0) return false;
    const int r = rowof[k];
    if (r < 0) {                      /* k no basica: c_B = 0 */
        ++n_k_no_basica;
        std::fill(y.begin(), y.end(), 0.0);
        return true;
    }
    ++n_k_basica;
    if (!asegurar_fact()) return false;
    fact.btran_row(r - 1, &lam_buf[0]);      /* fila r-1 de B^-1 */
    const double ck = maximize ? -1.0 : 1.0;
    for (int j = 0; j < m; ++j) y[j] = ck * lam_buf[j];
    return true;
}

bool DFBSimplex::refrescar_desde_factorizacion(int k, bool maximize,
                                               const IntervalVector& z,
                                               bool reubicar) {
    perf::Cron _c(&perf::t_refresh);
    xb_sucio = true;   /* lo recalcula el mismo mas abajo */
    if (m <= 0 || !asegurar_fact()) return false;

    /* SONDA DE DERIVA: lo que daria el tableau solo, antes de pisarlo. Solo
     * sin reubicar, para que `x_B` se compare con la misma ubicacion. */
    const bool sonda = sonda_deriva && !reubicar && fact_ok && !fact_sucia;
    std::vector<double> d_tab, xb_tab;
    std::vector<char> atup_antes;
    if (sonda) {
        atup_antes = atupper;
        set_objective(k, maximize, false);
        compute_basics(z);
        d_tab = d; xb_tab = xB;
        atupper = atup_antes;
        xb_sucio = true;   /* compute_basics lo limpio: se deja como estaba */
    }

    /* y exacto */
    if (!y_desde_factorizacion(k, maximize)) return false;

    /* d = c - Abar^T y, recorriendo solo los y no nulos */
    const double ck = maximize ? -1.0 : 1.0;
    std::fill(d.begin(), d.end(), 0.0);
    d[k] = ck;
    for (int j = 0; j < m; ++j) {
        const double yj = y[j];
        if (yj == 0.0) continue;
        const std::vector<std::pair<int,double> >& a = arow[j];
        for (size_t t = 0; t < a.size(); ++t) d[a[t].first] -= yj * a[t].second;
    }
    for (int j = 0; j < m; ++j) ylam[j] = -y[j];
    /* Las basicas tienen costo reducido exactamente 0. OJO: forzarlo
     * ENMASCARA una `y` inconsistente con la base —si la factorizacion y
     * `basic[]` se desincronizan, el residuo no es 0 y el test de factibilidad
     * dual pasa espuriamente—. Con DFB_SX_CHECK=1 se mide ese residuo. */
    const bool chk = false;
    if (chk) {
        double peor = 0.0;
        for (int r = 1; r <= m; ++r) {
            const int b = basic[r];
            if (b >= na) continue;
            const double e = std::fabs(d[b]);
            if (e > peor) peor = e;
        }
        static double acum = 0.0; static long n = 0;
        if (peor > acum) acum = peor;
        if (++n % 20 == 0)
            fprintf(stderr, "[res] residuo d[basica] max=%.3g (peor acumulado %.3g)\n",
                    peor, acum);
    }
    for (int r = 1; r <= m; ++r) { const int b = basic[r]; if (b < na) d[b] = 0.0; }

    /* Restauracion MINIMA de la factibilidad dual: se mueven solo las columnas
     * cuyo signo esta violado y las degeneradas se dejan donde estan. Ver la
     * documentacion de reubicar en la cabecera.
     *
     * LA COLUMNA `k` DEL OBJETIVO NO SE REUBICA, y eso no es una variante sino
     * la correccion de un defecto. Con `c = +-e_k` y `k` no basica, mandarla a
     * la cota que le da el signo correcto deja la base dual-factible con `k`
     * FUERA de la base, o sea `c_B = 0`, `y = 0` y `gamma = 0`: el LP termina
     * declarando que el optimo es la cota actual de `z_k` —ninguna
     * contraccion— y ademas sin certificado que aplicar. Reubicar `k` nunca
     * puede servir, porque tanto minimizando como maximizando equivale a
     * afirmar que el optimo es la cota de entrada. Dejandola donde esta, su
     * costo reducido queda violado y algun paso tiene que meterla en la base,
     * que es de donde sale la cota.
     *
     * Ademas, cuando `k` es no basica es la UNICA columna con costo reducido no
     * nulo (`y = 0` implica `d = c`), asi que en ese caso —el 72-82 % de las
     * resoluciones— saltearla equivale a no reubicar nada.
     *
     * Medido en el banco (137 instancias, 4 semillas) contra reubicarla:
     * celdas x0.986/0.975/0.972 en los cortes > 0.1 s / > 1 s / > 5 s, con cpu
     * igual, y desaparecen las dos regresiones que tenia el camino dual
     * (`bearing` y `ship-1`). `DFB_SX_REUBK=1` restaura el comportamiento
     * anterior, solo para reproducir mediciones viejas. */
    static const bool cuenta_reub = (getenv("DFB_SX_CUENTAREUB") != NULL);
    int n_mov = 0, n_mov_k = 0;
    if (reubicar) n_mov = reubicar_minima(k, n_mov_k);
    if (cuenta_reub) {
        /* Cuando `k` es NO BASICA, `c_B = 0` y por lo tanto `y = 0`, asi que
         * `d = c = +-e_k`: la unica columna con costo reducido no nulo es `k`.
         * La reubicacion minima no toca ninguna otra. Se verifica en vez de
         * deducirse. */
        static long s_kb = 0, s_knb = 0, mov_kb = 0, mov_knb = 0, mov_k_tot = 0, peor_knb = 0;
        if (rowof[k] < 0) { ++s_knb; mov_knb += n_mov; if (n_mov - n_mov_k > peor_knb) peor_knb = n_mov - n_mov_k; }
        else              { ++s_kb;  mov_kb  += n_mov; }
        mov_k_tot += n_mov_k;
        if (((s_kb + s_knb) % 20000) == 0)
            fprintf(stderr, "[reub] solves=%ld | k NO basica: %ld solves, %.2f columnas movidas por solve "
                    "(de las cuales k: %.2f; maximo de OTRAS en un solve: %ld) | "
                    "k basica: %ld solves, %.2f movidas por solve\n",
                    s_kb + s_knb, s_knb, s_knb ? (double)mov_knb/s_knb : 0.0,
                    s_knb ? (double)mov_k_tot/s_knb : 0.0, peor_knb,
                    s_kb, s_kb ? (double)mov_kb/s_kb : 0.0);
    }

    /* x_B = -B^-1 (Abar x_N) */
    for (int c = 0; c < na; ++c) {
        if (rowof[c] >= 0) { xhat[c] = 0.0; continue; }
        double v = atupper[c] ? z[c].ub() : z[c].lb();
        if (v == POS_INFINITY || v == NEG_INFINITY) {
            const double o = atupper[c] ? z[c].lb() : z[c].ub();
            if (o == POS_INFINITY || o == NEG_INFINITY) v = 0.0;
            else { v = o; atupper[c] = atupper[c] ? 0 : 1; }
        }
        xhat[c] = v;
    }
    /* Despues de reparar las cotas infinitas, para que el signo de la
     * perturbacion coincida con la ubicacion definitiva. */
    if (reubicar) perturbar_costos();
    fact.A_times(&xhat[0], &vbuf[0]);
    for (int j = 0; j < m; ++j) vbuf[j] = -vbuf[j];
    fact.ftran_vec(&vbuf[0], &lam_buf[0]);
    for (int r = 1; r <= m; ++r) xB[r] = lam_buf[r-1];

    if (sonda) {
        ++dv_refrescos;
        bool distinto = false;
        for (int c = 0; c < na; ++c) {
            if (rowof[c] >= 0) continue;
            ++dv_d_nobas;
            const double e = std::fabs(d_tab[c] - d[c]);
            if (e > dv_d_max) dv_d_max = e;
            const bool v_lu  = atupper[c] ? d[c]     > DUAL_TOL : d[c]     < -DUAL_TOL;
            const bool v_tab = atupper[c] ? d_tab[c] > DUAL_TOL : d_tab[c] < -DUAL_TOL;
            if (v_lu != v_tab) { ++dv_d_signo; distinto = true; }
        }
        for (int r = 1; r <= m; ++r) {
            const int c = basic[r];
            const double lo = z[c].lb(), hi = z[c].ub();
            const double esc = 1.0 + ((std::fabs(lo) > std::fabs(hi)) ? std::fabs(lo) : std::fabs(hi));
            const double tol = FEAS_TOL * esc;
            ++dv_xb_bas;
            const double e = std::fabs(xb_tab[r] - xB[r]) / esc;
            if (e > dv_xb_max) dv_xb_max = e;
            const bool f_lu  = xB[r]     >= lo - tol && xB[r]     <= hi + tol;
            const bool f_tab = xb_tab[r] >= lo - tol && xb_tab[r] <= hi + tol;
            if (f_lu != f_tab) { ++dv_xb_fact; distinto = true; }
        }
        if (distinto) ++dv_refr_dist;
    }
    return true;
}

int DFBSimplex::reubicar_minima(int k, int& n_k) {
    /* El dual puro reubica tambien `k`: la base queda dual-factible con `k`
     * no basica en su cota, y los pasos duales la meten si hace falta. La
     * mixta la excluye (ver `refrescar_desde_factorizacion`). */
    static const bool reub_k_viejo = (getenv("DFB_SX_REUBK") != NULL);
    const bool no_reub_k = !(dual_puro || reub_k_viejo);
    int n_mov = 0; n_k = 0;
    for (int c = 0; c < na; ++c) {
        if (rowof[c] >= 0) continue;
        if (no_reub_k && c == k) continue;
        const char antes = atupper[c];
        if (d[c] < -DUAL_TOL)      atupper[c] = 1;
        else if (d[c] > DUAL_TOL)  atupper[c] = 0;
        if (atupper[c] != antes) { ++n_mov; if (c == k) ++n_k; }
    }
    return n_mov;
}

/* GUARDIA DE LA FILA OBJETIVO. En un tableau sin deriva, `d = c - Abar^T y`
 * con `y = -ylam`, y `d = 0` sobre las basicas. Se rehace el producto sobre
 * las filas escaladas y se exige: residuo relativo bajo `GUARDIA_TOL` y ningun
 * veredicto de signo distinto en las no basicas, que es lo que decide pivotes
 * y optimalidad. */
bool DFBSimplex::fila_objetivo_consistente(int k, bool maximize) {
    if ((int)dchk.size() != na) dchk.assign(na, 0.0);
    std::fill(dchk.begin(), dchk.end(), 0.0);
    dchk[k] = maximize ? -1.0 : 1.0;
    for (int j = 0; j < m; ++j) {
        const double yj = -ylam[j];
        if (yj == 0.0) continue;
        const std::vector<std::pair<int,double> >& a = arow[j];
        for (size_t t = 0; t < a.size(); ++t) dchk[a[t].first] -= yj * a[t].second;
    }
    double escala = 1.0, peor = 0.0;
    for (int c = 0; c < na; ++c) {
        const double v = std::fabs(dchk[c]);
        if (v > escala) escala = v;
        const double dc = (rowof[c] >= 0) ? 0.0 : d[c];
        const double e = std::fabs(dchk[c] - dc);
        if (e > peor) peor = e;
        if (rowof[c] < 0) {
            const bool v1 = atupper[c] ? d[c]    > DUAL_TOL : d[c]    < -DUAL_TOL;
            const bool v2 = atupper[c] ? dchk[c] > DUAL_TOL : dchk[c] < -DUAL_TOL;
            if (v1 != v2) return false;
        }
    }
    return peor <= GUARDIA_TOL * escala;
}

/* GUARDIA DE `x_B`: residuo de `Abar z = 0` con `z = (x_N, x_B)`, relativo a
 * la suma de los modulos de los terminos de cada fila. */
bool DFBSimplex::basicas_consistentes() {
    for (int j = 0; j < m; ++j) {
        double acc = 0.0, mag = 0.0;
        const std::vector<std::pair<int,double> >& a = arow[j];
        for (size_t t = 0; t < a.size(); ++t) {
            const int c = a[t].first;
            const double v = (rowof[c] >= 0) ? xB[rowof[c]] : xN[c];
            const double term = a[t].second * v;
            acc += term; mag += std::fabs(term);
        }
        if (std::fabs(acc) > GUARDIA_TOL * (1.0 + mag)) return false;
    }
    return true;
}

bool DFBSimplex::y_desde_tableau(int k, bool maximize) {
    if (m <= 0) return false;
    if (rowof[k] < 0) {               /* k no basica: c_B = 0, y = 0 EXACTO */
        ++n_k_no_basica;
        std::fill(y.begin(), y.end(), 0.0);
        return true;
    }
    if (perturbado || !fila_objetivo_consistente(k, maximize)) { ++n_guardia_fallo_fin; return false; }
    ++n_k_basica;
    for (int j = 0; j < m; ++j) y[j] = -ylam[j];
    return true;
}

bool DFBSimplex::refrescar_desde_tableau(int k, bool maximize,
                                         const IntervalVector& z, bool reubicar) {
    set_objective(k, maximize, false);     /* d e ylam desde las filas; x_B sucio */
    if (!fila_objetivo_consistente(k, maximize)) { ++n_guardia_fallo_d; return false; }
    if (reubicar) {
        int n_k;
        reubicar_minima(k, n_k);
        reparar_cotas_infinitas(z);
        perturbar_costos();
    }
    compute_basics(z);
    if (!basicas_consistentes()) { ++n_guardia_fallo_xb; return false; }
    ++n_guardia_ok;
    return true;
}

bool DFBSimplex::refrescar(int k, bool maximize, const IntervalVector& z, bool reubicar) {
    if (lu_tableau) {
        /* La ubicacion se restaura si la guardia falla despues de reubicar,
         * para que el camino de la LU parta del mismo estado. */
        const std::vector<char> atup = atupper;
        if (refrescar_desde_tableau(k, maximize, z, reubicar)) return true;
        atupper = atup;
        perturbado = false;
    }
    return refrescar_desde_factorizacion(k, maximize, z, reubicar);
}

/* PERTURBACION DE COSTOS AL REUBICAR (`DFB_SX_PERTURB=eps`).
 *
 * Con objetivo +-e_j y x_j basica en la fila r, tras reubicar queda
 * `d = -abar_r` sobre las no basicas: exactamente 0 en toda columna que no
 * toca la fila r. El ratio test dual elige la entrante por |d_c / a| y esas
 * columnas dan cociente 0, asi que entran primero: el objetivo dual no se
 * mueve y el simplex camina por bases degeneradas. Es el 78-88 % de pivotes
 * degenerados del S25, y la razon de los 30 pivotes por LP contra 7 del
 * camino primal, que no sufre este empate porque su ratio test es primal.
 *
 * El remedio estandar (SoPlex, HiGHS, Koberstein) es dar a los d_c nulos un
 * valor chico, aleatorio en [eps, 2 eps), con el signo que su ubicacion ya
 * cumple: la base sigue dual-factible para el costo perturbado y los empates
 * desaparecen. Como `d` se actualiza por operaciones de fila en `pivot`, la
 * perturbacion viaja con los pivotes igual que si el costo fuera otro.
 *
 * No toca la validez de la cota: `lambda` sale de la factorizacion con el
 * costo verdadero. Al llegar al optimo perturbado, `solve` rehace `d` exacto
 * y, si quedo alguna columna dual-infactible, la arreglan pasos primales; esos
 * pivotes se cuentan en `n_pert_limpieza`. */
void DFBSimplex::perturbar_costos() {
    static const double eps = getenv("DFB_SX_PERTURB")
                            ? atof(getenv("DFB_SX_PERTURB")) : 0.0;
    if (eps <= 0.0) return;
    int n = 0;
    for (int c = 0; c < na; ++c) {
        if (rowof[c] >= 0) continue;
        if (std::fabs(d[c]) >= eps) continue;
        pert_estado = pert_estado * 1664525u + 1013904223u;
        const double u = 1.0 + (double)(pert_estado >> 8) / 16777216.0;   /* [1,2) */
        d[c] = atupper[c] ? -eps*u : eps*u;
        ++n;
    }
    if (n > 0) { perturbado = true; pert_eps = eps; ++n_pert_lps; }
}



void DFBSimplex::primal_solution(const IntervalVector& z,
                                 std::vector<double>& out) {
    out.assign(na, 0.0);
    if (m <= 0) return;
    for (int c = 0; c < na; ++c) zs[c] = z[c] / cscale[c];
    compute_basics(zs);
    for (int c = 0; c < na; ++c) {
        if (rowof[c] >= 0) out[c] = xB[rowof[c]] * cscale[c];
        else               out[c] = atupper[c] ? z[c].ub() : z[c].lb();
    }
}

int DFBSimplex::dual_infeasibilities(int k, bool maximize) const {
    if (m <= 0 || k < 0 || k >= na) return 0;
    const double ck = maximize ? -1.0 : 1.0;
    const int r = rowof[k];
    if (r < 0)   /* k no basica: d = c, la unica columna con costo reducido es k */
        return ((!atupper[k] && ck < 0.0) || (atupper[k] && ck > 0.0)) ? 1 : 0;
    int n = 0;
    const std::vector<int>&    ix = tidx[r];
    const std::vector<double>& vl = tval[r];
    for (size_t t = 0; t < ix.size(); ++t) {
        const int c = ix[t];
        if (c >= na) break;              /* columnas lambda: sin ubicacion */
        if (rowof[c] >= 0) continue;
        const double dc = -ck * vl[t];
        if (std::fabs(dc) <= DUAL_TOL) continue;
        if ((!atupper[c] && dc < 0.0) || (atupper[c] && dc > 0.0)) ++n;
    }
    return n;
}

int DFBSimplex::primal_infeasibilities(int k, bool maximize,
                                       const IntervalVector& z) const {
    if (m <= 0 || k < 0 || k >= na) return 0;
    const int r = rowof[k];
    if (r < 0) return 0;      /* k no basica: no se reubica nada (ver reubicar) */
    const double ck = maximize ? -1.0 : 1.0;

    /* F = columnas que se reubicarian, y cuanto se mueve cada una. */
    std::vector<double> dx(na, 0.0);
    const std::vector<int>&    ix = tidx[r];
    const std::vector<double>& vl = tval[r];
    bool alguna = false;
    for (size_t t = 0; t < ix.size(); ++t) {
        const int c = ix[t];
        if (c >= na) break;
        if (rowof[c] >= 0 || c == k) continue;
        const double dc = -ck * vl[t];
        if (std::fabs(dc) <= DUAL_TOL) continue;
        const bool a_arriba = (dc < 0.0);
        if (a_arriba == (atupper[c] != 0)) continue;   /* ya esta donde toca */
        const double lo = z[c].lb()/cscale[c], hi = z[c].ub()/cscale[c];
        if (lo <= NEG_INFINITY || hi >= POS_INFINITY) continue;
        dx[c] = a_arriba ? (hi - lo) : (lo - hi);
        alguna = true;
    }
    if (!alguna) return 0;

    /* dx_B(i) = -sum_c T(i,c) dx_c, y se cuentan las filas que se salen. */
    int n = 0;
    for (int i = 1; i <= m; ++i) {
        const std::vector<int>&    jx = tidx[i];
        const std::vector<double>& jv = tval[i];
        double d = 0.0;
        for (size_t t = 0; t < jx.size(); ++t) {
            const int c = jx[t];
            if (c >= na) break;
            if (dx[c] != 0.0) d -= jv[t] * dx[c];
        }
        if (d == 0.0) continue;
        const int b = basic[i];
        const double lo = z[b].lb()/cscale[b], hi = z[b].ub()/cscale[b];
        const double mag = (std::fabs(lo) > std::fabs(hi)) ? std::fabs(lo) : std::fabs(hi);
        const double tol = FEAS_TOL * (1.0 + (mag < 1e300 ? mag : 0.0));
        const double v = xB[i] + d;
        if (v < lo - tol || v > hi + tol) ++n;
    }
    return n;
}

/* Copia de `y` tal como quedaria si se detuviera AHORA: sale de la
 * factorizacion y se devuelve a la escala original, igual que al terminar. Solo
 * se paga cuando hay traza pedida. */
void DFBSimplex::registrar_traza(int k, bool maximize) {
    if (!traza_y) return;
    if (!y_desde_factorizacion(k, maximize))
        for (int r = 0; r < m; ++r) y[r] = -ylam[r];
    std::vector<double> c(m);
    for (int r = 0; r < m; ++r) c[r] = y[r] * rscale[r];
    traza_y->push_back(c);
}

DFBSimplex::Status DFBSimplex::solve(int k, bool maximize,
                                     const IntervalVector& z, int max_iter,
                                     bool warm) {
    perf::Cron _c(&perf::t_total);
    if (perf::on()) ++perf::n_solve;
    ++n_solves;
    perturbado = false;
    /* Cada resolucion arranca con `x_B` sucio: el objetivo cambio y la
     * ubicacion de las no basicas pudo moverse, asi que el `x_B` que quedo del
     * solve anterior no sirve. Sin esto la resolucion devolvia la respuesta de
     * la anterior —los pares (min, max) de una misma variable salian
     * intercambiados—, que es la firma de un desfase de un paso. */
    xb_sucio = true;
    pivote_fallido = false;
    iters = 0;
    bnd = 0.0;
    if (traza_y) traza_y->clear();
    if (m <= 0 || k < 0 || k >= na) return SINGULAR;

    /* Caja escalada: z = C z_s, o sea z_s = z / C. Con C potencia de 2 la
     * division es exacta y no ensancha nada. */
    {   perf::Cron _e(&perf::t_escala);
        for (int c = 0; c < na; ++c) zs[c] = z[c] / cscale[c];
    }
    const IntervalVector& zz = zs;

    const bool base_lista_previa = base_lista;
    perf::Cron* _frio = (perf::on() && (!warm || !base_lista)) ? new perf::Cron(&perf::t_frio) : NULL;
    if (!warm || !base_lista) {
        /* --- base inicial: las m columnas de b. B = -I, asi que el tableau es
         * -[Abar | I] y la columna basica queda con coeficiente +1. --- */
        std::fill(rowof.begin(), rowof.end(), -1);
        /* TERCER lugar donde se borraba la ubicacion de las no basicas. Con
         * `DFB_SX_HEREDA=1` se conserva la heredada del nodo anterior para las
         * columnas de `x`; las de `b` pasan a ser basicas, asi que su valor es
         * irrelevante. */
        {   static const bool hereda_sv = (getenv("DFB_SX_HEREDA") != NULL);
            if (hereda_sv) { for (int c = nx; c < nc; ++c) atupper[c] = 0; }
            else           std::fill(atupper.begin(), atupper.end(), 0);
        }
        for (int r = 1; r <= m; ++r) {
            tidx[r].clear(); tval[r].clear();
            const std::vector<std::pair<int,double> >& a = arow[r-1];
            const int lamcol = na + (r-1);
            for (size_t t = 0; t < a.size(); ++t) {
                tidx[r].push_back(a[t].first);
                tval[r].push_back(-a[t].second);
            }
            tidx[r].push_back(lamcol); tval[r].push_back(-1.0);
            basic[r] = nx + (r-1);
            rowof[nx + (r-1)] = r;
        }
        if (lu_tableau) fact_sucia = true;     /* se factoriza solo si la guardia lo pide */
        else {
            for (int r = 1; r <= m; ++r) cols_buf[r-1] = basic[r];
            fact_ok = fact.set_basis(cols_buf);
            fact_sucia = false;
        }
        dw.assign(m+1, 1.0);          /* reinicio del marco de referencia devex */
        if (USA_DSE) dse.assign(m+1, 1.0);   /* B = -I: cada fila de B^-1 tiene norma 1 */
        base_lista = true;
        ++n_cold;
    } else {
        ++n_warm;
    }
    delete _frio;

    const bool tibio = warm && base_lista_previa;
    /* En tibio se recalcula todo desde la factorizacion: el tableau deriva y
     * `d` es el test de optimalidad. En frio el tableau esta recien construido
     * y no hace falta. */
    /* Al cambiar de objetivo NO se reubican las no basicas.
     *
     * Reubicarlas restaura la factibilidad DUAL gratis, que es lo que permite
     * correr el simplex dual. Pero la factibilidad PRIMAL de la base anterior se
     * conserva sola —el poliedro no se movio, solo el objetivo— y reubicar la
     * destruye. Conservandola, el punto queda primal-factible y dual-infactible,
     * que es exactamente el arranque de un simplex PRIMAL: los `primal_step` del
     * bucle compuesto deberian re-optimizar en pocas iteraciones, que es lo que
     * le da a SoPlex su x5.7 de arranque tibio (§24).
     *
     * Medido en el arbol: pivotes por LP de 29.91 a 7.32 en `ex5_3_2` y de 32.47
     * a 8.39 en `launch`, con la degeneracion cayendo de 78-88 % a 14-30 % —o
     * sea que la degeneracion no era propiedad de las relajaciones sino
     * consecuencia de arrancar primal-infactible—. En el banco entero, celdas
     * x0.925 contra produccion (antes x1.007) y, en las 67 instancias que
     * cuestan mas de 0.1 s, la ventaja pasa de -0.013 a +0.036. */
    /* Reubicar las no basicas restaura la factibilidad DUAL en cada cambio de
     * objetivo, de modo que las `2n` cotas son simplex dual de verdad y sus
     * rayos son certificados de Farkas legitimos. Por omision NO se reubica
     * —cuesta mas pivotes y el barrido lo midio peor—, pero entonces solo la
     * cota 0 es dual y en las demas el estado `INFEASIBLE` no esta respaldado.
     * `DFB_SX_PIVOTEO=dual` para medir ese eje. */
    const bool reubicar = pivoteo_dual;
    if (!tibio || !refrescar(k, maximize, zz, reubicar)) {
        set_objective(k, maximize, !tibio);
        /* Respaldo en tibio: si la estrategia reubica, se reubica aca tambien,
         * o el LP arrancaria dual-infactible sin que nadie lo sepa. */
        if (tibio && reubicar) { int nk; reubicar_minima(k, nk); }
    }
    { perf::Cron _r(&perf::t_reparar); reparar_cotas_infinitas(zz); }

    /* ENTRADA EXPLICITA DE `k` (`DFB_SX_ENTRAK=1`, solo con el dual puro).
     * Si `k` no es basica se la mete con un pivote en la fila de mayor
     * `|a_rk|` y se reubica todo: la base queda dual-factible CON `k` basica,
     * o sea que el certificado ya es una cota util desde el primer paso. Es el
     * `interchange` del DFB original. */
    static const bool entrak = (getenv("DFB_SX_ENTRAK") != NULL);
    if (dual_puro && entrak && rowof[k] < 0) {
        tab_column(k, colq);
        int r_k = -1; double amax = 0.0;
        for (int r = 1; r <= m; ++r)
            if (std::fabs(colq[r]) > amax) { amax = std::fabs(colq[r]); r_k = r; }
        if (r_k > 0 && amax >= PIV_TOL && pivot(r_k, k)) {
            ++iters; ++n_entrak;
            int nk; reubicar_minima(k, nk);
            reparar_cotas_infinitas(zz);
            xb_sucio = true;
        }
    }

    Status st = ITER_LIMIT;
    int iters_limpieza = -1;   /* pivote en que se quito la perturbacion */

    /* El estado inicial, con 0 pivotes, ya es una parada valida: da la cota
     * trivial. Es el primer punto de la curva anytime. */
    registrar_traza(k, maximize);

    /* Bucle COMPUESTO sobre el tableau: paso dual mientras haya
     * infactibilidad primal, paso primal mientras haya dual. */
        /* Respaldo sin factorizacion: bucle compuesto sobre el tableau. */
        perf::Cron _b(&perf::t_bucle);
        const bool tope_activo = pivoteo_dual && tope_anytime == tope_anytime;
        const bool obs_activo  = pivoteo_dual && obs_tope == obs_tope;
        obs_pivote = -1;
        /* SONDA DE MONOTONIA (`DFB_SX_MONO=1`): bajo la estrategia dual, cuantos
         * LPs dan algun paso primal, y si el valor basico de `z_k` —la cota
         * dual cuando la base es dual-factible— se mueve siempre hacia el
         * optimo mientras `k` es basica. Solo lee. */
        static const bool sonda_mono = (getenv("DFB_SX_MONO") != NULL);
        struct Mono {
            static long& lps()      { static long v = 0; return v; }
            static long& lps_kb0()  { static long v = 0; return v; }   /* k basica al empezar */
            static long& lps_pri()  { static long v = 0; return v; }   /* con algun paso primal */
            static long& lps_pri_kb0() { static long v = 0; return v; }
            static long& pasos_d()  { static long v = 0; return v; }
            static long& pasos_p()  { static long v = 0; return v; }
            static long& lps_nomono() { static long v = 0; return v; }
            static long& lps_mono_ok() { static long v = 0; return v; } /* todo dual, k basica siempre, monotono */
        };
        bool m_pri = false, m_nomono = false, m_kb_siempre = true;
        const bool m_kb0 = rowof[k] >= 0;
        double m_prev = std::numeric_limits<double>::quiet_NaN();
        if (sonda_mono && m_kb0) { if (xb_sucio) compute_basics(zz); m_prev = xB[rowof[k]]; }
        while (iters < max_iter) {
            ++n_vueltas;   /* vueltas del bucle, para separarlas de los pivotes */
            /* OBSERVACION: en que pivote cruza, sin cortar. */
            if (obs_activo && obs_pivote < 0 && rowof[k] >= 0) {
                if (xb_sucio) compute_basics(zz);
                const double v = xB[rowof[k]] * cscale[k];
                if (maximize ? (v <= obs_tope) : (v >= obs_tope)) obs_pivote = iters;
            }
            /* CORTE ANYTIME POR BRECHA. Ver `tope_anytime`. */
            if (tope_activo && rowof[k] >= 0) {
                if (xb_sucio) compute_basics(zz);
                const double v = xB[rowof[k]] * cscale[k];
                if (maximize ? (v <= tope_anytime) : (v >= tope_anytime)) {
                    ++n_anytime; st = ITER_LIMIT; break;
                }
            }
            bool infeasible = false, unbounded = false;
            if (dual_step(zz, infeasible)) {
                ++iters; registrar_traza(k, maximize);
                if (sonda_mono) {
                    ++Mono::pasos_d();
                    if (rowof[k] >= 0) {
                        if (xb_sucio) compute_basics(zz);
                        const double v = xB[rowof[k]];
                        const double tol = 1e-9 * (1.0 + std::fabs(v));
                        /* min z_k (c=+e_k): la cota sube; max: baja. En el
                         * espacio del simplex siempre se minimiza c.z, y el
                         * valor basico de z_k es la cota en unidades de z_k. */
                        if (m_prev == m_prev && (maximize ? (v > m_prev + tol) : (v < m_prev - tol))) m_nomono = true;
                        m_prev = v;
                    } else { m_kb_siempre = false; m_prev = std::numeric_limits<double>::quiet_NaN(); }
                }
                continue;
            }
            if (infeasible) { st = INFEASIBLE; break; }
            /* Diagnostico del dual puro: por que hay un paso primal. Violacion
             * de signo del costo reducido mas grande entre las no basicas,
             * y si ocurre antes o despues del primer paso dual del LP. */
            double m_viol = 0.0; int m_col = -1;
            if (sonda_mono && dual_puro) {
                for (int c = 0; c < na; ++c) {
                    if (rowof[c] >= 0) continue;
                    const double v = atupper[c] ? d[c] : -d[c];
                    if (v > m_viol) { m_viol = v; m_col = c; }
                }
            }
            const int iters_antes = iters;
            if (primal_step(zz, unbounded)) {
                if (sonda_mono && dual_puro) {
                    static long n_ini = 0, n_desp = 0, b[5] = {0,0,0,0,0}, n_k = 0, n_inf = 0;
                    if (iters_antes == 0) ++n_ini; else ++n_desp;
                    const int q = m_viol < 1e-8 ? 0 : m_viol < 1e-6 ? 1 : m_viol < 1e-4 ? 2 : m_viol < 1e-2 ? 3 : 4;
                    ++b[q];
                    if (m_col == k) ++n_k;
                    if (m_col >= 0 && (zz[m_col].lb() == NEG_INFINITY || zz[m_col].ub() == POS_INFINITY)) ++n_inf;
                    if (((n_ini + n_desp) % 50) == 0)
                        fprintf(stderr, "[porque] pasos primales en el dual puro=%ld | al empezar el LP: %ld  despues de pasos duales: %ld | "
                                "violacion de signo |d|: <1e-8 %ld  <1e-6 %ld  <1e-4 %ld  <1e-2 %ld  >=1e-2 %ld | la columna es k: %ld | con cota infinita: %ld\n",
                                n_ini + n_desp, n_ini, n_desp, b[0], b[1], b[2], b[3], b[4], n_k, n_inf);
                }
                ++iters; registrar_traza(k, maximize);
                if (sonda_mono) {
                    ++Mono::pasos_p(); m_pri = true;
                    if (rowof[k] >= 0) { if (xb_sucio) compute_basics(zz); m_prev = xB[rowof[k]]; }
                    else { m_kb_siempre = false; m_prev = std::numeric_limits<double>::quiet_NaN(); }
                }
                continue;
            }
            if (unbounded) { st = UNBOUNDED; break; }
            /* Un pivote fallido NO es optimalidad. Va tras la misma bandera
             * que el filtro de estabilidad para que el brazo base del barrido
             * sea el default real del repositorio. */
            if (pivote_fallido) { st = SINGULAR; break; }
            /* LIMPIEZA DE LA PERTURBACION: la base es optima para el costo
             * perturbado. Se rehace `d` exacto desde la factorizacion, SIN
             * reubicar, y si quedo alguna columna dual-infactible la arreglan
             * los `primal_step` de este mismo bucle. */
            if (perturbado) {
                perturbado = false;
                iters_limpieza = iters;
                if (refrescar(k, maximize, zz, false)) continue;
            }
            st = OPTIMAL; break;
        }
    if (sonda_mono && iters > 0) {
        ++Mono::lps();
        if (m_kb0) ++Mono::lps_kb0();
        if (m_pri) { ++Mono::lps_pri(); if (m_kb0) ++Mono::lps_pri_kb0(); }
        if (m_nomono) ++Mono::lps_nomono();
        if (!m_pri && m_kb0 && m_kb_siempre && !m_nomono) ++Mono::lps_mono_ok();
        if ((Mono::lps() % 20000) == 0)
            fprintf(stderr, "[mono] LPs con pivotes=%ld | k basica al empezar: %.1f%% | con algun paso primal: %.1f%% "
                    "(de los que empiezan con k basica: %.1f%%) | pasos primales %.1f%% de los pasos | "
                    "no monotonos en un tramo dual: %ld | enteramente duales, k basica siempre y monotonos: %.1f%%\n",
                    Mono::lps(), 100.0*Mono::lps_kb0()/Mono::lps(), 100.0*Mono::lps_pri()/Mono::lps(),
                    Mono::lps_kb0() ? 100.0*Mono::lps_pri_kb0()/Mono::lps_kb0() : 0.0,
                    100.0*Mono::pasos_p()/std::max(1L, Mono::pasos_p()+Mono::pasos_d()),
                    Mono::lps_nomono(), 100.0*Mono::lps_mono_ok()/Mono::lps());
    }
    if (iters_limpieza >= 0) n_pert_limpieza += iters - iters_limpieza;

    /* ITER_LIMIT NO invalida la base: se agoto el presupuesto, pero la base es
     * perfectamente valida —el tableau y la factorizacion se mantuvieron
     * consistentes en cada pivote— solo que no es la optima. Invalidarla hacia
     * que cortar por presupuesto obligara a la cota siguiente a arrancar en
     * frio, que es justo lo que arruina el uso anytime: detenerse antes salia
     * MAS caro que llegar al optimo. Las otras salidas si dejan la base dudosa:
     * INFEASIBLE y UNBOUNDED cortan sin pivotear y SINGULAR viene de un pivote
     * que no se pudo aplicar. */
    if (st != OPTIMAL && st != ITER_LIMIT) base_lista = false;
    if (st == OPTIMAL)          ++n_optimal;
    else if (st == INFEASIBLE)  ++n_infeasible;
    else if (st == ITER_LIMIT)  ++n_iterlimit;

    /* VERIFICACION DEL OPTIMO (`DFB_SX_VERIF=1`). Se rehacen `x_B` y `d`
     * exactos desde la factorizacion, SIN reubicar, y se mide si la base que
     * el bucle declaro optima es de verdad primal- y dual-factible. Sirve para
     * distinguir un optimo real de una salida por estado desfasado. */
    static const bool verif = (getenv("DFB_SX_VERIF") != NULL);
    if (verif && st == OPTIMAL && asegurar_fact()
        && refrescar_desde_factorizacion(k, maximize, zz, false)) {
        double pv = 0.0; int dinf = 0;
        for (int r = 1; r <= m; ++r) {
            const int c = basic[r];
            const double lo = zz[c].lb(), hi = zz[c].ub();
            double mag = 0.0;
            if (lo > NEG_INFINITY && std::fabs(lo) > mag) mag = std::fabs(lo);
            if (hi < POS_INFINITY && std::fabs(hi) > mag) mag = std::fabs(hi);
            const double esc = 1.0 + mag;
            double v = 0.0;
            if (xB[r] < lo)      v = (lo - xB[r]) / esc;
            else if (xB[r] > hi) v = (xB[r] - hi) / esc;
            if (v > pv) pv = v;
        }
        for (int c = 0; c < na; ++c) {
            if (rowof[c] >= 0) continue;
            if ((!atupper[c] && d[c] < -1e-7) || (atupper[c] && d[c] > 1e-7)) ++dinf;
        }
        ++n_verif;
        if (pv > 1e-6) ++n_verif_pinf;
        if (dinf > 0)  ++n_verif_dinf;
    }

    /* lambda: de la FACTORIZACION si se puede, y del tableau como respaldo.
     * Ver y_desde_factorizacion: el tableau deriva y eso debilita la cota
     * certificada, que es de donde se saca la contraccion. */
    if (sonda_deriva) {
        y_tab.resize(m);
        for (int r = 0; r < m; ++r) y_tab[r] = -ylam[r] * rscale[r];
    }
    if (!(lu_tableau && y_desde_tableau(k, maximize))
        && !y_desde_factorizacion(k, maximize)) {
        ++n_fact_fallo;
        for (int r = 0; r < m; ++r) y[r] = -ylam[r];
    }
    /* Vuelta a los multiplicadores del problema ORIGINAL: y = R y_s. Es lo que
     * necesita la certificacion, que trabaja sobre refA sin escalar. */
    for (int r = 0; r < m; ++r) y[r] *= rscale[r];

    /* cota flotante: el valor de z_k en la solucion basica.
     *
     * Si `k` quedo NO basica, la cota sale de `atupper[k]`, asi que la
     * ubicacion tiene que estar reparada: una cota infinita daria +-inf. En el
     * camino que recalcula, `compute_basics` reparaba de paso en cada
     * iteracion; en el incremental hay que hacerlo explicito aca. */
    perf::Cron* _cola = perf::on() ? new perf::Cron(&perf::t_cola) : NULL;
    reparar_cotas_infinitas(zz);
    if (rowof[k] >= 0) bnd = xB[rowof[k]];
    else               bnd = atupper[k] ? zz[k].ub() : zz[k].lb();
    bnd *= cscale[k];
    if (maximize) { /* se minimizo -z_k; la cota sigue siendo el valor de z_k */ }
    delete _cola;

    return st;
}

void perf_volcar() { perf::volcar(); }

} /* namespace ibex */
