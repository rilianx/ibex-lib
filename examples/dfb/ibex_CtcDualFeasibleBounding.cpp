#include "ibex_CtcDualFeasibleBounding.h"
#include <cmath>
#include <tuple>
#include <set>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <algorithm>

/* Reloj para la instrumentacion por fase. */
namespace {
inline double dfb_now() {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}
}

using namespace std;
using namespace ibex;

/* Los interruptores leen el entorno en su inicializador para que valgan en
 * cualquier binario (bench_dfb, diag_dfb e ibexopt), no solo donde haya codigo
 * que los consulte. Los valores por defecto son las decisiones tomadas. */
static bool env_on(const char* v)  { return getenv(v) != NULL; }
static int  env_int(const char* v, int d) { const char* e = getenv(v); return e ? atoi(e) : d; }

bool   CtcDFB::harris_tie_break = !env_on("DFB_NO_HARRIS");   /* M3(a) aceptado */
double CtcDFB::harris_rel_band  = 1e-3;
long   CtcDFB::n_no_candidate   = 0;
long   CtcDFB::n_inconclusive   = 0;
long   CtcDFB::n_init_failed    = 0;
long   CtcDFB::n_empty_by_bound = 0;
bool   CtcDFB::prove_empty_by_bound = !env_on("DFB_NO_BOUND_EMPTY");
bool   CtcDFB::pricing_fused        = !env_on("DFB_NO_FUSED");
bool   CtcDFB::pair_init            = !env_on("DFB_NO_PAIR_INIT");
/* Por omision: resultados IDENTICOS al tableau denso (0 diferencias en las 148
 * instancias del banco y celdas exactamente iguales en las 190 de ibexopt) y
 * estrictamente mas rapido: x2.21 en el banco de una caja, x3.12 en n>=80, y
 * Eiger-1000 pasa de std::bad_alloc a terminar. Con DFB_DENSE=1 se vuelve al
 * tableau denso para comparar. Ver MEDICIONES_FACTORIZACION.md. */
bool   CtcDFB::sparse_tableau       = !env_on("DFB_DENSE");
double CtcDFB::stop_tol             = (getenv("DFB_STOP_TOL") ?
                                       atof(getenv("DFB_STOP_TOL")) : 1e-6);
bool   CtcDFB::dual_pricing        = env_on("DFB_DUAL");
/* DFB_SIMPLEX_RULE (antes DFB_SIMPLEX, que se acepta como alias). El nombre
 * importa: no son dos metodos pegados. DFB siempre tuvo la maquinaria de un
 * simplex —tableau, pivotes, ratio test— pero no su criterio: elegia fila por
 * un "impacto" calculado con evaluacion de intervalos y paraba cuando ninguna
 * mostraba mejora estimada, que es ascenso por coordenadas y se traba en puntos
 * no optimos. Esto reemplaza ESA REGLA por la de un simplex dual de variables
 * acotadas, con costos reducidos como test de optimalidad. */
bool   CtcDFB::simplex_mode        = (env_on("DFB_SIMPLEX_RULE") ||
                                      env_on("DFB_SIMPLEX"));
bool   CtcDFB::sparse_check         = env_on("DFB_SPARSE_CHECK");
static const bool dfb_trace         = (getenv("DFB_TRACE") != NULL);
static const bool price_stats       = (getenv("DFB_PRICE_STATS") != NULL);
long   CtcDFB::n_check_difs         = 0;
int    CtcDFB::chk_last_j           = -1;
int    CtcDFB::chk_last_i           = -1;
double CtcDFB::check_worst          = 0.0;
bool   CtcDFB::float_pivoting       = !env_on("DFB_INTERVAL_PIVOT");  /* decidido: flotantes */
bool   CtcDFB::light_mode           = env_on("DFB_LIGHT");
int    CtcDFB::pricing_norm         = env_int("DFB_PRICING_NORM", 1);  /* M2(b) aceptado */
bool   CtcDFB::no_revive            = env_on("DFB_NO_REVIVE");
long   CtcDFB::n_certifications     = 0;
long   CtcDFB::n_price_gtotal       = 0;
long   CtcDFB::n_price_gnz          = 0;
long   CtcDFB::n_price_anz          = 0;
long   CtcDFB::n_sx_fallback        = 0;
long   CtcDFB::n_lam_nz             = 0;
long   CtcDFB::n_lam_cert           = 0;
long   CtcDFB::n_hib_barato         = 0;
long   CtcDFB::n_hib_simplex        = 0;
long   CtcDFB::n_hib_memo           = 0;
long   CtcDFB::n_elim_rows          = 0;
long   CtcDFB::n_elim_calls         = 0;
double CtcDFB::t_certify            = 0;
bool   CtcDFB::partial_pricing      = env_on("DFB_PARTIAL");
int    CtcDFB::partial_pricing_size = env_int("DFB_PARTIAL_SIZE", 0);
int    CtcDFB::partial_pricing_refresh = env_int("DFB_PARTIAL_REFRESH", 0);
bool   CtcDFB::partial_pricing_adaptive = env_on("DFB_PARTIAL_ADAPT");
long   CtcDFB::n_regenerations  = 0;
double CtcDFB::t_pricing        = 0;
double CtcDFB::t_ratio_test     = 0;
double CtcDFB::t_pivot          = 0;
double CtcDFB::t_bound          = 0;
double CtcDFB::t_regen          = 0;
double CtcDFB::t_init           = 0;
double CtcDFB::t_init_copy      = 0;
double CtcDFB::t_init_elim      = 0;
long   CtcDFB::n_inits          = 0;
long   CtcDFB::n_applied        = 0;

void CtcDFB::init(IntervalMatrix& A, IntervalVector& x_ref){
    /* El despacho va aqui y no solo en la propagacion: hay llamadores que
     * usan init() directamente (por ejemplo los tests). */
    if (float_pivoting) { init_float(A, x_ref); return; }

    double t_ini_init = dfb_now();
    n_inits++;
    //cout << "[CtcDFB] Initializing with k=" << k << ", A dimensions: " 
    //     << A.nb_rows() << "x" << A.nb_cols() << endl;
    double t_c = dfb_now();
    this->A.resize(A.nb_rows(), A.nb_cols());
    this->A = A;
    t_init_copy += dfb_now() - t_c;
    this->x_ref.resize(x_ref.size());
    this->x_ref = IntervalVector(x_ref);
    state = INITIAL;

    if (upper_contract) 
        for (int i = 0; i < A.nb_rows(); ++i) this->A[i][k] = -this->A[i][k];

    /* makeColumnIdentity lanza si no hay fila pivoteable en la columna k.
     * Antes esa excepcion escapaba hasta abortar el programa (por ejemplo en
     * ibexopt --filtering=dfb); ahora el contractor queda simplemente
     * inutilizable y la propagacion lo ignora. */
    cand_rows.clear(); cand_age = 0;
    init_ok = true;
    double t_e = dfb_now();
    try {
        makeColumnIdentity(this->A, k, true, 0);
    } catch (std::exception& e) {
        init_ok = false;
        n_init_failed++;
    }
    t_init_elim += dfb_now() - t_e;
    identity_rows.clear();
 

    //cout << "[CtcDFB] Initialization complete for k=" << k << endl;
    t_init += dfb_now() - t_ini_init;
}

void CtcDFB::init_from_lower(const CtcDFB& lower, IntervalVector& x_ref){
    if (float_pivoting) {   /* no hay matriz de intervalos que copiar */
        IntervalMatrix& R = const_cast<IntervalMatrix&>(*lower.refA_ref);
        init_float(R, x_ref);
        return;
    }
    double t_ini_init = dfb_now();
    n_inits++;

    this->A.resize(lower.A.nb_rows(), lower.A.nb_cols());
    this->A = lower.A;
    this->x_ref.resize(x_ref.size());
    this->x_ref = IntervalVector(x_ref);
    state = INITIAL;
    init_ok = lower.init_ok;
    identity_rows.clear();
    cand_rows.clear(); cand_age = 0;

    /* Ver la deduccion en el header: solo cambia el signo de la fila 0 fuera
     * de la columna k. */
    if (init_ok) {
        const int n = this->A[0].size();
        for (int i = 0; i < n; ++i)
            if (i != k) this->A[0][i] = -this->A[0][i];
    }

    t_init += dfb_now() - t_ini_init;
}

double CtcDFB::get_virtual_bound(){
    return (upper_contract)? -virtual_x.lb():virtual_x.lb();
}

// impacto de var en la contracción de k (solo si bound de var k es activo)
/* Impacto estimado de la variable var sobre la cota de k, usado por
 * CtcDFBPropag para REACTIVAR este contractor cuando la propagacion HC4
 * aprieta var. Es la via por la que DFB interactua con la propagacion, que es
 * el nucleo del diseno incremental.
 *
 * Ojo: en el camino flotante la fila gamma vive en Af[0], no en A[0]. Leer
 * A[0] ahi devolvia siempre 0 —init_float no llena A— y la reactivacion
 * quedaba MUERTA: 0 reactivaciones contra 875..3380 en el camino de
 * intervalos, sobre decenas de miles de llamadas. */
double CtcDFB::real_impact(Interval& x_k, int var, double eps){
    if (var==k) return 0.0;
    if (!init_ok || var < 0 || no_revive) return 0.0;

    const bool activa =
        (upper_contract  && (virtual_x.lb()+eps >= -x_k.ub())) ||
        (!upper_contract && (virtual_x.lb()+eps >=  x_k.lb()));
    if (!activa) return 0.0;

    if (float_pivoting) {
        if (sparse_tableau) {
            if (var >= St.nb_cols_A()) return 0.0;
            return std::fabs(St.g[var]);
        }
        if (Af.nb_rows() < 1 || var >= Af.nb_cols()) return 0.0;
        return std::fabs(Af[0][var]);
    }

    if (var >= A[0].size()) return 0.0;
    return std::abs(A[0][var].mid());
}

double CtcDFB::get_Aerror(){
    /* En modo flotante no hay ancho acumulado que limpiar: la certificacion se
     * hace desde refA en cada llamada. */
    if (float_pivoting) return 0.0;
    double error = 0;
    //A[0][0] + A[0][1] + A[0][2] + A[0][3]...
    for (int i=0; i<A[0].size(); i++)
        error += A[0][i].diam();
    
    return error;
}

void CtcDFB::regenerateA(IntervalMatrix& Aref){
    //cout << "[CtcDFB] Regenerating A" << endl;
    double t_ini_regen = dfb_now();
    n_regenerations++;
    cand_rows.clear(); cand_age = 0;   /* la lista de candidatos queda obsoleta */

    list<int> columns;
    for (int i = 0; i < A[0].size(); ++i)
        if (A[0][i]==Interval(0)) 
            columns.push_back(i);

    A=Aref;
    if (upper_contract) 
        for (int i = 0; i < A.nb_rows(); ++i) A[i][k] = -A[i][k];

    init_ok = true;
    try {
        makeColumnIdentity(A, k, true, 0);
        identity_rows.clear();
        makeColumnsIdentity(columns);
    } catch (std::exception& e) {
        init_ok = false;
        n_init_failed++;
    }

    t_regen += dfb_now() - t_ini_regen;
}

void CtcDFB::makeColumnsIdentity(std::list<int> columns) {
    std::set<int> identity_rows_local;

    columns.push_front(k);
    for (int col : columns) {
        // Buscar la mejor fila para hacer la columna col identidad
        int best_row = -1;
        double max_abs_val = 0.0;

        for (int i = 0; i < A.nb_rows(); ++i) {
            if (identity_rows_local.count(i)) continue;
            double val = A[i][col].mid();  // Puedes usar diam() si prefieres
            if (std::abs(val) > max_abs_val) {
                max_abs_val = std::abs(val);
                best_row = i;
            }
        }

        if (best_row == -1) {
            cout << "[CtcDFB::makeColumnsIdentity] No available row to pivot column " << col << endl;
            continue;
        }

        // Hacer la columna col una identidad en la fila best_row
        makeColumnIdentity(A, col, false, best_row);

        // Guardar fila usada
        identity_rows_local.insert(best_row);
        identity_rows[best_row] = col;
        
        //cout << "[CtcDFB::makeColumnsIdentity] Made column " << col << " identity at row " << best_row << endl;
    }
}


//get delta_impr
double CtcDFB::get_perc_impr(int k){
    //sum last k perc_imprs
    double sum = 0;
    int count = 0;
    for (auto it = perc_imprs.rbegin(); it != perc_imprs.rend() && count < k; ++it, ++count) {
        sum += *it;
    }
    if (count < k) 
        return 1;
    else 
        return sum;
}



/* ===================== Modo liviano: sin tableau ========================= *
 * Ver light_mode en el header. Las direcciones son las filas originales de
 * refA (compartidas, de solo lectura) y el estado por cota es gam y lam.
 * ======================================================================== */

/* Elige la fila de arranque para acotar k: la de mayor |R[r][k]|, y normaliza
 * para que gam_k = 1. O(m + n+m). Es el "cambio de objetivo" barato. */
void CtcDFB::set_bound_light(int k_new){
    double t0 = dfb_now();
    k = k_new;
    const Matrix& R = *refAf;
    const int m = R.nb_rows(), na = R.nb_cols();

    int best = -1; double bestv = 0.0;
    for (int r = 0; r < m; ++r) {
        const double v = std::fabs(R[r][k]);
        if (v > bestv) { bestv = v; best = r; }
    }
    if (best < 0 || bestv < 1e-12) { init_ok = false; t_init += dfb_now()-t0; return; }

    const double inv = 1.0 / R[best][k];
    for (int i = 0; i < na; ++i) gam[i] = R[best][i] * inv;
    gam[k] = 1.0;
    for (int j = 0; j < m; ++j) lam[j] = 0.0;
    lam[best] = inv;

    lam_nz.clear();
    state = INITIAL;
    init_ok = true;
    t_init += dfb_now() - t0;
}

void CtcDFB::init_light(const Matrix& R, IntervalMatrix& A, IntervalVector& x_ref){
    n_inits++;
    refA_ref = &A;
    refAf = &R;
    if (gam.size() != R.nb_cols()) gam.resize(R.nb_cols());
    if (lam.size() != R.nb_rows()) lam.resize(R.nb_rows());
    this->x_ref.resize(x_ref.size());
    this->x_ref = IntervalVector(x_ref);
    set_bound_light(k);
}

bool CtcDFB::certify_gamma_light(IntervalVector& gamma_out){
    double t0 = dfb_now();
    n_certifications++;

    const IntervalMatrix& Ri = *refA_ref;
    const int m = Ri.nb_rows(), na = Ri.nb_cols();

    lam_nz.clear();
    for (int j = 0; j < m; ++j) if (lam[j] != 0.0) lam_nz.push_back(j);

    gamma_out.resize(na);
    for (int i = 0; i < na; ++i) gamma_out[i] = Interval(0.0);
    for (size_t t = 0; t < lam_nz.size(); ++t) {
        const int j = lam_nz[t];
        const Interval l(lam[j]);
        for (int i = 0; i < na; ++i) gamma_out[i] += l * Ri[j][i];
    }
    t_certify += dfb_now() - t0;
    return !gamma_out[k].contains(0.0);
}

void CtcDFB::contract_light(IntervalVector& x_new){
    if (!init_ok || !refAf) { iters = 0; state = FINAL; return; }
    n_applied++;

    const Matrix& R = *refAf;
    const int m = R.nb_rows(), na = R.nb_cols();

    const int next = x_ref.size();
    if (work.size() != next) work.resize(next);
    for (int i = 0; i < nb_var && i < x_new.size(); ++i) work[i] = x_new[i];
    for (int i = nb_var; i < next; ++i) work[i] = x_ref[i];

    if (state == INITIAL) { perc_imprs.clear(); virtual_x = work[k]; }
    state = CONTRACTING;

    IntervalVector gc(1);
    iters = 0;

    while (max_iters == -1 || iters < max_iters) {
        /* --- pricing: mejor fila entre las originales, misma regla que el
         *     camino con tableau, con la norma de M2(b) si esta activa --- */
        double t0 = dfb_now();
        double best = 0.0, best_dir = 1.0; int bj = -1;
        for (int j = 0; j < m; ++j) {
            double acc_incr = 0.0, acc_decr = 0.0, nrm = 0.0;
            for (int i = 0; i < na; ++i) {
                const double lo = work[i].lb(), hi = work[i].ub();
                const double g = gam[i], a = R[j][i];
                double xi_incr, xi_decr;
                if (g >= 1e-5)       { xi_incr = hi; xi_decr = hi; }
                else if (g <= -1e-5) { xi_incr = lo; xi_decr = lo; }
                else if (a >= 1e-5)  { xi_incr = hi; xi_decr = lo; }
                else if (a <= -1e-5) { xi_incr = lo; xi_decr = hi; }
                else                 { xi_incr = 0.0; xi_decr = 0.0; }
                acc_incr += a * xi_incr;
                acc_decr += a * xi_decr;
                if (pricing_norm) nrm += a * a;
            }
            double gi = -acc_incr, gd = acc_decr;
            if (pricing_norm && nrm > 1e-300) {
                const double inv = 1.0/std::sqrt(nrm); gi *= inv; gd *= inv;
            }
            if (gi > best) { best = gi; bj = j; best_dir =  1.0; }
            if (gd > best) { best = gd; bj = j; best_dir = -1.0; }
        }
        t_pricing += dfb_now() - t0;
        if (bj < 0 || best <= 1e-6) { state = FINAL; break; }

        /* --- ratio test sobre la fila elegida --- */
        t0 = dfb_now();
        double min_mag = POS_INFINITY, best_piv = 0.0; int mi = -1;
        for (int i = 0; i < na; ++i) {
            const double a = R[bj][i];
            if (std::fabs(a) < 1e-9) continue;
            const double alpha = (gam[i] / a) * best_dir;
            if (alpha >= 0.0) continue;
            const double mag = std::fabs(alpha), piv = std::fabs(a);
            if (mi == -1 || mag < min_mag * (1.0 - harris_rel_band)) {
                min_mag = mag; best_piv = piv; mi = i;
            } else if (harris_tie_break && mag <= min_mag*(1.0+harris_rel_band)
                       && piv > best_piv) { best_piv = piv; mi = i; }
        }
        t_ratio_test += dfb_now() - t0;
        if (mi == -1) { state = FINAL; break; }   /* sin fila de bloqueo: no se concluye */

        /* --- pivote: gam += alpha.R[bj], lam += alpha.e_bj --- */
        t0 = dfb_now();
        const double alpha = -(gam[mi] / R[bj][mi]);
        for (int i = 0; i < na; ++i) gam[i] += alpha * R[bj][i];
        gam[mi] = 0.0;                     /* solo en gam; lam queda veraz */
        lam[bj] += alpha;
        /* gamma_k pudo cambiar: se renormaliza para que la cota no degenere.
         * Escalar gam y lam por el mismo factor no cambia la relacion. */
        if (std::fabs(gam[k]) < 1e-9) { state = FINAL; t_pivot += dfb_now()-t0; break; }
        if (std::fabs(gam[k] - 1.0) > 1e-12) {
            const double sc = 1.0 / gam[k];
            for (int i = 0; i < na; ++i) gam[i] *= sc;
            for (int j = 0; j < m;  ++j) lam[j] *= sc;
            gam[k] = 1.0;
        }
        t_pivot += dfb_now() - t0;
        ++iters;
    }

    if (certify_gamma_light(gc)) {
        double t0 = dfb_now();
        Interval nb = gaussSeidel(work, k, gc);
        t_bound += dfb_now() - t0;
        if (work.is_empty()) { x_new.set_empty(); state = FINAL; return; }
        virtual_x = Interval(nb.lb(), work[k].ub());
    }
    write_back(x_new);
}

/* ===================== Pivoteo en punto flotante ======================= *
 *
 * Ver la explicacion de float_pivoting en el header. El tableau Af es
 * [A | I] con m filas y (nb_var+m + m) columnas; las ultimas m columnas de la
 * fila j son el lambda de esa fila.
 * ======================================================================= */

static const double DFB_F_TOL = 1e-9;   /* tolerancia relativa del ratio test */

void CtcDFB::init_float(IntervalMatrix& A, IntervalVector& x_ref){
    if (sparse_tableau) { init_sparse(A, x_ref); return; }
    double t_ini_init = dfb_now();
    n_inits++;
    refA_ref = &A;

    const int m  = A.nb_rows();
    const int na = A.nb_cols();          /* nb_var + m */

    this->x_ref.resize(x_ref.size());
    this->x_ref = IntervalVector(x_ref);

    Af.resize(m, na + m);
    for (int i = 0; i < m; ++i) {
        for (int c = 0; c < na; ++c) Af[i][c] = A[i][c].mid();
        for (int c = 0; c < m;  ++c) Af[i][na + c] = (i == c) ? 1.0 : 0.0;
        /* refA_k: columna k negada cuando se contrae la cota superior. */
        if (upper_contract) Af[i][k] = -Af[i][k];
    }

    lam_nz.clear();
    identity_rows.clear();
    cand_rows.clear(); cand_age = 0;
    state = INITIAL;
    init_ok = make_column_identity_f(k, 0, true);
    if (!init_ok) n_init_failed++;

    t_init += dfb_now() - t_ini_init;
}

/* ===================== camino de tableau disperso =======================
 *
 * Todas estas funciones son la traduccion literal de sus equivalentes densas,
 * recorriendo solo las entradas no nulas. Donde el camino denso multiplica por
 * un cero, aca simplemente no se visita la entrada: el resultado es el mismo
 * bit a bit, porque sumar 0.0 no cambia el acumulador y el orden de columnas
 * crecientes se conserva. test_dfbsparse.cpp lo verifica con tolerancia 0.
 */

/* Compara el tableau disperso con el denso, entrada por entrada. */
void CtcDFB::compare_tableaus(const char* donde){
    static long llamadas = 0;
    if (Af.nb_rows() != St.nb_rows()) {
        printf("[chk] %s: SIN COMPARAR (Af %dx%d vs St %dx%d)\n", donde,
               Af.nb_rows(), Af.nb_cols(), St.nb_rows(), St.nb_cols());
        return;
    }
    /* Invariantes de la representacion dispersa: indices estrictamente
     * crecientes y sin repetir. Un indice repetido haria que entry() devuelva
     * solo el primero y explicaria una divergencia de pocos ulps. */
    for (int j = 1; j < St.nb_rows(); ++j) {
        const std::vector<int>& ix = St.idx(j);
        for (size_t t = 1; t < ix.size(); ++t)
            if (ix[t] <= ix[t-1]) {
                printf("[chk] %s: fila %d indices no crecientes en t=%zu: %d luego %d\n",
                       donde, j, t, ix[t-1], ix[t]);
                ++n_check_difs;
                break;
            }
    }
    if (++llamadas % 500 == 0)
        printf("[chk] %ld comparaciones de tableau, difs=%ld peor=%.3g\n",
               llamadas, n_check_difs, check_worst);
    const int m = St.nb_rows(), nc = St.nb_cols();
    int reportados = 0;
    for (int j = 0; j < m; ++j)
        for (int c = 0; c < nc; ++c) {
            const double d = std::fabs(St.entry(j, c) - Af[j][c]);
            if (d > check_worst) check_worst = d;
            if (d != 0.0) {
                ++n_check_difs;
                if (n_check_difs == 1 && getenv("DFB_CHECK_ABORT")) {
                    printf("[chk] PRIMERA DIFERENCIA en %s: (%d,%d) "
                           "disperso=%.17g denso=%.17g\n",
                           donde, j, c, St.entry(j, c), Af[j][c]);
                    printf("[chk] contexto: k=%d up=%d m=%d na=%d nc=%d "
                           "ultimo pivote j=%d i=%d\n", k, (int)upper_contract,
                           m, St.nb_cols_A(), nc, chk_last_j, chk_last_i);
                    printf("[chk] fila %d completa (disperso | denso):\n", j);
                    for (int cc = 0; cc < nc; ++cc) {
                        const double sv = St.entry(j, cc), dv = Af[j][cc];
                        if (sv != 0.0 || dv != 0.0)
                            printf("   c=%-4d %+.17g   %+.17g%s\n", cc, sv, dv,
                                   (sv != dv) ? "   <---" : "");
                    }
                    fflush(stdout);
                    exit(7);
                }
                if (reportados < 3 && n_check_difs < 40) {
                    printf("[chk] %s: (%d,%d) disperso=%.17g denso=%.17g\n",
                           donde, j, c, St.entry(j, c), Af[j][c]);
                    ++reportados;
                }
            }
        }
}

void CtcDFB::set_basic(int row, int col){
    if (brow.empty()) return;
    const int viejo = brow[row];
    if (viejo >= 0 && viejo < (int)colbasic.size()) colbasic[viejo] = 0;
    brow[row] = col;
    if (col >= 0 && col < (int)colbasic.size()) colbasic[col] = 1;
}

/* Base inicial: la identidad artificial (fila j tiene basica la columna na+j),
 * con el intercambio de filas del init aplicado, y luego la columna k entrando
 * en la fila 0. */
void CtcDFB::init_basis_bookkeeping(int na_cols, int m_rows, int k_col,
                                    int swapped){
    brow.assign(m_rows, 0);
    colbasic.assign(na_cols + m_rows, 0);
    for (int j = 0; j < m_rows; ++j) {
        brow[j] = na_cols + j;
        colbasic[na_cols + j] = 1;
    }
    if (swapped > 0) std::swap(brow[0], brow[swapped]);
    set_basic(0, k_col);
    zN.assign(na_cols, 0.0);
}

/* Pricing dual: se elige la fila cuya variable basica mas viola sus cotas.
 * Ver la documentacion de dual_pricing en la cabecera. */
double CtcDFB::bound_proxy() const {
    const int na = St.nb_cols_A();
    double S_ub = 0.0, S_lb = 0.0;
    for (int i = 0; i < na; ++i) {
        if (i == k) continue;
        const double g = St.g[i];
        if (g == 0.0) continue;
        const double a = g * work[i].lb(), b = g * work[i].ub();
        if (a < b) { S_lb += a; S_ub += b; }
        else       { S_lb += b; S_ub += a; }
    }
    return (gamma_k_target > 0) ? -S_ub : S_lb;
}

int CtcDFB::largest_impact_dual(int& j_out, double& delta_out, double& dir_out){
    const int m  = St.nb_rows();
    const int na = St.nb_cols_A();
    if ((int)brow.size() != m) { j_out = -1; delta_out = 0.0; dir_out = 1.0; return -1; }

    /* Asignacion UNICA de las no basicas: la que alcanza la cota, o sea la que
     * dicta el signo de gamma. Es la misma que usa gaussSeidel para evaluar la
     * cota, y por eso el valor basico que sale de aca es el de la solucion
     * basica que define esa cota. */
    for (int i = 0; i < na; ++i) {
        if (colbasic[i]) { zN[i] = 0.0; continue; }
        zN[i] = (St.g[i] > 0.0) ? work[i].ub() : work[i].lb();
    }

    double best = 0.0; int jb = -1; double dirb = 1.0;
    for (int j = 1; j < m; ++j) {
        const std::vector<int>&    ix = St.idx(j);
        const std::vector<double>& vl = St.val(j);
        double acc = 0.0, nrm = 0.0;
        for (size_t t = 0; t < ix.size(); ++t) {
            const int i = ix[t];
            if (i >= na) break;
            const double a = vl[t];
            if (pricing_norm) nrm += a * a;
            if (colbasic[i]) continue;
            acc += a * zN[i];
        }
        const double zB = -acc;

        /* Cotas de la variable basica de la fila. Una columna artificial es una
         * variable que debe valer 0, asi que su violacion es |zB|. */
        const int c = brow[j];
        double lo, hi;
        if (c < na) { lo = work[c].lb(); hi = work[c].ub(); }
        else        { lo = 0.0; hi = 0.0; }

        double viol; double dir;
        if (zB > hi)      { viol = zB - hi; dir =  1.0; }
        else if (zB < lo) { viol = lo - zB; dir = -1.0; }
        else continue;

        double score = viol;
        if (pricing_norm && nrm > 1e-300) score *= 1.0 / std::sqrt(nrm);
        if (score > best) { best = score; jb = j; dirb = dir; }
    }

    if (jb < 0) { j_out = -1; delta_out = 0.0; dir_out = 1.0; return -1; }
    j_out = jb; dir_out = dirb; delta_out = best;
    if (delta_out <= stop_tol) j_out = -1;
    return j_out;
}

void CtcDFB::init_sparse(IntervalMatrix& A, IntervalVector& x_ref){
    double t_ini_init = dfb_now();
    n_inits++;
    refA_ref = &A;

    const int m  = A.nb_rows();
    const int na = A.nb_cols();

    this->x_ref.resize(x_ref.size());
    this->x_ref = IntervalVector(x_ref);

    /* Abar en flotantes, con la columna k negada si se contrae la cota
     * superior (refA_k). Se construye una vez y se descarta: el tableau
     * disperso se queda con sus filas. */
    Matrix Abar(m, na);
    for (int i = 0; i < m; ++i) {
        for (int c = 0; c < na; ++c) Abar[i][c] = A[i][c].mid();
        if (upper_contract) Abar[i][k] = -Abar[i][k];
    }
    St.build(Abar);

    if (sparse_check) {
        Af.resize(m, na + m);
        for (int i = 0; i < m; ++i) {
            for (int c = 0; c < na; ++c) Af[i][c] = Abar[i][c];
            for (int c = 0; c < m;  ++c) Af[i][na + c] = (i == c) ? 1.0 : 0.0;
        }
    }

    lam_nz.clear();
    identity_rows.clear();
    cand_rows.clear(); cand_age = 0;
    state = INITIAL;
    int swapped = -1;
    init_ok = St.make_identity(k, 0, true, DFB_F_TOL, &swapped);
    if (init_ok) init_basis_bookkeeping(na, m, k, swapped);
    if (sparse_check) {
        const bool okd = make_column_identity_f(k, 0, true);
        if (okd != init_ok) { ++n_check_difs; }
        compare_tableaus("init");
    }
    if (!init_ok) n_init_failed++;

    t_init += dfb_now() - t_ini_init;
}

void CtcDFB::init_from_lower_sparse(const CtcDFB& lower, IntervalVector& x_ref){
    double t_ini_init = dfb_now();
    n_inits++;
    refA_ref = lower.refA_ref;

    St.copy_from(lower.St);

    this->x_ref.resize(x_ref.size());
    this->x_ref = IntervalVector(x_ref);

    lam_nz.clear();
    identity_rows.clear();
    cand_rows.clear(); cand_age = 0;
    state = INITIAL;
    init_ok = lower.init_ok;

    if (sparse_check && lower.Af.nb_rows() == St.nb_rows()) {
        Af.resize(lower.Af.nb_rows(), lower.Af.nb_cols());
        Af = lower.Af;
    }

    if (init_ok) {
        brow = lower.brow;
        colbasic = lower.colbasic;
        zN.assign(St.nb_cols_A(), 0.0);
        St.negate_row0();
        St.g[k] = 1.0;
        if (sparse_check && Af.nb_rows() == St.nb_rows()) {
            const int nc = Af.nb_cols();
            for (int c = 0; c < nc; ++c) Af[0][c] = -Af[0][c];
            Af[0][k] = 1.0;
            compare_tableaus("init_from_lower");
        }
    }

    t_init += dfb_now() - t_ini_init;
}

/* Pricing disperso: identico a largest_impact_f, saltando los a == 0. */
int CtcDFB::largest_impact_s(int& j_out, double& delta_out, double& dir_out){
    const int m  = St.nb_rows();
    const int na = St.nb_cols_A();

    double best_incr = 0.0, best_decr = 0.0;
    int    j_incr = 0, j_decr = 0;

    for (int j = 1; j < m; ++j) {
        const std::vector<int>&    ix = St.idx(j);
        const std::vector<double>& vl = St.val(j);
        double acc_incr = 0.0, acc_decr = 0.0, nrm = 0.0;
        for (size_t t = 0; t < ix.size(); ++t) {
            const int i = ix[t];
            if (i >= na) break;          /* la parte lambda va al final */
            const double a = vl[t];
            const double lo = work[i].lb(), hi = work[i].ub();
            const double g = St.g[i];
            double xi_incr, xi_decr;
            if (g >= 1e-5)       { xi_incr = hi; xi_decr = hi; }
            else if (g <= -1e-5) { xi_incr = lo; xi_decr = lo; }
            else if (a >= 1e-5)  { xi_incr = hi; xi_decr = lo; }
            else if (a <= -1e-5) { xi_incr = lo; xi_decr = hi; }
            else                 { xi_incr = 0.0; xi_decr = 0.0; }
            acc_incr += a * xi_incr;
            acc_decr += a * xi_decr;
            if (pricing_norm) nrm += a * a;
        }
        double gi = -acc_incr, gd = acc_decr;
        if (pricing_norm && nrm > 1e-300) {
            const double inv = 1.0 / std::sqrt(nrm);
            gi *= inv; gd *= inv;
        }
        if (gi > best_incr)      { best_incr = gi;      j_incr = j; }
        if (gd > best_decr)      { best_decr = gd;      j_decr = j; }
    }
    if (best_incr > best_decr) { delta_out = best_incr; j_out = j_incr; dir_out =  1.0; }
    else                       { delta_out = best_decr; j_out = j_decr; dir_out = -1.0; }

    if (pricing_norm && j_out > 0) {
        const std::vector<int>&    ix = St.idx(j_out);
        const std::vector<double>& vl = St.val(j_out);
        double acc = 0.0;
        for (size_t t = 0; t < ix.size(); ++t) {
            const int i = ix[t];
            if (i >= na) break;
            const double a = vl[t];
            const double lo = work[i].lb(), hi = work[i].ub();
            const double g = St.g[i];
            double xi;
            if (g >= 1e-5)       xi = hi;
            else if (g <= -1e-5) xi = lo;
            else if (a >= 1e-5)  xi = (dir_out > 0) ? hi : lo;
            else if (a <= -1e-5) xi = (dir_out > 0) ? lo : hi;
            else                 xi = 0.0;
            acc += a * xi;
        }
        delta_out = (dir_out > 0) ? -acc : acc;
    }
    if (delta_out <= stop_tol) j_out = -1;
    return j_out;
}

/* Ratio test disperso: solo las entradas no nulas pueden ser candidatas,
 * porque el denso descarta |a| < DFB_F_TOL. */
int CtcDFB::ratio_test_s(int j, double dir, double& alpha_out){
    const int na = St.nb_cols_A();
    const std::vector<int>&    ix = St.idx(j);
    const std::vector<double>& vl = St.val(j);

    double min_mag = POS_INFINITY, best_piv = 0.0;
    int    min_i = -1;
    double a_min = 0.0;

    for (size_t t = 0; t < ix.size(); ++t) {
        const int i = ix[t];
        if (i >= na) break;
        const double a = vl[t];
        if (std::fabs(a) < DFB_F_TOL) continue;
        const double alpha = (St.g[i] / a) * dir;
        if (alpha >= 0.0) continue;
        const double mag = std::fabs(alpha), piv = std::fabs(a);
        if (min_i == -1 || mag < min_mag * (1.0 - harris_rel_band)) {
            min_mag = mag; best_piv = piv; min_i = i; a_min = a;
        } else if (harris_tie_break && mag <= min_mag * (1.0 + harris_rel_band)
                   && piv > best_piv) {
            best_piv = piv; min_i = i; a_min = a;
        }
    }
    if (min_i == -1) return -1;
    alpha_out = -(St.g[min_i] / a_min);
    return min_i;
}

/* gamma = lambda^T . refA, con lambda leido de la fila 0 densa. */
bool CtcDFB::certify_gamma_s(IntervalVector& gamma_out,
                             const std::vector<double>* lam_use){
    double t0 = dfb_now();
    n_certifications++;

    const IntervalMatrix& R = *refA_ref;
    const int m = R.nb_rows(), na = R.nb_cols();
    const std::vector<double>& L = lam_use ? *lam_use : St.lam;

    lam_nz.clear();
    for (int j = 0; j < m; ++j)
        if (L[j] != 0.0) lam_nz.push_back(j);

    gamma_out.resize(na);
    for (int i = 0; i < na; ++i) gamma_out[i] = Interval(0.0);
    for (size_t t = 0; t < lam_nz.size(); ++t) {
        const int j = lam_nz[t];
        const Interval lam_j(L[j]);
        for (int i = 0; i < na; ++i) gamma_out[i] += lam_j * R[j][i];
    }
    if (upper_contract && gamma_k_target > 0) gamma_out[k] = -gamma_out[k];

    t_certify += dfb_now() - t0;
    return !gamma_out[k].contains(0.0);
}

bool CtcDFB::certify_from_lambda(const std::vector<double>& lam,
                                 IntervalVector& gamma_out){
    double t0 = dfb_now();
    n_certifications++;

    const IntervalMatrix& R = *refA_ref;
    const int mm = R.nb_rows(), na = R.nb_cols();

    lam_nz.clear();
    for (int j = 0; j < mm && j < (int)lam.size(); ++j)
        if (lam[j] != 0.0) lam_nz.push_back(j);

    gamma_out.resize(na);
    for (int i = 0; i < na; ++i) gamma_out[i] = Interval(0.0);
    for (size_t t = 0; t < lam_nz.size(); ++t) {
        const int j = lam_nz[t];
        const Interval lj(lam[j]);
        for (int i = 0; i < na; ++i) gamma_out[i] += lj * R[j][i];
    }

    n_lam_nz += (long)lam_nz.size();
    ++n_lam_cert;

    t_certify += dfb_now() - t0;
    return !gamma_out[k].contains(0.0);
}

void CtcDFB::contract_simplex(IntervalVector& x_new){
    if (!init_ok || sx == NULL) { iters = 0; state = FINAL; return; }

    /* HIBRIDO (DFB_SX_HIBRIDO=1): primero la regla de impacto, que es ~21 veces
     * mas barata en pivotes, y el simplex SOLO donde aquella no rindio.
     *
     * La logica: la regla de impacto se traba en un punto no optimo, pero
     * cuando NO se traba llega a la misma cota que el simplex —el grupo A son
     * 64 de 147 instancias— asi que pagar el simplex ahi es puro costo. Y el
     * estancamiento deja una firma barata de detectar: contrajo poco o nada
     * teniendo la variable todavia ancha. */
    static const bool hibrido = (getenv("DFB_SX_HIBRIDO") != NULL);
    static const double hib_tol = (getenv("DFB_SX_HIB_TOL") ?
                                   atof(getenv("DFB_SX_HIB_TOL")) : 0.01);
    if (hibrido) {
        const Interval antes_h = (k < x_new.size()) ? x_new[k] : Interval::all_reals();
        contract_sparse(x_new);
        if (x_new.is_empty()) return;
        const Interval desp_h = (k < x_new.size()) ? x_new[k] : Interval::all_reals();
        /* Si la regla barata ya contrajo bien, no se paga el simplex. */
        if (antes_h.ratiodelta(desp_h) >= hib_tol) { ++n_hib_barato; return; }
        /* Si el simplex ya se probo en esta caja y no dio nada, no se repite:
         * es el grupo C, donde no hay nada que sacar. */
        if (sx_skip > 0) { --sx_skip; ++n_hib_memo; return; }
        ++n_hib_simplex;
        /* si no, se sigue con el simplex sobre la caja ya contraida */
    }

    n_applied++;

    const int next = x_ref.size();
    if (work.size() != next) work.resize(next);
    for (int i = 0; i < nb_var && i < x_new.size(); ++i) work[i] = x_new[i];
    for (int i = nb_var; i < next; ++i) work[i] = x_ref[i];

    if (state == INITIAL) { perc_imprs.clear(); virtual_x = work[k]; }
    const Interval antes = work[k];

    double t0 = dfb_now();
    const int cap = 4 * refA_ref->nb_rows() + 100;
    /* Arranque en FRIO por omision, y esto es contraintuitivo pero esta medido:
     * reusar la base anterior con el objetivo nuevo cuesta MAS pivotes (8 864
     * contra 1 925 en DiscreteBoundary-0040, 286 contra 156 en Katsura-12).
     * El motivo es que reubicar las no basicas por el signo del costo reducido
     * nuevo restaura la factibilidad DUAL pero destruye la PRIMAL, y la base
     * optima para otro objetivo queda lejos de ser primal-factible; la base de
     * las b, en cambio, tiene una infactibilidad chica y natural.
     *
     * El warm start que si sirve —el que le da a PolytopeHull sus 1.08
     * iteraciones por cota— conserva la base Y la ubicacion de las no basicas,
     * o sea que se queda primal-factible, y corre un simplex PRIMAL para el
     * objetivo nuevo. Eso pide implementar el primal, que es el paso siguiente.
     * Con DFB_SX_WARM=1 se puede volver a probar el dual tibio. */
    static const bool sx_warm = (getenv("DFB_SX_WARM") != NULL);
    /* DFB_SX_WARM_MIN=1: tibio solo para las cotas inferiores. Sirve para
     * aislar el defecto del camino tibio, que sobre la MISMA caja devuelve para
     * `maximizar` la respuesta de `minimizar` (Brown-15: 16 en vez de
     * 6.67e6), con base dual-factible y primal-factible verificadas. */
    static const bool sx_warm_min = (getenv("DFB_SX_WARM_MIN") != NULL);
    /* Con un simplex POR CONTRACTOR el tibio es lo natural: la base heredada
     * es la de la llamada anterior de ESTA misma cota, con el mismo objetivo y
     * la caja apenas cambiada. */
    static const bool sx_propio = (getenv("DFB_SX_PROPIO") != NULL);
    const bool tibio = sx_warm || sx_propio || (sx_warm_min && !upper_contract);
    DFBSimplex::Status st = sx->solve(k, upper_contract, work, cap, tibio);
    iters = sx->iterations();
    t_pivot += dfb_now() - t0;

    /* Si el simplex no llega al optimo se cae a la regla vieja, de modo que el
     * modo no pueda ser peor que el actual. INFEASIBLE viene de un simplex
     * flotante, asi que NO prueba vacio: tampoco se usa para vaciar (eso
     * pediria un certificado de Farkas, que es otro asunto). */
    if (st != DFBSimplex::OPTIMAL) {
        n_sx_fallback++;
        contract_sparse(x_new);
        return;
    }

    IntervalVector gamma(1);
    state = FINAL;
    const Interval antes_sx = (k < x_new.size()) ? x_new[k] : Interval::all_reals();
    static const bool sx_trace = (getenv("DFB_SX_TRACE") != NULL);
    const bool cert_ok = certify_from_lambda(sx->lambda(), gamma);
    if (sx_trace) {
        static int n = 0;
        if (n++ < 400)
            fprintf(stderr, "[sxt] k=%d up=%d piv=%d cota_float=%.10g gamma_k=%s cert=%d "
                   "work_k=[%.6g,%.6g]\n", k, (int)upper_contract, iters,
                   sx->bound(),
                   gamma.size() > k ? std::to_string(gamma[k].mid()).c_str() : "?",
                   (int)cert_ok, work[k].lb(), work[k].ub());
    }
    if (cert_ok) {
        double t1 = dfb_now();
        Interval nb = gaussSeidel(work, k, gamma);
        t_bound += dfb_now() - t1;
        if (work.is_empty()) { x_new.set_empty(); return; }
        virtual_x = Interval(nb.lb(), work[k].ub());

        /* Reencolarse mientras siga rindiendo. El optimo del LP es el mejor
         * posible **para esta caja**, pero la propagacion sigue apretando otras
         * variables, asi que volver a resolver sobre la caja nueva puede dar una
         * cota mejor. Sin esto el simplex hace una sola pasada donde la regla
         * vieja hace varias, y en las Brown-* eso le costaba mas de lo que
         * ganaba por llegar al optimo. */
        /* Umbral de reencolado. Con el simplex tibio cada resolucion cuesta
         * unos pocos pivotes, asi que conviene reencolarse ante mejoras mas
         * chicas: lo que pierde el tibio no es el optimo del LP —el grupo B en
         * 1 lo confirma— sino la contraccion EXTRA de pasar repetidas veces
         * sobre la caja ya apretada. */
        static const double req_tol = (getenv("DFB_SX_REQUEUE") ?
                                       atof(getenv("DFB_SX_REQUEUE")) : 0.01);
        if (!work.is_empty() && antes.ratiodelta(work[k]) >= req_tol)
            state = CONTRACTING;
    }
    write_back(x_new);
    /* Si el simplex no dio nada, se saltea en los proximos nodos: es la firma
     * del grupo C, donde no hay nada que sacar. Se reintenta despues, por si la
     * caja cambio lo suficiente. */
    /* APAGADA por omision (0). Medida y es claramente peor: saltear el simplex
     * en los nodos siguientes destruye la poda —`alkyl` pasa de 82 celdas y
     * 1.44 s a 1 712 celdas y timeout—. En retrospectiva es razonable: «no dio
     * nada en este nodo» no predice el siguiente, porque la biseccion parte una
     * variable al medio y lo que no tenia nada que sacar puede tenerlo de
     * golpe. */
    static const int memo_n = (getenv("DFB_SX_MEMO_N") ?
                               atoi(getenv("DFB_SX_MEMO_N")) : 0);
    if (!x_new.is_empty() && k < x_new.size() &&
        antes_sx.ratiodelta(x_new[k]) < 1e-9)
        sx_skip = memo_n;
}

void CtcDFB::contract_sparse(IntervalVector& x_new){
    if (!init_ok) { iters = 0; state = FINAL; return; }
    n_applied++;

    const int next = x_ref.size();
    if (work.size() != next) work.resize(next);
    for (int i = 0; i < nb_var && i < x_new.size(); ++i) work[i] = x_new[i];
    for (int i = nb_var; i < next; ++i) work[i] = x_ref[i];
    const bool flip_caja = box_sign_flip();
    if (flip_caja) work[k] = -work[k];

    if (state == INITIAL) { perc_imprs.clear(); virtual_x = work[k]; }
    state = CONTRACTING;

    IntervalVector gamma(1);
    iters = 0;

    /* Punto de partida del seguimiento del mejor lambda: la fila 0 con la que
     * se entra ya da una cota, y hay que poder volver a ella. */
    best_lam = St.lam;
    best_bound = bound_proxy();
    have_best = true;

    while (max_iters == -1 || iters < max_iters) {
        int j; double delta, dir;
        double t0 = dfb_now();
        if (dual_pricing) largest_impact_dual(j, delta, dir);
        else              largest_impact_s(j, delta, dir);
        t_pricing += dfb_now() - t0;

        if (check_ready()) {
            int jd; double dd, dird;
            largest_impact_f(jd, dd, dird);
            if (jd != j || dird != dir || dd != delta) {
                ++n_check_difs;
                if (n_check_difs < 20)
                    printf("[chk] pricing difiere: disperso j=%d dir=%g delta=%.17g "
                           "| denso j=%d dir=%g delta=%.17g\n", j, dir, delta, jd, dird, dd);
            }
        }
        if (j == -1) { state = FINAL; break; }

        double alpha;
        t0 = dfb_now();
        const int i = ratio_test_s(j, dir, alpha);
        t_ratio_test += dfb_now() - t0;

        if (check_ready()) {
            double ad; const int id = ratio_test_f(j, dir, ad);
            if (id != i || (id != -1 && ad != alpha)) {
                ++n_check_difs;
                if (n_check_difs < 20)
                    printf("[chk] ratio test difiere: disperso i=%d alpha=%.17g "
                           "| denso i=%d alpha=%.17g\n", i, alpha, id, ad);
            }
        }
        if (i == -1) { state = FINAL; break; }

        t0 = dfb_now();
        St.add_to_row0(alpha, j);
        St.g[i] = 0.0;

        if (check_ready()) {
            /* Paso 1 del pivote, comparado por separado para poder atribuir
             * cualquier divergencia a add_to_row0 o a la eliminacion. */
            chk_last_j = j; chk_last_i = i;
            const int nc = Af.nb_cols();
            for (int c = 0; c < nc; ++c) Af[0][c] += alpha * Af[j][c];
            Af[0][i] = 0.0;
            compare_tableaus("add_to_row0");
        }

        if (!St.make_identity(i, j, false, DFB_F_TOL)) { state = FINAL; break; }
        St.g[k] = gamma_k_target;
        set_basic(j, i);
        t_pivot += dfb_now() - t0;

        const double b = bound_proxy();
        if (b > best_bound) { best_bound = b; best_lam = St.lam; }

        if (check_ready()) {
            make_column_identity_f(i, j, false);
            Af[0][k] = gamma_k_target;
            compare_tableaus("eliminacion");
        }

        if (dfb_trace)
            printf("T k=%d up=%d j=%d i=%d alpha=%.17g\n", k, (int)upper_contract, j, i, alpha);
        identity_rows[j] = i;
        ++iters;
    }

    if (certify_gamma_s(gamma, have_best ? &best_lam : NULL)) {
        double t0 = dfb_now();
        Interval nb = gaussSeidel(work, k, gamma);
        t_bound += dfb_now() - t0;
        if (work.is_empty()) { x_new.set_empty(); state = FINAL; return; }
        virtual_x = Interval(nb.lb(), work[k].ub());
    }

    if (flip_caja && !work.is_empty()) work[k] = -work[k];
    write_back(x_new);
}

void CtcDFB::flip_side(){
    if (sparse_tableau) {
        St.negate_row0();
        gamma_k_target = -gamma_k_target;
        St.g[k] = gamma_k_target;
        upper_contract = !upper_contract;
        lam_nz.clear();
        cand_rows.clear(); cand_age = 0;
        return;
    }
    const int nc = Af.nb_cols();
    for (int c = 0; c < nc; ++c) Af[0][c] = -Af[0][c];
    gamma_k_target = -gamma_k_target;
    Af[0][k] = gamma_k_target;
    upper_contract = !upper_contract;
    lam_nz.clear();
    cand_rows.clear(); cand_age = 0;
}

bool CtcDFB::renormalize_to(int k_new){
    double t0 = dfb_now();
    k = k_new;
    identity_rows.clear();
    cand_rows.clear(); cand_age = 0;
    lam_nz.clear();
    gamma_k_target = 1.0;
    upper_contract = false;
    const bool ok = make_column_identity_f(k_new, 0, true);
    init_ok = ok;
    t_init += dfb_now() - t0;
    return ok;
}

void CtcDFB::init_from_lower_float(const CtcDFB& lower, IntervalVector& x_ref){
    if (sparse_tableau) { init_from_lower_sparse(lower, x_ref); return; }
    double t_ini_init = dfb_now();
    n_inits++;
    refA_ref = lower.refA_ref;

    const int na = (*refA_ref).nb_cols();
    Af.resize(lower.Af.nb_rows(), lower.Af.nb_cols());
    Af = lower.Af;

    this->x_ref.resize(x_ref.size());
    this->x_ref = IntervalVector(x_ref);

    lam_nz.clear();
    identity_rows.clear();
    cand_rows.clear(); cand_age = 0;
    state = INITIAL;
    init_ok = lower.init_ok;

    if (init_ok) {
        /* Fila 0 completa negada (parte A y parte lambda), salvo la entrada k
         * de la parte A, que vale 1 en los dos casos. */
        const int nc = Af.nb_cols();
        for (int c = 0; c < nc; ++c) Af[0][c] = -Af[0][c];
        Af[0][k] = 1.0;
        (void)na;
    }

    t_init += dfb_now() - t_ini_init;
}

/* Hace identidad la columna col en la fila row, sobre TODO el tableau
 * aumentado (para que lambda siga siendo veraz). */
bool CtcDFB::make_column_identity_f(int col, int row, bool interchange){
    const int m = Af.nb_rows();
    const int nc = Af.nb_cols();

    if (interchange && std::fabs(Af[row][col]) < DFB_F_TOL) {
        int found = -1;
        for (int i = 0; i < m; ++i)
            if (i != row && std::fabs(Af[i][col]) >= DFB_F_TOL) { found = i; break; }
        if (found < 0) return false;
        for (int c = 0; c < nc; ++c) std::swap(Af[row][c], Af[found][c]);
    }
    const double piv = Af[row][col];
    if (std::fabs(piv) < DFB_F_TOL) return false;

    const double inv = 1.0 / piv;
    for (int c = 0; c < nc; ++c) Af[row][c] *= inv;
    Af[row][col] = 1.0;                       /* solo la parte A se fuerza */

    ++n_elim_calls;
    for (int i = 0; i < m; ++i) {
        if (i == row) continue;
        const double f = Af[i][col];
        if (f == 0.0) continue;
        ++n_elim_rows;
        for (int c = 0; c < nc; ++c) Af[i][c] -= f * Af[row][c];
        Af[i][col] = 0.0;
    }
    return true;
}

/* Pricing: misma regla que la version de intervalos, en doubles. Devuelve la
 * fila entrante o -1. */
int CtcDFB::largest_impact_f(int& j_out, double& delta_out, double& dir_out){
    const int m  = Af.nb_rows();
    const int na = (*refA_ref).nb_cols();

    double best_incr = 0.0, best_decr = 0.0;
    int    j_incr = 0, j_decr = 0;

    /* Diagnostico de dispersion (DFB_PRICE_STATS=1). Apagado no cuesta una
     * rama por pivote; encendido recorre na columnas de mas, asi que no hay que
     * medir tiempos con esto puesto. */
    if (price_stats) {
        for (int i = 0; i < na; ++i) {
            ++n_price_gtotal;
            if (std::fabs(Af[0][i]) >= 1e-5) ++n_price_gnz;
        }
    }

    for (int j = 1; j < m; ++j) {
        double acc_incr = 0.0, acc_decr = 0.0, nrm = 0.0;
        for (int i = 0; i < na; ++i) {
            const double lo = work[i].lb(), hi = work[i].ub();
            const double g = Af[0][i], a = Af[j][i];
            if (price_stats && a != 0.0) ++n_price_anz;
            double xi_incr, xi_decr;
            if (g >= 1e-5)       { xi_incr = hi; xi_decr = hi; }
            else if (g <= -1e-5) { xi_incr = lo; xi_decr = lo; }
            else if (a >= 1e-5)  { xi_incr = hi; xi_decr = lo; }
            else if (a <= -1e-5) { xi_incr = lo; xi_decr = hi; }
            else                 { xi_incr = 0.0; xi_decr = 0.0; }
            acc_incr += a * xi_incr;
            acc_decr += a * xi_decr;
            if (pricing_norm) nrm += a * a;   /* misma pasada: sale gratis */
        }
        double gi = -acc_incr, gd = acc_decr;
        if (pricing_norm && nrm > 1e-300) {
            /* M2(b): delta/||A_j|| aproxima la mejora real delta*alpha. Solo
             * cambia el ORDEN de eleccion, no el criterio de parada, que sigue
             * comparandose sobre la mejora sin normalizar. */
            const double inv = 1.0 / std::sqrt(nrm);
            gi *= inv; gd *= inv;
        }
        if (gi > best_incr)      { best_incr = gi;      j_incr = j; }
        if (gd > best_decr)      { best_decr = gd;      j_decr = j; }
    }
    if (best_incr > best_decr) { delta_out = best_incr; j_out = j_incr; dir_out =  1.0; }
    else                       { delta_out = best_decr; j_out = j_decr; dir_out = -1.0; }

    if (pricing_norm && j_out > 0) {
        /* Con normalizacion, delta_out esta escalado y no sirve como criterio
         * de parada: se recalcula la mejora sin normalizar de la fila elegida. */
        double acc = 0.0;
        for (int i = 0; i < na; ++i) {
            const double lo = work[i].lb(), hi = work[i].ub();
            const double g = Af[0][i], a = Af[j_out][i];
            double xi;
            if (g >= 1e-5)       xi = (dir_out > 0) ? hi : hi;
            else if (g <= -1e-5) xi = lo;
            else if (a >= 1e-5)  xi = (dir_out > 0) ? hi : lo;
            else if (a <= -1e-5) xi = (dir_out > 0) ? lo : hi;
            else                 xi = 0.0;
            acc += a * xi;
        }
        delta_out = (dir_out > 0) ? -acc : acc;
    }
    if (delta_out <= stop_tol) j_out = -1;
    return j_out;
}

/* Ratio test en doubles, con desempate por magnitud de pivote (Harris). */
int CtcDFB::ratio_test_f(int j, double dir, double& alpha_out){
    const int na = (*refA_ref).nb_cols();
    double min_mag = POS_INFINITY, best_piv = 0.0;
    int    min_i = -1;

    for (int i = 0; i < na; ++i) {
        const double a = Af[j][i];
        if (std::fabs(a) < DFB_F_TOL) continue;
        const double alpha = (Af[0][i] / a) * dir;
        if (alpha >= 0.0) continue;
        const double mag = std::fabs(alpha), piv = std::fabs(a);
        if (min_i == -1 || mag < min_mag * (1.0 - harris_rel_band)) {
            min_mag = mag; best_piv = piv; min_i = i;
        } else if (harris_tie_break && mag <= min_mag * (1.0 + harris_rel_band)
                   && piv > best_piv) {
            best_piv = piv; min_i = i;               /* no se excede el minimo */
        }
    }
    if (min_i == -1) return -1;
    alpha_out = -(Af[0][min_i] / Af[j][min_i]);
    return min_i;
}

/* gamma = lambda^T . refA, en intervalos, con la columna k negada si se
 * contrae la cota superior (equivale a usar refA_k). */
bool CtcDFB::certify_gamma(IntervalVector& gamma_out){
    double t0 = dfb_now();
    n_certifications++;

    const IntervalMatrix& R = *refA_ref;
    const int m = R.nb_rows(), na = R.nb_cols();

    lam_nz.clear();
    for (int j = 0; j < m; ++j)
        if (Af[0][na + j] != 0.0) lam_nz.push_back(j);

    gamma_out.resize(na);
    for (int i = 0; i < na; ++i) gamma_out[i] = Interval(0.0);
    for (size_t t = 0; t < lam_nz.size(); ++t) {
        const int j = lam_nz[t];
        const Interval lam(Af[0][na + j]);
        for (int i = 0; i < na; ++i) gamma_out[i] += lam * R[j][i];
    }
    /* Con gamma_k_target = -1 el signo ya esta en la fila 0 (y por lo tanto en
     * lambda), asi que no hay que negar la entrada k: eso es solo para el
     * camino de la cola, donde el signo esta en la columna k de refA_k. */
    if (upper_contract && gamma_k_target > 0) gamma_out[k] = -gamma_out[k];

    t_certify += dfb_now() - t0;

    /* Sin un gamma_k acotado lejos de 0 la cota no sirve (dividir por un
     * intervalo que contiene el 0 no acota nada). No es un problema de
     * solidez: simplemente no se usa. */
    return !gamma_out[k].contains(0.0);
}

void CtcDFB::contract_float(IntervalVector& x_new){
    if (simplex_mode)   { contract_simplex(x_new); return; }
    if (sparse_tableau) { contract_sparse(x_new); return; }
    if (!init_ok) { iters = 0; state = FINAL; return; }
    n_applied++;

    const int next = x_ref.size();
    if (work.size() != next) work.resize(next);
    for (int i = 0; i < nb_var && i < x_new.size(); ++i) work[i] = x_new[i];
    for (int i = nb_var; i < next; ++i) work[i] = x_ref[i];
    const bool flip_caja = box_sign_flip();
    if (flip_caja) work[k] = -work[k];

    if (state == INITIAL) { perc_imprs.clear(); virtual_x = work[k]; }
    state = CONTRACTING;

    IntervalVector gamma(1);
    iters = 0;

    while (max_iters == -1 || iters < max_iters) {
        int j; double delta, dir;
        double t0 = dfb_now();
        largest_impact_f(j, delta, dir);
        t_pricing += dfb_now() - t0;
        if (j == -1) { state = FINAL; break; }

        double alpha;
        t0 = dfb_now();
        const int i = ratio_test_f(j, dir, alpha);
        t_ratio_test += dfb_now() - t0;
        if (i == -1) {
            /* Sin fila de bloqueo no se puede concluir infactibilidad: el
             * pivoteo es flotante y no constituye prueba. Se detiene. */
            state = FINAL; break;
        }

        t0 = dfb_now();
        const int nc = Af.nb_cols();
        for (int c = 0; c < nc; ++c) Af[0][c] += alpha * Af[j][c];
        Af[0][i] = 0.0;
        if (!make_column_identity_f(i, j, false)) { state = FINAL; break; }
        Af[0][k] = gamma_k_target;
        t_pivot += dfb_now() - t0;

        if (dfb_trace)
            printf("T k=%d up=%d j=%d i=%d alpha=%.17g\n", k, (int)upper_contract, j, i, alpha);
        identity_rows[j] = i;
        ++iters;
    }

    /* Certificacion: una sola vez, con el lambda acumulado. */
    if (certify_gamma(gamma)) {
        double t0 = dfb_now();
        Interval nb = gaussSeidel(work, k, gamma);
        t_bound += dfb_now() - t0;
        if (work.is_empty()) { x_new.set_empty(); state = FINAL; return; }
        virtual_x = Interval(nb.lb(), work[k].ub());
    }

    if (flip_caja && !work.is_empty()) work[k] = -work[k];
    write_back(x_new);
}

void CtcDFB::write_back(IntervalVector& x_new) const {
    if (work.is_empty()) { x_new.set_empty(); return; }
    const int n = (nb_var < x_new.size()) ? nb_var : x_new.size();
    for (int i = 0; i < n; ++i) x_new[i] = work[i];
}

void CtcDFB::contract(IntervalVector& x_new) {
    if (light_mode)     { contract_light(x_new); return; }
    if (float_pivoting) { contract_float(x_new); return; }

    /* Si init()/regenerateA() no pudo construir la base, este contractor no
     * puede hacer nada; la caja se devuelve intacta. */
    if (!init_ok) { iters = 0; state = FINAL; return; }
    n_applied++;

    /* Se trabaja sobre un buffer propio con la dimension extendida (x U b).
     * NO se redimensiona la caja del llamador: ver el comentario de `work`. */
    const int next = x_ref.size();
    if (work.size() != next) work.resize(next);
    for (int i = 0; i < nb_var && i < x_new.size(); ++i) work[i] = x_new[i];
    for (int i = nb_var; i < next; ++i) work[i] = x_ref[i];

    if(state==INITIAL){
        perc_imprs.clear();
        if (upper_contract)  x_ref[k] = -x_ref[k];
        virtual_x = gaussSeidel(x_ref, k, A[0]);
        /* Se comprueba antes de leer x_ref[k]: si quedo vacia, ese componente
         * ya no tiene cotas utilizables. x_ref contiene a la caja actual, asi
         * que su vacuidad la implica. */
        if (x_ref.is_empty()) {
            x_new.set_empty();
            iters = 0; state = FINAL; return;
        }
        virtual_x = Interval(virtual_x.lb(), x_ref[k].ub());
        if (upper_contract) x_ref[k] = -x_ref[k];
    }
    
    int i;
    int j;
    Interval alpha;
    Interval delta, direction;
    state = CONTRACTING;

    if (upper_contract) work[k] = -work[k];

    iters = 0;

    while (max_iters == -1 || iters < max_iters) {
        double t0 = dfb_now();
        tie(j, delta, direction) = largestImpact(A, work, A[0]);
        t_pricing += dfb_now() - t0;
       
        if (j == -1 || delta == Interval(0)) {
            if (upper_contract) work[k] = -work[k];
            state = FINAL;
            write_back(x_new);
            return;
        } 

        t0 = dfb_now();
        tie(alpha, i) = calculateAlpha(A[j], A[0], direction);
        t_ratio_test += dfb_now() - t0;

        if (i == -1) {
            if (upper_contract) work[k] = -work[k];

            if (ratio_test_status == NO_CANDIDATE) {
                /* Ninguna fila puede bloquear el paso dual: la relajacion
                 * lineal es infactible y la caja no contiene soluciones. */
                x_new.set_empty();
            } else {
                /* INCONCLUSIVE: no se puede concluir nada. Se detiene el
                 * pivoteo conservando la caja. Declararla vacia aqui era el
                 * origen de las pruebas de infactibilidad falsas. */
                state = FINAL;
                write_back(x_new);
            }
            return;
        }

        t0 = dfb_now();
        A[0] = A[0] + alpha * A[j];
        
        A[0][i] = Interval(0);

        makeColumnIdentity(A, i, false, j);
        A[0][k] = Interval(1);
        t_pivot += dfb_now() - t0;

        double x_lb = virtual_x.lb();

        double old_size = virtual_x.diam();
        t0 = dfb_now();
        virtual_x = gaussSeidel(work, k, A[0]);
        t_bound += dfb_now() - t0;
        /* La cota deducida de gamma.x = 0 puede no intersectar la caja, lo que
         * prueba que no hay solucion. Hay que salir aqui: seguir pivoteando
         * sobre una caja vacia leia cotas inexistentes. */
        if (work.is_empty()) {
            x_new.set_empty();
            state = FINAL; return;
        }
        virtual_x = Interval(virtual_x.lb(), work[k].ub());
        double impr = (virtual_x.lb() - x_lb)/old_size;
        perc_imprs.push_back(impr);

        /* Refresco adaptativo: si este pivote no mejoro la cota, la lista de
         * candidatos se quedo sin filas utiles y hay que rehacerla. */
        if (partial_pricing && partial_pricing_adaptive &&
            (!std::isfinite(impr) || impr <= 1e-12)) {
            cand_rows.clear(); cand_age = 0;
        }

        identity_rows[j] = i;
        if (contract_all) {
            for (const auto& r : identity_rows) {
                gaussSeidel(work, r.second, A[r.first]);
                if (work.is_empty()) {
                    x_new.set_empty();
                    state = FINAL; return;
                }
            }
        }
        ++iters;
    }

    if (upper_contract && !work.is_empty()) work[k] = -work[k];

    write_back(x_new);
}

std::pair<IntervalVector, IntervalVector> CtcDFB::calculateImpacts(
    const IntervalMatrix& A, const IntervalVector& x_new, const IntervalVector& gamma) {
   
            //cout << "[CtcDFB] Calculating impacts. A dimensions: " << A.nb_rows() << "x" << A.nb_cols()
            //     << ", x_new size: " << x_new.size() << ", gamma size: " << gamma.size() << endl;
    
  //  cout << "x_new" << x_new << endl;
    assert(A.nb_cols() == x_new.size() && "Matrix column count must match x_new size");
            
    int m = A.nb_rows();
    int n = A.nb_cols();

    IntervalVector grad_incr = IntervalVector(m, Interval(0));
    IntervalVector grad_decr = IntervalVector(m, Interval(0));

    for (int j = 1; j < m; ++j) {
        IntervalVector x_prime_incr = IntervalVector(n, Interval(0));
        IntervalVector x_prime_decr = IntervalVector(n, Interval(0));

        for (int i = 0; i < n; ++i) {
            //Interval signOfGammai = sign(gamma[i]);
       //     cout << "gamma[" << i << "] = " << gamma[i] << endl;
            if (gamma[i].lb() >= 0.00001) {
                x_prime_incr[i] = Interval(x_new[i].ub());
                x_prime_decr[i] = Interval(x_new[i].ub());
            } else if (gamma[i].ub() <= -0.00001) {
                x_prime_incr[i] = Interval(x_new[i].lb());
                x_prime_decr[i] = Interval(x_new[i].lb());
            }
            else{
                //Interval signOfAji = sign(A[j][i]);
       //         cout << "A[" << j << "][" << i << "] = " << A[j][i] << endl;
                if (A[j][i].lb() >= 0.00001) {
                    x_prime_incr[i] = Interval(x_new[i].ub());
                    x_prime_decr[i] = Interval(x_new[i].lb());
     //               cout << "case 1" << endl;
                } else if (A[j][i].ub() <= -0.00001) {
                    x_prime_incr[i] = Interval(x_new[i].lb());
                    x_prime_decr[i] = Interval(x_new[i].ub());
       //             cout << "case 2" << endl;
                }
         //       cout << "x_prime_incr[" << i << "] = " << x_prime_incr[i] << endl;
           //     cout << "x_prime_decr[" << i << "] = " << x_prime_decr[i] << endl;
            }
        }

        grad_incr[j] = Interval(A[j] * x_prime_incr);
        grad_incr[j] = - grad_incr[j];
        grad_decr[j] = Interval(A[j] * x_prime_decr);
    }

    return {grad_incr, grad_decr};
}

void CtcDFB::prepare_columns(const IntervalVector& x_new, const IntervalVector& gamma, int n) {
    if ((int)col_case.size() != n) {
        col_case.resize(n); col_lb.resize(n); col_ub.resize(n);
    }
    for (int i = 0; i < n; ++i) {
        col_ub[i] = Interval(x_new[i].ub());
        col_lb[i] = Interval(x_new[i].lb());
        col_case[i] = (gamma[i].lb() >= 0.00001) ? 0
                    : (gamma[i].ub() <= -0.00001) ? 1 : 2;
    }
}

/* Impacto de una fila, acumulando en el mismo orden de indices que mulVV
 * (y += v1[i]*v2[i] desde y=0), de modo que el resultado es identico al de la
 * version original con IntervalVector intermedios. */
void CtcDFB::row_impact(const IntervalVector& Aj, int n,
                        Interval& g_incr, Interval& g_decr) const {
    Interval acc_incr(0.0), acc_decr(0.0);
    const Interval zero(0.0);
    for (int i = 0; i < n; ++i) {
        const Interval& a = Aj[i];
        switch (col_case[i]) {
        case 0:
            acc_incr += a * col_ub[i]; acc_decr += a * col_ub[i]; break;
        case 1:
            acc_incr += a * col_lb[i]; acc_decr += a * col_lb[i]; break;
        default:
            if (a.lb() >= 0.00001) {
                acc_incr += a * col_ub[i]; acc_decr += a * col_lb[i];
            } else if (a.ub() <= -0.00001) {
                acc_incr += a * col_lb[i]; acc_decr += a * col_ub[i];
            } else {
                acc_incr += a * zero;     acc_decr += a * zero;
            }
        }
    }
    g_incr = -acc_incr;
    g_decr = acc_decr;
}

std::tuple<int, Interval, Interval> CtcDFB::largestImpact(
    const IntervalMatrix& A, const IntervalVector& x_new, const IntervalVector& gamma) {

    if (pricing_fused || partial_pricing) {
        const int m = A.nb_rows();
        const int n = A.nb_cols();
        prepare_columns(x_new, gamma, n);

        /* Equivale a getMaxValue sobre grad, cuyo elemento 0 es Interval(0). */
        Interval best_incr(0.0), best_decr(0.0);
        int j_incr = 0, j_decr = 0;
        Interval gi, gd;

        int L = partial_pricing_size > 0 ? partial_pricing_size
                                         : (int)std::ceil(std::sqrt((double)m));
        int R = partial_pricing_refresh > 0 ? partial_pricing_refresh : L;

        bool use_list = partial_pricing && !cand_rows.empty() && cand_age < R;
        bool did_full = false;

        for (int attempt = 0; attempt < 2; ++attempt) {
            best_incr = Interval(0.0); best_decr = Interval(0.0);
            j_incr = 0; j_decr = 0;

            if (use_list) {
                for (size_t c = 0; c < cand_rows.size(); ++c) {
                    int j = cand_rows[c];
                    if (j < 1 || j >= m) continue;
                    row_impact(A[j], n, gi, gd);
                    if (gi.lb() > best_incr.lb()) { best_incr = gi; j_incr = j; }
                    if (gd.lb() > best_decr.lb()) { best_decr = gd; j_decr = j; }
                }
                cand_age++;
            } else {
                /* Pasada completa: elige la mejor fila y, si hay partial
                 * pricing, rearma la lista con las L mejores. */
                did_full = true;
                std::vector< std::pair<double,int> > score;
                if (partial_pricing) score.reserve(m);
                for (int j = 1; j < m; ++j) {
                    row_impact(A[j], n, gi, gd);
                    if (gi.lb() > best_incr.lb()) { best_incr = gi; j_incr = j; }
                    if (gd.lb() > best_decr.lb()) { best_decr = gd; j_decr = j; }
                    if (partial_pricing)
                        score.push_back(std::make_pair(std::max(gi.lb(), gd.lb()), j));
                }
                if (partial_pricing) {
                    int keep = std::min((size_t)L, score.size());
                    std::partial_sort(score.begin(), score.begin() + keep, score.end(),
                        [](const std::pair<double,int>& a, const std::pair<double,int>& b) {
                            return a.first > b.first;
                        });
                    cand_rows.clear();
                    for (int c = 0; c < keep; ++c) cand_rows.push_back(score[c].second);
                    cand_age = 0;
                }
            }

            Interval delta; int j; Interval direction;
            if (best_incr.lb() > best_decr.lb()) {
                delta = best_incr; j = j_incr; direction = Interval(1);
            } else {
                delta = best_decr; j = j_decr; direction = Interval(-1);
            }
            if (delta.ub() <= 0.000001) j = -1;

            /* Si la lista no ofrecio ninguna fila que mejore, se rehace la
             * pasada completa antes de concluir que no hay ninguna. */
            if (j == -1 && !did_full) { use_list = false; continue; }
            return {j, delta, direction};
        }
    }

    IntervalVector grad_incr, grad_decr;
    tie(grad_incr, grad_decr) = calculateImpacts(A, x_new, gamma);
    //cout << "grad_incr: " << grad_incr << endl;
    //cout << "grad_decr: " << grad_decr << endl;

    Interval delta_incr, delta_decr;
    int j_incr, j_decr;
    tie(delta_incr, j_incr) = getMaxValue(grad_incr);
    tie(delta_decr, j_decr) = getMaxValue(grad_decr);

    Interval delta = Interval(0);
    int j;
    Interval direction;

  //  cout << "delta_incr: " << delta_incr << " (j=" << j_incr << "), delta_decr: " 
    if (delta_incr.lb() > delta_decr.lb()) {
        delta = delta_incr;
        j = j_incr;
        direction = Interval(1);
    } else {
        delta = delta_decr;
        j = j_decr;
        direction = Interval(-1);
    }

    //Interval signOfDelta = sign(delta);
    
    
    if (delta.ub() <= 0.000001) {
        j = -1;
    }

    return {j, delta, direction};
}

/* Ratio test.
 *
 * Elige la fila de bloqueo i: aquella cuyo cociente gamma_i/A[j][i] limita el
 * paso alpha. Devuelve i==-1 cuando no hay fila elegible, pero deja en
 * ratio_test_status *por que*, porque las consecuencias son opuestas:
 *
 *   NO_CANDIDATE  ninguna fila puede bloquear el paso. Toda fila o no afecta a
 *                 gamma_i (A[j][i] es exactamente 0) o lo mueve en el sentido
 *                 que no bloquea (cociente estrictamente positivo). Recien
 *                 entonces el paso dual es no acotado y la relajacion lineal
 *                 es infactible.
 *
 *   INCONCLUSIVE  hay filas cuyo signo no se puede determinar en aritmetica de
 *                 intervalos: A[j][i] contiene el 0 (y entonces el cociente no
 *                 tiene signo ni magnitud acotada) o el cociente mismo contiene
 *                 el 0. No se puede descartar que alguna bloquee, asi que no se
 *                 puede concluir nada.
 *
 * La version anterior no distinguia estos casos y ademas descartaba candidatos
 * legitimos con una tolerancia ABSOLUTA (alpha.ub() <= -1e-5) sobre una
 * cantidad cuya escala depende del problema; el resultado era declarar vacias
 * cajas que contenian soluciones. Aqui la condicion de candidatura es el signo
 * riguroso del cociente (alpha.ub() < 0), sin umbral de magnitud: un paso
 * pequeno es legitimo, y si no aporta contraccion la propagacion ya no vuelve a
 * encolar al contractor.
 */
std::pair<Interval, int> CtcDFB::calculateAlpha(
    const IntervalVector& Aj, const IntervalVector& gamma, const Interval& direction) {
    Interval min_alpha(1e20);
    int min_index = -1;
    double min_mag = POS_INFINITY;
    bool inconclusive = false;
    int n = Aj.size();

    /* Candidatos, solo si hace falta el desempate de Harris (segunda pasada). */
    std::vector< std::pair<double, int> > cand;   // (|alpha.lb()|, i)
    std::vector<Interval> cand_alpha;

    for (int i = 0; i < n; ++i) {
        const Interval& a_ji = Aj[i];

        /* La fila j no puede modificar gamma_i: no es candidata y no aporta duda. */
        if (a_ji.lb() == 0.0 && a_ji.ub() == 0.0) continue;

        /* A[j][i] contiene el 0: el cociente no tiene signo determinado y su
         * magnitud no esta acotada. No se puede descartar que bloquee. */
        if (a_ji.lb() <= 0.0 && a_ji.ub() >= 0.0) { inconclusive = true; continue; }

        Interval alpha = (gamma[i] / a_ji) * direction;

        if (alpha.ub() < 0.0) {
            /* Candidata rigurosa: el cociente es estrictamente negativo. */
            double mag = std::abs(alpha.lb());
            if (harris_tie_break) {
                cand.push_back(std::make_pair(mag, i));
                cand_alpha.push_back(alpha);
            }
            if (mag < min_mag) { min_mag = mag; min_alpha = alpha; min_index = i; }
        } else if (alpha.lb() <= 0.0) {
            /* El cociente contiene el 0 (o es 0): podria bloquear con paso nulo. */
            inconclusive = true;
        }
        /* alpha estrictamente positiva: no bloquea y no aporta duda. */
    }

    /* Desempate tipo Harris: entre los candidatos con cociente dentro de la
     * banda relativa del minimo, se toma el de mayor |A[j][i]|. El minimo nunca
     * se excede, de modo que la eleccion sigue siendo valida. */
    if (harris_tie_break && min_index != -1) {
        double best_pivot = std::abs(Aj[min_index].mid());
        double bound = min_mag * (1.0 + harris_rel_band);
        for (size_t c = 0; c < cand.size(); ++c) {
            if (cand[c].first > bound) continue;
            double piv = std::abs(Aj[cand[c].second].mid());
            if (piv > best_pivot) {
                best_pivot = piv;
                min_index  = cand[c].second;
                min_alpha  = cand_alpha[c];
            }
        }
    }

    if (min_index != -1) {
        ratio_test_status = BLOCKING_ROW;
    } else if (inconclusive) {
        ratio_test_status = INCONCLUSIVE;
        n_inconclusive++;
    } else {
        ratio_test_status = NO_CANDIDATE;
        n_no_candidate++;
    }

    min_alpha = min_alpha * direction * Interval(-1.0);
    return {min_alpha, min_index};
}

void CtcDFB::changeSigns(IntervalMatrix& A, IntervalVector& x_new) {
    
    //A[0] = -A[0];
    //A[0][this->k] = -A[0][this->k];

    for (int i = 0; i < A.nb_rows(); ++i) 
        A[i][k] = -A[i][k];
    

    x_new[k] = -x_new[k];
}

                                // AUX FUNCTIONS
int CtcDFB::makeColumnIdentity(IntervalMatrix& A, const int k, bool interchange, 
                               int j){
    int m = A.nb_rows();
    int n = A.nb_cols();

    // Paso 1: Encontrar la fila adecuada para el pivoteo
    if (interchange && (A[j][k].lb() == 0 || A[j][k].ub() == 0)) {
        bool found = false;
        for (int jAux = 0; jAux < m; ++jAux) {
            if (jAux == j) continue; // Ya la estamos evaluando
            if (!(A[jAux][k].lb() == 0 || A[jAux][k].ub() == 0)) {
                // Intercambiar las filas j y jAux
                for (int i = 0; i < n; ++i) {
                    std::swap(A[j][i], A[jAux][i]);
                }
                //j = jAux; // Actualizar j con la fila válida
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("No hay ninguna fila con valor no nulo en la columna k");
        }
    }
    

    // Paso 2: Normalizar la fila j respecto del valor en posición k
    if (A[j][k].lb() == 0 || A[j][k].ub() == 0) {
        throw std::invalid_argument("The divider has a bound equal to 0. Which is not allowed.");
    }
    
    //if (interchange){
        A[j] = (Interval(1) / A[j][k]) * A[j];

        A[j][k] = Interval(1);
    //}

    // Paso 3: Hacer ceros los demás elementos en la columna k
    for (int jj = 0; jj < m; ++jj) {
        if (jj == j) continue;
        Interval factor = A[jj][k];///A[j][k];
        for (int i = 0; i < n; ++i) {
            A[jj][i] -= factor * A[j][i];
        }
        A[jj][k] = Interval(0);
    }

    return j;
}

// Aplica la fórmula de Gauss-Seidel sobre el vector x
Interval CtcDFB::gaussSeidel(IntervalVector& x, int k, IntervalVector& gamma){
    double epsilon = 1e-6;
    int n = gamma.size();
    

    if (k == -1 || k >= n){
        throw std::invalid_argument("Invalid k.");
    }

//    if (gamma[k] != Interval(1)){
//        cout << "gamma[k] = " << gamma[k] << endl;
//        throw std::invalid_argument("Gamma[k] is not 1.");
//    }

    Interval tmp = Interval(gamma[k]);
    gamma[k] = Interval(0);

    Interval xContract = -(IntervalVector(gamma) * IntervalVector(x));
    gamma[k] = Interval(tmp);
    //cout << -gamma << "*" << x << " = " << xContract << endl;
    if (gamma[k] != Interval(1))
        xContract = xContract / gamma[k];

    if (x[k].intersects(xContract))
    {
        if (contract_all){
    //        std::cout << "x[k] = " << x[k] << std::endl;
    //        std::cout << "xContract = " << xContract << std::endl;
        }
        x[k] = x[k] & xContract;
        if (contract_all){
    //        std::cout << "x[k] = " << x[k] << std::endl;
        }
    } else {
        /* gamma es una combinacion lineal de las filas de A, de modo que
         * gamma.x = 0 vale para toda solucion de la relajacion. Si la cota que
         * se deduce de esa ecuacion no intersecta a x[k], la relajacion no
         * tiene solucion en la caja: es una prueba de vacio rigurosa.
         *
         * Antes esta rama no hacia nada y la deteccion se perdia; el vaciado
         * ocurria en cambio por la via del ratio test, donde no era valido. */
        if (prove_empty_by_bound) {
            x[k].set_empty();
            n_empty_by_bound++;
        }
    }
    if (x[k].is_empty()) x.set_empty();

    return xContract;
}

// Encuentra el índice del valor máximo en un vector
std::pair<Interval, int> CtcDFB::getMaxValue(const IntervalVector& vector) {
    Interval max_value = Interval(vector[0]);
    int max_index = 0;
    int n = vector.size();

    for (int i = 1; i < n; ++i) {
        if (vector[i].lb() > max_value.lb()) {
            max_value = Interval(vector[i]);
            max_index = i;
        }
    }

    return {max_value, max_index};
}
