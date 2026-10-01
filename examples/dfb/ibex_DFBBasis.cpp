#include "ibex_DFBBasis.h"

#include <cmath>
#include <cstdlib>

using namespace soplex;

namespace ibex {

long DFBBasis::n_factorizations = 0;
long DFBBasis::n_updates       = 0;
long DFBBasis::n_btran         = 0;
long DFBBasis::n_ftran         = 0;

bool DFBBasis::refactor_always  = (getenv("DFB_REFACTOR_ALWAYS") != NULL);
int  DFBBasis::refactor_period  = (getenv("DFB_REFACTOR_PERIOD") ?
                                   atoi(getenv("DFB_REFACTOR_PERIOD")) : 50);

void DFBBasis::reset_counters() {
    n_factorizations = n_updates = n_btran = n_ftran = 0;
}

DFBBasis::DFBBasis() : m(0), na(0), obj_row(0), rhs(1), sol(1), eta(NULL),
                       updates_since_factor(0) {
}

DFBBasis::~DFBBasis() {
    delete eta;
}

void DFBBasis::load_matrix(const Matrix& Abar) {
    m  = Abar.nb_rows();
    na = Abar.nb_cols();

    row_nz.assign(m, std::vector<std::pair<int,double> >());
    col_nz.clear();  col_nz.reserve(na);
    unit_nz.clear(); unit_nz.reserve(m);

    for (int c = 0; c < na; ++c) col_nz.push_back(DSVector(m));
    for (int j = 0; j < m; ++j) {
        unit_nz.push_back(DSVector(1));
        unit_nz[j].add(j, 1.0);
    }

    for (int j = 0; j < m; ++j) {
        for (int c = 0; c < na; ++c) {
            const double v = Abar[j][c];
            if (v == 0.0) continue;
            row_nz[j].push_back(std::make_pair(c, v));
            col_nz[c].add(j, v);
        }
    }

    rhs.reDim(m); sol.reDim(m);
    delete eta; eta = new soplex::SSVector(m);
}


bool DFBBasis::init_basis(int k, double tol) {
    if (m <= 0 || k < 0 || k >= na) return false;

    /* Base de holguras: la fila j tiene basica la columna na+j, o sea B = I. */
    basic.assign(m, 0);
    where_basic.assign(na + m, -1);
    for (int j = 0; j < m; ++j) {
        basic[j] = na + j;
        where_basic[na + j] = j;
    }
    updates_since_factor = 0;
    if (!factorize()) return false;

    /* La columna k entra en la base. Va a la fila 0 salvo que su entrada ahi
     * sea demasiado chica, en cuyo caso se toma la primera fila usable: es el
     * mismo criterio que el `interchange` del camino denso, pero en vez de
     * permutar filas se mueve la fila objetivo, con lo que los indices de
     * lambda siguen refiriendose a las filas originales de refA. */
    obj_row = -1;
    for (int j = 0; j < m; ++j) {
        double v = 0.0;
        for (size_t t = 0; t < row_nz[j].size(); ++t)
            if (row_nz[j][t].first == k) { v = std::fabs(row_nz[j][t].second); break; }
        if (v >= tol) { obj_row = j; break; }
    }
    if (obj_row < 0) return false;

    return change_basis(obj_row, k);
}

bool DFBBasis::set_basis(const std::vector<int>& cols) {
    if ((int)cols.size() != m) return false;
    basic.assign(m, 0);
    where_basic.assign(na + m, -1);
    for (int j = 0; j < m; ++j) {
        const int c = cols[j];
        if (c < 0 || c >= na + m) return false;
        basic[j] = c;
        where_basic[c] = j;
    }
    updates_since_factor = 0;
    return factorize();
}

/* OJO: las excepciones de SoPlex (SPxException y derivadas) NO derivan de
 * std::exception, asi que se escapan de cualquier `catch (std::exception&)` y
 * llegan como `catch (...)` al llamador. Aca se atajan y se degrada a "la
 * factorizacion no sirve", que el simplex maneja cayendo a su camino de
 * respaldo. Sin esto, 6 de 153 instancias del banco terminaban en
 * `unknown_error`. */
bool DFBBasis::factorize() {
    std::vector<const SVector*> cols(m);
    for (int j = 0; j < m; ++j) cols[j] = &column(basic[j]);
    ++n_factorizations;
    updates_since_factor = 0;
    try {
        return lu.load(&cols[0], m) == SLinSolver::OK;
    } catch (...) {
        return false;
    }
}

bool DFBBasis::change_basis(int j, int col) {
    if (j < 0 || j >= m || col < 0 || col >= na + m) return false;

    const int out = basic[j];
    if (out == col) return true;
    where_basic[out] = -1;
    basic[j] = col;
    where_basic[col] = j;

    if (refactor_always || updates_since_factor >= refactor_period)
        return factorize();

    /* Actualizacion tipo Forrest-Tomlin: mucho mas barata que refactorizar,
     * pero acumula error, de ahi el refactor_period. */
    const SVector& newcol = column(col);
    try {
        lu.solveRight4update(*eta, newcol);
        SLinSolver::Status st = lu.change(j, newcol, eta);
        ++n_updates; ++updates_since_factor;
        if (st != SLinSolver::OK) return factorize();
    } catch (...) {
        return factorize();
    }
    return true;
}

void DFBBasis::btran_row(int j, double* out) {
    try {
        for (int r = 0; r < m; ++r) rhs[r] = 0.0;
        rhs[j] = 1.0;
        lu.solveLeft(sol, rhs);
        ++n_btran;
        for (int r = 0; r < m; ++r) out[r] = sol[r];
    } catch (...) {
        for (int r = 0; r < m; ++r) out[r] = 0.0;   /* SoPlex tiro: se degrada */
    }
}

void DFBBasis::ftran_col(int col, double* out) {
    try {
        for (int r = 0; r < m; ++r) rhs[r] = 0.0;
        const SVector& c = column(col);
        for (int t = 0; t < c.size(); ++t) rhs[c.index(t)] = c.value(t);
        lu.solveRight(sol, rhs);
        ++n_ftran;
        for (int r = 0; r < m; ++r) out[r] = sol[r];
    } catch (...) {
        for (int r = 0; r < m; ++r) out[r] = 0.0;   /* SoPlex tiro: se degrada */
    }
}

void DFBBasis::ftran_vec(const double* v, double* out) {
    try {
        for (int r = 0; r < m; ++r) rhs[r] = v[r];
        lu.solveRight(sol, rhs);
        ++n_ftran;
        for (int r = 0; r < m; ++r) out[r] = sol[r];
    } catch (...) {
        for (int r = 0; r < m; ++r) out[r] = 0.0;   /* SoPlex tiro: se degrada */
    }
}

void DFBBasis::row_of_A(const double* lambda, double* out) const {
    for (int c = 0; c < na; ++c) out[c] = 0.0;
    for (int j = 0; j < m; ++j) {
        const double l = lambda[j];
        if (l == 0.0) continue;
        const std::vector<std::pair<int,double> >& r = row_nz[j];
        for (size_t t = 0; t < r.size(); ++t) out[r[t].first] += l * r[t].second;
    }
}

void DFBBasis::A_times(const double* u, double* out) const {
    for (int j = 0; j < m; ++j) {
        double acc = 0.0;
        const std::vector<std::pair<int,double> >& r = row_nz[j];
        for (size_t t = 0; t < r.size(); ++t) acc += r[t].second * u[r[t].first];
        out[j] = acc;
    }
}

} /* namespace ibex */
