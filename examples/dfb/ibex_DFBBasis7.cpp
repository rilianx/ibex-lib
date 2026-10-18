/* Base del motor revisado de SC con la LU de SoPlex 7 (variante `-DDFB_LU7`,
 * compilada con -Dsoplex=soplex7 contra los encabezados de SoPlex 7). Mismo
 * algoritmo que ibex_DFBBasis.cpp; solo cambia la biblioteca de la LU. */
#include "soplex.h"
#include "soplex/slufactor.h"
#include "ibex_DFBBasis.h"
#include <cmath>
#include <cstdlib>
#include <memory>

using namespace soplex;

namespace ibex {

struct DFBBasis::Impl7 {
    std::shared_ptr<Tolerances> tol;
    std::vector<DSVector> col_nz, unit_nz;
    SLUFactor<Real> lu;
    DVector rhs, sol;
    SSVector* eta;
    Impl7() : tol(std::make_shared<Tolerances>()), rhs(1), sol(1), eta(NULL) { lu.setTolerances(tol); }
    ~Impl7() { delete eta; }
    const SVector& column(int c, int na) const {
        return (c < na) ? (const SVector&)col_nz[c] : (const SVector&)unit_nz[c - na];
    }
};

long DFBBasis::n_factorizations = 0;
long DFBBasis::n_updates       = 0;
long DFBBasis::n_btran         = 0;
long DFBBasis::n_ftran         = 0;
bool DFBBasis::refactor_always  = (getenv("DFB_REFACTOR_ALWAYS") != NULL);
int  DFBBasis::refactor_period  = (getenv("DFB_REFACTOR_PERIOD") ? atoi(getenv("DFB_REFACTOR_PERIOD")) : 50);

void DFBBasis::reset_counters() { n_factorizations = n_updates = n_btran = n_ftran = 0; }

DFBBasis::DFBBasis() : m(0), na(0), obj_row(0), updates_since_factor(0), col_preparada(-1), p7(new Impl7()) {}
DFBBasis::~DFBBasis() { delete p7; }

void DFBBasis::load_matrix(const Matrix& Abar) {
    m  = Abar.nb_rows();
    na = Abar.nb_cols();
    row_nz.assign(m, std::vector<std::pair<int,double> >());
    p7->col_nz.clear();  p7->col_nz.reserve(na);
    p7->unit_nz.clear(); p7->unit_nz.reserve(m);
    for (int c = 0; c < na; ++c) p7->col_nz.push_back(DSVector(m));
    for (int j = 0; j < m; ++j) { p7->unit_nz.push_back(DSVector(1)); p7->unit_nz[j].add(j, 1.0); }
    for (int j = 0; j < m; ++j)
        for (int c = 0; c < na; ++c) {
            const double v = Abar[j][c];
            if (v == 0.0) continue;
            row_nz[j].push_back(std::make_pair(c, v));
            p7->col_nz[c].add(j, v);
        }
    p7->rhs.reDim(m); p7->sol.reDim(m);
    delete p7->eta; p7->eta = new SSVector(m, p7->tol);
}

bool DFBBasis::init_basis(int k, double tol) {
    if (m <= 0 || k < 0 || k >= na) return false;
    basic.assign(m, 0);
    where_basic.assign(na + m, -1);
    for (int j = 0; j < m; ++j) { basic[j] = na + j; where_basic[na + j] = j; }
    updates_since_factor = 0;
    if (!factorize()) return false;
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
        basic[j] = c; where_basic[c] = j;
    }
    updates_since_factor = 0;
    return factorize();
}

bool DFBBasis::factorize() {
    std::vector<const SVector*> cols(m);
    for (int j = 0; j < m; ++j) cols[j] = &p7->column(basic[j], na);
    ++n_factorizations;
    updates_since_factor = 0;
    col_preparada = -1;
    try { return p7->lu.load(&cols[0], m) == SLinSolver<Real>::OK; }
    catch (...) { return false; }
}

bool DFBBasis::change_basis(int j, int col) {
    if (j < 0 || j >= m || col < 0 || col >= na + m) return false;
    const int out = basic[j];
    if (out == col) return true;
    where_basic[out] = -1; basic[j] = col; where_basic[col] = j;
    if (refactor_always || updates_since_factor >= refactor_period) return factorize();
    const SVector& newcol = p7->column(col, na);
    try {
        p7->lu.solveRight4update(*p7->eta, newcol);
        typename SLinSolver<Real>::Status st = p7->lu.change(j, newcol, p7->eta);
        ++n_updates; ++updates_since_factor;
        if (st != SLinSolver<Real>::OK) return factorize();
    } catch (...) { return factorize(); }
    return true;
}

bool DFBBasis::preparar_cambio(int col, double* out) {
    col_preparada = -1;
    if (col < 0 || col >= na + m) return false;
    try {
        p7->lu.solveRight4update(*p7->eta, p7->column(col, na));
        ++n_ftran;
        for (int r = 0; r < m; ++r) out[r] = (*p7->eta)[r];
        col_preparada = col;
        return true;
    } catch (...) { for (int r = 0; r < m; ++r) out[r] = 0.0; return false; }
}

bool DFBBasis::aplicar_cambio(int j, int col) {
    if (col_preparada != col) { col_preparada = -1; return change_basis(j, col); }
    col_preparada = -1;
    if (j < 0 || j >= m || col < 0 || col >= na + m) return false;
    const int out = basic[j];
    if (out == col) return true;
    where_basic[out] = -1; basic[j] = col; where_basic[col] = j;
    if (refactor_always || updates_since_factor >= refactor_period) return factorize();
    try {
        typename SLinSolver<Real>::Status st = p7->lu.change(j, p7->column(col, na), p7->eta);
        ++n_updates; ++updates_since_factor;
        if (st != SLinSolver<Real>::OK) return factorize();
    } catch (...) { return factorize(); }
    return true;
}

void DFBBasis::btran_row(int j, double* out) {
    try {
        for (int r = 0; r < m; ++r) p7->rhs[r] = 0.0;
        p7->rhs[j] = 1.0;
        p7->lu.solveLeft(p7->sol, p7->rhs);
        ++n_btran;
        for (int r = 0; r < m; ++r) out[r] = p7->sol[r];
    } catch (...) { for (int r = 0; r < m; ++r) out[r] = 0.0; }
}

void DFBBasis::ftran_col(int col, double* out) {
    try {
        for (int r = 0; r < m; ++r) p7->rhs[r] = 0.0;
        const SVector& c = p7->column(col, na);
        for (int t = 0; t < c.size(); ++t) p7->rhs[c.index(t)] = c.value(t);
        p7->lu.solveRight(p7->sol, p7->rhs);
        ++n_ftran;
        for (int r = 0; r < m; ++r) out[r] = p7->sol[r];
    } catch (...) { for (int r = 0; r < m; ++r) out[r] = 0.0; }
}

void DFBBasis::ftran_vec(const double* v, double* out) {
    try {
        for (int r = 0; r < m; ++r) p7->rhs[r] = v[r];
        p7->lu.solveRight(p7->sol, p7->rhs);
        ++n_ftran;
        for (int r = 0; r < m; ++r) out[r] = p7->sol[r];
    } catch (...) { for (int r = 0; r < m; ++r) out[r] = 0.0; }
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
