#include "ibex_DFBSparse.h"

#include <algorithm>
#include <cmath>

namespace ibex {

DFBSparseTableau::DFBSparseTableau() : m(0), na(0), nc(0) {
}

void DFBSparseTableau::clear() {
    sidx.clear(); sval.clear(); g.clear(); lam.clear();
    m = na = nc = 0;
}

void DFBSparseTableau::build(const Matrix& Abar) {
    m  = Abar.nb_rows();
    na = Abar.nb_cols();
    nc = na + m;

    g.assign(na, 0.0);
    lam.assign(m, 0.0);
    sidx.assign(m, std::vector<int>());
    sval.assign(m, std::vector<double>());

    /* fila 0 densa */
    for (int c = 0; c < na; ++c) g[c] = Abar[0][c];
    lam[0] = 1.0;

    /* filas 1..m-1 dispersas: [Abar_j | e_j] */
    for (int j = 1; j < m; ++j) {
        std::vector<int>&    ix = sidx[j];
        std::vector<double>& vl = sval[j];
        for (int c = 0; c < na; ++c) {
            const double v = Abar[j][c];
            if (v == 0.0) continue;
            ix.push_back(c); vl.push_back(v);
        }
        ix.push_back(na + j); vl.push_back(1.0);
    }

    bidx.reserve(nc); bval.reserve(nc);
}

void DFBSparseTableau::copy_from(const DFBSparseTableau& o) {
    m = o.m; na = o.na; nc = o.nc;
    g = o.g; lam = o.lam;
    sidx = o.sidx; sval = o.sval;
    bidx.reserve(nc); bval.reserve(nc);
}

long DFBSparseTableau::total_nnz() const {
    long t = (long)na + (long)m;      /* la fila 0 es densa */
    for (int j = 1; j < m; ++j) t += (long)sidx[j].size();
    return t;
}

double DFBSparseTableau::entry(int j, int col) const {
    if (j == 0) return (col < na) ? g[col] : lam[col - na];
    const std::vector<int>& ix = sidx[j];
    const std::vector<int>::const_iterator it =
        std::lower_bound(ix.begin(), ix.end(), col);
    if (it == ix.end() || *it != col) return 0.0;
    return sval[j][it - ix.begin()];
}

void DFBSparseTableau::negate_row0() {
    for (int c = 0; c < na; ++c) g[c] = -g[c];
    for (int r = 0; r < m; ++r) lam[r] = -lam[r];
}

void DFBSparseTableau::add_to_row0(double alpha, int j) {
    const std::vector<int>&    ix = sidx[j];
    const std::vector<double>& vl = sval[j];
    for (size_t t = 0; t < ix.size(); ++t) {
        const int c = ix[t];
        if (c < na) g[c]      += alpha * vl[t];
        else        lam[c-na] += alpha * vl[t];
    }
}

/* fila_r -= f * fila_j, las dos dispersas y ordenadas. */
void DFBSparseTableau::axpy_sparse(int r, int j, double f, int drop_col) {
    const std::vector<int>&    ai = sidx[r]; const std::vector<double>& av = sval[r];
    const std::vector<int>&    bi = sidx[j]; const std::vector<double>& bv = sval[j];

    bidx.clear(); bval.clear();
    size_t p = 0, q = 0;
    while (p < ai.size() || q < bi.size()) {
        int c; double v;
        if (q == bi.size() || (p < ai.size() && ai[p] < bi[q])) {
            c = ai[p]; v = av[p]; ++p;
        } else if (p == ai.size() || bi[q] < ai[p]) {
            /* OJO: `0.0 - f*b` y no `-f*b`. Ibex corre con redondeo dirigido
             * (-frounding-math, y gaol deja el modo hacia arriba), y con
             * redondeo dirigido el redondeo NO es simetrico respecto del signo:
             * round_up(-(f*b)) != -round_up(f*b). El camino denso calcula
             * `Af[i][c] -= f*Af[row][c]` con Af[i][c] == 0, o sea `0.0 - f*b`,
             * asi que hay que escribir la misma expresion para obtener el mismo
             * resultado bit a bit. Con `-f*b` la diferencia es de 1 ulp, se
             * amplifica en los pivotes siguientes y termina cambiando los
             * desempates del ratio test. */
            c = bi[q]; v = 0.0 - f * bv[q]; ++q;
        } else {
            c = ai[p]; v = av[p] - f * bv[q]; ++p; ++q;
        }
        if (c == drop_col) continue;          /* se fuerza a 0, como el denso */
        if (v == 0.0) continue;               /* cero exacto: no se almacena */
        bidx.push_back(c); bval.push_back(v);
    }
    sidx[r].swap(bidx); sval[r].swap(bval);
}

/* fila_r -= f * fila0, con la fila 0 densa. Se recorre la fila 0 salteando sus
 * ceros, de modo que el resultado tiene solo entradas realmente no nulas. */
void DFBSparseTableau::axpy_from_row0(int r, double f, int drop_col) {
    const std::vector<int>&    ai = sidx[r]; const std::vector<double>& av = sval[r];

    bidx.clear(); bval.clear();
    size_t p = 0;
    for (int c = 0; c < nc; ++c) {
        const double b = (c < na) ? g[c] : lam[c - na];
        const bool tiene_a = (p < ai.size() && ai[p] == c);
        if (b == 0.0 && !tiene_a) continue;
        double v = tiene_a ? av[p] : 0.0;
        if (tiene_a) ++p;
        v -= f * b;
        if (c == drop_col) continue;
        if (v == 0.0) continue;
        bidx.push_back(c); bval.push_back(v);
    }
    sidx[r].swap(bidx); sval[r].swap(bval);
}

bool DFBSparseTableau::make_identity(int col, int row, bool interchange,
                                     double tol, int* swapped_with) {
    if (swapped_with) *swapped_with = -1;
    if (m <= 0 || col < 0 || col >= nc) return false;

    if (interchange && std::fabs(entry(row, col)) < tol) {
        int found = -1;
        for (int i = 0; i < m; ++i)
            if (i != row && std::fabs(entry(i, col)) >= tol) { found = i; break; }
        if (found < 0) return false;
        if (swapped_with) *swapped_with = found;
        /* Intercambio de filas. La fila 0 es densa, asi que si esta en juego se
         * la convierte a dispersa y viceversa. */
        if (row == 0 || found == 0) {
            const int other = (row == 0) ? found : row;
            std::vector<int> ni; std::vector<double> nv;
            for (int c = 0; c < nc; ++c) {
                const double v = (c < na) ? g[c] : lam[c - na];
                if (v != 0.0) { ni.push_back(c); nv.push_back(v); }
            }
            for (int c = 0; c < na; ++c) g[c] = 0.0;
            for (int r = 0; r < m; ++r) lam[r] = 0.0;
            const std::vector<int>& oi = sidx[other];
            const std::vector<double>& ov = sval[other];
            for (size_t t = 0; t < oi.size(); ++t) {
                const int c = oi[t];
                if (c < na) g[c] = ov[t]; else lam[c - na] = ov[t];
            }
            sidx[other].swap(ni); sval[other].swap(nv);
        } else {
            sidx[row].swap(sidx[found]); sval[row].swap(sval[found]);
        }
    }

    const double piv = entry(row, col);
    if (std::fabs(piv) < tol) return false;
    const double inv = 1.0 / piv;

    if (row == 0) {
        for (int c = 0; c < na; ++c) g[c] *= inv;
        for (int r = 0; r < m; ++r) lam[r] *= inv;
        if (col < na) g[col] = 1.0; else lam[col - na] = 1.0;
    } else {
        std::vector<double>& vl = sval[row];
        for (size_t t = 0; t < vl.size(); ++t) vl[t] *= inv;
        const std::vector<int>& ix = sidx[row];
        const std::vector<int>::const_iterator it =
            std::lower_bound(ix.begin(), ix.end(), col);
        if (it != ix.end() && *it == col) vl[it - ix.begin()] = 1.0;
    }

    /* Eliminacion: solo las filas que tienen entrada no nula en col. */
    for (int i = 0; i < m; ++i) {
        if (i == row) continue;
        const double f = entry(i, col);
        if (f == 0.0) continue;
        if (i == 0) {
            /* fila0 -= f * fila_row, con fila_row dispersa */
            const std::vector<int>&    ix = sidx[row];
            const std::vector<double>& vl = sval[row];
            for (size_t t = 0; t < ix.size(); ++t) {
                const int c = ix[t];
                if (c < na) g[c]      -= f * vl[t];
                else        lam[c-na] -= f * vl[t];
            }
            if (col < na) g[col] = 0.0; else lam[col - na] = 0.0;
        } else if (row == 0) {
            axpy_from_row0(i, f, col);
        } else {
            axpy_sparse(i, row, f, col);
        }
    }
    return true;
}

} /* namespace ibex */
