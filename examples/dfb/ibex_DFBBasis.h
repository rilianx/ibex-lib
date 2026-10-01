//============================================================================
//                                  I B E X
// File        : ibex_DFBBasis.h
// Author      : (DFB) motor de base factorizada
//============================================================================

#ifndef __IBEX_DFB_BASIS_H__
#define __IBEX_DFB_BASIS_H__

#include "ibex_Matrix.h"
#include "soplex.h"
#include "soplex/slufactor.h"

#include <vector>

namespace ibex {

/**
 * \brief Base factorizada del tableau aumentado de DFB.
 *
 * DFB mantiene hoy el tableau aumentado `Af = [Abar | I]` materializado y lo
 * transforma con eliminaciones completas. Esta clase reemplaza esa
 * representacion por la que usa un simplex moderno: **la base y su
 * factorizacion LU**, calculando las filas del tableau a demanda.
 *
 * ## Por que es equivalente
 *
 * El pivoteo de DFB es, mirado de cerca, un cambio de base. Sea `Abar` la
 * matriz flotante `m x na` (con la columna `k` negada si se contrae la cota
 * superior) y `[Abar | I]` el tableau aumentado. En todo momento el tableau
 * transformado es
 *
 *     Af = B^-1 [Abar | I] = [ B^-1 Abar | B^-1 ]
 *
 * donde `B` es la matriz de columnas basicas de `[Abar | I]`, una por fila. De
 * ahi salen las dos cosas que DFB necesita:
 *
 *  - la **parte lambda** de una fila es una fila de `B^-1`, o sea `e_j^T B^-1`,
 *    que es exactamente un **BTRAN**;
 *  - la **parte A** de esa fila es `lambda^T Abar`.
 *
 * La fila objetivo (la que DFB llama fila 0) no es una excepcion: su columna
 * basica es `k`, y por eso tiene 1 en `k` y 0 en las demas columnas basicas,
 * que es la condicion de «costo reducido nulo sobre las basicas». El paso
 * `Af[0] += alpha*Af[j]` seguido de la eliminacion de la columna entrante es
 * precisamente la matriz elemental del cambio de base, con
 * `alpha = -t_0i/t_ji`.
 *
 * ## Estado
 *
 * `basic[j]` es el indice de columna basica en la fila `j`; las columnas
 * `0..na-1` son las de `Abar` y las `na..na+m-1` las de la identidad
 * (holguras). `obj_row` es la fila objetivo: normalmente 0, pero si
 * `|Abar[0][k]|` es demasiado chico el pivote inicial va a otra fila, igual que
 * el `interchange` del camino denso.
 *
 * La memoria pasa de `O(m*(na+m))` del tableau denso a `O(nnz(Abar))`, que es
 * lo que impide hoy terminar las instancias grandes.
 */
class DFBBasis {
public:
    DFBBasis();
    ~DFBBasis();

    /** \brief Guarda `Abar` en formato disperso por filas y por columnas.
     *
     * Se copia una vez por linealizacion. Las dos copias cuestan `2*nnz`, muy
     * por debajo del tableau denso que reemplazan. */
    void load_matrix(const Matrix& Abar);


    /** \brief Base de holguras y pivote de la columna `k` en la fila objetivo.
     *
     * Devuelve false si ninguna fila tiene una entrada usable en la columna
     * `k`, que es el mismo criterio de fallo que `make_column_identity_f`. */
    bool init_basis(int k, double tol);

    /** \brief Cambia la columna basica de la fila `j` por `col`. */
    bool change_basis(int j, int col);

    /** \brief Fija la base a `cols` (una columna por fila) y factoriza. */
    bool set_basis(const std::vector<int>& cols);

    /** \brief `out = e_j^T B^-1` (BTRAN). `out` tiene dimension `m`. */
    void btran_row(int j, double* out);

    /** \brief `out = B^-1 * columna(col)` (FTRAN). `out` tiene dimension `m`. */
    void ftran_col(int col, double* out);

    /** \brief `out = B^-1 * v` con `v` denso de dimension `m` (FTRAN). */
    void ftran_vec(const double* v, double* out);

    /** \brief `out = lambda^T Abar`, recorriendo solo los `lambda` no nulos.
     *
     * Misma tecnica que la certificacion: el dual de una base es disperso, asi
     * que el trabajo es `nnz(lambda) * nnz(fila)` y no `m*na`. */
    void row_of_A(const double* lambda, double* out) const;

    /** \brief `out = Abar * u` con `u` denso de dimension `na`. */
    void A_times(const double* u, double* out) const;

    int nb_rows() const { return m; }
    int nb_cols_A() const { return na; }
    int nb_cols() const { return na + m; }

    int basic_col(int j) const { return basic[j]; }
    int objective_row() const { return obj_row; }

    /** \brief true si la columna `c` es basica en alguna fila. */
    bool is_basic(int c) const { return where_basic[c] >= 0; }
    int  row_of_basic(int c) const { return where_basic[c]; }

    /** \brief Refactorizaciones y actualizaciones acumuladas (diagnostico). */
    static long n_factorizations;
    static long n_updates;
    static long n_btran;
    static long n_ftran;
    static void reset_counters();

    /** \brief Si es true se refactoriza en cada cambio de base en vez de usar
     *  `change()`. Es la etapa 1: mas lenta, pero sin la deriva de los updates,
     *  asi que sirve de referencia de correccion. */
    static bool refactor_always;

    /** \brief Cada cuantos updates se refactoriza de todos modos. */
    static int  refactor_period;

private:
    int m, na;
    int obj_row;

    /* Abar disperso, por filas y por columnas. */
    std::vector<std::vector<std::pair<int,double> > > row_nz;
    std::vector<soplex::DSVector> col_nz;      //!< columnas de Abar
    std::vector<soplex::DSVector> unit_nz;     //!< columnas de la identidad

    std::vector<int> basic;        //!< columna basica de cada fila
    std::vector<int> where_basic;  //!< fila de cada columna, o -1

    soplex::SLUFactor lu;
    soplex::DVector   rhs, sol;
    soplex::SSVector* eta;   //!< se recrea en load_matrix con la dimension m
    int updates_since_factor;

    const soplex::SVector& column(int c) const {
        return (c < na) ? (const soplex::SVector&)col_nz[c]
                        : (const soplex::SVector&)unit_nz[c - na];
    }

    bool factorize();
};

} /* namespace ibex */

#endif /* __IBEX_DFB_BASIS_H__ */
