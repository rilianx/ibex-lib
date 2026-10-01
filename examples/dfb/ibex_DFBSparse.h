//============================================================================
//                                  I B E X
// File        : ibex_DFBSparse.h
// Author      : (DFB) tableau aumentado disperso
//============================================================================

#ifndef __IBEX_DFB_SPARSE_H__
#define __IBEX_DFB_SPARSE_H__

#include "ibex_Matrix.h"

#include <vector>

namespace ibex {

/**
 * \brief Tableau aumentado `[Abar | I]` con las filas dispersas.
 *
 * Reemplaza la matriz densa `Af` de `CtcDFB` sin cambiar el algoritmo: las
 * mismas operaciones de fila, sobre la misma secuencia de pivotes, pero
 * recorriendo solo las entradas no nulas. Medido en
 * [MEDICIONES_FACTORIZACION.md](MEDICIONES_FACTORIZACION.md):
 *
 *  - el pricing recorre `(m-1)*na` entradas de las que solo el **2.3 %** son no
 *    nulas en las instancias con `n >= 40`;
 *  - un pivote modifica `nnz(columna entrante)` filas, que es un numero
 *    **constante** —16.8 con `m = 400` y 16.7 con `m = 280`— y no una fraccion
 *    de `m`, mientras la eliminacion densa recorre las `m`.
 *
 * Sobre BroydenBanded-200 eso son 604 400 operaciones por pivote contra ~25 000
 * utiles en la eliminacion, y 283 689 contra ~6 500 en el pricing.
 *
 * ## Representacion
 *
 * El espacio de columnas es unico, `0..na+m-1`: las primeras `na` son la parte
 * A y las ultimas `m` la parte lambda. Que sea unico simplifica las
 * combinaciones de filas, que tienen que tocar las dos partes para que lambda
 * siga siendo veraz.
 *
 *  - **Las filas `1..m-1` son dispersas**, con los indices ordenados.
 *  - **La fila 0 es densa** (`g` para la parte A, `lam` para la parte lambda).
 *    Es la fila objetivo: el pricing y el ratio test la leen en orden aleatorio
 *    por columna, asi que densa es la representacion correcta, y ademas es la
 *    que acumula fill-in. Cuesta `na + m` doubles, no `m*(na+m)`.
 *
 * Las entradas que quedan exactamente en cero se descartan; las que quedan
 * chicas pero no nulas **se conservan**, porque el camino denso tambien las
 * conserva y el objetivo es que las dos representaciones den la misma secuencia
 * de pivotes.
 */
class DFBSparseTableau {
public:
    DFBSparseTableau();

    /** \brief Construye `[Abar | I]`: fila j = `[Abar_j | e_j]`. */
    void build(const Matrix& Abar);

    /** \brief Vacia el tableau (para reconstruirlo). */
    void clear();

    /**
     * \brief Hace identidad la columna `col` en la fila `row`.
     *
     * Equivalente a `CtcDFB::make_column_identity_f`: normaliza la fila `row`
     * para que su entrada en `col` valga 1, y elimina `col` de todas las demas
     * filas. Con `interchange` y `row == 0`, si la entrada de la fila 0 es
     * demasiado chica se intercambia con la primera fila que sirva.
     */
    bool make_identity(int col, int row, bool interchange, double tol,
                       int* swapped_with = NULL);

    /** \brief `fila0 += alpha * fila_j` (sobre las dos partes). */
    void add_to_row0(double alpha, int j);

    /** \brief Niega la fila 0 completa (para `flip_side`). */
    void negate_row0();

    /** \brief Entrada (j, col); busqueda binaria en las filas dispersas. */
    double entry(int j, int col) const;

    int nb_rows() const { return m; }
    int nb_cols_A() const { return na; }
    int nb_cols() const { return nc; }

    /* Fila 0, densa. */
    std::vector<double> g;     //!< parte A de la fila 0 (tamano na)
    std::vector<double> lam;   //!< parte lambda de la fila 0 (tamano m)

    /* Filas 1..m-1, dispersas y ordenadas por indice. */
    const std::vector<int>&    idx(int j) const { return sidx[j]; }
    const std::vector<double>& val(int j) const { return sval[j]; }
    int nnz(int j) const { return (int)sidx[j].size(); }

    /** \brief Entradas totales almacenadas (diagnostico de fill-in). */
    long total_nnz() const;

    /** \brief Copia el estado de otro tableau (para `init_from_lower`). */
    void copy_from(const DFBSparseTableau& o);

private:
    int m, na, nc;
    std::vector<std::vector<int> >    sidx;
    std::vector<std::vector<double> > sval;

    /* Buffers reusados por las combinaciones de filas. */
    std::vector<int>    bidx;
    std::vector<double> bval;

    /** \brief `fila_r -= f * fila_j`, las dos dispersas. Descarta el cero
     *  exacto y fuerza a 0 la entrada de `drop_col`. */
    void axpy_sparse(int r, int j, double f, int drop_col);

    /** \brief `fila_r -= f * fila0`, con la fila 0 densa. */
    void axpy_from_row0(int r, double f, int drop_col);
};

} /* namespace ibex */

#endif /* __IBEX_DFB_SPARSE_H__ */
