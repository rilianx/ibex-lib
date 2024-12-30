#include "ibex.h"
#include "IntervalStructuresOperations.h"

using namespace std;
using namespace ibex;

// Producto punto entre columnas de una matriz y un vector
IntervalVector IntervalStructuresOperations::matrixVectorColumnProduct(const IntervalVector& vector, const IntervalMatrix& matrix) {
    int m = matrix.nb_rows();
    int n = matrix.nb_cols();
    IntervalVector result(n, Interval(0));
    for (int i = 0; i < n; ++i) {
        Interval sum(0);
        for (int j = 0; j < m; ++j) {
            sum += matrix[j][i] * vector[j];
        }
        result[i] = sum;
    }
    return result;
}

// Divide un vector por un escalar
IntervalVector IntervalStructuresOperations::divideVectorByScalar(const IntervalVector& row, const Interval& divisor) {
    if (divisor.lb() == 0 || divisor.ub() == 0) {
        throw std::invalid_argument("El divisor tiene un borde igual a 0, lo que no es permitido.");
    }
    IntervalVector result(row.size());
    for (int i = 0; i < row.size(); ++i) {
        result[i] = row[i] / divisor;
    }
    return result;
}

// Intercambia filas de una matriz
void IntervalStructuresOperations::swapRows(IntervalMatrix& matrix, int row1, int row2) {
    if (row1 < 0 || row1 >= matrix.nb_rows() || row2 < 0 || row2 >= matrix.nb_rows()) {
        throw std::out_of_range("Índice de fila fuera de rango.");
    }
    for (int j = 0; j < matrix.nb_cols(); ++j) {
        std::swap(matrix[row1][j], matrix[row2][j]);
    }
}

// Hace identidad una columna de la matriz mediante Operaciones Elementales Fila
size_t IntervalStructuresOperations::makeColumnIdentity(
    IntervalMatrix& A, 
    size_t k, 
    bool interchange, 
    size_t j, 
    const set<size_t> identity_rows
) {
    size_t m = A.nb_rows();
    size_t n = A.nb_cols();

    // Paso 1: Encontrar la fila adecuada para el pivoteo
    if (!identity_rows.empty()) {
        for (size_t jAux = 0; jAux < m; ++jAux) {
            if (identity_rows.count(jAux) > 0){
                j = jAux;

                if (interchange){
                    swapRows(A, j, jAux);
                    break;
                }
            }
        }
    }

    // Paso 2: Normalizar la fila j respecto del valor en posición k
    A[j] = divideVectorByScalar(A[j], A[j][k]);
    A[j][k] = Interval(1.0);

    // Paso 3: Hacer ceros los demás elementos en la columna k
    for (size_t jj = 0; jj < m; ++jj) {
        if (jj == j) continue;
        Interval factor = A[jj][k];
        for (size_t i = 0; i < n; ++i) {
            A[jj][i] -= factor * A[j][i];
        }
        A[jj][k] = Interval(0.0);
    }

    return j;
}

// Aplica la fórmula de Gauss-Seidel sobre el vector x
Interval IntervalStructuresOperations::gaussSeidel(
    size_t k, 
    IntervalVector& x, 
    IntervalVector& gamma
){
    double epsilon = 1e-3;
    size_t n = x.size();

    if (gamma[k].lb() < 1 - epsilon || gamma[k].ub() > 1 + epsilon) {
        throw std::invalid_argument("gamma[k] debía ser 1");
    }

    if (k == -1) return Interval(0);

    gamma[k] = Interval(1);

    Interval xContract = Interval(0);
    for (size_t i = 0; i < n; ++i) {
        if (i != k) {
            xContract = xContract - gamma[i] * x[i];
        }
    }
    
    
    if (x[k].intersects(xContract))
    {
        x[k] = x[k] & xContract;
    }

    return x[k];
}

// Encuentra el índice del valor máximo en un vector
std::pair<Interval, size_t> IntervalStructuresOperations::getMaxValue(const IntervalVector& vector) {
    Interval max_value = vector[0];
    size_t max_index = 0;
    for (size_t i = 1; i < vector.size(); ++i) {
        if (vector[i].lb() > max_value.lb()) {
            max_value = vector[i];
            max_index = i;
        }
    }
    return {max_value, max_index};
}

// Multiplica un vector por un escalar
IntervalVector IntervalStructuresOperations::scalarMultiply(const IntervalVector& vector, const Interval& scalar) {
    IntervalVector result(vector.size());
    for (size_t i = 0; i < vector.size(); ++i) {
        result[i] = vector[i] * scalar;
    }
    return result;
}
