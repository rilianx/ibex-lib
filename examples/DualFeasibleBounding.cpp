#include "DualFeasibleBounding.h"
#include <cmath>
#include <tuple>

using namespace std;
using namespace ibex;

int DualFeasibleBounding::contract(IntervalVector& x_new, int maxIters, bool improveUpper) {
    IntervalMatrix A = getMatrix();
    Interval alpha;
    int i;
    int j;
    Interval delta, direction;
    int k = getKIndex();


    if (improveUpper) {
        changeSigns(A, x_new);
    }


    int iters = 0;
    
    gaussSeidel(x_new, k, A[0]);
    while (maxIters == -1 || iters < maxIters) {
        tie(j, delta, direction) = largestImpact(A, x_new, A[0]);
       
        if (j == -1) {
            if (improveUpper) {
                changeSigns(A, x_new);
            }
            setMatrix(A);
            return iters;
        } 

        tie(alpha, i) = calculateAlpha(A[j], A[0], direction);

        if (i == -1) {
            x_new.set_empty();
            continue;
        }

        A[0] = A[0] + alpha * A[j];
        
        A[0][i] = Interval(0);


        makeColumnIdentity(A, i, false, j);
        A[0][k] = Interval(1);
        gaussSeidel(x_new, k, A[0]);


        identity_rows[j] = i;
        if (contract_all) {
            for (const auto& r : identity_rows) {
            gaussSeidel(x_new, r.second, A[r.first]);
        }
}
        ++iters;
    }

    if (improveUpper) {
        changeSigns(A, x_new);
    }
    setMatrix(A);

    return iters;
}

int DualFeasibleBounding::totalContraction(IntervalVector& x_new) {
    int total_iters = 0;

    int itersLB = contract(x_new, -1);

    int itersUB = contract(x_new, -1, true);

    total_iters = itersLB + itersUB;
    return total_iters;
}

std::pair<IntervalVector, IntervalVector> DualFeasibleBounding::calculateImpacts(
    const IntervalMatrix& A, const IntervalVector& x_new, const IntervalVector& gamma) {
    int m = A.nb_rows();
    int n = A.nb_cols();

    IntervalVector grad_incr = IntervalVector(m, Interval(0));
    IntervalVector grad_decr = IntervalVector(m, Interval(0));

    for (int j = 1; j < m; ++j) {
        IntervalVector x_prime_incr = IntervalVector(n, Interval(0));
        IntervalVector x_prime_decr = IntervalVector(n, Interval(0));

        for (int i = 0; i < n; ++i) {
            Interval signOfGammai = sign(gamma[i]);
            if (signOfGammai == Interval(1)) {
                x_prime_incr[i] = Interval(x_new[i].ub());
                x_prime_decr[i] = Interval(x_new[i].ub());
            } else if (signOfGammai == Interval(-1)) {
                x_prime_incr[i] = Interval(x_new[i].lb());
                x_prime_decr[i] = Interval(x_new[i].lb());
            }
            else{
                Interval signOfAji = sign(A[j][i]);
                if (signOfAji == Interval(1)){
                    x_prime_incr[i] = Interval(x_new[i].ub());
                    x_prime_decr[i] = Interval(x_new[i].lb());
                } else if (signOfAji == Interval(-1)) {
                    x_prime_incr[i] = Interval(x_new[i].lb());
                    x_prime_decr[i] = Interval(x_new[i].ub());
                }
            }
        }

        grad_incr[j] = Interval(A[j] * x_prime_incr);
        grad_incr[j] = - grad_incr[j];
        grad_decr[j] = Interval(A[j] * x_prime_decr);
    }

    return {grad_incr, grad_decr};
}

std::tuple<int, Interval, Interval> DualFeasibleBounding::largestImpact(
    const IntervalMatrix& A, const IntervalVector& x_new, const IntervalVector& gamma) {
    IntervalVector grad_incr, grad_decr;
    tie(grad_incr, grad_decr) = calculateImpacts(A, x_new, gamma);

    Interval delta_incr, delta_decr;
    int j_incr, j_decr;
    tie(delta_incr, j_incr) = getMaxValue(grad_incr);
    tie(delta_decr, j_decr) = getMaxValue(grad_decr);

    Interval delta = Interval(0);
    int j;
    Interval direction;

    if (delta_incr.lb() > delta_decr.lb()) {
        delta = delta_incr;
        j = j_incr;
        direction = Interval(1);
    } else {
        delta = delta_decr;
        j = j_decr;
        direction = Interval(-1);
    }

    Interval signOfDelta = sign(delta);
    if (signOfDelta != Interval(1)) {
        j = -1;
    }

    return {j, delta, direction};
}

std::pair<Interval, int> DualFeasibleBounding::calculateAlpha(
    const IntervalVector& Aj, const IntervalVector& gamma, const Interval& direction) {
    Interval min_alpha(1e20);
    int min_index = -1;
    int n = Aj.size();

    for (int i = 0; i < n; ++i) {
        if (Aj[i].lb() != 0 && Aj[i].ub() != 0) {
            Interval alpha = (gamma[i] / Aj[i]) * direction;
            Interval signOfAlpha = sign(alpha);
            if (signOfAlpha == Interval(-1)) {
                if (abs(alpha.lb()) < abs(min_alpha.lb())) {
                    min_alpha = alpha;
                    min_index = i;
                }
            }
        }
    }

    min_alpha = min_alpha * direction * Interval(-1.0);
    return {min_alpha, min_index};
}

void DualFeasibleBounding::changeSigns(IntervalMatrix& A, IntervalVector& x_new) {
    int k = getKIndex();
    A[0] = -A[0];
    A[0][k] = -A[0][k];
    x_new[k] = -x_new[k];
}

                                // AUX FUNCTIONS
int DualFeasibleBounding::makeColumnIdentity(IntervalMatrix& A, const int k, bool interchange, int j, 
                                             map<int, int> identity_rows){
    int m = A.nb_rows();
    int n = A.nb_cols();

    // Paso 1: Encontrar la fila adecuada para el pivoteo
    if (interchange){
        if (!identity_rows.empty()) {
            for (int jAux = 0; jAux < m; ++jAux) {
                if (identity_rows.count(jAux) > 0){
                    j = jAux;

                    for (int i = 0; i < n; ++j) {
                        std::swap(A[j][i], A[jAux][i]);
                    }
                    break;
                }
            }
        }
    }

    

    // Paso 2: Normalizar la fila j respecto del valor en posición k
    if (A[j][k].lb() == 0 || A[j][k].ub() == 0) {
        throw std::invalid_argument("The divider has a bound equal to 0. Which is not allowed.");
    }
    
    A[j] = (Interval(1) / A[j][k]) * A[j];

    A[j][k] = Interval(1);

    // Paso 3: Hacer ceros los demás elementos en la columna k
    for (int jj = 0; jj < m; ++jj) {
        if (jj == j) continue;
        Interval factor = A[jj][k];
        for (int i = 0; i < n; ++i) {
            A[jj][i] -= factor * A[j][i];
        }
        A[jj][k] = Interval(0);
    }

    return j;
}

// Aplica la fórmula de Gauss-Seidel sobre el vector x
Interval DualFeasibleBounding::gaussSeidel(IntervalVector& x, int k, IntervalVector& gamma){
    double epsilon = 1e-6;
    int n = gamma.size();

    if (k == -1 || k >= n){
        throw std::invalid_argument("Invalid k.");
    }

    if (gamma[k] != Interval(1)){
        throw std::invalid_argument("Gamma in k is not 1.");
    }

    gamma[k] = Interval(0);

    Interval xContract = -(gamma * x);

    gamma[k] = Interval(1);
    
    if (x[k].intersects(xContract))
    {
        if (contract_all){
            std::cout << "x[k] = " << x[k] << std::endl;
            std::cout << "xContract = " << xContract << std::endl;
        }
        x[k] = x[k] & xContract;
        if (contract_all){
            std::cout << "x[k] = " << x[k] << std::endl;
        }
    }

    return x[k];
}

// Encuentra el índice del valor máximo en un vector
std::pair<Interval, int> DualFeasibleBounding::getMaxValue(const IntervalVector& vector) {
    Interval max_value = vector[0];
    int max_index = 0;
    int n = vector.size();

    for (int i = 1; i < n; ++i) {
        if (vector[i].lb() > max_value.lb()) {
            max_value = vector[i];
            max_index = i;
        }
    }

    return {max_value, max_index};
}
