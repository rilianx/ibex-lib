#include "DualFeasibleBounding.h"
#include "IntervalStructuresOperations.h"
#include <cmath>
#include <tuple>

using namespace std;
using namespace ibex;

int DualFeasibleBounding::contract(IntervalVector& x_new, int maxIters, bool improveUpper) {
    if (improveUpper) {
        adjustSignsOfSystem(A, x_new);
    }

    int iters = 0;
    while (maxIters == -1 || iters < maxIters) {
        size_t j;
        Interval delta, direction;
        tie(j, delta, direction) = largestImpact(A, x_new, A[0], k);
       
        if (j == -1) {
            if (improveUpper) {
                adjustSignsOfSystem(A, x_new);
            }
            return iters;
        }

        Interval alpha;
        size_t i;
        tie(alpha, i) = calculateAlpha(A[j], A[0], direction);

        if (i == -1) {
            x_new.set_empty();
            continue;
        }

        A[0] = A[0] + Operator.scalarMultiply(A[j], alpha);
        A[0][i] = Interval(0.0);

        identity_rows.insert(j);
        Operator.makeColumnIdentity(A, i, false, j);
        A[0][k] = Interval(1.0);
        Operator.gaussSeidel(k, x_new, A[0]);

        if (contract_all) {
            for (const auto& r : identity_rows) {
                x_new[r] = Operator.gaussSeidel(-1, x_new, A[r]);
            }
        }

        ++iters;
    }

    if (improveUpper) {
        adjustSignsOfSystem(A, x_new);
    }

    return iters;
}

std::pair<IntervalVector, IntervalVector> DualFeasibleBounding::calculateImpacts(
    const IntervalMatrix& A, const IntervalVector& x_new, const IntervalVector& gamma, size_t k) {
    int m = A.nb_rows();
    int n = A.nb_cols();

    IntervalVector grad_incr = IntervalVector(m, Interval(0));
    IntervalVector grad_decr = IntervalVector(m, Interval(0));

    for (int j = 1; j < m; ++j) {
        IntervalVector x_prime_incr = IntervalVector(n, Interval(0));
        IntervalVector x_prime_decr = IntervalVector(n, Interval(0));

        for (int i = 0; i < n; ++i) {
            if (gamma[i].lb() > 0) {
                x_prime_incr[i] = Interval(x_new[i].ub());
                x_prime_decr[i] = Interval(x_new[i].ub());
            } else if (gamma[i].ub() < 0) {
                x_prime_incr[i] = Interval(x_new[i].lb());
                x_prime_decr[i] = Interval(x_new[i].lb());
            }
        }
        for (int i = 0; i < n; ++i) {
            if (gamma[i].lb() == 0 && gamma[i].ub() == 0) {
                if (A[j][i].lb() > 0){
                    x_prime_incr[i] = Interval(x_new[i].ub());
                    x_prime_decr[i] = Interval(x_new[i].lb());
                } else if (A[j][i].ub() < 0) {
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

std::tuple<size_t, Interval, Interval> DualFeasibleBounding::largestImpact(
    const IntervalMatrix& A, const IntervalVector& x_new, const IntervalVector& gamma, size_t k) {
    IntervalVector grad_incr, grad_decr;
    tie(grad_incr, grad_decr) = calculateImpacts(A, x_new, gamma, k);

    Interval delta_incr, delta_decr;
    size_t j_incr, j_decr;
    tie(delta_incr, j_incr) = Operator.getMaxValue(grad_incr);
    tie(delta_decr, j_decr) = Operator.getMaxValue(grad_decr);

    Interval delta = Interval(0.0);
    size_t j;
    Interval direction;

    if (delta_incr.lb() > delta_decr.lb()) {
        delta = delta_incr;
        j = j_incr;
        direction = Interval(1.0);
    } else {
        delta = delta_decr;
        j = j_decr;
        direction = Interval(-1.0);
    }

    if (delta.lb() <= 0.0) {
        j = -1;
    }

    return {j, delta, direction};
}

std::pair<Interval, size_t> DualFeasibleBounding::calculateAlpha(
    const IntervalVector& Aj, const IntervalVector& gamma, const Interval& direction) {
    Interval min_alpha(1e20);
    size_t min_idx = -1;

    for (size_t i = 0; i < Aj.size(); ++i) {
        if (Aj[i].lb() != 0 && Aj[i].ub() != 0) {
            Interval possible_alpha = (gamma[i] / Aj[i]) * direction;
            if (possible_alpha.ub() < 0) {
                if (abs(possible_alpha.lb()) < abs(min_alpha.lb())) {
                    min_alpha = possible_alpha;
                    min_idx = i;
                }
            }
        }
    }

    min_alpha = min_alpha * direction * Interval(-1.0);
    return {min_alpha, min_idx};
}

void DualFeasibleBounding::adjustSignsOfSystem(IntervalMatrix& A, IntervalVector& x_new) {
    A[0] = -A[0];
    A[0][k] = -A[0][k];
    x_new[k] = -x_new[k];
}

int DualFeasibleBounding::totalContraction(IntervalVector& x_new) {
    int total_iters = 0;

    int itersLB = contract(x_new, -1);

    int itersUB = contract(x_new, -1, true);

    total_iters = itersLB + itersUB;
    return total_iters;
}
