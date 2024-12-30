#include "ibex.h"
#include "DualFeasibleBounding.h"
#include "TestingDualFeasibleBounding.h"

int main() {
    standar_test();

    standar_test2();

    no_solution_test();

    no_solution_test2();

    unbounded_test();

    ill_conditioned_matrix_test();

    //dimensions_15x20_test();
    return 0;
}
