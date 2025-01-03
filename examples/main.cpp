#include "ibex.h"
#include "DualFeasibleBounding.h"
#include "TestingDualFeasibleBounding.h"

int main() {
    //standar_test();

 //   standar_test2();

    //no_solution_test();

    //no_solution_test2();

  //  unbounded_test();

//    ill_conditioned_matrix_test();

   // dimensions_15x20_test();

    /*int n = 0;
    int m = 0;
    cout << "Ingrese el valor de n: ";
    std::cin >> n;

    cout << "Ingrese el valor de m: ";
    std::cin >> m;

    random_test(m, n);*/
    special_test_contraction_box();
    
    return 0;
}
