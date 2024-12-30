#include "ibex.h"
#include "DualFeasibleBounding.h"
#include "TestingDualFeasibleBounding.h"

void test_case(IntervalMatrix& A, IntervalVector& x, bool visualize_iters){
    int m = A.nb_rows();
    int n = A.nb_cols();
    IntervalVector xOG = IntervalVector(x);
    std::vector<int> iters;

    for (int k = 0; k < n; ++k) {
        DualFeasibleBounding DFB(A, k);
        int iter = DFB.totalContraction(x);
        iters.push_back(iter);
    }
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "xOG: " << xOG << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "xDFB: " << x << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;


    if (visualize_iters) {
        std::cout << "ITERACIONES: ";
        for (const auto& iter : iters) {
            std::cout << iter << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "\n\n";
}

void standar_test(){
    std::cout << "                                   TESTING: STANDAR TEST 1" << std::endl;

    // Matriz A (3x5) de intervalos
    IntervalMatrix A(3, 5);
    A[0][0] = Interval(-7.31);
    A[0][1] = Interval(6.95);
    A[0][2] = Interval(5.28);
    A[0][3] = Interval(-4.90);
    A[0][4] = Interval(-0.08);

    A[1][0] = Interval(-1.01);
    A[1][1] = Interval(3.03);
    A[1][2] = Interval(5.77);
    A[1][3] = Interval(-8.12);
    A[1][4] = Interval(-9.43);

    A[2][0] = Interval(6.72);
    A[2][1] = Interval(-1.34);
    A[2][2] = Interval(5.25);
    A[2][3] = Interval(-9.96);
    A[2][4] = Interval(-1.09);

    // Inicializar el vector x con intervalos
    IntervalVector x = IntervalVector(5, Interval(0));
    x[0] = Interval(-1.565, 2.880);
    x[1] = Interval(-0.478, 4.463);
    x[2] = Interval(-1.038, 6.032);
    x[3] = Interval(-0.048, 3.615);
    x[4] = Interval(-1.076, 2.647);

    test_case(A, x, true);
}


void standar_test2(){
    std::cout << "                                   TESTING: STANDAR TEST 2" << std::endl;

    IntervalMatrix A(3, 5);
    A[0][0] = Interval(6.78);
    A[0][1] = Interval(0.45);
    A[0][2] = Interval(9.06);
    A[0][3] = Interval(7.21);
    A[0][4] = Interval(-5.20);

    A[1][0] = Interval(-1.54);
    A[1][1] = Interval(-4.08);
    A[1][2] = Interval(-3.41);
    A[1][3] = Interval(-3.07);
    A[1][4] = Interval(4.74);

    A[2][0] = Interval(-2.45);
    A[2][1] = Interval(-2.63);
    A[2][2] = Interval(0.45);
    A[2][3] = Interval(-2.58);
    A[2][4] = Interval(6.10);

    // Inicializar el vector x con intervalos
    IntervalVector x = IntervalVector(5, Interval(0));
    x[0] = Interval(-9.15, 6.97);
    x[1] = Interval(-9.72, 0.92);
    x[2] = Interval(-9.33, 5.51);
    x[3] = Interval(-0.06, 5.53);
    x[4] = Interval(-3.66, 6.73);

    // Ejecutar el caso de prueba
    test_case(A, x, true);
}


void no_solution_test(){
    std::cout << "                                   TESTING: NO SOLUTION TEST 1" << std::endl;

    IntervalMatrix A(3, 5);
    A[0][0] = Interval(-7.31);
    A[0][1] = Interval(6.95);
    A[0][2] = Interval(5.28);
    A[0][3] = Interval(-4.90);
    A[0][4] = Interval(-0.08);

    A[1][0] = Interval(-1.01);
    A[1][1] = Interval(3.03);
    A[1][2] = Interval(5.77);
    A[1][3] = Interval(-8.12);
    A[1][4] = Interval(-9.43);

    A[2][0] = Interval(6.72);
    A[2][1] = Interval(-1.34);
    A[2][2] = Interval(5.25);
    A[2][3] = Interval(-9.96);
    A[2][4] = Interval(-1.09);

    // Inicializar el vector x con intervalos
    IntervalVector x = IntervalVector(5, Interval(0));
    x[0] = Interval(1.565, 2.880);
    x[1] = Interval(0.478, 4.463);
    x[2] = Interval(1.038, 6.032);
    x[3] = Interval(0.048, 3.615);
    x[4] = Interval(1.076, 2.647);

    // Ejecutar el caso de prueba
    test_case(A, x, true);
}

void no_solution_test2(){
    std::cout << "                                   TESTING: NO SOLUTION TEST 2" << std::endl;

    // Matriz A (3x5) de intervalos 
    IntervalMatrix A(3, 5);
    A[0][0] = Interval(-7.31);
    A[0][1] = Interval(6.95);
    A[0][2] = Interval(5.28);
    A[0][3] = Interval(-4.90);
    A[0][4] = Interval(-0.08);

    A[1][0] = Interval(-1.01);
    A[1][1] = Interval(3.03);
    A[1][2] = Interval(5.77);
    A[1][3] = Interval(-8.12);
    A[1][4] = Interval(-9.43);

    A[2][0] = Interval(6.72);
    A[2][1] = Interval(-1.34);
    A[2][2] = Interval(5.25);
    A[2][3] = Interval(-9.96);
    A[2][4] = Interval(-1.09);

    // Inicializar el vector x con intervalos
    IntervalVector x = IntervalVector(5, Interval(0));
    x[0] = Interval(-1.565, -0.30);
    x[1] = Interval(-0.478, 4.463);
    x[2] = Interval(-1.038, 6.032);
    x[3] = Interval(-0.048, 3.615);
    x[4] = Interval(-1.076, 2.647);

    // Ejecutar el caso de prueba
    test_case(A, x, true);
}

// Caso de prueba: UNBOUNDED TEST
void unbounded_test() {
    std::cout << "                                   TESTING: UNBOUNDED TEST" << std::endl;

    // Matriz A (3x5) de intervalos 
    IntervalMatrix A(3, 5);
    A[0][0] = Interval(-7.31);
    A[0][1] = Interval(6.95);
    A[0][2] = Interval(5.28);
    A[0][3] = Interval(-4.90);
    A[0][4] = Interval(-0.08);

    A[1][0] = Interval(-1.01);
    A[1][1] = Interval(3.03);
    A[1][2] = Interval(5.77);
    A[1][3] = Interval(-8.12);
    A[1][4] = Interval(-9.43);

    A[2][0] = Interval(6.72);
    A[2][1] = Interval(-1.34);
    A[2][2] = Interval(5.25);
    A[2][3] = Interval(-9.96);
    A[2][4] = Interval(-1.09);

    // Inicializar el vector x con intervalos
    IntervalVector x = IntervalVector(5, Interval(0));
    x[0] = Interval(-1e20, 2.880);
    x[1] = Interval(-0.478, 4.463);
    x[2] = Interval(-1e20, 6.032);
    x[3] = Interval(-0.048, 3.615);
    x[4] = Interval(-1.076, 2.647);

    // Ejecutar el caso de prueba
    test_case(A, x, true);
}

void ill_conditioned_matrix_test() {
    // Matriz A (5x8) de intervalos con elementos escalados a 1e-6
    IntervalMatrix A(5, 8);
    A[0][0] = Interval(9.99750e-1);
    A[0][1] = Interval(1.00090e0);
    A[0][2] = Interval(1.00046e0);
    A[0][3] = Interval(1.00020e0);
    A[0][4] = Interval(9.99310e-1);
    A[0][5] = Interval(9.99310e-1);
    A[0][6] = Interval(9.99120e-1);
    A[0][7] = Interval(1.00073e0);

    A[1][0] = Interval(1.00020e0);
    A[1][1] = Interval(1.00042e0);
    A[1][2] = Interval(9.99040e-1);
    A[1][3] = Interval(1.00094e0);
    A[1][4] = Interval(1.00066e0);
    A[1][5] = Interval(9.99420e-1);
    A[1][6] = Interval(9.99360e-1);
    A[1][7] = Interval(9.99370e-1);

    A[2][0] = Interval(9.99610e-1);
    A[2][1] = Interval(1.00005e0);
    A[2][2] = Interval(9.99860e-1);
    A[2][3] = Interval(9.99580e-1);
    A[2][4] = Interval(1.00022e0);
    A[2][5] = Interval(9.99280e-1);
    A[2][6] = Interval(9.99580e-1);
    A[2][7] = Interval(9.99730e-1);

    A[3][0] = Interval(9.99910e-1);
    A[3][1] = Interval(1.00057e0);
    A[3][2] = Interval(9.99400e-1);
    A[3][3] = Interval(1.00003e0);
    A[3][4] = Interval(1.00018e0);
    A[3][5] = Interval(9.99090e-1);
    A[3][6] = Interval(1.00022e0);
    A[3][7] = Interval(9.99340e-1);

    A[4][0] = Interval(9.99130e-1);
    A[4][1] = Interval(1.00090e0);
    A[4][2] = Interval(1.00093e0);
    A[4][3] = Interval(1.00062e0);
    A[4][4] = Interval(9.99610e-1);
    A[4][5] = Interval(9.99200e-1);
    A[4][6] = Interval(1.00037e0);
    A[4][7] = Interval(9.99880e-1);

    // Inicializar el vector x con intervalos
    IntervalVector x(8, Interval(0));
    x[0] = Interval(-1.565, 2.880);
    x[1] = Interval(-0.478, 4.463);
    x[2] = Interval(-1.038, 6.032);
    x[3] = Interval(-0.048, 3.615);
    x[4] = Interval(-1.076, 2.647);
    x[5] = Interval(-1.565, 2.880);
    x[6] = Interval(-0.478, 4.463);
    x[7] = Interval(-1.038, 6.032);

    // Llamada al caso de prueba
    test_case(A, x, true);
}


void dimensions_15x20_test() {
    // Matriz A (15x20) de intervalos
    IntervalMatrix A(15, 20);
    A[0][0] = Interval(-3.57);
    A[0][1] = Interval(-5.2);
    A[0][2] = Interval(-1.7);
    A[0][3] = Interval(-2.86);
    A[0][4] = Interval(-9.57);
    A[0][5] = Interval(4.26);
    A[0][6] = Interval(8.96);
    A[0][7] = Interval(8.48);
    A[0][8] = Interval(3.59);
    A[0][9] = Interval(-6.75);
    A[0][10] = Interval(-3.2);
    A[0][11] = Interval(2.74);
    A[0][12] = Interval(-7.45);
    A[0][13] = Interval(6.38);
    A[0][14] = Interval(2.54);
    A[0][15] = Interval(-5.98);
    A[0][16] = Interval(3.65);
    A[0][17] = Interval(-4.87);
    A[0][18] = Interval(8.21);
    A[0][19] = Interval(-1.32);

    // (Añade las demás filas de la matriz A aquí, siguiendo el patrón de arriba...)

    // Inicializar el vector x con intervalos
    IntervalVector x = IntervalVector(20, Interval(0));
    x[0] = Interval(-1.11, 6.71);
    x[1] = Interval(-8.7, 3.65);
    x[2] = Interval(-3.7, 2.92);
    x[3] = Interval(-3.2, 5.1);
    x[4] = Interval(-7.56, 3.15);
    x[5] = Interval(-9.74, 8.38);
    x[6] = Interval(-2.34, 7.37);
    x[7] = Interval(-8.61, 7.17);
    x[8] = Interval(-6.86, 7.69);
    x[9] = Interval(-6.99, 1.1);
    x[10] = Interval(-4.73, 2.86);
    x[11] = Interval(-3.65, 5.27);
    x[12] = Interval(-5.49, 3.85);
    x[13] = Interval(-7.14, 8.13);
    x[14] = Interval(-2.32, 6.48);
    x[15] = Interval(-6.57, 4.31);
    x[16] = Interval(-8.12, 7.34);
    x[17] = Interval(-4.26, 5.68);
    x[18] = Interval(-5.87, 6.91);
    x[19] = Interval(-9.31, 2.76);

    // Ejecutar el caso de prueba
    test_case(A, x, true);
}