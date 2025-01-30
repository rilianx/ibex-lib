#include <random>
#include <cstdlib>
#include <vector> 
#include "ibex.h"
#include "ibex_CtcDualFeasibleBounding.h"
#include "TestingDualFeasibleBounding.h"

void test_case(IntervalMatrix& A, IntervalVector& x, bool visualize_iters){
    int m = A.nb_rows();
    int n = A.nb_cols();
    IntervalVector xOG = IntervalVector(x);
    std::vector<int> iters;

    for (int k = 0; k < n; ++k) {
        CtcDFB DFB(n, k);
        CtcDFB DFB_upper(n, k, true);

        DFB.init(IntervalMatrix(A));
        DFB.contract(x);

        DFB_upper.init(IntervalMatrix(A));
        DFB_upper.contract(x);
    }
    
    
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "xOG: " << xOG << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "xDFB: " << x << std::endl;
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
    std::cout << "                                   TESTING: ILL-CONDITIONED-MATRIX-TEST" << std::endl;

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
    std::cout << "                                   TESTING: 15X20 DIMENSIONS TEST" << std::endl;

  // Matriz A (15x20) de intervalos
    IntervalMatrix A(15, 20);
    A[0][0] = Interval(-2.96292886);
    A[0][1] = Interval(3.16497243);
    A[0][2] = Interval(5.25512492);
    A[0][3] = Interval(-3.06567412);
    A[0][4] = Interval(-0.773952131);
    A[0][5] = Interval(1.71046441);
    A[0][6] = Interval(-6.52370722);
    A[0][7] = Interval(0.784382043);
    A[0][8] = Interval(-8.63718489);
    A[0][9] = Interval(9.76196200);
    A[0][10] = Interval(2.76863150);
    A[0][11] = Interval(4.34192947);
    A[0][12] = Interval(1.81791509);
    A[0][13] = Interval(8.88620986);
    A[0][14] = Interval(7.69284717);
    A[0][15] = Interval(-5.87094454);
    A[0][16] = Interval(-3.73419425);
    A[0][17] = Interval(-4.53074796);
    A[0][18] = Interval(-5.13155560);
    A[0][19] = Interval(-6.84112827);

    A[1][0] = Interval(7.50698372);
    A[1][1] = Interval(2.65942777);
    A[1][2] = Interval(0.194975342);
    A[1][3] = Interval(7.46691182);
    A[1][4] = Interval(-9.93152881);
    A[1][5] = Interval(-5.34092981);
    A[1][6] = Interval(4.18822810);
    A[1][7] = Interval(-9.20999448);
    A[1][8] = Interval(-7.98060028);
    A[1][9] = Interval(1.25607812);
    A[1][10] = Interval(3.82887859);
    A[1][11] = Interval(-2.60828162);
    A[1][12] = Interval(-0.235515110);
    A[1][13] = Interval(-1.18677816);
    A[1][14] = Interval(-2.35092640);
    A[1][15] = Interval(7.59992020);
    A[1][16] = Interval(7.20333674);
    A[1][17] = Interval(-5.82820795);
    A[1][18] = Interval(5.05968507);
    A[1][19] = Interval(-0.947119568);

    A[2][0] = Interval(-8.33587012);
    A[2][1] = Interval(0.734911570);
    A[2][2] = Interval(4.35577785);
    A[2][3] = Interval(1.94107811);
    A[2][4] = Interval(3.22574873);
    A[2][5] = Interval(-3.19129536);
    A[2][6] = Interval(8.38687758);
    A[2][7] = Interval(-5.43097397);
    A[2][8] = Interval(9.12216834);
    A[2][9] = Interval(1.31358988);
    A[2][10] = Interval(3.72293250);
    A[2][11] = Interval(-6.54344326);
    A[2][12] = Interval(8.49657154);
    A[2][13] = Interval(-8.60585334);
    A[2][14] = Interval(-3.97240766);
    A[2][15] = Interval(-3.27605941);
    A[2][16] = Interval(9.54440624);
    A[2][17] = Interval(-9.75848682);
    A[2][18] = Interval(5.85335765);
    A[2][19] = Interval(-7.83895676);

    A[3][0] = Interval(-3.25589114);
    A[3][1] = Interval(5.58980772);
    A[3][2] = Interval(3.42271315);
    A[3][3] = Interval(0.649510131);
    A[3][4] = Interval(-5.46790895);
    A[3][5] = Interval(6.39404629);
    A[3][6] = Interval(-1.76039534);
    A[3][7] = Interval(2.41389755);
    A[3][8] = Interval(-3.60675089);
    A[3][9] = Interval(3.17050003);
    A[3][10] = Interval(-5.46225469);
    A[3][11] = Interval(5.96680153);
    A[3][12] = Interval(-7.30551401);
    A[3][13] = Interval(-0.446746377);
    A[3][14] = Interval(-4.45955774);
    A[3][15] = Interval(-0.371929158);
    A[3][16] = Interval(9.43244852);
    A[3][17] = Interval(-9.91649136);
    A[3][18] = Interval(-6.91720120);
    A[3][19] = Interval(3.46017912);

    A[4][0] = Interval(-9.23889348);
    A[4][1] = Interval(-1.20591512);
    A[4][2] = Interval(-2.43457251);
    A[4][3] = Interval(-5.94779301);
    A[4][4] = Interval(-5.46563896);
    A[4][5] = Interval(4.21882970);
    A[4][6] = Interval(4.25798728);
    A[4][7] = Interval(-0.696477254);
    A[4][8] = Interval(-0.217176906);
    A[4][9] = Interval(-5.92986713);
    A[4][10] = Interval(-4.27140293);
    A[4][11] = Interval(9.20233750);
    A[4][12] = Interval(3.03109645);
    A[4][13] = Interval(-7.84524343);
    A[4][14] = Interval(7.29927434);
    A[4][15] = Interval(7.37937525);
    A[4][16] = Interval(3.99079651);
    A[4][17] = Interval(8.32552826);
    A[4][18] = Interval(-3.17767506);
    A[4][19] = Interval(-1.67081256);
    A[5][0] = Interval( 9.53017145e+00);
    A[5][1] = Interval(-9.23602844e+00);
    A[5][2] = Interval( 3.76430759e+00);
    A[5][3] = Interval( 4.38021382e+00);
    A[5][4] = Interval(-9.32140722e+00);
    A[5][5] = Interval( 3.55071933e+00);
    A[5][6] = Interval( 5.27240087e+00);
    A[5][7] = Interval(-1.16961986e+00);
    A[5][8] = Interval(-3.02318388e+00);
    A[5][9] = Interval( 1.13668142e+00);
    A[5][10] = Interval(-3.06443623e+00);
    A[5][11] = Interval( 3.15591467e+00);
    A[5][12] = Interval( 3.37728716e+00);
    A[5][13] = Interval( 7.83822353e+00);
    A[5][14] = Interval( 6.13825333e+00);
    A[5][15] = Interval(-1.84730421e+00);
    A[5][16] = Interval( 7.77327673e+00);
    A[5][17] = Interval( 4.32468059e+00);
    A[5][18] = Interval( 2.55227209e+00);
    A[5][19] = Interval(-1.47966281e+00);

    A[6][0] = Interval( 1.88722089e+00);
    A[6][1] = Interval( 4.21792949e+00);
    A[6][2] = Interval(-7.01523889e+00);
    A[6][3] = Interval( 5.29471912e+00);
    A[6][4] = Interval(-9.83440967e+00);
    A[6][5] = Interval(-1.01334474e+00);
    A[6][6] = Interval( 3.11258273e+00);
    A[6][7] = Interval(-9.38126499e+00);
    A[6][8] = Interval( 5.04298485e+00);
    A[6][9] = Interval(-7.48076672e+00);
    A[6][10] = Interval(-1.07539097e+00);
    A[6][11] = Interval( 4.55784043e+00);
    A[6][12] = Interval(-5.27633107e+00);
    A[6][13] = Interval( 2.36518387e+00);
    A[6][14] = Interval(-3.50411675e+00);
    A[6][15] = Interval( 6.46363987e+00);
    A[6][16] = Interval( 9.25822375e+00);
    A[6][17] = Interval( 7.62203882e+00);
    A[6][18] = Interval( 1.19147015e+00);
    A[6][19] = Interval( 7.55728711e+00);

    A[7][0] = Interval( 8.26878739e+00);
    A[7][1] = Interval(-7.62848114e+00);
    A[7][2] = Interval( 1.49270743e+00);
    A[7][3] = Interval( 3.71529571e+00);
    A[7][4] = Interval( 2.59504271e+00);
    A[7][5] = Interval(-2.40792867e+00);
    A[7][6] = Interval( 2.34701762e+00);
    A[7][7] = Interval(-6.68945644e+00);
    A[7][8] = Interval( 7.72475743e+00);
    A[7][9] = Interval( 8.45213371e+00);
    A[7][10] = Interval( 7.05522976e+00);
    A[7][11] = Interval( 4.34538797e+00);
    A[7][12] = Interval( 4.16769602e+00);
    A[7][13] = Interval(-6.52870361e+00);
    A[7][14] = Interval( 3.46684799e+00);
    A[7][15] = Interval( 7.87745791e+00);
    A[7][16] = Interval( 4.30383281e+00);
    A[7][17] = Interval(-3.91171625e+00);
    A[7][18] = Interval(-3.05442661e+00);
    A[7][19] = Interval( 6.52492666e+00);
    A[8][0] = Interval( 4.53561669e+00);
    A[8][1] = Interval(-5.85061848e+00);
    A[8][2] = Interval(-1.02928562e+00);
    A[8][3] = Interval(-5.32525871e+00);
    A[8][4] = Interval( 3.85813972e+00);
    A[8][5] = Interval( 4.17958093e+00);
    A[8][6] = Interval(-2.21926435e+00);
    A[8][7] = Interval( 6.25654043e+00);
    A[8][8] = Interval( 1.19605358e+00);
    A[8][9] = Interval(-1.73851852e+00);
    A[8][10] = Interval( 5.48903167e+00);
    A[8][11] = Interval( 5.26942050e+00);
    A[8][12] = Interval(-3.65430356e+00);
    A[8][13] = Interval( 1.99382882e+00);
    A[8][14] = Interval( 1.29270065e+00);
    A[8][15] = Interval( 7.22017598e+00);
    A[8][16] = Interval(-9.45115998e+00);
    A[8][17] = Interval( 9.23144094e+00);
    A[8][18] = Interval( 6.09157622e+00);
    A[8][19] = Interval(-1.27938786e+00);

    A[9][0] = Interval(-9.01357490e+00);
    A[9][1] = Interval( 9.07196574e+00);
    A[9][2] = Interval( 2.29093348e+00);
    A[9][3] = Interval( 4.61739978e+00);
    A[9][4] = Interval(-9.26517322e+00);
    A[9][5] = Interval( 2.07811485e+00);
    A[9][6] = Interval( 1.89865738e+00);
    A[9][7] = Interval( 1.79731245e+00);
    A[9][8] = Interval( 4.81840818e+00);
    A[9][9] = Interval(-5.45626314e+00);
    A[9][10] = Interval( 1.83280183e+00);
    A[9][11] = Interval( 1.69668019e+00);
    A[9][12] = Interval( 3.94691328e+00);
    A[9][13] = Interval( 1.22610116e+00);
    A[9][14] = Interval( 9.86614096e+00);
    A[9][15] = Interval( 7.48662523e+00);
    A[9][16] = Interval( 2.34657479e+00);
    A[9][17] = Interval( 1.47297897e+00);
    A[9][18] = Interval(-2.21942096e+00);
    A[9][19] = Interval( 4.58690427e+00);

    A[10][0] = Interval( 2.34792656e+00);
    A[10][1] = Interval(-5.19759877e+00);
    A[10][2] = Interval(-9.22941972e+00);
    A[10][3] = Interval( 4.26359771e+00);
    A[10][4] = Interval( 3.17272072e+00);
    A[10][5] = Interval( 9.07936849e+00);
    A[10][6] = Interval(-9.82669592e+00);
    A[10][7] = Interval(-9.38078151e+00);
    A[10][8] = Interval(-3.51823461e+00);
    A[10][9] = Interval(-1.67942304e+00);
    A[10][10] = Interval(-9.85464958e+00);
    A[10][11] = Interval( 2.42801945e+00);
    A[10][12] = Interval( 2.55432597e+00);
    A[10][13] = Interval( 7.52305522e+00);
    A[10][14] = Interval( 7.45401424e+00);
    A[10][15] = Interval( 1.97343589e+00);
    A[10][16] = Interval(-7.25548215e+00);
    A[10][17] = Interval(-9.97884009e+00);
    A[10][18] = Interval( 1.09357487e+00);
    A[10][19] = Interval( 1.84079152e+00);

    A[11][0] = Interval( 6.30994283e+00);
    A[11][1] = Interval(-4.14903796e+00);
    A[11][2] = Interval( 6.53768581e+00);
    A[11][3] = Interval( 2.49471682e+00);
    A[11][4] = Interval(-6.71677488e+00);
    A[11][5] = Interval( 6.03019429e+00);
    A[11][6] = Interval( 7.22728758e+00);
    A[11][7] = Interval( 8.16350784e+00);
    A[11][8] = Interval( 9.89303669e+00);
    A[11][9] = Interval( 3.41557304e+00);
    A[11][10] = Interval( 5.34456716e+00);
    A[11][11] = Interval(-7.09845306e+00);
    A[11][12] = Interval( 3.47811339e+00);
    A[11][13] = Interval( 9.33501967e+00);
    A[11][14] = Interval( 9.17918928e+00);
    A[11][15] = Interval(-2.27410442e+00);
    A[11][16] = Interval(-3.18427857e+00);
    A[11][17] = Interval(-7.11579752e+00);
    A[11][18] = Interval( 4.00455239e+00);
    A[11][19] = Interval( 6.61665817e+00);

    A[12][0] = Interval( 4.44893849e+00);
    A[12][1] = Interval( 3.22296558e+00);
    A[12][2] = Interval( 6.37345702e+00);
    A[12][3] = Interval( 3.28169712e+00);
    A[12][4] = Interval( 4.39313496e+00);
    A[12][5] = Interval(-2.36633423e+00);
    A[12][6] = Interval( 1.49116022e+00);
    A[12][7] = Interval( 7.09702355e+00);
    A[12][8] = Interval( 2.72623376e+00);
    A[12][9] = Interval( 2.47334338e+00);
    A[12][10] = Interval(-1.45507482e+00);
    A[12][11] = Interval(-9.44112360e+00);
    A[12][12] = Interval( 1.53562458e+00);
    A[12][13] = Interval(-1.70422635e+00);
    A[12][14] = Interval(-5.63440290e+00);
    A[12][15] = Interval( 8.35000814e+00);
    A[12][16] = Interval( 4.74885866e+00);
    A[12][17] = Interval( 6.45926884e+00);
    A[12][18] = Interval(-3.69502777e+00);
    A[12][19] = Interval( 1.18139198e+00);

    A[13][0] = Interval( 7.66052972e+00);
    A[13][1] = Interval(-2.83811783e+00);
    A[13][2] = Interval( 9.44434873e+00);
    A[13][3] = Interval( 1.01058188e+00);
    A[13][4] = Interval( 4.07177567e+00);
    A[13][5] = Interval( 3.43257084e+00);
    A[13][6] = Interval( 4.49158986e+00);
    A[13][7] = Interval( 5.25824793e+00);
    A[13][8] = Interval( 2.71721144e+00);
    A[13][9] = Interval( 3.85547622e+00);
    A[13][10] = Interval(-2.55951882e+00);
    A[13][11] = Interval( 9.75466337e+00);
    A[13][12] = Interval( 1.84760367e+00);
    A[13][13] = Interval(-2.15323233e+00);
    A[13][14] = Interval( 9.36832894e+00);
    A[13][15] = Interval( 9.06361186e+00);
    A[13][16] = Interval( 3.75442138e+00);
    A[13][17] = Interval( 5.49045521e+00);
    A[13][18] = Interval( 2.58830339e+00);
    A[13][19] = Interval( 4.05656485e+00);

    A[14][0] = Interval(-5.60015447e+00);
    A[14][1] = Interval(-9.30790735e+00);
    A[14][2] = Interval( 3.69344167e+00);
    A[14][3] = Interval( 2.65001028e+00);
    A[14][4] = Interval( 3.72440643e+00);
    A[14][5] = Interval( 1.20737692e+00);
    A[14][6] = Interval(-2.70047056e+00);
    A[14][7] = Interval( 1.02197091e+00);
    A[14][8] = Interval( 2.83953375e+00);
    A[14][9] = Interval( 2.06144193e+00);
    A[14][10] = Interval(-4.73214576e+00);
    A[14][11] = Interval( 1.27633696e+00);
    A[14][12] = Interval( 3.07287417e+00);
    A[14][13] = Interval(-6.07094390e+00);
    A[14][14] = Interval( 1.09177381e+00);
    A[14][15] = Interval( 6.76290271e+00);
    A[14][16] = Interval(-3.51907432e+00);
    A[14][17] = Interval(-1.78945858e+00);
    A[14][18] = Interval( 4.74151625e+00);
    A[14][19] = Interval( 7.10645780e+00);


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

pair<IntervalMatrix, IntervalVector> generate_test_case(int m, int n){
    IntervalMatrix A = IntervalMatrix(m, n, Interval(0));
    IntervalVector x = IntervalVector(n, Interval(0));

    // Crear el generador y la distribución con valores por defecto
    mt19937 gen(random_device{}()); // Generador con semilla automática
    uniform_real_distribution<> dis(0.5, 10.0); // Distribución [0.5, 10.0]
    mt19937 gen2(random_device{}()); // Generador con semilla automática
    uniform_real_distribution<> dis2(-10.0, 10.0); // Distribución [0.5, 10.0]
    // Llenar la matriz A con intervalos aleatorios
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            double value = dis2(gen2);
            A[j][i] = Interval(value);
        }
    }

    // Llenar el vector x con intervalos aleatorios
    for (int i = 0; i < n; ++i) {
        double lower = dis(gen);
        double upper = dis(gen);
        x[i] = Interval(-lower, upper);
    }

    return {A, x};
}

void random_test(int m, int n, int quant_cases){
    if (m >= n){
        throw invalid_argument("Number of restrictions are greater or equal than the number of variables. Which is not allowed..");
    }

    cout << "Random Tests " << endl;
    cout << "Dimensions: " << "m = " << m << " n = " << n << endl;

    for (int i = 0; i < quant_cases; ++i){
        IntervalMatrix A = IntervalMatrix(m, n, Interval(0));
        IntervalVector x = IntervalVector(n, Interval(0));
        cout << "Caso de Prueba " << (i + 1) << endl;
        tie(A, x) = generate_test_case(m, n);
        cout << "A = " << A  << "\n" << endl;

        test_case(A, x, true);
    }
}

/*
pair<CtcDFB, int> incremental_test_case(IntervalMatrix& A, IntervalVector& x, int k, CtcDFB* DFB){
    int m = A.nb_rows();
    int n = x.size();
    IntervalVector xOG = IntervalVector(x);
    std::vector<int> iters;
    if (DFB == nullptr){
        DFB = new CtcDFB(A, k);
    }
    int iter = DFB->totalContraction(x);

    return {*DFB, iter};
}


tuple<vector<int>, vector<int>, vector<double>> special_test_contraction_box(int m, int n, int k,  int percentage){
    IntervalMatrix A = IntervalMatrix(m, n);
    IntervalVector x = IntervalVector(n);    

    tie(A, x) = generate_test_case(m, n);
    cout << "A = " << A  << "\n" << endl;

    IntervalMatrix AOG = IntervalMatrix(A);
    IntervalVector xOG = IntervalVector(x);

    vector<pair<CtcDFB, int>> DFBsOG;
    for (int i = 0; i < n; ++i){
        // Llamar a la función y almacenar el resultado directamente en un pair
        pair<CtcDFB, int> result = incremental_test_case(A, x, i);
        
        // Acceder a los elementos del pair
        CtcDFB DFB = result.first;  // Objeto DFB
        int iter = result.second;  // Valor de iteraciones
        DFBsOG.push_back({DFB, iter});
    }

    cout << "ITERACIONES: " << endl;
    for (const auto& pair : DFBsOG) {
        CtcDFB DFBOGAux = pair.first;
        cout << "Para k = : " << DFBOGAux.getKIndex();
        cout << "~ Iteraciones: " << pair.second << endl;
    }


      
    cout << "AProcessed = " << A  << "\n" << endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "xOG: " << xOG << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "xDFB: " << x << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;


    if (k > n){
        k = n; 
    }

    if (percentage > 100){
        percentage = 100;
    }

    if (k <= 0){
        k = 1 + rand() % n;
    }

    if (percentage <= 0){
        k = 1 + rand() % 100;
    }
    double weight =  static_cast<double>(percentage) / 100;

    cout << "Número de variables a modificar de la caja: " << k <<" de un total de " << n << " variables." << endl;
    cout << "Porcentaje a modificar de las " << k << " variables: " <<  percentage << "%" << endl;
    cout << "Peso a obtener de las " << k << " variables: "<<  weight << endl;

    for (int i = 0; i < k; ++i) {
        int index = rand() % n;
        x[index] = weight * x[index];
    }


    xOG = IntervalVector(x);
    IntervalVector xModified = IntervalVector(x);
    
    cout << "Special Test Case ~ A modified contracted vector x and ORIGINAL matrix A: \n" << endl;
    test_case(AOG, xModified, true);




    cout << "Special Test Case ~ A modified contracted vector x and PROCESSED matrix A: \n" << endl;

    vector<pair<CtcDFB, int>> DFBsFinal;
    for (int i = 0; i < n; ++i){
        CtcDFB DFBTemp = DFBsOG[i].first;
        pair<CtcDFB, int> result = incremental_test_case(A, x, i, &DFBTemp);
        
        // Acceder a los elementos del pair
        CtcDFB DFB = result.first;  // Objeto DFB
        int iter = result.second;  // Valor de iteraciones
        DFBsFinal.push_back({DFB, iter});
    }

    cout << "ITERACIONES: " << endl;
    for (const auto& pair : DFBsFinal) {
        CtcDFB DFBAux = pair.first;
        cout << "Para k = : " << DFBAux.getKIndex();
        cout << " ~ Iteraciones: " << pair.second << endl;
    }

    vector<int> itersAOG;
    vector<int> itersAProcessed;

    cout << "ITERS DIFF: " << endl;
    for (const auto& finalPair : DFBsFinal) {
        CtcDFB DFBAuxFinal = finalPair.first;
        int finalValue = finalPair.second;
        int kIndexFinal = DFBAuxFinal.getKIndex();

        // Buscar el kIndex correspondiente en DFBsOG
        for (const auto& ogPair : DFBsOG) {
            CtcDFB DFBAuxOG = ogPair.first;
            int ogValue = ogPair.second;
            int kIndexOG = DFBAuxOG.getKIndex();

            if (kIndexFinal == kIndexOG) {
                // Resta de los valores
                int iterAOG = ogValue;
                int iterAProcessed = finalValue;
                // Mostrar resultados
                cout << "Para k = " << kIndexFinal << " ~ Iters AOG: " << iterAOG << endl;
                cout << "Para k = " << kIndexFinal << " ~ Iters AProcessed: " << iterAProcessed << endl;
                itersAOG.push_back(iterAOG);
                itersAProcessed.push_back(iterAProcessed);
                break;
            }
        }
    }


    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "xOG: " << xOG << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "xDFB: " << x << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;

    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
    vector<double> errores;
    for (int i = 0; i < n; ++i){
        double error = abs(xModified[i].lb() - x[i].lb()) + abs(xModified[i].ub() - x[i].ub());
        errores.push_back(error);
    }
    
    cout << "Vector de errores: " << endl;
    for (int i = 0; i < n; ++i){
        cout << errores[i] << " ";
    }
    cout << "\n" << endl;

    return {itersAOG, itersAProcessed, errores};
}



void special_test(){
    int n = 0;
    int m = 0;
    int quant_cases = 0;
    
    cout << "Ingrese la cantidad de casos de prueba: ";
    std::cin >> quant_cases;


    while (n <= 0 || m <= 0){
        cout << "Ingrese el valor de m: ";
        std::cin >> m;

        cout << "Ingrese el valor de n: ";
        std::cin >> n;
    }

    vector<int> itersAOG;
    vector<int> itersAProcessed;
    vector<double> errors;

    vector<tuple<vector<int>, vector<int>, vector<double>>> testInfo;
    for (int i = 0; i < quant_cases; ++i){
        tie(itersAOG, itersAProcessed, errors) = special_test_contraction_box(m, n);

        testInfo.push_back({itersAOG, itersAProcessed, errors});
    }

    int avgItersAOG = 0;
    int avgItersAProcessed = 0;
    double avgError = 0.0;
    for (int i = 0; i < quant_cases; ++i){
        itersAOG = get<0>(testInfo[i]);
        itersAProcessed = get<1>(testInfo[i]);
        errors = get<2>(testInfo[i]);
        for (int k = 0; k < n; ++k){
            avgItersAOG += itersAOG[k];
            avgItersAProcessed += itersAProcessed[k];
            avgError += errors[k];
        }
    }

    avgItersAOG /= (quant_cases);
    avgItersAProcessed /= (quant_cases);
    avgError /= (quant_cases);

    cout << "Average Iters A OG " << avgItersAOG << endl;
    cout << "Average Iters A PROCESSED " << avgItersAProcessed << endl;
    cout << "Average Total Error " << avgError << endl;
}*/
