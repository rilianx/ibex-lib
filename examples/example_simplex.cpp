#include "ibex.h"
#include "affine_arithm.h"
#include "ibex_AffineEval.h"
#include <iostream>
#include <vector>

using namespace std;
using namespace ibex;

#include <iostream>
#include <vector>
#include <cstdlib> // Para std::atoi y std::atof

int main(int argc, char *argv[]) {
    // Verificar que se han proporcionado suficientes argumentos
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " n m timeout matrix_elements... b1 b2 ... bm" << std::endl;
        return 1;
    }

    // Leer los parámetros desde la línea de comandos
    int nb_vars = std::atoi(argv[1]);
    int nb_ctrs = std::atoi(argv[2]);
    IntervalVector bounds(nb_vars, Interval(std::atof(argv[3]), std::atof(argv[4])));
    double timeout = std::atof(argv[5]);

    int matrix_size = nb_vars * nb_ctrs;

    // Leer los elementos de la matriz
    std::vector<double> matrix(matrix_size);
    for (int i = 0; i < matrix_size; ++i) {
        matrix[i] = std::atof(argv[6 + i]);
    }


    double tolerance = 1e-6;
    int max_iter = 1000;

    LPSolver lp_solver(nb_vars, LPSolver::Mode::NotCertified, tolerance, timeout, max_iter);


    Matrix rows(nb_ctrs, nb_vars, matrix.data());
    cout << "A:" << rows << endl;
    Vector rhs(nb_ctrs);
    for (int i = 0; i < nb_ctrs; ++i) {
        rhs[i] = std::atof(argv[6 + matrix_size + i]);
    }

    lp_solver.add_constraints(rows, CmpOp::LEQ, rhs);

    Vector obj(nb_vars, 0.0);

    lp_solver.set_bounds(bounds);

    IntervalVector new_box(nb_vars);
    for (int i = 0; i < nb_vars; ++i) {
        if (i>0) obj[i-1] = 0;
        obj[i] = 1;
        lp_solver.set_cost(obj);

        //Status minimize();
        LPSolver::Status status = lp_solver.minimize();
        cout << "status:" << status << endl;
        cout << "minimum:" << lp_solver.minimum() << endl;
        cout << "iterations:" << lp_solver.mysoplex->numIterations() << endl;

        new_box[i] = lp_solver.minimum();

        obj[i] = -1;
        lp_solver.set_cost(obj);
        status = lp_solver.minimize();
        cout << "status:" << status << endl;
        cout << "minimum:" << -lp_solver.minimum() << endl;
        cout << "iterations:" << lp_solver.mysoplex->numIterations() << endl;

        new_box[i] = Interval(new_box[i].lb(), -lp_solver.minimum().ub());

    }


    cout << "new_box:" << new_box << endl;

    lp_solver.write_to_file("output.txt");
}