//============================================================================
//                                  I B E X
// File        : ibexsolve_dfb.cpp
// Author      : (DFB) ibexsolve con DFB en el hueco de CtcPolytopeHull
//============================================================================
//
// Misma composicion que DefaultSolver (ibex_DefaultSolver.cpp):
//   HC4, ACID(HC4), Newton si el sistema es cuadrado, y
//   CtcFixPoint(CtcCompo(<contractor lineal>(LinearizerXTaylor), HC4)),
// con el mismo bisector (SmearSumRelative si es cuadrado, RoundRobin si no) y
// el mismo buffer (pila: profundidad primero). Lo unico que cambia es el
// contractor lineal: CtcPolytopeHull (--lr=xn, produccion) o CtcDFBHull
// (--lr=dfbhull), los dos sobre el mismo LinearizerXTaylor.
//
// Uso: ibexsolve_dfb archivo.bch [--lr=xn|dfbhull] [-e eps_min] [-t timeout]
//                                [--seed=s] [--bfs]
// Ultima linea: celdas tiempo estado soluciones frontera desconocidas pendientes

#include "ibex.h"
#include "ibex_CtcDFBHull.h"
#include <cstring>
#include <cstdlib>

using namespace std;
using namespace ibex;

int main(int argc, char** argv) {
    if (argc < 2) { cerr << "uso: ibexsolve_dfb archivo.bch [--lr=xn|dfbhull] [-e eps] [-t s] [--seed=s] [--bfs]" << endl; return 1; }
    string lr = "xn"; double eps = 1e-3, timeout = -1, seed = 1.0; bool bfs = false;
    for (int i = 2; i < argc; ++i) {
        string a = argv[i];
        if (a.rfind("--lr=", 0) == 0) lr = a.substr(5);
        else if (a == "-e" && i+1 < argc) eps = atof(argv[++i]);
        else if (a == "-t" && i+1 < argc) timeout = atof(argv[++i]);
        else if (a.rfind("--seed=", 0) == 0) seed = atof(a.substr(7).c_str());
        else if (a == "--bfs") bfs = true;
        else { cerr << "argumento desconocido: " << a << endl; return 1; }
    }
    if (lr != "xn" && lr != "dfbhull") { cerr << "--lr debe ser xn o dfbhull" << endl; return 1; }

    System sys(argv[1]);
    RNG::srand(seed);

    // Sistema cuadrado de igualdades, como en DefaultSolver
    int nb_eq = 0;
    for (int i = 0; i < sys.nb_ctr; i++) if (sys.ctrs[i].op == EQ) nb_eq += sys.ctrs[i].f.image_dim();
    System* eqs = (sys.nb_var == nb_eq) ? new System(sys, System::EQ_ONLY) : NULL;

    CtcHC4 hc4(sys.ctrs, 0.01);
    CtcHC4 hc4acid(sys.ctrs, 0.1, true);
    CtcAcid acid(sys, hc4acid);
    CtcNewton* newton = eqs ? new CtcNewton(eqs->f_ctrs, 5e8, eps, 1.e-4) : NULL;
    LinearizerXTaylor lin(sys);
    Ctc* clin = (lr == "dfbhull") ? (Ctc*) new CtcDFBHull(lin) : (Ctc*) new CtcPolytopeHull(lin);
    CtcHC4 hc4xn(sys.ctrs, 0.01);
    CtcCompo compo_xn(*clin, hc4xn);
    CtcFixPoint fix(compo_xn);

    Array<Ctc> lista(4); int k = 0;
    lista.set_ref(k++, hc4);
    lista.set_ref(k++, acid);
    if (newton) lista.set_ref(k++, *newton);
    lista.set_ref(k++, fix);
    lista.resize(k);
    CtcCompo ctc(lista);

    Bsc* bsc = eqs ? (Bsc*) new SmearSumRelative(*eqs, eps) : (Bsc*) new RoundRobin(eps);
    CellBuffer* buf = bfs ? (CellBuffer*) new CellList() : (CellBuffer*) new CellStack();

    Solver s(sys, ctc, *bsc, *buf, Vector(sys.nb_var, eps), Vector(sys.nb_var, POS_INFINITY));
    if (timeout > 0) s.time_limit = timeout;
    Solver::Status st = s.solve(sys.box);
    const CovSolverData& d = s.get_data();
    s.report();
    cout << s.get_nb_cells() << " " << s.get_time() << " " << (int) st << " "
         << d.nb_solution() << " " << d.nb_boundary() << " " << d.nb_unknown() << " " << d.nb_pending() << endl;
    return 0;
}
