//============================================================================
//                                  I B E X
// File        : ibexsolve_dfb.cpp
// Author      : (DFB) ibexsolve con DFB en el hueco de CtcPolytopeHull
//============================================================================
//
// Misma composicion que DefaultSolver (ibex_DefaultSolver.cpp):
//   HC4, ACID(HC4), Newton si el sistema es cuadrado, y
//   CtcFixPoint(CtcCompo(<contractor lineal>(LinearizerXTaylor), HC4)),
// y el mismo buffer (pila: profundidad primero). Cambia el contractor lineal:
// CtcPolytopeHull (--lr=xn, produccion) o CtcDFBHull (--lr=dfbhull), y el
// linearizador (--lin=art por omision, --lin=xt el de DefaultSolver), y el
// bisector: -b lsmear (omision desde 2026-10-06), -b smearsumrel, el de
// DefaultSolver (SmearSumRelative si es cuadrado, RoundRobin si no), o
// -b lsmeardfb, la regla de LSmear con las sensibilidades de DFB (requiere
// --lr=dfbhull). --lin=compo escribe las filas del afin y de XTaylor.
//
// Uso: ibexsolve_dfb archivo.bch [--lr=xn|dfbhull] [--lin=art|xt|compo] [-b lsmear|smearsumrel|lsmeardfb] [-e eps_min] [-t timeout]
//                                [--seed=s] [--bfs]
// Ultima linea: celdas tiempo estado soluciones frontera desconocidas pendientes

#include "ibex.h"
#include "ibex_CtcDFBHull.h"
#include "affine/ibex_LinearizerAffine2.h"
#include "ibex_LSmearDFB.h"
#include "ibex_LinearizerCompo.h"
#include <cstring>
#include <cstdlib>

using namespace std;
using namespace ibex;

int main(int argc, char** argv) {
    if (argc < 2) { cerr << "uso: ibexsolve_dfb archivo.bch [--lr=xn|dfbhull] [-e eps] [-t s] [--seed=s] [--bfs]" << endl; return 1; }
    string lr = "xn", lin = "art", bis = "lsmear"; bool sols = false; double eps = 1e-3, timeout = -1, seed = 1.0; bool bfs = false;
    for (int i = 2; i < argc; ++i) {
        string a = argv[i];
        if (a.rfind("--lr=", 0) == 0) lr = a.substr(5);
        else if (a == "-b" && i+1 < argc) bis = argv[++i];
        else if (a == "-e" && i+1 < argc) eps = atof(argv[++i]);
        else if (a == "-t" && i+1 < argc) timeout = atof(argv[++i]);
        else if (a.rfind("--lin=", 0) == 0) lin = a.substr(6);
        else if (a.rfind("--seed=", 0) == 0) seed = atof(a.substr(7).c_str());
        else if (a == "--bfs") bfs = true;
        else if (a == "--sols") sols = true;
        else { cerr << "argumento desconocido: " << a << endl; return 1; }
    }
    if (lr != "xn" && lr != "dfbhull") { cerr << "--lr debe ser xn o dfbhull" << endl; return 1; }
    if (lin != "xt" && lin != "art" && lin != "compo") { cerr << "--lin debe ser art, xt o compo" << endl; return 1; }
    if (bis != "lsmear" && bis != "smearsumrel" && bis != "lsmeardfb") { cerr << "-b debe ser lsmear, smearsumrel o lsmeardfb" << endl; return 1; }
    if (bis == "lsmeardfb" && lr != "dfbhull") { cerr << "-b lsmeardfb requiere --lr=dfbhull" << endl; return 1; }

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
    /* --lin=art (omision desde 2026-10-05: el linearizador afin, una fila por
     * lado de cada restriccion; mas rapido con los dos contractores) o --lin=xt
     * (el de DefaultSolver: XTaylor con dos esquinas). Para reproducir al
     * ibexsolve oficial: --lr=xn --lin=xt. */
    Linearizer* linp;
    if (lin == "art") linp = new LinearizerAffine2(sys);
    else if (lin == "xt") linp = new LinearizerXTaylor(sys);
    else linp = new LinearizerCompo(*new LinearizerAffine2(sys), *new LinearizerXTaylor(sys));
    Ctc* clin = (lr == "dfbhull") ? (Ctc*) new CtcDFBHull(*linp) : (Ctc*) new CtcPolytopeHull(*linp);
    if (bis == "lsmeardfb") ((CtcDFBHull*) clin)->activar_sensibilidades();
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

    /* -b lsmear: sobre las igualdades si el sistema es cuadrado, sobre todo el
     * sistema si no. Sin objetivo, cada LP de LSmear minimiza o maximiza una
     * variable al azar. -b smearsumrel reproduce a DefaultSolver. */
    Bsc* bsc;
    if (bis == "lsmear") bsc = new LSmear(eqs ? *eqs : sys, eps);
    else if (bis == "lsmeardfb") bsc = new LSmearDFB(sys, eps);
    else bsc = eqs ? (Bsc*) new SmearSumRelative(*eqs, eps) : (Bsc*) new RoundRobin(eps);
    CellBuffer* buf = bfs ? (CellBuffer*) new CellList() : (CellBuffer*) new CellStack();

    Solver s(sys, ctc, *bsc, *buf, Vector(sys.nb_var, eps), Vector(sys.nb_var, POS_INFINITY));
    if (timeout > 0) s.time_limit = timeout;
    Solver::Status st = s.solve(sys.box);
    const CovSolverData& d = s.get_data();
    s.report();
    if (bis == "lsmeardfb")
        cerr << "[lsmeardfb] llamadas=" << LSmearDFB::n_llamadas << " sin_sens=" << LSmearDFB::n_sin_sens
             << " fuera=" << LSmearDFB::n_fuera << " eligio_dfb=" << LSmearDFB::n_dfb << endl;
    if (sols) cout << d << endl;
    cout << s.get_nb_cells() << " " << s.get_time() << " " << (int) st << " "
         << d.nb_solution() << " " << d.nb_boundary() << " " << d.nb_unknown() << " " << d.nb_pending() << endl;
    if (lr == "xn") delete clin;   /* imprime las sondas de CtcPolytopeHull */
    return 0;
}
