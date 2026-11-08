//============================================================================
// sc_solve: system solver with CtcSC in the slot of CtcPolytopeHull.
//
// Same composition as DefaultSolver and as examples/dfb/ibexsolve_dfb
// (--lr=dfbhull), so that both can be compared:
//     HC4, ACID, interval Newton (square systems), FixPoint(CtcSC ; HC4),
//     LSmear on the equations, depth-first search.
//
// Usage: sc_solve file.bch [--lin=art|xt|compo] [-e eps] [-t timeout] [--seed=s] [--bfs] [--sols]
// --lin=art (default): affine relaxation; dual warm start, tolerance 1e-9.
// --lin=compo: affine + XTaylor rows; ADAPTIVE warm start, tolerance 1e-6.
// Last line: cells time status solutions boundary unknown pending
//============================================================================

#include "ibex.h"
#include "ibex_CtcSC.h"
#include "ibex_LinearizerAffine2.h"
#include <cstdlib>

using namespace std;
using namespace ibex;

int main(int argc, char** argv) {
	if (argc < 2) { cerr << "usage: sc_solve file.bch [--lin=art|xt|compo] [-e eps] [-t s] [--seed=s] [--bfs] [--sols]" << endl; return 1; }
	string lin = "art"; bool sols = false, bfs = false; double eps = 1e-3, timeout = -1, seed = 1.0;
	for (int i = 2; i < argc; ++i) {
		string a = argv[i];
		if (a == "-e" && i+1 < argc) eps = atof(argv[++i]);
		else if (a == "-t" && i+1 < argc) timeout = atof(argv[++i]);
		else if (a.rfind("--lin=", 0) == 0) lin = a.substr(6);
		else if (a.rfind("--seed=", 0) == 0) seed = atof(a.substr(7).c_str());
		else if (a == "--lr=sc" || a == "--lr=dfbhull") { }
		else if (a == "--bfs") bfs = true;
		else if (a == "--sols") sols = true;
		else { cerr << "unknown argument: " << a << endl; return 1; }
	}
	if (lin != "xt" && lin != "art" && lin != "compo") { cerr << "--lin must be art, xt or compo" << endl; return 1; }

	System sys(argv[1]);
	RNG::srand(seed);

	// square system of equations, as in DefaultSolver
	int nb_eq = 0;
	for (int i = 0; i < sys.nb_ctr; i++) if (sys.ctrs[i].op == EQ) nb_eq += sys.ctrs[i].f.image_dim();
	System* eqs = (sys.nb_var == nb_eq) ? new System(sys, System::EQ_ONLY) : NULL;

	CtcHC4 hc4(sys.ctrs, 0.01);
	CtcHC4 hc4acid(sys.ctrs, 0.1, true);
	CtcAcid acid(sys, hc4acid);
	CtcNewton* newton = eqs ? new CtcNewton(eqs->f_ctrs, 5e8, eps, 1.e-4) : NULL;

	Linearizer* linp;
	if (lin == "art") linp = new LinearizerAffine2(sys);
	else if (lin == "xt") linp = new LinearizerXTaylor(sys);
	else linp = new LinearizerCompo(*new LinearizerAffine2(sys), *new LinearizerXTaylor(sys));
	CtcSC* sc = (lin == "compo") ? new CtcSC(*linp, CtcSC::ADAPTIVE, 1e-6)
	                             : new CtcSC(*linp, CtcSC::DUAL, 1e-9);
	CtcHC4 hc4xn(sys.ctrs, 0.01);
	CtcCompo compo_xn(*sc, hc4xn);
	CtcFixPoint fix(compo_xn);

	Array<Ctc> list(4); int k = 0;
	list.set_ref(k++, hc4);
	list.set_ref(k++, acid);
	if (newton) list.set_ref(k++, *newton);
	list.set_ref(k++, fix);
	list.resize(k);
	CtcCompo ctc(list);

	Bsc* bsc = new LSmear(eqs ? *eqs : sys, eps);
	CellBuffer* buf = bfs ? (CellBuffer*) new CellList() : (CellBuffer*) new CellStack();

	Solver s(sys, ctc, *bsc, *buf, Vector(sys.nb_var, eps), Vector(sys.nb_var, POS_INFINITY));
	if (timeout > 0) s.time_limit = timeout;
	Solver::Status st = s.solve(sys.box);
	const CovSolverData& d = s.get_data();
	s.report();
	cerr << "[sc] passes=" << sc->nb_passes() << " LPs=" << sc->nb_lps() << " pivots=" << sc->nb_pivots() << endl;
	if (sols) cout << d << endl;
	cout << s.get_nb_cells() << " " << s.get_time() << " " << (int) st << " "
	     << d.nb_solution() << " " << d.nb_boundary() << " " << d.nb_unknown() << " " << d.nb_pending() << endl;
	return 0;
}
