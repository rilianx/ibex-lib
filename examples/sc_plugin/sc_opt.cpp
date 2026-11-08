//============================================================================
// sc_opt: optimizer with CtcSC in the slot of CtcPolytopeHull.
//
// Same strategy as examples/dfb/ibexopt (Optimizer05Config) with
// --lr=dfbhull, so that both can be compared: HC4 + ACID (or HC4, or
// HC4 + 3BCID), then FixPoint(CtcSC ; HC4) with ratio 0.2, KKT for
// unconstrained problems, LSmear (lsmearmg) bisection, beam search,
// LoupFinderDefault.
//
// Usage: sc_opt file.bch [--filtering=acidhc4|hc4|3bcidhc4] [--lin=art|compo|xn]
//                        [-b bisector] [-t timeout] [--seed=s] [--eps_x=..] [--eps=..]
// --lin=art (default): affine relaxation. CtcSC with the PRIMAL warm start
// and tolerance 1e-9.
// Last line: cells time loup loup
//============================================================================

#include "ibex.h"
#include "parse_args.h"
#include "ibex_CtcSC.h"
#include "ibex_LinearizerAffine2.h"
#include "ibex_OptimMemory.h"

#include <sstream>

using namespace std;
using namespace ibex;

namespace {

class SCOptimizerConfig : public OptimizerConfig, protected Memory {
public:
	SCOptimizerConfig(const string& filename, const string& filtering, const string& lin,
	                  const string& bisection, int beamsize, double prec, double goalprec,
	                  double timelimit, int seed) :
		filtering(filtering), lin(lin), bisection(bisection), beamsize(beamsize) {
		sys = &rec(new System(filename.c_str()));
		RNG::srand(seed);
		double eqeps = 1e-8;
		ext_sys = &rec(new ExtendedSystem(*sys, eqeps));
		norm_sys = &rec(new NormalizedSystem(*sys, eqeps));
		set_eps_x(Vector(sys->nb_var, prec));
		set_rel_eps_f(goalprec);
		set_abs_eps_f(goalprec);
		set_timeout(timelimit);
	}

	System* sys;
	CtcSC* sc = NULL;

protected:
	virtual unsigned int nb_var() { return sys->nb_var; }

	virtual int goal_var() { return ext_sys->goal_var(); }

	Linearizer* get_linear_relax() {
		if (lin == "art")
			return &rec(new LinearizerAffine2(*ext_sys));
		else if (lin == "compo")
			return &rec(new LinearizerCompo(rec(new LinearizerAffine2(*ext_sys)),
			                                rec(new LinearizerXTaylor(*ext_sys))));
		else
			return &rec(new LinearizerXTaylor(*ext_sys));
	}

	virtual Ctc& get_ctc() {
		CtcHC4* hc4 = &rec(new CtcHC4(ext_sys->ctrs, 0.01, true));
		CtcHC4* hc44cid = &rec(new CtcHC4(ext_sys->ctrs, 0.1, true));
		CtcHC4* hc44xn = &rec(new CtcHC4(ext_sys->ctrs, 0.01, false));

		Ctc* ctc;
		if (filtering == "hc4")
			ctc = hc4;
		else if (filtering == "acidhc4")
			ctc = &rec(new CtcCompo(*hc4, rec(new CtcAcid(*ext_sys, *hc44cid, true))));
		else if (filtering == "3bcidhc4")
			ctc = &rec(new CtcCompo(*hc4, rec(new Ctc3BCid(*hc44cid))));
		else {
			stringstream ss;
			ss << "[sc_opt] " << filtering << " is not an implemented contraction mode";
			ibex_error(ss.str().c_str());
		}

		Linearizer* lr = get_linear_relax();
		sc = &rec(new CtcSC(*lr, CtcSC::PRIMAL, 1e-9));
		Ctc* cxn;
		if (sys->nb_ctr == 0)
			cxn = sc;
		else
			cxn = &rec(new CtcFixPoint(rec(new CtcCompo(*sc, *hc44xn)), 0.2));

		Ctc* ctcxn = &rec(new CtcCompo(*ctc, *cxn));
		if (sys->nb_ctr == 0)
			ctcxn = &rec(new CtcCompo(*ctcxn, rec(new CtcKuhnTucker(*norm_sys, true))));
		return *ctcxn;
	}

	virtual Bsc& get_bsc() {
		const Vector& eps_x = get_eps_x();
		Vector prec(ext_sys->nb_var);
		if (eps_x.size() == 1)
			ext_sys->write_ext_vec(Vector(ext_sys->nb_var, eps_x[0]), prec);
		else
			ext_sys->write_ext_vec(eps_x, prec);
		prec[ext_sys->goal_var()] = OptimizerConfig::default_eps_x;
		OptimLargestFirst& lf = rec(new OptimLargestFirst(ext_sys->goal_var(), true, prec, 0.5));
		if (bisection == "lsmearmg")
			return rec(new LSmear(*ext_sys, prec, lf));
		else if (bisection == "lsmear")
			return rec(new LSmear(*ext_sys, prec, lf, LSMEAR));
		else if (bisection == "smearsumrel")
			return rec(new SmearSumRelative(*ext_sys, prec, lf));
		else if (bisection == "largestfirst")
			return lf;
		else if (bisection == "roundrobin")
			return rec(new RoundRobin(prec, 0.5));
		stringstream ss;
		ss << "[sc_opt] " << bisection << " is not an implemented bisection mode";
		ibex_error(ss.str().c_str());
		return lf;
	}

	virtual LoupFinder& get_loup_finder() {
		return rec(new LoupFinderDefault(*norm_sys, true));
	}

	virtual CellBufferOptim& get_cell_buffer() {
		CellHeap* future = &rec(new CellHeap(*ext_sys));
		CellHeap* current = &rec(new CellHeap(*ext_sys));
		return rec(new CellBeamSearch(*current, *future, *ext_sys, beamsize));
	}

	NormalizedSystem* norm_sys;
	ExtendedSystem* ext_sys;
	string filtering, lin, bisection;
	int beamsize;
};

} // anonymous namespace

int main(int argc, char** argv) {
	args::ArgumentParser parser("sc_opt", "Ibex optimizer with the CtcSC contractor");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	args::ValueFlag<string> _filtering(parser, "filtering", "hc4|acidhc4*|3bcidhc4", {"filtering"});
	args::ValueFlag<string> _lin(parser, "lin", "linear relaxation: art*|compo|xn", {"lin"});
	args::ValueFlag<string> _lr(parser, "lr", "ignored (sc|dfbhull), for compatibility", {"lr"});
	args::ValueFlag<string> _bisection(parser, "bisection", "lsmearmg*|lsmear|smearsumrel|largestfirst|roundrobin", {'b'});
	args::ValueFlag<double> _prec(parser, "prec", "eps_x (default 1e-7)", {"eps_x"});
	args::ValueFlag<double> _goalPrec(parser, "goal precision", "eps (default 1e-6)", {"eps"});
	args::ValueFlag<double> _timeLimit(parser, "time limit", "time limit (default 1000)", {'t'});
	args::ValueFlag<int> _randomSeed(parser, "random seed", "random seed (default 42)", {"seed"});
	args::Positional<string> filename(parser, "filename", "The name of the MINIBEX file.");

	try {
		parser.ParseCLI(argc, argv);
	} catch (const args::Help&) {
		cout << parser;
		return 0;
	} catch (const args::ParseError& e) {
		cerr << e.what() << endl;
		return 1;
	} catch (const args::ValidationError& e) {
		cerr << e.what() << endl;
		return 1;
	}
	if (!filename) { cerr << "Filename is required." << endl; return 1; }

	SCOptimizerConfig config(filename.Get(),
		_filtering ? _filtering.Get() : "acidhc4",
		_lin ? _lin.Get() : "art",
		_bisection ? _bisection.Get() : "lsmearmg",
		1,
		_prec ? _prec.Get() : 1e-7,
		_goalPrec ? _goalPrec.Get() : 1e-6,
		_timeLimit ? _timeLimit.Get() : 1000,
		_randomSeed ? _randomSeed.Get() : 42);

	Optimizer o(config);
	cout.precision(12);
	o.optimize(config.sys->box);
	o.report();
	if (config.sc)
		cerr << "[sc] passes=" << config.sc->nb_passes() << " LPs=" << config.sc->nb_lps()
		     << " pivots=" << config.sc->nb_pivots() << endl;
	cout << o.get_nb_cells() << " " << o.get_time() << " " << o.get_loup() << " " << o.get_loup() << endl;
	return 0;
}
