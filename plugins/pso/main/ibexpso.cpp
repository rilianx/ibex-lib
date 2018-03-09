/*
 * optimsolver.cpp
 *
 *  Created on: 21 dic. 2017
 *      Author: iaraya
 */

#include "ibex.h"
#include "ibex_System.h"
#include "args.hxx"

using namespace std;
using namespace ibex;

double get_obj_abs_prec(double uplo, double loup);
double get_obj_rel_prec(double uplo, double loup);

int main(int argc, char** argv){

	try {

	args::ArgumentParser parser("********* IbexOpt (defaultoptimizer) *********.", "Solve a Minibex file.");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	args::ValueFlag<std::string> _filtering(parser, "string", "the filtering method", {'f', "filt"});
	args::ValueFlag<std::string> _linear_relax(parser, "string", "the linear relaxation method", {"linear_relax"});
	args::ValueFlag<std::string> _bisector(parser, "string", "the bisection method", {'b', "bis"});
	args::ValueFlag<std::string> _strategy(parser, "string", "the search strategy (ibex_pso, ibex_opt, ibex_rootpso)", {'s', "search"});
	args::ValueFlag<double> _eps_x(parser, "float", "eps_x (the precision of the boxes)", {"exp_x"});
	args::ValueFlag<double> _eps(parser, "float", "eps (the precision of the objective)", {"eps"});
	args::ValueFlag<double> _timelimit(parser, "float", "timelimit", {'t',"timelimit"});
	args::ValueFlag<int> _seed(parser, "int", "seed", {"seed"});
	args::Flag _trace(parser, "trace", "Activate trace. Updates of loup/uplo are printed while minimizing.", {"trace"});
	args::ValueFlag<double> _c1(parser, "double", "cognitive parameter", {"c1"});
	args::ValueFlag<double> _c2(parser, "double", "social parameter", {"c2"});
	args::ValueFlag<double> _x(parser, "double", "constriction factor", {"x"});
	args::ValueFlag<int> _np(parser, "int", "number of particles (inner pso)", {"np"});
	args::ValueFlag<int> _iter(parser, "int", "number of iterations per node (inner pso)", {"iter"});
	args::ValueFlag<int> _np_node(parser, "int", "number of particles (node pso)", {"np0"});
	args::ValueFlag<int> _iter_node(parser, "int", "number of iterations per node (node pso)", {"iter0"});

	args::Positional<std::string> filename(parser, "filename", "The name of the MINIBEX file.");

	try
	{
		parser.ParseCLI(argc, argv);
	}
	catch (args::Help&)
	{
		std::cout << parser;
		return 0;
	}
	catch (args::ParseError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}
	catch (args::ValidationError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}

	string filtering = (_filtering)? _filtering.Get() : "acidhc4";
	cout << "Filtering: " << filtering << endl;
	string linearrelaxation= (_linear_relax)? _linear_relax.Get() : "compo";
	cout << "Relax method: " << linearrelaxation << endl;
	string bisection= (_bisector)? _bisector.Get() : "lsmear";
	cout << "Bisector: " << bisection << endl;
	string strategy= (_strategy)? _strategy.Get() : "ibex_pso";
	cout << "Search strategy: " << strategy << endl;
	double prec= (_eps_x)? _eps_x.Get() : 1e-7 ;
	cout << "prec_x: " << prec << endl;
	double goalprec= (_eps)? _eps.Get() : 1e-6 ;
	cout << "prec_goal: " << goalprec << endl;
	double timelimit = (_timelimit)? _timelimit.Get() : 3600 ;
	cout << "timelimit: " << timelimit << endl;
	int nseed= (_seed)? _seed.Get() : 1 ;
	cout << "seed: " << nseed << endl;

 	double eqeps= 1.e-8;
 	double default_relax_ratio = 0.2;
	srand(nseed);

    System* orig_sys,*sys;
    LoupFinderDefault* loupfinder;
    std::size_t found = string(filename.Get().c_str()).find(".nl");
	if (found!=std::string::npos){
	       AmplInterface interface (argv[1]);
	       orig_sys= new System(interface);
     }else
           orig_sys = new System(argv[1]);

	// the extended system
	if(strategy=="solver") sys=orig_sys;
	else {
		sys=new ExtendedSystem(*orig_sys,eqeps);
		NormalizedSystem* norm_sys = new NormalizedSystem(*orig_sys,eqeps);
		loupfinder = new  LoupFinderDefault(*norm_sys,true);
	}


	// Build the bisection heuristic
	// --------------------------
	Bsc * bs;

	if (bisection=="roundrobin")
	  bs = new RoundRobin (prec);
	else if (bisection== "largestfirst")
          bs= new LargestFirst();
	else if (bisection=="smearsum")
	  bs = new SmearSum(*sys,prec);
	else if (bisection=="smearmax")
	  bs = new SmearMax(*sys,prec);
	else if (bisection=="smearsumrel")
	  bs = new SmearSumRelative(*sys,prec);
	else if (bisection=="smearmaxrel")
	  bs = new SmearMaxRelative(*sys,prec);
	else if (bisection=="lsmear")
	  bs = new LSmear(*sys,prec);

	else {cout << bisection << " is not an implemented  bisection mode "  << endl; return -1;}

	// Build the contractors

	// the first contractor called
	CtcHC4 hc4(sys->ctrs,0.01,true);
	// hc4 inside acid and 3bcid : incremental propagation beginning with the shaved variable
	CtcHC4 hc44cid(sys->ctrs,0.1,true);
	// hc4 inside xnewton loop
	CtcHC4 hc44xn (sys->ctrs,0.01,false);
	// The 3BCID contractor on all variables (component of the contractor when filtering == "3bcidhc4")
	Ctc3BCid c3bcidhc4(hc44cid);
	// hc4 followed by 3bcidhc4 : the actual contractor used when filtering == "3bcidhc4"
	CtcCompo hc43bcidhc4 (hc4, c3bcidhc4);
	// The ACID contractor (component of the contractor  when filtering == "acidhc4")
	CtcAcid acidhc4(*sys,hc44cid,1);
	// hc4 followed by acidhc4 : the actual contractor used when filtering == "acidhc4"
	CtcCompo hc4acidhc4 (hc4, acidhc4);


	string filtering2="acidhc4";

	Ctc* ctc;
	if (filtering == "hc4")
	  ctc= &hc4;
	else if (filtering =="acidhc4")
	  ctc= &hc4acidhc4;
	else if (filtering =="3bcidhc4")
	  ctc= &hc43bcidhc4;
	else{cout << filtering <<  " is not an implemented  contraction  mode "  << endl; return -1;}

    CtcNewton* ctcnewton=NULL;
    if(strategy=="solver"){
	    ctcnewton=new CtcNewton(sys->f_ctrs,5e8,prec,1.e-4);
		ctc =  new CtcCompo (*ctc, *ctcnewton);
	}


	Linearizer* lr;
	if (linearrelaxation=="art")
		  lr= new LinearizerCombo(*sys,LinearizerCombo::ART);
	else if  (linearrelaxation=="compo")
		  lr= new LinearizerCombo(*sys,LinearizerCombo::COMPO);
	else if (linearrelaxation=="xn")
	  lr= new LinearizerXTaylor (*sys, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM_OPP);
	//	else {cout << linearrelaxation  <<  " is not an implemented  linear relaxation mode "  << endl; return -1;}
	// fixpoint linear relaxation , hc4  with default fix point ratio 0.2
	CtcFixPoint* cxn;
	CtcPolytopeHull* cxn_poly;
	CtcCompo* cxn_compo;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn")
          {
		cxn_poly = new CtcPolytopeHull(*lr);
		cxn_compo =new CtcCompo(*cxn_poly, hc44xn);
		cxn = new CtcFixPoint (*cxn_compo, default_relax_ratio);
	  }
	//  the actual contractor  ctc + linear relaxation
	Ctc* ctcxn;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn")
          ctcxn= new CtcCompo  (*ctc, *cxn);
	else
	  ctcxn = ctc;

	// one point probed when looking for a new feasible point (updating the loup)
	int samplesize=1;

	if(strategy=="solver"){
		// A "CellStack" means a depth-first search.
		CellStack buff;


		Solver s(*sys, *ctcxn, *bs, buff, Vector(sys->nb_var,prec), Vector(sys->nb_var,POS_INFINITY));
		s.time_limit = timelimit;;
		s.trace=_trace;  // solutions are printed as soon as found when trace=1

		// Solve the system and get the solutions
		Solver::Status state=s.solve(sys->box);


		cout << state << endl;
		s.report();


		// Display the number of boxes (called "cells")
		// generated during the search
		cout << "number of cells=" << s.get_nb_cells() << endl;
		// Display the cpu time used
		cout << "cpu time used=" << s.get_time() << "s."<< endl;

		cout << argv[1] << " " << s.get_manifold().size() << " " << s.get_time() << " " <<
		s.get_nb_cells() << " " << (s.get_time()>timelimit) << endl;

		return 1;

	}


	// the optimizer : the same precision goalprec is used as relative and absolute precision
	OptimizerPSO* o=NULL;
	Optimizer* oo=NULL;
	//double c1, double c2, int particles, int iterations, double p) :
	// ** parameters for PSO algorithm **
	double c1 = 2.837;				// cognitive parameter (def. 2)
	double c2 = 1.597;				// social parameter (def. 2)
	double x=0.754;
	int particles = 10;		// amount of particles in pso algorithm
	int iterations = 5;		// number of iterations for pso algorithm

	int particles_node = 10;
	int iterations_node = 5;

	if(_c1) c1 =_c1.Get();
	cout << "cognitive parameter (c1):" << c1 << endl;
	if(_c2) c2 =_c2.Get();
	cout << "social parameter (c2):" << c2 << endl;
	if(_x) x =_x.Get();
	cout << "constrinction factor (x):" << x << endl;
	if(_np) particles =_np.Get();
	cout << "number of particles (inner pso):" << particles << endl;
	if(_iter) iterations =_iter.Get();
	cout << "number of iterations per node (inner pso):" << iterations << endl;
	if(_np_node) particles_node =_np_node.Get();
	cout << "number of particles (node pso):" << particles_node << endl;
	if(_iter_node) iterations_node =_iter_node.Get();
	cout << "number of iterations per node (node pso):" << iterations_node << endl;



	CellBufferOptim* buffer;
	BufferPSO* bufferPSO;
	TreeCellOpt* tree = new TreeCellOpt(orig_sys);
	TreeCellOpt* tree_single_node = new TreeCellOpt(orig_sys);

	PSOSwarm* swarm = new PSOSwarm(tree,orig_sys,c1,c2,x,particles,iterations);
	PSOSwarm* swarm_node = new PSOSwarm(tree_single_node,orig_sys,c1,c2,x,particles_node,iterations_node);

	if(strategy=="ibex_pso"){
		bufferPSO = new BufferPSO(swarm);
	}else if(strategy=="ibex_psodiv"){
		bufferPSO = new BufferPSO(swarm);
		buffer = new CellBeamSearch(*new CellHeap (*dynamic_cast<ExtendedSystem*>(sys)),
						*new CellHeap (*dynamic_cast<ExtendedSystem*>(sys)), *dynamic_cast<ExtendedSystem*>(sys));
	}else if(strategy == "ibex_opt" || strategy=="ibex_nodepso" || strategy=="ibex_rootpso"){
		buffer = new CellBeamSearch(*new CellHeap (*dynamic_cast<ExtendedSystem*>(sys)),
				*new CellHeap (*dynamic_cast<ExtendedSystem*>(sys)), *dynamic_cast<ExtendedSystem*>(sys));
	}



	/*
	if(strategy=="lbub")
		buffer = new CellDoubleHeap(*dynamic_cast<ExtendedSystem*>(sys));
	else if(strategy=="feasdiv")
		buffer = new CellBeamSearch(*new CellHeap (*dynamic_cast<ExtendedSystem*>(sys)),
								       *new CellHeap (*dynamic_cast<ExtendedSystem*>(sys)), *dynamic_cast<ExtendedSystem*>(sys));
	else {cout << strategy <<  " is not an implemented search strategy "  << endl; return -1;}*/

	cout << sys->nb_var << endl;
	cout << orig_sys->box << endl;
	cout << sys->box << endl;

	if(strategy=="ibex_pso"){
		o=new OptimizerPSO(orig_sys->nb_var, *ctcxn,*bs, *loupfinder, *bufferPSO, dynamic_cast<ExtendedSystem*>(sys)->goal_var(),
	    		prec,goalprec,goalprec, bufferPSO, swarm_node,false);
	}else if(strategy=="ibex_psodiv"){
		o=new OptimizerPSO(orig_sys->nb_var, *ctcxn,*bs, *loupfinder, *buffer, dynamic_cast<ExtendedSystem*>(sys)->goal_var(),
	    		prec,goalprec,goalprec, bufferPSO, swarm_node,false);
	}else if(strategy=="ibex_nodepso"){
		o=new OptimizerPSO(orig_sys->nb_var, *ctcxn,*bs, *loupfinder, *buffer, dynamic_cast<ExtendedSystem*>(sys)->goal_var(),
	    		prec,goalprec,goalprec, NULL, swarm_node,true);
	}else if(strategy=="ibex_rootpso"){
		o=new OptimizerPSO(orig_sys->nb_var, *ctcxn,*bs, *loupfinder, *buffer, dynamic_cast<ExtendedSystem*>(sys)->goal_var(),
	    		prec,goalprec,goalprec, NULL, swarm_node,false);
	}else if(strategy=="ibex_opt"){
		oo=new Optimizer(orig_sys->nb_var, *ctcxn,*bs, *loupfinder, *buffer, dynamic_cast<ExtendedSystem*>(sys)->goal_var(),
			    		prec,goalprec,goalprec);
	}

	// the trace
	if(o) o->trace=_trace;
	else oo->trace=_trace;

	if (o && o->trace)	cout << " sys.box " << sys->box << endl;

	// the allowed time for search
	if(o) o->timeout=timelimit;
	else oo->timeout=timelimit;

    vector<IntervalVector> sols;
    cout.precision(10);
	// the search itself
    std::cout << "optimize" << endl;
    if(o) o->optimize(orig_sys->box);
    else oo->optimize(orig_sys->box);

	// printing the results
	if (o && o->trace)
	  o->report();
	else if(oo->trace)
	  oo->report();


	delete bs;
	if (linearrelaxation=="xn" ||linearrelaxation=="compo" || linearrelaxation=="art" )
	  {
		//delete lr;
	    delete ctcxn;
	  //  delete cxn;
	  }


  list<pair<double,pair<double,int>>> loup_time;
	double uplo;
  if(o){
    cout << argv[1] << " " << o->get_uplo() << " " << o->get_loup() << " " << double(o->get_time()) << " " <<
             double(o->get_nb_cells()) << " " << (o->get_time()>timelimit) << " ";

		loup_time=o->loup_time;
		uplo=o->get_uplo();


  }else{
	    cout << argv[1] << " " << oo->get_uplo() << " " << oo->get_loup() << " " << double(oo->get_time()) << " " <<
	             double(oo->get_nb_cells()) << " " << (oo->get_time()>timelimit) << " ";

			loup_time=oo->loup_time;
			uplo= oo->get_uplo();
	}

	bool flag=true;
	for(auto pair:loup_time){

		if(flag && (get_obj_rel_prec(uplo,pair.first) <= 1e-4 ||
								get_obj_abs_prec(uplo,pair.first) <= 1e-4 )) {flag=false;
									cout << pair.second.first << " " << pair.second.second << " ";}

		if(get_obj_rel_prec(uplo,pair.first) <= 1e-6 ||
								get_obj_abs_prec(uplo,pair.first) <= 1e-6 ) {
									cout << pair.second.first << " " << pair.second.second << endl; break;}
	}

	return 0;

	}


	catch(ibex::SyntaxError& e) {
	  cout << e << endl;
	}
}

double get_obj_rel_prec(double uplo, double loup) {
	if (loup==POS_INFINITY)
		return POS_INFINITY;
	else if (loup==0)
		if (uplo<0) return POS_INFINITY;
		else return 0;
	else
		return (loup-uplo)/(fabs(loup));
}

double get_obj_abs_prec(double uplo, double loup) {
	return loup-uplo;
}
