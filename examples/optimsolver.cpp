/*
 * optimsolver.cpp
 *
 *  Created on: 05-10-2017
 *      Author: victor
 */


using namespace std;
using namespace ibex;

int main(int argc, char** argv){

	try {

	if (argc<9) {
		cerr << "usage: optimsolver filename filtering linear_relaxation strategy bisection prec goal_prec timelimit nseed lsmear_mode"  << endl;
		exit(1);
	}

	string filtering = argv[2];
	string linearrelaxation= argv[3];
    string strategy = argv[4];
	string bisection= argv[5];
	double prec= atof(argv[6]);
	double goalprec= atof (argv[7]);
	double timelimit = atof(argv[8]);
	int nseed= atoi(argv[9]);
	bool extended_embedded_system = atoi(argv[10]);

	/*
    ibex::lsmear_mode ls_mode;
    switch(atoi(argv[10])){
		case 0: ls_mode=BASIC; break;
		case 1: ls_mode=BASIC_JMID; break;
	}
	*/

 	double eqeps= 1.e-8;
 	double default_relax_ratio = 0.2;
	srand(nseed);

    System* orig_sys,*sys;
    std::size_t found = string(argv[1]).find(".nl");
	if (found!=std::string::npos){
	       AmplInterface interface (argv[1]);
	       orig_sys= new System(interface);
     }else
           orig_sys = new System(argv[1]);

	// the extended system
	if(strategy=="solver") sys=orig_sys;
	else sys=new ExtendedSystem(*orig_sys,eqeps);


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
	/*
	else if (bisection=="lsmear")
	  bs = new LSmear(*sys,prec, 0.45, ls_mode);
	*/
	else {cout << bisection << " is not an implemented  bisection mode "  << endl; return -1;}

	// Build the contractors

	// the first contractor called
	CtcHC4 hc4(sys->ctrs,0.01,true);
	// hc4 inside acid and 3bcid : incremental propagation beginning with the shaved variable
	CtcHC4 hc44cid(sys->ctrs,0.1,true);
	// hc4 inside xnewton loop
	CtcHC4 hc44xn (sys->ctrs,0.01,false);
	// The 3BCID contractor on all variables (component of the contractor when filtering == "3bcidhc4")
	Ctc3BCid c3bcidhc4(sys->nb_var,hc44cid);
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
	else if (filtering =="cse+hc4" || filtering =="cse+hc4+simplex"
	|| filtering =="cse+hc4+pseudoinverse" || filtering =="cse+hc4+gauss"
	|| filtering =="cse+hc4+gauss_pseudoinv" || filtering =="cse+hc4+gauss_mult"){
		int ctc_type=-1;
//		if(filtering =="cse+hc4+simplex") ctc_type=LinearSystem::SIMPLEX;
//		else if(filtering =="cse+hc4+pseudoinverse") ctc_type=LinearSystem::PSEUDOINVERSE;
//		else if(filtering =="cse+hc4+gauss") ctc_type=LinearSystem::GAUSS_JORDAN;
//		else if(filtering =="cse+hc4+gauss_mult") ctc_type=LinearSystem::MULT_GAUSS_JORDAN;
//		else if(filtering =="cse+hc4+gauss_pseudoinv") ctc_type=LinearSystem::GAUSS_PSEUDOINV;

		Array<const ExprSymbol> vars_dag(sys->args.size());
        const ExprNode& dag_root=ExprNode2Dag::generate(*sys, vars_dag);
	    Function *f = new Function(vars_dag, dag_root, "mi_dag");

	    CtcDag* dagctc=new CtcDag(*f, sys->ctrs);
	    CtcFixPoint* dagctc_fp = new CtcFixPoint(*dagctc,0.01);
	    if(filtering =="cse+hc4") ctc = dagctc_fp;

//	    else{
//
//			Array<Ctc> a;
//			a.add(getEmbeddedLinearSystems<ExprAdd>(*dagctc, false, ctc_type, extended_embedded_system/*, sys->nb_ctr*/));
//			a.add(getEmbeddedLinearSystems<ExprMul>(*dagctc, true, ctc_type, extended_embedded_system));
//
//			CtcCompo* linear_systems = (a.size()>0)? new CtcCompo(a) : NULL;
//
//			ctc = (linear_systems)? new CtcFixPoint(*new CtcCompo(*dagctc_fp,*linear_systems),0.01) : dagctc_fp;
//
//	    }
		if(filtering2 == "acidhc4"){
			CtcFixPoint* dagctc_fp4cid = new CtcFixPoint(*dagctc,0.1);
			CtcAcid* acid_dagctc = new CtcAcid(*sys,*dagctc_fp4cid,1);
			ctc = new CtcCompo(*ctc, *acid_dagctc);
		}

	}else{cout << filtering <<  " is not an implemented  contraction  mode "  << endl; return -1;}

    CtcNewton* ctcnewton=NULL;
    if(strategy=="solver"){
	    ctcnewton=new CtcNewton(sys->f,5e8,prec,1.e-4);
		ctc =  new CtcCompo (*ctc, *ctcnewton);
	}

	// The CtcXNewton contractor
	// corner selection for linearizations : two corners are selected, a random one and its opposite
//	vector<LinearRelaxXTaylor::corner_point> cpoints;
//	cpoints.push_back(LinearRelaxXTaylor::RANDOM);
//	cpoints.push_back(LinearRelaxXTaylor::RANDOM_INV);

//	LinearRelax* lr=NULL;
//	if (linearrelaxation=="art")
//	  lr= new LinearRelaxCombo(*sys,LinearRelaxCombo::ART);
//	else if  (linearrelaxation=="compo")
//	  lr= new LinearRelaxCombo(*sys,LinearRelaxCombo::COMPO);
//	else if (linearrelaxation=="xn" || linearrelaxation=="xn_onlyy")
//	  lr= new LinearRelaxXTaylor (*sys,cpoints);





	//	else {cout << linearrelaxation  <<  " is not an implemented  linear relaxation mode "  << endl; return -1;}
	// fixpoint linear relaxation , hc4  with default fix point ratio 0.2
//	CtcFixPoint* cxn;
//	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn" )
//	  //cxn = new CtcLinearRelaxation (*lr, hc44xn);
//		cxn = new CtcFixPoint (*new CtcCompo(*new CtcPolytopeHull(*lr, CtcPolytopeHull::ALL_BOX), hc44xn), default_relax_ratio);
//	//  the actual contractor  ctc + linear relaxation
//	else if (linearrelaxation=="xn_onlyy")
//        cxn = new CtcFixPoint (*new CtcCompo(*new CtcPolytopeHull(*lr, CtcPolytopeHull::ONLY_Y), hc44xn), default_relax_ratio);
// 	else if (linearrelaxation=="xn_art"){
//		cxn = new CtcFixPoint (*new CtcCompo(
//		*new CtcPolytopeHull(*new LinearRelaxXTaylor (*sys,cpoints), CtcPolytopeHull::ALL_BOX),
//		*new CtcPolytopeHull(*new LinearRelaxCombo(*sys,LinearRelaxCombo::ART), CtcPolytopeHull::ALL_BOX), hc44xn), default_relax_ratio);
//    }

//	if (linearrelaxation=="taylor")
//		cxn = new CtcFixPoint (*new CtcCompo(*new TaylorLinear(*sys), hc44xn), default_relax_ratio);


	Ctc* ctcxn;
//	if (linearrelaxation=="xn_art" || linearrelaxation=="compo" ||
//			linearrelaxation=="art"|| linearrelaxation=="xn" || linearrelaxation=="xn_onlyy" || linearrelaxation=="taylor")
//          ctcxn= new CtcCompo  (*ctc, *cxn);
//	else
	  ctcxn = ctc;
	// one point probed when looking for a new feasible point (updating the loup)
	int samplesize=1;

	if(strategy=="solver"){
		// A "CellStack" means a depth-first search.
		CellStack buff;
		Solver s(*ctcxn,*bs,buff);
		s.time_limit = timelimit;;
		s.trace=1;  // solutions are printed as soon as found when trace=1

		// Solve the system and get the solutions
		vector<IntervalVector> sols=s.solve(sys->box);

		// Display the number of solutions
		cout << "number of solutions=" << sols.size() << endl;

		// Display the number of boxes (called "cells")
		// generated during the search
		cout << "number of cells=" << s.nb_cells << endl;
		// Display the cpu time used
		cout << "cpu time used=" << s.time << "s."<< endl;

		cout << argv[1] << " " << sols.size() << " " << s.time << " " <<
				s.nb_cells << " " << (s.time>timelimit) << endl;

		return 1;

	}


	// the optimizer : the same precision goalprec is used as relative and absolute precision
	Optimizer* o=NULL;

	if(strategy=="lbub")
	    o=new Optimizer(*orig_sys,*ctcxn,*bs, prec,goalprec,goalprec,samplesize,Optimizer::default_equ_eps, false);
//	else if(strategy=="feasdiv")
//		o=new FeasibleDiving(*orig_sys,*ctcxn,*bs, prec,goalprec,goalprec,samplesize,Optimizer::default_equ_eps, false);
	else {cout << strategy <<  " is not an implemented search strategy "  << endl; return -1;}

	// the trace
	o->trace=1;

	if (o && o->trace)	cout << " sys.box " << sys->box << endl;

	// the allowed time for search
	o->timeout=timelimit;

    vector<IntervalVector> sols;
    cout.precision(10);
	// the search itself
	o->optimize(orig_sys->box);

	// printing the results
	if (o->trace)
	  o->report();


	delete bs;
	if (linearrelaxation=="xn" ||linearrelaxation=="compo" || linearrelaxation=="art" )
	  {
		//delete lr;
	    delete ctcxn;
	  //  delete cxn;
	  }




    cout << argv[1] << " " << o->uplo << "," << o->loup << " " << double(o->time) << " " <<
             double(o->nb_cells) << " " << (o->time>timelimit) << endl;


	return 0;

	}


	catch(ibex::SyntaxError& e) {
	  cout << e << endl;
	}
}

