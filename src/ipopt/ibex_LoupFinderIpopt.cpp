//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderIpopt.cpp
// Author      : Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Aug 03, 2023
// Last update : Mar 07, 2024
//============================================================================





#include <stdio.h>
#include <limits>
#include "ibex_LoupFinderIpopt.h"
#include "ibex_Timer.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <ctime>
#include "ibex_BscHijackGuard.h"



#include <cassert>

using namespace Ipopt;
using namespace std;


namespace ibex {
  double expansion_precisionB=1.e-6;
  extern double optimizer_node_lb;
  extern int optimizer_node_depth;

  namespace {
    // IBEX_IPOPT_LOG=<file>: one JSON line per Ipopt call (what makes it pay off)
    FILE* ipopt_log() {
      static FILE* f = NULL; static bool init = false;
      if (!init) { init = true; const char* p = getenv("IBEX_IPOPT_LOG"); if (p) f = fopen(p, "a"); }
      return f;
    }
    std::vector<Vector> ipopt_sols;   // previous Ipopt solutions (for the distance)
  }
  double ipopt_diam=1.e8;
  // half-width (relative when |v|>1) of the box around Ipopt's point that the
  // correction search explores; <= 0: expansion_precisionB (1e-6)
  double ipopt_expansion=-1;
  // with ipopt_diam < 0: Ipopt only starts from points with |x_i| <= this
  double ipopt_start_max=1e8;
  // CPU limit of one Ipopt call (Ipopt's max_cpu_time); <= 0: none
  double ipopt_max_time=-1;
  // at the first Ipopt call, this many extra Ipopt runs from random points of
  // the root box (each certified like the main one); counts of the extra runs
  // that improved the loup
  // CPU limit of one Ipopt call enforced from the evaluation callbacks: past
  // it they report an evaluation error and Ipopt gives up. max_cpu_time is
  // only checked between iterations, and Ipopt 3.11 can loop inside one line
  // search forever (matrix2 and ex9_1_5 from far starting points). <= 0: none
  double ipopt_call_limit=5;
  namespace { std::clock_t ipopt_deadline=0; bool ipopt_deadline_on=false;
    struct IpoptDeadline {};
    // throws out of Ipopt (RethrowNonIpoptException): returning false from a
    // callback inside the restoration phase crashes Ipopt 3.11
    inline bool past_deadline() { if (ipopt_deadline_on && std::clock() > ipopt_deadline) throw IpoptDeadline(); return false; } }
  int ipopt_multistart=0; int ipopt_ms_improved=0; int ipopt_ms_seed=1;
  namespace { int ms_index=0; std::mt19937 ms_rng; }
  // Guided restarts: at a node where Ipopt is not scheduled, an extra Ipopt
  // run from a random point of the node's box, as long as these runs have
  // used less than ipopt_guided (fraction) of the CPU time since the first
  // call, and the box contains no previous Ipopt solution (a basin already
  // visited) -- replaced by: of 10 random points of the box, the farthest
  // from every previous Ipopt solution (maximin). The B&B order (lowest lower
  // bound first) makes the node a promising region not yet discarded. ipopt_guided_local: Ipopt bounded by
  // the node's box instead of the root box. <= 0: off.
  double ipopt_guided=-1; bool ipopt_guided_local=false;
  int ipopt_guided_calls=0, ipopt_guided_improved=0; double ipopt_guided_time=0;
  // patience: after this many guided runs in a row that do not improve the
  // loup, they stop until the loup improves by any means. <= 0: never stop
  int ipopt_guided_patience=0;
  namespace { int guided_streak=0; double guided_last_loup=POS_INFINITY; }
  namespace { bool guided_now=false; std::clock_t guided_t0=0; bool guided_started=false; std::mt19937 g_rng(7); }
  // with ipopt_expansion: first try the 1e-6 box (with this CPU timeout) and
  // widen to ipopt_expansion only when it gives no loup; <= 0: no first try
  double ipopt_expansion_first=-1;
  namespace { bool corr_first_try=false; }
  bool ipopt_defer=false;
  
  LoupFinderIpopt::LoupFinderIpopt(const System& sys,const System& normsys, const ExtendedSystem& extsys) : sys(sys), normsys(normsys), extsys(extsys), solution(sys.nb_var), the_box(sys.box), ipopt_box(sys.box) {
	try {
		df = new Function(*sys.goal,Function::DIFF);
		//                cout << " nb ctr " << sys.nb_ctr << endl;
		if (sys.nb_ctr>0) {
			dg = new Function*[sys.f_ctrs.image_dim()];

			for (int i=0; i<sys.f_ctrs.image_dim(); i++) {
				dg[i] = new Function(sys.f_ctrs[i], Function::DIFF);
			}
		} else {
			dg = NULL;
		}

	}
	catch(Exception&) {
		//TODO: replace with ExprDiffException.
		// Currently, DimException is also sometimes raised.
		cerr << "Warning: symbolic differentiation has failed ==> ipopt hessian disabled" << endl;
		df = NULL;
		dg = NULL;
	}
	//	cout << " fin initialisation gradients " << endl;
      //The IpOPT initialization should be here
      app = IpoptApplicationFactory();
      app->RethrowNonIpoptException(true);
     

      // Change some options
      // Note: The following choices are only examples, they might not be
      // suitable for your optimization problem.
      
      //      app->Options()->SetNumericValue("tol", 1e-8);
      //      app->Options()->SetNumericValue("constr_viol_tol", 1e-10);
      // app->Options()->SetStringValue("mu_strategy", "adaptive");
      app->Options()->SetStringValue("output_file", "ipopt.out");
      // Ipopt prints its license banner on stdout at the first solve, which
      // print_level does not cover. ibexopt-ml speaks a line-of-JSON protocol on
      // stdout, so one banner is enough to break every client: "sb" (suppress
      // banner) is Ipopt's own option for exactly this.
      app->Options()->SetStringValue("sb", "yes");
      app->Options()->SetIntegerValue("print_level",0);
      app->Options()->SetIntegerValue("max_iter",1000);
      app->Options()->SetStringValue("linear_solver","mumps");
      // Ipopt checks f and g for NaN/Inf but not their derivatives by default;
      // a NaN in the Jacobian or Hessian goes to MUMPS, which may then crash
      // (water, starting near a division by 0 once Ipopt runs on wide boxes)
      app->Options()->SetStringValue("check_derivatives_for_naninf","yes");
      // Initialize the IpoptApplication and process the options
      
      ApplicationReturnStatus status;
      status = app->Initialize();
      //      cout << "initialisation ipopt " << endl;
      if (status != Solve_Succeeded) {
      	std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
      	exit(0);
      }

      // This object is its own Ipopt::TNLP and is passed to OptimizeTNLP() as a
      // raw "this": Ipopt wraps it in SmartPtrs (TNLPAdapter keeps one), and when
      // the application is torn down the last of them would "delete this" -- on an
      // object that is not owned by Ipopt and, worse, that is already inside its
      // own destructor. That is the "error with this delete" that made the
      // original optimizer04-ipopt.cpp leak the finder rather than free it.
      // One self-reference, dropped by hand in the destructor once the
      // application is gone, keeps the count from ever reaching zero.
      AddRef(NULL);
      }

    LoupFinderIpopt::~LoupFinderIpopt()
    {
      // Order matters: the application (and with it every SmartPtr Ipopt holds on
      // this TNLP) goes first, then the self-reference taken in the constructor,
      // so that ~ReferencedObject sees a count of zero.
      app = NULL;
      ReleaseRef(NULL);
    }


    std::pair<IntervalVector, double> LoupFinderIpopt::find(const IntervalVector& box, const IntervalVector& loup_point, double loup) {
      double loup0=loup;
      if (ipopt_calls== -1) ipopt_box=box;  // initialization of ipopt_box before first call.
      // MINLP is not supported in this port: Ipopt always runs on the root box.
 
      the_box  = box;
      double newloup=false;
      IntervalVector loup_point0=loup_point;
      if (recursive_call){
	ipopt_calls++;  // at first call, ipopt_calls=0
	bool scheduled = (ipopt_calls%ipopt_frequency==0 || ipopt_calls<ipopt_first || ipopt_calls==10 || ipopt_calls==20  || ipopt_calls==50);
	bool narrow;
	if (ipopt_diam >= 0) narrow = box.max_diam() <= ipopt_diam;
	else {
	  // no width threshold: what matters is where Ipopt starts. A box "at
	  // infinity" (x in [9e307, 1.8e308]) has no reasonable starting point:
	  // the evaluations overflow and Ipopt fails at once
	  // A forced call starts from the point another loup finder just found,
	  // which can be just as far (matrix2: 1e20, where Ipopt hangs in a line
	  // search that max_cpu_time does not interrupt).
	  narrow = true;
	  ibex::Vector v = force ? solution : start_point(box);
	  for (int i=0; i<v.size(); i++) if (!(std::fabs(v[i]) <= ipopt_start_max)) { narrow=false; break; }
	}
	guided_now = false;
	if (loup < guided_last_loup) { guided_last_loup = loup; guided_streak = 0; }
	if (ipopt_guided > 0 && !scheduled && !force && narrow
	    && (ipopt_guided_patience <= 0 || guided_streak < ipopt_guided_patience)) {
	  if (!guided_started) { guided_started = true; guided_t0 = std::clock(); }
	  double elapsed = double(std::clock() - guided_t0) / CLOCKS_PER_SEC;
	  if (ipopt_guided_time < ipopt_guided * elapsed) {
	    guided_now = true;
	  }
	}
	// ipopt_defer: a scheduled call that falls on a box too wide is not lost,
	// it waits for the next box narrow enough
	if (ipopt_defer && scheduled && !narrow) ipopt_pending = true;
	if (
	    narrow
	    &&
	    (scheduled || (ipopt_defer && ipopt_pending) || force || guided_now)
	    ){
	  std::clock_t guided_c0 = std::clock();
	  ipopt_pending = false;
	
	  //	  cout << "nb_cells " <<  optimizer->get_nb_cells() << endl;
	  bool was_forced = force;
	  // the node being bounded: read now, the correction search below runs an
	  // Optimizer of its own that overwrites these
	  log_node_lb = optimizer_node_lb; log_depth = optimizer_node_depth;
	  int nstarts = (ipopt_calls==0 && ipopt_multistart>0) ? ipopt_multistart : 0;
	  if (nstarts>0) ms_rng.seed(ipopt_ms_seed);
	  for (int ms=0; ms<=nstarts; ms++) {
	  ms_index = ms;
	  double loup_ms = loup;
	  double loup_in = loup;
	  Timer ipopt_timer; ipopt_timer.start();
	  if (ipopt_max_time > 0) app->Options()->SetNumericValue("max_cpu_time", ipopt_max_time);
	  if (getenv("IBEX_IPOPT_TRACE")) { std::cerr.precision(17); std::cerr << "call " << ipopt_calls << " force " << force << " start " << (force ? solution : start_point(box)) << " box " << box << std::endl; }
	  ipopt_deadline_on = ipopt_call_limit > 0;
	  ipopt_deadline = std::clock() + (std::clock_t) (ipopt_call_limit * CLOCKS_PER_SEC);
	  ApplicationReturnStatus status;
	  try {
	    status = app->OptimizeTNLP(this);
	  } catch (IpoptDeadline&) {
	    status = Maximum_CpuTime_Exceeded;
	    optimalValue = POS_INFINITY;   // no usable result
	  }
	  ipopt_deadline_on = false;
	  force=0;  // after call , force is reset to 0
	  ipopt_timer.stop();
	  log_pending = ipopt_log()!=NULL;
	  log_forced = was_forced; log_loup_in = loup_in; log_time = ipopt_timer.get_time();
	  log_value = optimalValue; log_status = (int) status; log_width = box.max_diam();
	  log_dist = POS_INFINITY;
	  log_inner=-1; log_corr_status=-1; log_corr_nodes=0; log_corr_time=0; log_corr_loup=POS_INFINITY; log_goal=POS_INFINITY; log_viol=POS_INFINITY;
	  if (ipopt_sols.size()<2000 && (int) status > -10 && std::isfinite(optimalValue)) {
	    bool fin = true; for (int i=0; i<solution.size(); i++) fin = fin && std::isfinite(solution[i]);
	    if (fin) ipopt_sols.push_back(solution);
	  }
	  if (log_pending) {
	    for (size_t k=0; k+1<ipopt_sols.size(); k++) {
	      double d = 0;
	      for (int i=0; i<solution.size(); i++) {
	        double w = sys.box[i].diam(); if (!(w>0) || w>1e12) w = 1+std::fabs(solution[i]);
	        d = std::max(d, std::fabs(solution[i]-ipopt_sols[k][i])/w);
	      }
	      log_dist = std::min(log_dist, d);
	    }
	  }
	  //	  cout << " status " << status << endl;
	  //	  if (status == Solve_Succeeded) {
	  //        	std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
	  //	  cout << " optimal_value ipopt "<< optimalValue << endl;
	  //	  cout << " solution ipopt " << solution << endl;
	  // An Ipopt failure (status <= -10: invalid number, problem definition,
	  // option...) leaves optimalValue and solution meaningless (value 0 at the
	  // starting point): nothing to certify, and the correction search around
	  // that point would burn its whole timeout for nothing.
	  bool usable = (int) status > -10 && (int) status != 2 && (int) status != 4
	                && std::isfinite(optimalValue);
	  for (int i=0; usable && i<solution.size(); i++) usable = std::isfinite(solution[i]);
	  if( usable && optimalValue< loup){
	    
	    for (int i=0; i< box.size(); i++) {
	      bound_check_i(sys,solution,i);
	    }
	    //	    cout << " pt avant check "  << solution << endl;
	    if (log_pending) {
	      log_goal = normsys.goal_ub(solution);
	      log_inner = normsys.is_inner(IntervalVector(solution));
	      if (normsys.nb_ctr>0) {
	        IntervalVector g = normsys.f_ctrs.eval_vector(IntervalVector(solution));
	        log_viol = NEG_INFINITY;
	        for (int k=0; k<g.size(); k++) log_viol = std::max(log_viol, g[k].ub());
	      } else log_viol = NEG_INFINITY;
	    }
	    if (check(normsys,solution,loup,false))
		  { 
		    loup_point0=solution;
		    if (optimizer->trace) 
		      cout << "*** ipopt      " ;
		  }

	    if (true){  // no integer variables in this port
	    //	    else if (optimalValue< loup && sys.get_integer_variables()->size() < sys.nb_var){
	      //	      cout << " solution " << solution << endl;
	      //	      cout << " optimal value " << optimalValue << " loup " << loup << endl;
		  double ipoptloup=POS_INFINITY;
		  if (ipopt_expansion > 0 && ipopt_expansion_first > 0) {
		    corr_first_try = true;
		    correct_ipopt_sol(solution, ipoptloup);
		    corr_first_try = false;
		    if (!(ipoptloup < loup)) correct_ipopt_sol(solution, ipoptloup);
		  } else
		    correct_ipopt_sol(solution, ipoptloup);
		  if (ipoptloup < loup)
		    {
		      loup_point0=solution;
		      loup=ipoptloup;
		      if (optimizer->trace)
			cout << "*** ipopt+corr " ;
		      newloup=true;
		    }
	    }
	  }
	  if (ms>0 && loup < loup_ms) ipopt_ms_improved++;
	  if (guided_now && loup < loup_ms) ipopt_guided_improved++;
	  } // multistart
	  ms_index = 0;
	  if (guided_now) {
	    if (loup < guided_last_loup) { guided_last_loup = loup; guided_streak = 0; } else guided_streak++;
	    ipopt_guided_calls++;
	    ipopt_guided_time += double(std::clock() - guided_c0) / CLOCKS_PER_SEC;
	    guided_now = false;
	  }
	}
	
      }
        
      if (log_pending && recursive_call) {   // not from inside the correction search
        log_pending = false;
        fprintf(ipopt_log(), "{\"call\": %d, \"depth\": %d, \"node_lb\": %.17g, \"loup_in\": %.17g, \"width\": %.6g, "
                "\"forced\": %d, \"status\": %d, \"value\": %.17g, \"loup_out\": %.17g, \"improved\": %d, "
                "\"dist\": %.6g, \"time\": %.6g, \"goal\": %.17g, \"inner\": %d, \"viol\": %.6g, "
                "\"corr_status\": %d, \"corr_nodes\": %d, \"corr_time\": %.6g, \"corr_loup\": %.17g}\n",
                ipopt_calls, log_depth, log_node_lb, log_loup_in, log_width, (int) log_forced,
                log_status, log_value, loup, (int) (loup < log_loup_in), log_dist, log_time,
                log_goal, log_inner, log_viol, log_corr_status, log_corr_nodes, log_corr_time, log_corr_loup);
        fflush(ipopt_log());
      }
      if (loup < loup0)
	return std::make_pair(loup_point0, loup);
      else  throw NotFound();


}



    bool LoupFinderIpopt::get_nlp_info(int& n, int& m, int& nnz_jac_g,
					int& nnz_h_lag, IndexStyleEnum& index_style)
    {
        n = sys.nb_var;

        m = sys.nb_ctr;

        nnz_jac_g = n * m;

        nnz_h_lag = (n * (n+1))/2;

        index_style = TNLP::C_STYLE;

        return true;
    }



    // returns the variable bounds
    bool LoupFinderIpopt::get_bounds_info(int n, Number* x_l, Number* x_u, int m, Number* g_l, Number* g_u)
    {
        // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
        // If desired, we could assert to make sure they are what we think they are.


        assert(n == sys.nb_var);
        assert(m == sys.nb_ctr);
	//        cout << "ipopt box " << ipopt_box << endl;
        Array<NumConstraint> contrains = sys.ctrs;
        for (int i = 0; i < n; i++) {
	  x_l[i] = (guided_now && ipopt_guided_local) ? the_box[i].lb() : ipopt_box[i].lb();   // Ipopt has for domain the initial box after contraction

        }


        for (int i = 0; i < n; i++) {
	  x_u[i] = (guided_now && ipopt_guided_local) ? the_box[i].ub() : ipopt_box[i].ub();
        }
	
	
        for (int i = 0; i < m; i++)
        {

            if (contrains[i].op == EQ)
            {
	      g_u[i] = 0;
	      g_l[i] = 0;
	      // a "thick" equality f0(x) - [a,b] = 0 means f0(x) in [a,b]; eval_g
	      // returns the midpoint f0(x) - mid([a,b]), so the constraint is
	      // g in [-rad, rad] (rad: the width the constant gives at a point)
	      ibex::Vector pt = start_point(ipopt_box);
	      Interval fi = contrains[i].f.eval(IntervalVector(pt));
	      if (!fi.is_empty() && !fi.is_unbounded() && fi.rad() > 1e-9) {
	        g_u[i] = fi.rad();
	        g_l[i] = -fi.rad();
	      }
	    }

            else if (contrains[i].op == GEQ)
            {
                g_u[i] = 2.e19;
                g_l[i] = 0;

            }

            else if (contrains[i].op == LEQ)
            {
                g_u[i] = 0;
                g_l[i] = -2.e19;

            }
        }

        return true;
    }


    // returns the initial point for the problem
    // the center of the box; an unbounded box has no finite center: there,
    // the point of the interval closest to 0 (a finite one), so that Ipopt
    // can run there too
    ibex::Vector LoupFinderIpopt::start_point(const IntervalVector& b) {
	ibex::Vector v = b.mid();
	for (int i=0; i<b.size(); i++) {
	  const Interval& xi = b[i];
	  if (!std::isfinite(v[i]) || std::fabs(v[i]) > 1e300 || xi.is_unbounded()) {
	    if (xi.contains(0)) v[i] = 0;
	    else if (xi.lb() > 0) v[i] = std::isfinite(xi.lb()) ? xi.lb() + std::max(1.0, std::fabs(xi.lb())) * 1e-3 : 1.0;
	    else v[i] = std::isfinite(xi.ub()) ? xi.ub() - std::max(1.0, std::fabs(xi.ub())) * 1e-3 : -1.0;
	    if (!xi.contains(v[i])) v[i] = xi.lb() > 0 ? xi.lb() : xi.ub();
	  }
	}
	return v;
    }

    bool LoupFinderIpopt::get_starting_point(int n, bool init_x, Number* x,
        bool init_z, Number* z_L, Number* z_U,
        int
        m, bool init_lambda,
        Number* lambda)
    {
        // the starting point is the center of the current box, and in case of a
        // new loup found by another loup finder (force=true) , the corresponding point (solution)

        assert(init_x == true);

        assert(init_z == false);
        assert(init_lambda == false);

        // initialize to the given starting point
	//	ibex::Vector v = the_box.random();
	ibex::Vector v = start_point(the_box);
	if (ms_index>0) {
	  // a random point of the root box; around the default start, at most
	  // max(10, 10|c|) away from it in each unbounded or very wide direction
	  for (int i=0; i<n; i++) {
	    const Interval& xi = ipopt_box[i];
	    double c = v[i], S = std::max(10.0, 10*std::fabs(c));
	    double lo = std::max(xi.lb(), c-S), hi = std::min(xi.ub(), c+S);
	    if (xi.diam() <= 2*S) { lo = xi.lb(); hi = xi.ub(); }
	    std::uniform_real_distribution<double> U(lo, hi);
	    v[i] = (hi > lo) ? U(ms_rng) : lo;
	  }
	}
	if (guided_now) {
	  // maximin: of 10 random points of the node's box (same clamp around its
	  // default start), the one farthest from every previous Ipopt solution
	  // (max-norm, each coordinate scaled by its sampling range)
	  ibex::Vector c0 = v, best = v, w(n);
	  double best_d = -1;
	  for (int t=0; t<10; t++) {
	    ibex::Vector cand(n);
	    for (int i=0; i<n; i++) {
	      const Interval& xi = the_box[i];
	      double c = c0[i], S = std::max(10.0, 10*std::fabs(c));
	      double lo = std::max(xi.lb(), c-S), hi = std::min(xi.ub(), c+S);
	      if (xi.diam() <= 2*S) { lo = xi.lb(); hi = xi.ub(); }
	      w[i] = (hi > lo) ? hi - lo : 1.0;
	      std::uniform_real_distribution<double> U(lo, hi);
	      cand[i] = (hi > lo) ? U(g_rng) : lo;
	    }
	    double d = POS_INFINITY;
	    for (size_t k=0; k<ipopt_sols.size(); k++) {
	      double dk = 0;
	      for (int i=0; i<n; i++) dk = std::max(dk, std::fabs(cand[i]-ipopt_sols[k][i])/w[i]);
	      d = std::min(d, dk);
	    }
	    if (d > best_d) { best_d = d; best = cand; }
	  }
	  v = best;
	}
	//	cout << " force " << force << endl;
	if (force && ms_index==0)  v=solution;

	for (int i=0; i<n ; i++)
	  x[i]=v[i];
        return true;
    }

    // returns the value of the objective function
    bool LoupFinderIpopt::eval_f(int n, const Number* x, bool new_x, Number& obj_value)
    {
        if (past_deadline()) return false;
        assert(n == the_box.size());
        ibex::Vector v(n);

        for (int i = 0; i < n; i++) {
            v[i] = x[i];
        }

        obj_value = sys.goal_ub(v);
     
        return true;
    }


    // return the gradient of the objective function grad_{x} f(x)

    bool LoupFinderIpopt::eval_grad_f(int n, const Number* x, bool new_x, Number* grad_f)
    {
        if (past_deadline()) return false;
      assert(n == the_box.size());

      
        ibex::IntervalVector v(n);

        Function* funcion = sys.goal;
	//	cout << " goal " << *funcion << endl;
	
        for (int i = 0; i < v.size(); i++)
	  {
            v[i] = x[i];
	  }
	ibex::IntervalVector gradiente(n);
  	gradiente= funcion->gradient(v);
        //  gradiente= df->eval_vector(v);
        //	  cout << " gradient " << gradiente << endl;
	//        cout << " grad size " << gradiente.size() << endl;
	//	cout << " v " << v << endl;
        //        cout << "df[0] " << df[0] << endl;
        //  cout << "df[0] val" << *df[0].eval(v)[0].mid() <<  " " << gradiente[0].mid() << endl;
        for (int i = 0; i < n; i++)
        {
	  grad_f[i] = gradiente[i].mid();
	  //	  double val=df[i].eval(v).mid();
	  //	  cout << i << " df[0] val" << val << " "  << gradiente[i].mid()        << endl;
	  //	  grad_f[i] = df[i].eval(v).mid();

        }


        return true;
    }

    // return the value of the constraints: g(x)
    bool LoupFinderIpopt::eval_g(int n, const Number* x, bool new_x, int m, Number* g)
    {
        if (past_deadline()) return false;
        assert(n == sys.nb_var);
        assert(m == sys.nb_ctr);

        ibex::IntervalVector v(n);

        for (int i = 0; i < v.size(); i++)
            v[i] = x[i];


        for (int i = 0; i < m; i++)
        {
	    g[i] = sys.ctrs[i].f.eval(v).mid();
        }
        return true;
    }

    // return the structure or values of the jacobian
    bool LoupFinderIpopt::eval_jac_g(int n, const Number* x, bool new_x, int m, int nele_jac, int* iRow, int* jCol, Number* values)
    {

        ibex::IntervalVector v(n);

        int cont = 0;
        if (values == NULL) {
            // return the structure of the jacobian
            int p = 0;
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++) {
                    iRow[p] = i;
                    jCol[p] = j;
                    p++;
                }


        }
        else {
            // return the values of the jacobian of the constraints
            for (int i = 0; i < v.size(); i++) {
                v[i] = x[i];
               // cout << v[i] << endl;
            }

            Array<NumConstraint> contrains = sys.ctrs;
       
            for (int c = 0; c < m; c++)
            {

	      IntervalVector grad(n);
	      grad= sys.ctrs[c].f.gradient(v);
	       //               cout << "restriccion:" << c << endl;
	       //               cout <<grad<< endl;

	       for (int k = 0; k < n; k++)
                    {

                            values[cont] = grad[k].mid();
			    //		    cout << cont << " " << values[cont] << endl;
			    cont++;

                    }

            }

        }

        return true;
    }
 

  
    bool LoupFinderIpopt::eval_h(int n, const Number* x, bool new_x, Number obj_factor, int m, const Number* lambda, bool new_lambda, int nele_hess, int* iRow, int* jCol, Number* values)
    {

      if( values == NULL )
	{
      // return the structure. This is a symmetric matrix, fill the lower left
      // triangle only.
 
      // the hessian  is actually dense
	  int idx = 0;
	  for( int row = 0; row < n; row++ )
	    {
	      for( int col = 0; col <= row; col++ )
		{
		  iRow[idx] = row;
		  jCol[idx] = col;
		  idx++;
		}
	    }
	  assert(idx == nele_hess);
        }
  
      else {
	ibex::IntervalVector v(n);
	for (int i = 0; i < n; i++) {
	  v[i] = x[i];

	}

	IntervalMatrix jacobianobj(n,n);
	df->jacobian(v,jacobianobj);

	int nn=(n*(n+1))/2;

	for (int i=0; i< nn ; i++)
	  values[i]=0;
	int l=0;
	for (int i = 0; i < n; i++){
          for (int j=0; j<=i; j++){
	    //	    cout << " jacobianobj " << i << "  " << j << " " << jacobianobj[i][j].mid() << endl;
	   
	    values[l]+=obj_factor*jacobianobj[i][j].mid();
	    l++;
	  }
	}
	for (int k=0 ; k< m; k++){
          l=0;
	  Function* f= dg[k];
	  //	  cout << "function " << *f << endl;
	  IntervalMatrix jacobianct(n,n);
	  f->jacobian(v, jacobianct);
	  for (int i = 0; i < n; i++){
	    for (int j=0; j<=i; j++){
	      // cout << k << " jacobianct " << i << "  " << j << " " << jacobianct[i][j].mid() << endl;
	      values[l]+=lambda[k]*jacobianct[i][j].mid();
	      l++;
	    }
	  }
	}
      }
      
      return true;
    }
    

    void LoupFinderIpopt::finalize_solution(SolverReturn status,
        int n, const Number* x, const Number* z_L, const Number* z_U,
        int m, const Number* g, const Number* lambda,
        Number obj_value,
        const IpoptData* ip_data,
        IpoptCalculatedQuantities* ip_cq)
    {
        // here is where we would store the solution to variables, or write to a file, etc
        // so we could use the solution.

        // For this example, we write the solution to the console
        /*
        std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
        for (int
            i = 0; i < n; i++) {
            std::cout << "x[" << i << "] = " << x[i] << std::endl;
        }

        std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
        for (int
            i = 0; i < n; i++) {
            std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
        }
        for (int
            i = 0; i < n; i++) {
            std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
        }

        std::cout << std::endl << std::endl << "Objective value" << std::endl;
        std::cout << "f(x*) = " << obj_value << std::endl;

        std::cout << std::endl << "Final value of the constraints:" << std::endl;
        for (int i = 0; i < m; i++) {
            std::cout << "g(" << i << ") = " << g[i] << std::endl;
        }
        */

        optimalValue = obj_value;

        for (int i = 0; i < n; i++) {
            solution[i] = x[i];
        }

    }


  // Without integer variables this reduces to clamping the point to the box,
  // which is the branch the original took for continuous problems.
  bool LoupFinderIpopt::integer_check(Vector& pt)
           {bound_check(normsys,pt); return true;}
  bool LoupFinderIpopt::is_inner(Vector& pt)
           {return is_inner0(normsys,pt);}
  double LoupFinderIpopt::goal_ub(Vector& pt)
           {return goal_ub0(normsys,pt);}
  void LoupFinderIpopt::sysbound(Vector& pt) {bound_check(normsys,pt);}
    void LoupFinderIpopt::sysbound(IntervalVector& vec) {bound_check(normsys,vec);}
  // OK only for QP problems (quadratic objective and linear constraints)
  void LoupFinderIpopt::set_quadratic(bool quadra){
    if (quadra) app->Options()->SetStringValue("hessian_constant", "yes");
  }
  
  void LoupFinderIpopt::correct_ipopt_sol (Vector&v, double& loup){
    if (recursive_call){
      
      recursive_call=false;
      IntervalVector box = sys.box;
      double eps=expansion_precisionB;
      if (ipopt_expansion > 0 && !corr_first_try) eps = ipopt_expansion;
      //      if (sys.minlp) eps=0.01;
      IntervalVector boxsol(v.size());
      for ( int i=0; i< v.size() ; i++){
	double epsi = eps;
	if (fabs(v[i])>1) epsi= eps*fabs(v[i]);
	boxsol[i]= sys.box[i] & Interval(v[i]- epsi, v[i]+ epsi);
      }
      CellHeap buffer(extsys);

      Optimizer opt(sys.nb_var,optimizer->ctc,optimizer->bsc,optimizer->loup_finder,buffer,extsys.goal_var(),optimizer->eps_x[0],optimizer->rel_eps_f, optimizer->abs_eps_f);


      opt.set_uplo(optimizer->get_uplo());
      opt.set_loup(optimizer->get_loup());
      opt.timeout = corr_first_try ? ipopt_expansion_first : 1;
      //      cout << " boxsol " << boxsol << endl;
      // The nested search shares the enclosing one's bisector. A stateless
      // bisector does not mind; lsmear-guard would count the nested decisions
      // in its window and horizon, and since this search is cut by a CPU
      // timeout, how many it adds depends on the machine. Its state belongs
      // to the enclosing search: put it back.
      BscHijackGuard* guard = dynamic_cast<BscHijackGuard*>(&optimizer->bsc);
      BscHijackGuard::State guard_state;
      if (guard) guard_state = guard->get_state();
      try {
        opt.optimize(boxsol);
      } catch (...) {
        if (guard) guard->set_state(guard_state);
        throw;
      }
      if (guard) guard->set_state(guard_state);
      recursive_call=true;
      log_corr_status = (int) opt.get_status(); log_corr_nodes = opt.get_nb_cells();
      log_corr_time = opt.get_time(); log_corr_loup = opt.get_loup();
      correction_nodes+=opt.get_nb_cells();
      correction_time+=opt.get_time();
      if (optimizer->trace)
	cout << " correction nodes " << correction_nodes << " correction_time " << correction_time << endl;
      if (opt.get_loup() < optimizer->get_loup()){
	if (optimizer->trace)
	  cout << "new loup after correction " << opt.get_loup() << endl;
	loup= opt.get_loup();
	v = opt.get_loup_point().mid();
      }
    
    }
}


}
