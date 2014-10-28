//                                  I B E X                                   
// File        : ibex_Optimizer.cpp
// Author      : Gilles Chabert, Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 14, 2012
// Last Update : December 24, 2012
//============================================================================

#include "ibex_Optimizer-BS.h"
#include "ibex_EmptyBoxException.h"
#include "ibex_EmptySystemException.h"
#include "ibex_Timer.h"
//#include "ibex_OptimProbing.cpp_"
//#include "ibex_OptimSimplex.cpp_"
#include "ibex_CtcFwdBwd.h"
#include "ibex_CtcOptimShaving.h"
#include "ibex_CtcHC4.h"

#include "ibex_ExprCopy.h"
#include "ibex_Function.h"
#include "ibex_NoBisectableVariableException.h"
//#include "ibex_Multipliers.h"
#include "ibex_PdcFirstOrder.h"

#include <float.h>
#include <stdlib.h>

#include <tr1/unordered_set>

using namespace std;
using namespace std::tr1;
namespace ibex {

//~ const double OptimizerBS::default_prec = 1e-07;
//~ const double OptimizerBS::default_goal_rel_prec = 1e-07;
//~ const double OptimizerBS::default_goal_abs_prec = 1e-07;
//~ const int    OptimizerBS::default_sample_size = 10;
//~ const double OptimizerBS::default_equ_eps = 1e-08;
//~ const double OptimizerBS::default_loup_tolerance = 0.1;





OptimizerBS::OptimizerBS(System& user_sys, Ctc& ctc, Bsc& bsc, double prec,
					 double goal_rel_prec, double goal_abs_prec, int sample_size, double equ_eps,
					 bool rigor, int N, bool exploitation) : Optimizer(user_sys, ctc, bsc, prec,
					 goal_rel_prec, goal_abs_prec, sample_size, equ_eps,
					 rigor,0), N(N), exploitation(exploitation)
				 {

}


int iter=0;
Optimizer::Status OptimizerBS::optimize(const IntervalVector& init_box, double obj_init_bound) {

	loup=obj_init_bound;
	double old_loup=loup;
	pseudo_loup=obj_init_bound;

	uplo=NEG_INFINITY;
	uplo_of_epsboxes=POS_INFINITY;

	nb_cells=0;
	nb_simplex=0;
	diam_simplex=0;
	nb_rand=0;
	diam_rand=0;

    buffer.flush();
	
	OptimCell* root=new OptimCell(IntervalVector(n+1));

	write_ext_box(init_box,root->box);

	// add data required by the bisector
	bsc.add_backtrackable(*root);

	// add data required by optimizer + Fritz John contractor
	root->add<EntailedCtr>();
	//root->add<Multipliers>();
	entailed=&root->get<EntailedCtr>();
	entailed->init_root(user_sys,sys);

	loup_changed=false;
	loup_point=init_box.mid();
	time=0;
	Timer::start();
	handle_cell(*root,init_box);

	
	double ymax= POS_INFINITY;
	
	try {	
		while (!buffer.empty()) {
						
		  iter++;

		  //~ update_uplo();

		  if (buffer.empty()) 
			  break;
		 
		  
			loup_changed=false;
			
			set<OptimCell*,minLB> SS;

		    SS.insert(buffer.pop());

	        int depth=0;

            set<OptimCell*,minLB> S;
 
        while(!SS.empty()){

         S=SS;
         SS.clear();
			depth++;
		 
		 
         while(!S.empty()){


			OptimCell *c=(*S.begin());
              //cout  << c->id << "," << c->lb << ";" << c->ub << " depth:" << depth << endl;
			try {

				pair<IntervalVector,IntervalVector> boxes=bsc.bisect(*c);

				pair<OptimCell*,OptimCell*> new_cells=c->bisect(boxes.first,boxes.second);
				
				int v=new_cells.first->get<BisectedVar>().var;

				//buffer.erase(S.front());
				S.erase(S.begin());

               // int loup_subtree=c->loup_subtree;
				delete c; 
				
			    

               for(int i=0;i<2;i++){
				 try{
				    c=(i==0)? new_cells.first:new_cells.second;
				    contract_and_bound(*c, init_box);
				    c->lb=c->box[c->box.size()-1].lb();
	                c->id=nb_cells; //mid_depth
	                
	                 //~ c->id=(rand()%100)*1000000 + nb_cells;
	                //c->id=nb_cells;
	                //~ cout << c->id << endl;
	               // c->loup_subtree=(loup_changed)? nb_cells+1:loup_subtree;
 
                    nb_cells++;
			     }catch(EmptyBoxException&) {
                    delete c;
                    if(i==0) new_cells.first=NULL; 
                    else new_cells.second=NULL;
                 }
			    }
			    
			    
			    if(!exploitation){
					if(new_cells.first){ buffer.push(new_cells.first);}
				    if(new_cells.second){  buffer.push(new_cells.second);}		
				}else{
					if(new_cells.first) SS.insert(new_cells.first);
				    if(new_cells.second) SS.insert(new_cells.second); 
				}
						  
				if (uplo_of_epsboxes == NEG_INFINITY) {
					cout << " possible infinite minimum " << endl;
					goto end;
					//~ break;
				}
				
				
				
				time_limit_check();

               
 
			}
			catch (NoBisectableVariableException& ) {
				update_uplo_of_epsboxes ((c->box)[ext_sys.goal_var()].lb());
                S.erase(S.begin());
				delete c;
			}

			
		}


        //cout << SS.size() << endl;
        if(SS.size()>N){
		  std::set<OptimCell *>::iterator it=SS.begin();
		  for(int i=0;i<N;i++,it++);
		  
          for (;it!=SS.end();){
			buffer.push(*it);
			SS.erase(it++);
		  }
	    }


		}
		
        if (loup<old_loup) {

					// In case of a new upper bound (loup_changed == true), all the boxes
					// with a lower bound greater than (loup - goal_prec) are removed and deleted.
					// Note: if contraction was before bisection, we could have the problem
					// that the current cell is removed by contract_heap. See comments in
					// older version of the code (before revision 284).
					ymax= compute_ymax();

                    buffer.contract_heap(ymax);

                    
					if (ymax <=NEG_INFINITY) {
					    cout << " infinite value for the minimum " << endl;
						goto end;
						break;
					}
					if (trace) cout << setprecision(12) << "ymax=" << ymax << " uplo= " <<  uplo<< endl;

					old_loup=loup; 
		}	

		
		if (uplo_of_epsboxes == NEG_INFINITY || ymax <=NEG_INFINITY) {cout << "uplo = neginfinity"<< endl; break;}
		update_uplo();

	  }
	}
	catch (TimeOutException& ) {
		return TIME_OUT;
	}
	
	end:

	Timer::stop();
	time+= Timer::VIRTUAL_TIMELAPSE();


	if (uplo_of_epsboxes == POS_INFINITY && (loup==POS_INFINITY || (loup==initial_loup && goal_abs_prec==0 && goal_rel_prec==0)))
		return INFEASIBLE;
	else if (loup==initial_loup)
		return NO_FEASIBLE_FOUND;
	else if (uplo_of_epsboxes == NEG_INFINITY)
		return UNBOUNDED_OBJ;
	else
		return SUCCESS;

}

		  


} // end namespace ibex

