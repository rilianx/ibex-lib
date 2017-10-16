/*
 * ibex_CtcSystems.cpp
 *
 *  Created on: 13-10-2017
 *      Author: victor
 */

#include "ibex_CtcSystems.h"

using namespace std;

namespace ibex {


	void LinearSystem::contract(IntervalVector& box){

		if(ctc_type==SIMPLEX){
			ctc.contract(box);
			if(box.is_empty()) return;
		}

		if(ctc_type==PSEUDOINVERSE || ctc_type==GAUSS_JORDAN || ctc_type == GAUSS_PSEUDOINV ){
			IntervalMatrix tmp(PA);
			if(!bwd_mul(Pb, tmp, box, 1e8)){
				box.set_empty();
				return;
			}
		}
	}
	virtual int LinearSystem::linearization(const IntervalVector& x, LinearSolver& lp_solver){
        //cout <<  LinearRelax::nb_var() << endl;
		int num=0;
		for (int i=0; i<A.nb_rows(); i++) {


			Interval error = (IntervalVector(A[i].diam())* IntervalVector(x.mag())).ub();
			Vector row=A[i].mid();
			try {
			    Interval eval=(row*x); // left side image
			   // cout << (bb - error).lb() << "<" << eval << "<" << (bb + error).ub() << endl;

				if(eval.lb() < (b[i] - error).lb()){
				   lp_solver.add_constraint(row,GEQ, (b[i] - error).lb());
				   num++;
				}
				if(eval.ub() > (b[i] + error).ub()){
					lp_solver.add_constraint(row,LEQ, (b[i] + error).ub());
					num++;

				}
			} catch (LPException&) { }
		}
		/*
		lp_solver.writeFile("out.txt");
		system("cat out.txt");

       cout << A << endl;
       cout << b << endl;
       cout << x << endl;
*/
		return num;
	}


}

