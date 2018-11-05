//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderDefault.cpp
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Jul 09, 2017
//============================================================================

#include "ibex_LoupFinderTrustRegion.h"

using namespace std;

namespace ibex {


LoupFinderTrustRegion::LoupFinderTrustRegion(const System& sys,const IntervalVector& initial_box,double alpha) :
	finder_abs_taylor(sys,true),finder_xt(sys,false),initial_box(initial_box),f_goal(sys.goal),alpha(alpha) {

}

std::pair<IntervalVector, double> LoupFinderTrustRegion::find(const IntervalVector& box, const IntervalVector& old_loup_point, double old_loup) {

	pair<IntervalVector,double> p=make_pair(old_loup_point, old_loup);

	bool found=false;
	pair<IntervalVector,double> new_ub = p;
	pair<IntervalVector,double> old_ub = p;
 	bool flag = true;
 	bool line_search;
 	double pp;

		try{
			pair<IntervalVector,double> new_ub=finder_abs_taylor.find(box,box.mid(),p.second);
			if(new_ub.second < p.second){
				found = true;
				p = new_ub;
			}
			else throw NotFound();
		} catch(NotFound&) { }
		try{
			pair<IntervalVector,double> new_ub=finder_xt.find(box,box.mid(),p.second);
			if(new_ub.second < p.second){
				found = true;
				p = new_ub;
			}
			else throw NotFound();
		} catch(NotFound&) { }

		if (!found) throw NotFound();

	 	IntervalVector box_aux(box.size());
	 	while ((old_ub.second-new_ub.second>1e-6) || (flag)){
	 		if (old_ub.second-new_ub.second < 1e-6) flag = false;
	 		else flag = true;
		 	line_search = true;
	 		Vector old_exp = new_ub.first.mid();
	 		for (int i = 0 ; i < box.size() ; i++){
	 			if ((std::abs(box[i].lb()-new_ub.first[i].lb())<1e-6) || (std::abs(box[i].ub()-new_ub.first[i].ub())<1e-6)){
	 				pp = 1;
	 			}
	 			else{
	 				pp = 1;
	 			}
				if ((old_exp[i]-(alpha/pp)*box[i].diam()>=initial_box[i].lb()) && (old_exp[i]+(alpha/pp)*box[i].diam()<=initial_box[i].ub()))
					box_aux[i] = Interval(old_exp[i]-(alpha/pp)*box[i].diam(),old_exp[i]+(alpha/pp)*box[i].diam());
				else if ((old_exp[i]-(alpha/pp)*box[i].diam()>=initial_box[i].lb()) && (old_exp[i]+(alpha/pp)*box[i].diam()>initial_box[i].ub()))
					box_aux[i] = Interval(old_exp[i]-(alpha/pp)*box[i].diam(),initial_box[i].ub());
				else if ((old_exp[i]-(alpha/pp)*box[i].diam()<initial_box[i].lb()) && (old_exp[i]+(alpha/pp)*box[i].diam()<=initial_box[i].ub())){
					box_aux[i] = Interval(initial_box[i].lb(),old_exp[i]+(alpha/pp)*box[i].diam());
					line_search=false;
				}
				else if ((old_exp[i]-alpha*box[i].diam()<initial_box[i].lb()) && (old_exp[i]+alpha*box[i].diam()>initial_box[i].ub())){
					box_aux[i] = Interval(initial_box[i].lb(),initial_box[i].ub());
					line_search=false;
				}
	 		}
	 		try {
	 			pair<IntervalVector,double> newnew_ub = make_pair(old_loup_point, old_loup);
 		 		old_ub = new_ub;

//	 			if (line_search){
//	 				newnew_ub=finder_abs_taylor.find(box_aux,box_aux.mid(),new_ub.second);
//		 			newnew_ub=finder_abs_taylor.find(box_aux,new_ub.first.mid(),new_ub.second);
//					pair<double, double> c1_t1 = make_pair(POS_INFINITY,0);
//					IntervalVector point_a = old_ub.first;
//					PFunction pf(*f_goal,  newnew_ub.first.mid(),new_ub.first.mid());
//					c1_t1 = pf.optimize(Interval(POS_INFINITY), PFunction::MIN, POS_INFINITY,  0.1);
//
//					if (f_goal->eval(point_a+c1_t1.second*(newnew_ub.first-point_a)).ub() < newnew_ub.second){
//						newnew_ub.first = point_a+c1_t1.second*(newnew_ub.first-point_a);
//						newnew_ub.second = f_goal->eval(old_loup_point).ub();
//					}
//	 			}
//	 			else
		 		newnew_ub=finder_abs_taylor.find(box_aux,box_aux.mid(),new_ub.second);

	 			if(newnew_ub.second < new_ub.second){
	 		 		new_ub = newnew_ub;
	 		 		p = new_ub;
	 		 	}
	 		}
	 		catch(NotFound&) {}
	 		try {
				pair<IntervalVector,double> newnew_ub = make_pair(old_loup_point, old_loup);
				old_ub = new_ub;

//	 			if (line_search){
//	 				newnew_ub=finder_abs_taylor.find(box_aux,box_aux.mid(),new_ub.second);
//		 			newnew_ub=finder_abs_taylor.find(box_aux,new_ub.first.mid(),new_ub.second);
//					pair<double, double> c1_t1 = make_pair(POS_INFINITY,0);
//					IntervalVector point_a = old_ub.first;
//					PFunction pf(*f_goal,  newnew_ub.first.mid(),new_ub.first.mid());
//					c1_t1 = pf.optimize(Interval(POS_INFINITY), PFunction::MIN, POS_INFINITY,  0.1);
//
//					if (f_goal->eval(point_a+c1_t1.second*(newnew_ub.first-point_a)).ub() < newnew_ub.second){
//						newnew_ub.first = point_a+c1_t1.second*(newnew_ub.first-point_a);
//						newnew_ub.second = f_goal->eval(old_loup_point).ub();
//					}
//	 			}
//	 			else
				newnew_ub=finder_xt.find(box_aux,box_aux.mid(),new_ub.second);

				if(newnew_ub.second < new_ub.second){
					new_ub = newnew_ub;
					p = new_ub;
				}
			}
			catch(NotFound&) {}
 			if (!flag)
 				break;
	 	}
	 	if (found)
	 		return p;
	 	else
	 		throw NotFound();
	 }

LoupFinderTrustRegion::~LoupFinderTrustRegion() {
	delete &finder_abs_taylor;
}

} /* namespace ibex */
