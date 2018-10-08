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


LoupFinderTrustRegion::LoupFinderTrustRegion(const System& sys,const IntervalVector& initial_box) :
	finder_abs_taylor(sys,true),initial_box(initial_box) {

}

std::pair<IntervalVector, double> LoupFinderTrustRegion::find(const IntervalVector& box, const IntervalVector& old_loup_point, double old_loup) {

	pair<IntervalVector,double> p=make_pair(old_loup_point, old_loup);

	bool found=false;
		try {

		pair<IntervalVector,double> new_ub=finder_abs_taylor.find(box,box.mid(),p.second);
		pair<IntervalVector,double> old_ub = p;
	 	if(new_ub.second < p.second){
	 		found = true;
	 		p = new_ub;
	 	}
	 	else throw NotFound();
	 	double current_alpha;
	 	bool flag = true;
	 	IntervalVector box_aux(box.size());
	 	while ((old_ub.second-new_ub.second>1e-5) || (flag)){
	 		if (old_ub.second-new_ub.second < 1e-5) flag = false;
	 		else flag = true;

	 		Vector old_exp = new_ub.first.mid();
	 		for (int i = 0 ; i < box.size() ; i++){
	 			if ((std::abs(box[i].lb()-new_ub.first[i].lb())<1e-5) || (std::abs(box[i].ub()-new_ub.first[i].ub())<1e-5))
	 				current_alpha = 1;
	 			else{ current_alpha = 0.2;
	 			}
				if ((old_exp[i]-current_alpha*box[i].diam()>=initial_box[i].lb()) && (old_exp[i]+current_alpha*box[i].diam()<=initial_box[i].ub()))
					box_aux[i] = Interval(old_exp[i]-current_alpha*box[i].diam(),old_exp[i]+current_alpha*box[i].diam());
				else if ((old_exp[i]-current_alpha*box[i].diam()>=initial_box[i].lb()) && (old_exp[i]+current_alpha*box[i].diam()>initial_box[i].ub()))
					box_aux[i] = Interval(old_exp[i]-current_alpha*box[i].diam(),initial_box[i].ub());
				else if ((old_exp[i]-current_alpha*box[i].diam()<initial_box[i].lb()) && (old_exp[i]+current_alpha*box[i].diam()<=initial_box[i].ub()))
					box_aux[i] = Interval(initial_box[i].lb(),old_exp[i]+current_alpha*box[i].diam());
				else if ((old_exp[i]-current_alpha*box[i].diam()<initial_box[i].lb()) && (old_exp[i]+current_alpha*box[i].diam()>initial_box[i].ub()))
					box_aux[i] = Interval(initial_box[i].lb(),initial_box[i].ub());
			}
	 		try {
	 			pair<IntervalVector,double> newnew_ub = make_pair(old_loup_point, old_loup);
 		 		old_ub = new_ub;
	 			newnew_ub=finder_abs_taylor.find(box_aux,box_aux.mid(),new_ub.second);
	 		 	if(newnew_ub.second < new_ub.second){
	 		 		new_ub = newnew_ub;
	 		 		found = true;
	 		 		p = new_ub;
	 		 	}
	 		}
	 		catch(NotFound&) {
	 			if (!flag)
	 				return p;
	 		}
	 	}

	 	} catch(NotFound&) { }

	 	if (found)
	 		return p;
	 	else
	 		throw NotFound();
	 }

LoupFinderTrustRegion::~LoupFinderTrustRegion() {
	delete &finder_abs_taylor;
}

} /* namespace ibex */
