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
			/*
			 * old_exp <-- mid(x)
			 * (exp,y) <-- abstaylor(old_exp)
			 * if y:
			 *   step <-- (old_exp,exp) #maximo step relativo
			 *   while step>0.1:
			 *      x <-- generate(exp,step*size(x),ratio)
			 *      old_exp <-- exp
			 *      (exp,y) <-- abstaylor(old_exp)
			 *      step <-- (old_exp,exp) #maximo step relativo
			 *
			 */

			Vector old_exp = box.mid();
			pair<IntervalVector,double> p2=finder_abs_taylor.find(box,old_exp,p.second);
			if (p2.second < p.second){ //si es mas pequeño o si encuentra cualquiera?
				double step = 0.0;
				for (int i = 0 ; i < old_exp.size() ; i++)
					if (std::abs((old_exp[i]-p2.first.mid()[i])/old_exp[i]) > step)
						step = std::abs((old_exp[i]-p2.first.mid()[i])/old_exp[i]);
				IntervalVector box_aux(box.size());
				while(step > 0.1){
					old_exp = p2.first.mid();
					for (int i = 0 ; i < box.size() ; i++){
						if ((old_exp[i]+step*box[i].lb()>=initial_box[i].lb()) && (old_exp[i]+step*box[i].ub()<=initial_box[i].ub()))
							box_aux[i] = Interval(old_exp[i]+step*box[i].lb(),old_exp[i]+step*box[i].ub());
						else if ((old_exp[i]+step*box[i].lb()>=initial_box[i].lb()) && (old_exp[i]+step*box[i].ub()>initial_box[i].ub()))
							box_aux[i] = Interval(old_exp[i]+step*box[i].lb(),initial_box[i].ub());
						else if ((old_exp[i]+step*box[i].lb()<initial_box[i].lb()) && (old_exp[i]+step*box[i].ub()<=initial_box[i].ub()))
							box_aux[i] = Interval(initial_box[i].lb(),old_exp[i]+step*box[i].ub());
						else if ((old_exp[i]+step*box[i].lb()<initial_box[i].lb()) && (old_exp[i]+step*box[i].ub()>initial_box[i].ub()))
							box_aux[i] = Interval(initial_box[i].lb(),initial_box[i].ub());
					}
					p2 = finder_abs_taylor.find(box_aux,old_exp,p.second);
					for (int i = 0 ; i < old_exp.size() ; i++)
						if (std::abs((old_exp[i]-p2.first.mid()[i])/old_exp[i]) > step)
							step = std::abs((old_exp[i]-p2.first.mid()[i])/old_exp[i]);
				}
			}
			if(p2.second < p.second){
				p=p2;
				found=true;
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
