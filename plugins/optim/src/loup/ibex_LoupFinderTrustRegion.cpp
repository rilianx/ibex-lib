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


LoupFinderTrustRegion::LoupFinderTrustRegion(const System& sys) :
	finder_abs_taylor(sys,true) {

}

std::pair<IntervalVector, double> LoupFinderTrustRegion::find(const IntervalVector& box, const IntervalVector& old_loup_point, double old_loup) {

	pair<IntervalVector,double> p=make_pair(old_loup_point, old_loup);

	bool found=false;

		try {



			pair<IntervalVector,double> p2=finder_abs_taylor.find(box,box.mid(),p.second);
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
