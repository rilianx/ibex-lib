//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderDefault.cpp
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Jul 09, 2017
//============================================================================

#include "ibex_LoupFinderDefault.h"
#include "ibex_LoupFinderInHC4.h"
#include "ibex_LoupFinderFwdBwd.h"

using namespace std;

namespace ibex {

string LoupFinderDefault::foundby="";

LoupFinderDefault::LoupFinderDefault(const System& sys,const IntervalVector& initial_box, bool inHC4, mode m) :
	finder_probing(inHC4? (LoupFinder&) *new LoupFinderInHC4(sys) : (LoupFinder&) *new LoupFinderFwdBwd(sys)),
	finder_abs_taylor(sys,true), finder_x_taylor(sys,false),initial_box(initial_box), finder_trustregion(sys,initial_box), m(m) {

}

std::pair<IntervalVector, double> LoupFinderDefault::find(const IntervalVector& box, const IntervalVector& old_loup_point, double old_loup) {

	pair<IntervalVector,double> p=make_pair(old_loup_point, old_loup);

	bool found=false;

	try {
		p=finder_probing.find(box,p.first,p.second);
		//cout << "inhc4:" << p.first.size() << endl;
		foundby="inhc4";
		found=true;
	} catch(NotFound&) { }

	if(m==xt || m==both)
		try {
			// TODO
			// in_x_taylor.set_inactive_ctr(entailed->norm_entailed);
			pair<IntervalVector,double> p2=finder_x_taylor.find(box,p.first,p.second);
			if(p2.second < p.second){
				p=p2;
				foundby="xtaylor";
				found=true;
			}
		} catch(NotFound&) { }

	if(m==abst || m==both)
		try {
			// TODO
			// in_x_taylor.set_inactive_ctr(entailed->norm_entailed);
			pair<IntervalVector,double> p2=finder_abs_taylor.find(box,box.mid(),p.second);
			if(p2.second < p.second){
				p=p2;
				foundby="abstaylor";
				found=true;
			}
		} catch(NotFound&) { }

		if(m==trustregion)
			try {
				pair<IntervalVector,double> p2=finder_trustregion.find(box,box.mid(),p.second);
				if(p2.second < p.second){
					p=p2;
					foundby="trust_region";
					found=true;
				}
			} catch(NotFound&) { }

	if (found)
		return p;
	else
		throw NotFound();
}

LoupFinderDefault::~LoupFinderDefault() {
	delete &finder_probing;
}

} /* namespace ibex */
