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
#include "ibex_PdcHansenFeasibility.h"
#include "ibex_FncActivation.h"
#include "ibex_NormalizedSystem.h"

using namespace std;

namespace ibex {

LoupFinderDefault::LoupFinderDefault(const System& sys, bool inHC4, bool safe_loup) :
		sys(sys), finder_probing(inHC4? (LoupFinder&) *new LoupFinderInHC4(sys) : (LoupFinder&) *new LoupFinderFwdBwd(sys)),
	finder_x_taylor(sys), safe_loup(safe_loup) {

}

namespace {

bool __is_unfeasible(const Interval& gx, CmpOp op) {
	bool unfeas;
	switch (op) {
	case LT:  unfeas=gx.lb()>=0; break;
	case LEQ: unfeas=gx.lb()>0; break;
	case EQ:  unfeas=(!gx.contains(0.0)); break;
	case GEQ: unfeas=gx.ub()<0; break;
	case GT:  unfeas=gx.ub()<=0; break;
	}
	return unfeas;
}

}

std::pair<IntervalVector, double> LoupFinderDefault::find(const IntervalVector& box, const IntervalVector& old_loup_point, double old_loup) {

	pair<IntervalVector,double> p=make_pair(old_loup_point, old_loup);

	bool found=false;

	try {
		p=finder_probing.find(box,p.first,p.second);
		found=true;
	} catch(NotFound&) { }

	try {
		// TODO
		// in_x_taylor.set_inactive_ctr(entailed->norm_entailed);

		//if the fourth parameter is false, the solution can be unfeasible
		p=finder_x_taylor.find(box, p.first, p.second, safe_loup);
		found=true;
	} catch(NotFound&) { }

	if (found){
		//we certify (or correct) the found solution
		/*if(real_sys){
			IntervalVector loup_point(p.first);
			if (p.second<old_loup) {

				IntervalVector epsbox(loup_point);
				epsbox.inflate(1e-6);

				IntervalVector res_ctr = sys.f_ctrs.eval_vector(epsbox);

				//TODO: restrict to the initial large box
				epsbox &= box;

				//cout << epsbox << endl;

				for (int c=0; c<sys.f_ctrs.image_dim(); c++)
					if (__is_unfeasible(res_ctr[c], sys.ops[c])) throw NotFound();

				//the inflated box epsbox cannot be discarded thus it may contain a real solution
				return certify(epsbox, old_loup);
			}
		}else*/
		  return p;
	}

	throw NotFound();
}

LoupFinderDefault::~LoupFinderDefault() {
	delete &finder_probing;
}

} /* namespace ibex */
