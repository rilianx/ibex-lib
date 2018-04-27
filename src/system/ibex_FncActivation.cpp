/* ============================================================================
 * I B E X - ibex_FncActivation.cpp
 * ============================================================================
 * Copyright   : IMT Atlantique (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Gilles Chabert
 * Created     : Jun, 16 2017
 * ---------------------------------------------------------------------------- */

#include "ibex_FncActivation.h"

using namespace std;

namespace ibex {

FncActivation::FncActivation(const System& sys, const Vector& pt, double activation_threshold, bool trace) : sys(sys), _activated(sys.f_ctrs.image_dim()) {

	BitSet ith=BitSet::empty(sys.f_ctrs.image_dim());

	for (int i=0; i<sys.f_ctrs.image_dim(); i++) {
		if (sys.ops[i]==EQ) _activated.add(i);
		else {
			ith.add(i);

			// call to eval_vector to avoid component generation
			if(sys.ops[i]==LEQ || sys.ops[i]==LT){
				if (sys.f_ctrs.eval_vector(pt,ith)[0].lb() >= -activation_threshold) {
					if (trace) cout << " activate inequality n°" << i << endl;
					_activated.add(i);
				}
			}else if(sys.ops[i]==GEQ || sys.ops[i]==GT){
				if (sys.f_ctrs.eval_vector(pt,ith)[0].ub() <= activation_threshold) {
					if (trace) cout << " activate inequality n°" << i << endl;
					_activated.add(i);
				}
			}

			ith.remove(i);
		}
	}

	(int &) _nb_var = sys.nb_var;
	(int &) _image_dim = _activated.size();
}


FncActivation::~FncActivation() {

}

IntervalVector FncActivation::eval_vector(const IntervalVector& x, const BitSet& components) const {
	return sys.f_ctrs.eval_vector(x, _activated.compose(components));
}

void FncActivation::jacobian(const IntervalVector& x, IntervalMatrix& J, const BitSet& components, int v) const {
	sys.f_ctrs.jacobian(x, J, _activated.compose(components),v);
}

} /* namespace ibex */
