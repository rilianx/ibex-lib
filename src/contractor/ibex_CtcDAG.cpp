/*
 * ibex_CtcDAG.cpp
 *
 *  Created on: 05-10-2017
 *      Author: victor
 */

#include "ibex_CtcDAG.h"

using namespace std;

namespace ibex {

	ExprDomain& Eval2::getDomains(){
		return d;
	};

	void Eval2::idx_cp_fwd(int x, int y) {
		assert(dynamic_cast<const ExprIndex*> (&f.node(y)));

		const ExprIndex& e = (const ExprIndex&) f.node(y);

		d[y] = d[x][e.index];
	}

	bool HC4Revise2::proj(const Domain& y, IntervalVector& x) {
		IntervalVector x2(x);

		d.read_arg_domains(x2);
		if(x.is_subset(x2))	eval2.eval(x);
		else eval.eval(x);



		//return false;
		//std::cout << "forward:" << std::endl; f.cf.print(d);

		bool is_inner=false;

		try {
			is_inner = backward(y);

			d.read_arg_domains(x);

			return is_inner;

		} catch(EmptyBoxException&) {
			x.set_empty();
			return false;
		}

		return false;
	}

	CtcDag::CtcDag(Function& f, Array<NumConstraint>& ctrs) : Ctc(f.nb_var()), f(f),
				y(NumConstraint(f,EQ).right_hand_side()), eval(f),
				hc4r(eval), d(eval.getDomains()){

			for(int i=0; i<y.dim.size(); i++) y[i]=ctrs[i].right_hand_side();
	}

	void CtcDag::contract(IntervalVector& box) {
		assert(box.size()==f.nb_var());
		hc4r.proj(y,box);
	}

	Function& CtcDag::get_f() {
		return f;;
	}

	ExprDomain& CtcDag::getDomains() {
		return d;
	}

}

