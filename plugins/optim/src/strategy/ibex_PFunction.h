/*
 * ibex_PFunction.h
 *
 *  Created on: 23 may. 2018
 *      Author: iaraya
 */

#include "ibex_Function.h"


#ifndef OPTIM_MOP_SRC_STRATEGY_IBEX_PFUNCTION_H_
#define OPTIM_MOP_SRC_STRATEGY_IBEX_PFUNCTION_H_


using namespace std;

namespace ibex {

/**
 * Parameterized function f(t) ← f1(xt) - m*f2(xt)
 * xt = xa + t*(xb-xa)
 */
class PFunction{

public:
	PFunction(const Function& f1, const IntervalVector& xa, const IntervalVector& xb);

	/**
	 * convert pf.t to t in inter
	 * t = inter.lb() + pf.t*(inter.ub() - inter.lb());
	 */

	static bool MIN;
	static bool MAX;


	void contract_curve(const Interval& t);

	Interval eval(const Interval& t, const Interval& m, bool minimize) const;

	Interval deriv(const Interval& t, const Interval& m, bool minimize) const;

	IntervalVector get_point(const Interval& t) const;
	static Interval eval_goal(const Function& goal, const IntervalVector& x, int n);
	static IntervalVector deriv_goal(const Function& goal, const IntervalVector& x, int n);
	/**
	 * \brief gets the image of the segment line xa-xb
	 */
	void get_curve_y(std::vector< pair <double, double> >& curve_y );

	/**
	 * \brief minimize/maximize the function pf: f1(t)+w*f2(t)
	 * returning the lb/ub of its evaluation (c) and the best solution found t and its
	 * input m, minimize, max_c=max_value
	 */
	pair<double,double> optimize(const Interval& m, bool minimize, double max_c, Interval init);

	bool newton_lcontract(const Interval& m, bool minimize, Interval& inter, const Interval& derivate, double lb);
	bool newton_rcontract(const Interval& m, bool minimize, Interval& inter, const Interval& derivate, double lb);

private:

	const Function& f1;
	IntervalVector xa;
	IntervalVector xb;

	static double _min_newton_step;
	static double _min_diam;
	static double _eps_opt;
};

} /* namespace ibex */

#endif /* OPTIM_MOP_SRC_STRATEGY_IBEX_PFUNCTION_H_ */
