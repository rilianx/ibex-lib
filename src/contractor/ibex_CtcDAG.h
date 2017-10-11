/*
 * ibex_CtcDAG.h
 *
 *  Created on: 05-10-2017
 *      Author: victor
 */

#ifndef SRC_CONTRACTOR_IBEX_CTCDAG_H_
#define SRC_CONTRACTOR_IBEX_CTCDAG_H_

#include "ibex_Ctc.h"
#include "ibex_Function.h"
#include "ibex_Eval.h"

using namespace std;

namespace ibex {

class Eval2 : public Eval{

public:
	Eval2(Function &f):Eval(f) {};

	Domain& eval(const IntervalVector& box) {

		d.write_arg_domains(box);

		try {
			f.forward<Eval2>(*this);
		} catch(EmptyBoxException&) {
			d.top->set_empty();
		}
		return *d.top;
	}

	ExprDomain& getDomains();
	void idx_cp_fwd(int x, int y);

	inline void chi_fwd(int x1, int x2, int x3, int y) { d[y].i() &= chi(d[x1].i(),d[x2].i(),d[x3].i()); }
	inline void add_fwd(int x1, int x2, int y)   {	d[y].i()&=d[x1].i()+d[x2].i(); }
	inline void mul_fwd(int x1, int x2, int y)   { d[y].i()&=d[x1].i()*d[x2].i(); }
	inline void sub_fwd(int x1, int x2, int y)   { d[y].i()&=d[x1].i()-d[x2].i(); }
	inline void div_fwd(int x1, int x2, int y)   { d[y].i()&=d[x1].i()/d[x2].i(); }
	inline void max_fwd(int x1, int x2, int y)   { d[y].i()&=max(d[x1].i(),d[x2].i()); }
	inline void min_fwd(int x1, int x2, int y)   { d[y].i()&=min(d[x1].i(),d[x2].i()); }
	inline void atan2_fwd(int x1, int x2, int y) { d[y].i()&=atan2(d[x1].i(),d[x2].i()); }

	inline void minus_fwd(int x, int y)          { d[y].i()&=-d[x].i(); }
	inline void minus_V_fwd(int x, int y)        { d[y].v()&=-d[x].v(); }
	inline void minus_M_fwd(int x, int y)        { d[y].m()&=-d[x].m(); }
	inline void sign_fwd(int x, int y)           { d[y].i()&=sign(d[x].i()); }
	inline void abs_fwd(int x, int y)            { d[y].i()&=abs(d[x].i()); }
	inline void power_fwd(int x, int y, int p)   { d[y].i()&=pow(d[x].i(),p); }
	inline void sqr_fwd(int x, int y)            { d[y].i()&=sqr(d[x].i()); }
	inline void sqrt_fwd(int x, int y)           { if ((d[y].i()&=sqrt(d[x].i())).is_empty()) throw EmptyBoxException(); }
	inline void exp_fwd(int x, int y)            { d[y].i()&=exp(d[x].i()); }
	inline void log_fwd(int x, int y)            { if ((d[y].i()&=log(d[x].i())).is_empty()) throw EmptyBoxException(); }
	inline void cos_fwd(int x, int y)            { d[y].i()&=cos(d[x].i()); }
	inline void sin_fwd(int x, int y)            { d[y].i()&=sin(d[x].i()); }
	inline void tan_fwd(int x, int y)            { if ((d[y].i()&=tan(d[x].i())).is_empty()) throw EmptyBoxException(); }
	inline void cosh_fwd(int x, int y)           { d[y].i()&=cosh(d[x].i()); }
	inline void sinh_fwd(int x, int y)           { d[y].i()&=sinh(d[x].i()); }
	inline void tanh_fwd(int x, int y)           { d[y].i()&=tanh(d[x].i()); }
	inline void acos_fwd(int x, int y)           { if ((d[y].i()&=acos(d[x].i())).is_empty()) throw EmptyBoxException(); }
	inline void asin_fwd(int x, int y)           { if ((d[y].i()&=asin(d[x].i())).is_empty()) throw EmptyBoxException(); }
	inline void atan_fwd(int x, int y)           { d[y].i()&=atan(d[x].i()); }
	inline void acosh_fwd(int x, int y)          { if ((d[y].i()&=acosh(d[x].i())).is_empty()) throw EmptyBoxException(); }
	inline void asinh_fwd(int x, int y)          { d[y].i()&=asinh(d[x].i()); }
	inline void atanh_fwd(int x, int y)          { if ((d[y].i()&=atanh(d[x].i())).is_empty()) throw EmptyBoxException(); }

	inline void trans_V_fwd(int x, int y)        { d[y].v()&=d[x].v(); }
	inline void trans_M_fwd(int x, int y)        { d[y].m()&=d[x].m().transpose(); }
	inline void add_V_fwd(int x1, int x2, int y) { d[y].v()&=d[x1].v()+d[x2].v(); }
	inline void add_M_fwd(int x1, int x2, int y) { d[y].m()&=d[x1].m()+d[x2].m(); }
	inline void mul_SV_fwd(int x1, int x2, int y){ d[y].v()&=d[x1].i()*d[x2].v(); }
	inline void mul_SM_fwd(int x1, int x2, int y){ d[y].m()&=d[x1].i()*d[x2].m(); }
	inline void mul_VV_fwd(int x1, int x2, int y){ d[y].i()&=d[x1].v()*d[x2].v(); }
	inline void mul_MV_fwd(int x1, int x2, int y){ d[y].v()&=d[x1].m()*d[x2].v(); }
	inline void mul_VM_fwd(int x1, int x2, int y){ d[y].v()&=d[x1].v()*d[x2].m(); }
	inline void mul_MM_fwd(int x1, int x2, int y){ d[y].m()&=d[x1].m()*d[x2].m(); }
	inline void sub_V_fwd(int x1, int x2, int y) { d[y].v()&=d[x1].v()-d[x2].v(); }
	inline void sub_M_fwd(int x1, int x2, int y) { d[y].m()&=d[x1].m()-d[x2].m(); }


};

class HC4Revise2 : public HC4Revise{

public:
	Eval2& eval2;

	HC4Revise2(Eval2& e) : eval2(e), HC4Revise(e) {

	}

	bool proj(const Domain& y, IntervalVector& x);

	bool backward(const Domain& y) {

		Domain& root=*d.top;

		if (root.is_empty())
			throw EmptyBoxException();

	    root &= y;

		if (root.is_empty())
			throw EmptyBoxException();

		// may throw an EmptyBoxException().
		eval.f.backward<HC4Revise2>(*this);

		return false;
		//std::cout << "backward:" << std::endl; f.cf.print();
	}



};

class CtcDag : public Ctc {
public:

	Eval2 eval;
	HC4Revise2 hc4r;

	Function& f;
	ExprDomain& d;

	/** The domain "y". */
	Domain y;


	CtcDag(Function& f, Array<NumConstraint>& ctrs);

	void contract(IntervalVector& box);

	Function& get_f();

	ExprDomain& getDomains();

};

} /* namespace ibex */



#endif /* SRC_CONTRACTOR_IBEX_CTCDAG_H_ */
