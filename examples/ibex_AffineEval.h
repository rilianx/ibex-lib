/* ============================================================================
 * I B E X - Function basic evaluation
 * ============================================================================
 * Copyright   : IMT Atlantique (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Ignacio Araya
 * Created     : 
 * Last update : 
 * ---------------------------------------------------------------------------- */

//#ifndef __IBEX_EVAL_H__
//#define __IBEX_EVAL_H__

#include <iostream>
#include "ibex.h"
#include "ibex_ExprDomain.h"
#include "ibex_Function.h"
#include "ibex_Interval.h"
#include "affine_arithm.h"
#include <vector>
#include <typeinfo>

using namespace std;

namespace ibex {

class Function;

/**
 * \ingroup symbolic
 *
 * \brief Function evaluator.
 *
 */
class AffineEval : public FwdAlgorithm {

public:
	/**
	 * \brief Build the evaluator for the function f.
	 */
	AffineEval(Function &f, int n);

	/**
	 * \brief Delete this.
	 */
	~AffineEval();

	/**
	 * \brief Run the forward algorithm with input domains.
	 */
	Domain& eval(const Array<const Domain>& d);

	/**
	 * \brief Run the forward algorithm with input domains.
	 */
	Domain& eval(const Array<Domain>& d);

	/**
	 * \brief Run the forward algorithm with an input box.
	 */
	Domain& eval(const IntervalVector& box);

protected:
	/**
	 * Class used internally to interrupt the forward procedure
	 * when an empty domain occurs (<=> the input box is outside
	 * the definition domain of the function).
	 */
	class EmptyBoxException { };

public: // because called from CompiledFunction

	       void vector_fwd (int* x, int y);
	       void apply_fwd  (int* x, int y);
	inline void idx_fwd    (int x, int y);
	inline void idx_cp_fwd (int x, int y);
	inline void symbol_fwd (int y);
	inline void cst_fwd    (int y);
	inline void chi_fwd    (int x1, int x2, int x3, int y);
	inline void gen2_fwd   (int x, int x2, int y);
	inline void add_fwd    (int x1, int x2, int y);
	inline void mul_fwd    (int x1, int x2, int y);
	inline void sub_fwd    (int x1, int x2, int y);
	inline void div_fwd    (int x1, int x2, int y);
	inline void max_fwd    (int x1, int x2, int y);
	inline void min_fwd    (int x1, int x2, int y);
	inline void atan2_fwd  (int x1, int x2, int y);
	inline void gen1_fwd   (int x, int y);
	inline void minus_fwd  (int x, int y);
	inline void minus_V_fwd(int x, int y);
	inline void minus_M_fwd(int x, int y);
	inline void trans_V_fwd(int x, int y);
	inline void trans_M_fwd(int x, int y);
	inline void sign_fwd   (int x, int y);
	inline void abs_fwd    (int x, int y);
	inline void power_fwd  (int x, int y, int p);
	inline void sqr_fwd    (int x, int y);
	inline void sqrt_fwd   (int x, int y);
	inline void exp_fwd    (int x, int y);
	inline void log_fwd    (int x, int y);
	inline void cos_fwd    (int x, int y);
	inline void sin_fwd    (int x, int y);
	inline void tan_fwd    (int x, int y);
	inline void cosh_fwd   (int x, int y);
	inline void sinh_fwd   (int x, int y);
	inline void tanh_fwd   (int x, int y);
	inline void acos_fwd   (int x, int y);
	inline void asin_fwd   (int x, int y);
	inline void atan_fwd   (int x, int y);
	inline void acosh_fwd  (int x, int y);
	inline void asinh_fwd  (int x, int y);
	inline void atanh_fwd  (int x, int y);
	inline void floor_fwd  (int x, int y);
	inline void ceil_fwd   (int x, int y);
	inline void saw_fwd    (int x, int y);
	inline void add_V_fwd  (int x1, int x2, int y);
	inline void add_M_fwd  (int x1, int x2, int y);
	inline void mul_SV_fwd (int x1, int x2, int y);
	inline void mul_SM_fwd (int x1, int x2, int y);
	inline void mul_VV_fwd (int x1, int x2, int y);
	inline void mul_MV_fwd (int x1, int x2, int y);
	inline void mul_VM_fwd (int x1, int x2, int y);
	inline void mul_MM_fwd (int x1, int x2, int y);
	inline void sub_V_fwd  (int x1, int x2, int y);
	inline void sub_M_fwd  (int x1, int x2, int y);

	Function& f;
	ExprDomain d;
	vector<Affine> af;
};

/* ============================================================================
 	 	 	 	 	 	 	 implementation
  ============================================================================*/

inline void AffineEval::idx_fwd(int, int) { /* nothing to do */ }

inline void AffineEval::symbol_fwd(int y) { 

	for (int i=0; i<d.args.size(); i++){
		if(&d[y] == &d.args[i]){
			af[y].a[i] = 1.0;
			af[y].ev = af[y].box->operator[](i);
			af[y].is_constant = false; 
		}else 
			af[y].a[i] = 0.0; 
	}
	af[y].err = Interval(0,0);

 }

inline void AffineEval::cst_fwd(int y) {
	const ExprConstant& c = (const ExprConstant&) f.node(y);
	switch (c.type()) {
	case Dim::SCALAR:       d[y].i() = c.get_value(); 
							af[y].err =  c.get_value(); 
							af[y].is_constant = true;       
							break;
	case Dim::ROW_VECTOR:
	case Dim::COL_VECTOR:   d[y].v() = c.get_vector_value();  break;
	case Dim::MATRIX:       d[y].m() = c.get_matrix_value();  break;
	}
}

inline void AffineEval::chi_fwd(int x1, int x2, int x3, int y) { d[y].i() = chi(d[x1].i(),d[x2].i(),d[x3].i()); }
inline void AffineEval::add_fwd(int x1, int x2, int y)   { 
	// interval arithmetic
	d[y].i()=d[x1].i()+d[x2].i(); 
	af[y] = af[x1] + af[x2];

}
inline void AffineEval::mul_fwd(int x1, int x2, int y)   { 
	d[y].i()=d[x1].i()*d[x2].i(); 
	af[y] = af[x1]*af[x2];
}

inline void AffineEval::sub_fwd(int x1, int x2, int y)   { 
	d[y].i()=d[x1].i()-d[x2].i();
	af[y] = af[x1] - af[x2];
}

inline void AffineEval::div_fwd(int x1, int x2, int y)   { d[y].i()=d[x1].i()/d[x2].i(); }
inline void AffineEval::max_fwd(int x1, int x2, int y)   { d[y].i()=max(d[x1].i(),d[x2].i()); }
inline void AffineEval::min_fwd(int x1, int x2, int y)   { d[y].i()=min(d[x1].i(),d[x2].i()); }
inline void AffineEval::atan2_fwd(int x1, int x2, int y) { d[y].i()=atan2(d[x1].i(),d[x2].i()); }

inline void AffineEval::minus_fwd(int x, int y)          { 
	d[y].i()=-d[x].i();
	 
	af[y] = -af[x];
	
}
inline void AffineEval::minus_V_fwd(int x, int y)        { d[y].v()=-d[x].v(); }
inline void AffineEval::minus_M_fwd(int x, int y)        { d[y].m()=-d[x].m(); }
inline void AffineEval::sign_fwd(int x, int y)           { d[y].i()=sign(d[x].i()); }
inline void AffineEval::abs_fwd(int x, int y)            { d[y].i()=abs(d[x].i()); }
inline void AffineEval::power_fwd(int x, int y, int p)   { d[y].i()=pow(d[x].i(),p); }
inline void AffineEval::sqr_fwd(int x, int y)            { 
	d[y].i()=sqr(d[x].i()); 	
	af[y] = sqr(af[x]);
}
inline void AffineEval::sqrt_fwd(int x, int y)           { if ((d[y].i()=sqrt(d[x].i())).is_empty()) throw EmptyBoxException(); }
inline void AffineEval::exp_fwd(int x, int y)            { d[y].i()=exp(d[x].i()); }
inline void AffineEval::log_fwd(int x, int y)            { if ((d[y].i()=log(d[x].i())).is_empty()) throw EmptyBoxException(); }
inline void AffineEval::cos_fwd(int x, int y)            { d[y].i()=cos(d[x].i()); }
inline void AffineEval::sin_fwd(int x, int y)            { d[y].i()=sin(d[x].i()); }
inline void AffineEval::tan_fwd(int x, int y)            { if ((d[y].i()=tan(d[x].i())).is_empty()) throw EmptyBoxException(); }
inline void AffineEval::cosh_fwd(int x, int y)           { d[y].i()=cosh(d[x].i()); }
inline void AffineEval::sinh_fwd(int x, int y)           { d[y].i()=sinh(d[x].i()); }
inline void AffineEval::tanh_fwd(int x, int y)           { d[y].i()=tanh(d[x].i()); }
inline void AffineEval::acos_fwd(int x, int y)           { if ((d[y].i()=acos(d[x].i())).is_empty()) throw EmptyBoxException(); }
inline void AffineEval::asin_fwd(int x, int y)           { if ((d[y].i()=asin(d[x].i())).is_empty()) throw EmptyBoxException(); }
inline void AffineEval::atan_fwd(int x, int y)           { d[y].i()=atan(d[x].i()); }
inline void AffineEval::acosh_fwd(int x, int y)          { if ((d[y].i()=acosh(d[x].i())).is_empty()) throw EmptyBoxException(); }
inline void AffineEval::asinh_fwd(int x, int y)          { d[y].i()=asinh(d[x].i()); }
inline void AffineEval::atanh_fwd(int x, int y)          { if ((d[y].i()=atanh(d[x].i())).is_empty()) throw EmptyBoxException(); }
inline void AffineEval::floor_fwd(int x, int y)          { if ((d[y].i()=floor(d[x].i())).is_empty()) throw EmptyBoxException(); }
inline void AffineEval::ceil_fwd(int x, int y)           { if ((d[y].i()=ceil(d[x].i())).is_empty()) throw EmptyBoxException(); }
inline void AffineEval::saw_fwd(int x, int y)            { if ((d[y].i()=saw(d[x].i())).is_empty()) throw EmptyBoxException(); }

inline void AffineEval::trans_V_fwd(int x, int y)        { d[y].v()=d[x].v(); }
inline void AffineEval::trans_M_fwd(int x, int y)        { d[y].m()=d[x].m().transpose(); }
inline void AffineEval::add_V_fwd(int x1, int x2, int y) { d[y].v()=d[x1].v()+d[x2].v(); }
inline void AffineEval::add_M_fwd(int x1, int x2, int y) { d[y].m()=d[x1].m()+d[x2].m(); }
inline void AffineEval::mul_SV_fwd(int x1, int x2, int y){ d[y].v()=d[x1].i()*d[x2].v(); }
inline void AffineEval::mul_SM_fwd(int x1, int x2, int y){ d[y].m()=d[x1].i()*d[x2].m(); }
inline void AffineEval::mul_VV_fwd(int x1, int x2, int y){ d[y].i()=d[x1].v()*d[x2].v(); }
inline void AffineEval::mul_MV_fwd(int x1, int x2, int y){ d[y].v()=d[x1].m()*d[x2].v(); }
inline void AffineEval::mul_VM_fwd(int x1, int x2, int y){ d[y].v()=d[x1].v()*d[x2].m(); }
inline void AffineEval::mul_MM_fwd(int x1, int x2, int y){ d[y].m()=d[x1].m()*d[x2].m(); }
inline void AffineEval::sub_V_fwd(int x1, int x2, int y) { d[y].v()=d[x1].v()-d[x2].v(); }
inline void AffineEval::sub_M_fwd(int x1, int x2, int y) { d[y].m()=d[x1].m()-d[x2].m(); }


AffineEval::AffineEval(Function& f, int n) : f(f), d(f) {
	af.clear();
	for (int i=0; i<d.data.size(); i++)
		af.push_back(Affine(n));
	
}

AffineEval::~AffineEval() {
}

Domain& AffineEval::eval(const Array<const Domain>& d2) {

	d.write_arg_domains(d2);

	//------------- for debug
	//	cout << "Function " << f.name << ", domains before eval:" << endl;
	//	for (int i=0; i<f.nb_arg(); i++) {
	//		cout << "arg[" << i << "]=" << f.arg_domains[i] << endl;
	//	}

	try {
		f.forward<AffineEval>(*this);
	} catch(EmptyBoxException&) {
		d.top->set_empty();
	}
	return *d.top;
}

Domain& AffineEval::eval(const Array<Domain>& d2) {

	d.write_arg_domains(d2);

	try {
		f.forward<AffineEval>(*this);
	} catch(EmptyBoxException&) {
		d.top->set_empty();
	}
	return *d.top;
}

Domain& AffineEval::eval(const IntervalVector& box) {

	d.write_arg_domains(box);

	for (int i=0; i<d.data.size(); i++)	af[i].box = &box;
	
	try {
		f.forward<AffineEval>(*this);
	} catch(EmptyBoxException&) {
		d.top->set_empty();
	}

	cout << af[0] << endl;
	d[0].i() &= af[0].ev;
	return *d.top;
}


void AffineEval::idx_cp_fwd(int x, int y) {
	assert(dynamic_cast<const ExprIndex*> (&f.node(y)));

	const ExprIndex& e = (const ExprIndex&) f.node(y);

	d[y] = d[x][e.index];

}

void AffineEval::apply_fwd(int* x, int y) {
	assert(dynamic_cast<const ExprApply*> (&f.node(y)));

	const ExprApply& a = (const ExprApply&) f.node(y);

	assert(&a.func!=&f); // recursive calls not allowed

	Array<const Domain> d2(a.func.nb_arg());

	for (int i=0; i<a.func.nb_arg(); i++) {
		d2.set_ref(i,d[x[i]]);
	}

	d[y] = a.func.basic_evaluator().eval(d2);
}

void AffineEval::vector_fwd(int* x, int y) {
	//not implemented
	throw std::invalid_argument("vector_fwd not implemented");

}

void AffineEval::gen1_fwd(int x, int y) {
	assert(dynamic_cast<const ExprGenericUnaryOp*>(&(f.node(y))));

	const ExprGenericUnaryOp& e = (const ExprGenericUnaryOp&) f.node(y);
	d[y]=e.eval(d[x]);
}

void AffineEval::gen2_fwd(int x1, int x2, int y) {
	assert(dynamic_cast<const ExprGenericBinaryOp*>(&(f.node(y))));

	const ExprGenericBinaryOp& e = (const ExprGenericBinaryOp&) f.node(y);
	d[y]=e.eval(d[x1],d[x2]);
}


} // namespace ibex

//#endif // __IBEX_EVAL_H__