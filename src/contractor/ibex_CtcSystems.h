/*
 * ibex_CtcSystems.h
 *
 *  Created on: 13-10-2017
 *      Author: victor
 */

#ifndef SRC_CONTRACTOR_IBEX_CTCSYSTEMS_H_
#define SRC_CONTRACTOR_IBEX_CTCSYSTEMS_H_
#include <list>




#define CAST(m,n) dynamic_cast<const m*>(n)

using namespace std;

namespace ibex {

/**
 * \brief A contractor using a linear system A.x=b
 *
 * Linear system of the form A.x=b.
 *
 * If P exists, the contractor use it for contracting x by using the system PA.x=Pb
 * and a simple projection based method
 * If P does not exist the contractor perform at most 2n simplex for contracting x
 */

class LinearSystem : public Ctc, Linearizer {
public:


	enum {SIMPLEX, PSEUDOINVERSE, GAUSS_JORDAN, GAUSS_PSEUDOINV, MULT_GAUSS_JORDAN};

	CtcPolytopeHull ctc;
	int ctc_type;


	/**
	 * \brief Create the linear inequalities Ax<=b.
	 */
	LinearSystem(int nb_var, int ctc_type = SIMPLEX) :
		ctc(*this, LinearSolver::default_max_iter,
		LinearSolver::default_max_time_out, LinearSolver::default_eps, Interval (1e-14, 1e10)),
		ctc_type(ctc_type), Ctc(nb_var), A(1,1), b(1), P(1,1), PA(1,1), Pb(1) {

	}

	/**
	 * \brief Create the linear inequalities Ax<=b.
	 */
	LinearSystem(const IntervalMatrix& A, const IntervalMatrix& P, const IntervalMatrix& PA,  int ctc_type = GAUSS_JORDAN) :
		ctc(*this, LinearSolver::default_max_iter,
		LinearSolver::default_max_time_out, LinearSolver::default_eps, Interval (1e-14, 1e10)),
		ctc_type(ctc_type), Ctc(A.nb_cols()), A(A), b(1), P(P), PA(PA), Pb(1) {

	}


	void contract(IntervalVector& box);

	/**
	 * \brief Add the inequalities in the solver (does not perform a linearization)
	 *
	 * This method is required by the PolytopeHull
	 */
	virtual int linearization(const IntervalVector& x, LinearSolver& lp_solver);




protected:

	IntervalMatrix A;
	IntervalVector b;

	IntervalMatrix P;
	IntervalMatrix PA;
	IntervalVector Pb;

};


class EmbeddedLinearSystem : public LinearSystem {
public:

	/**
	 * All but SIMPLEX
	 */
	EmbeddedLinearSystem(const IntervalMatrix& A, const IntervalMatrix& P, const IntervalMatrix& PA,
			Array<const ExprNode>& xn, Array<const ExprNode>& bn, map<const ExprNode*, int>& node2i,
			ExprDomain& d, bool is_mult=false, bool extended=false);

	/**
	 * \brief Create the linear inequalities Ax<=b.
	 */
	EmbeddedLinearSystem(const IntervalMatrix& AA, Array<const ExprNode>& xn,
			Array<const ExprNode>& bn, map<const ExprNode*, int>& node2i,
			ExprDomain& d, bool is_mult=false, int ctc_type = SIMPLEX, bool extended=false);


	virtual void contract(IntervalVector& box);
	virtual int linearization(const IntervalVector& x, LinearSolver& lp_solver);

	bool is_mult;
	bool extended;


private:

	Array<const ExprNode> xn; // size:nb_cols
	Array<const ExprNode> bn; // size:nb_rows
	IntervalVector x;

	map<const ExprNode*, int> node2i;
	ExprDomain& d;

};


class EmbeddedLinearSystemBuilder {

public:

	EmbeddedLinearSystemBuilder(const IntervalMatrix& A, Array<const ExprNode>& xn,
			Array<const ExprNode>& bn) : A(A), xn(xn), bn(bn), PA(NULL), P(NULL), node2i(NULL), d(NULL), is_mult(false), ctc_type(-1), extended(false) { }

	void set_P(Matrix* P2);
	void set_PA(IntervalMatrix* PA2);

	void set_node2i(map<const ExprNode*, int>* n) { node2i=n; }
	void set_domain(ExprDomain* dd) { d=dd; }
	void set_is_mult(bool m) { is_mult=m; }
	void set_ctc_type(int t) { ctc_type=t; }
	void set_extended(bool e) { extended=e; }


	EmbeddedLinearSystem* create(){
		if(PA && P && node2i && d && ctc_type!=-1) return new EmbeddedLinearSystem(A, *P, *PA, xn, bn, *node2i, *d, is_mult, extended);
		else if(node2i && d && ctc_type!=-1) return new EmbeddedLinearSystem(A, xn, bn, *node2i, *d, is_mult, ctc_type, extended);
		return NULL;
	}

	inline Array<const ExprNode>& get_bn() { return bn;}
	inline Array<const ExprNode>& get_xn() { return xn;}


    inline IntervalMatrix& getA() { return A; }


private:
	IntervalMatrix A;
	IntervalMatrix* P;
	IntervalMatrix* PA;

	Array<const ExprNode> xn; // size:nb_cols
	Array<const ExprNode> bn; // size:nb_rows

	map<const ExprNode*, int>* node2i;
	ExprDomain* d;

	bool is_mult;
	int ctc_type;
	bool extended;
};

	void create_subsystems(list<EmbeddedLinearSystemBuilder *> &ls_list, IntervalMatrix& A, Array<const ExprNode> &x, Array<const ExprNode> &b,
		vector<pair <set <int>,set <int> > >& subsets);


	void find_subsystems(list<EmbeddedLinearSystemBuilder *> &ls_list, IntervalMatrix& A, Array<const ExprNode> &x,
		Array<const ExprNode> &b, int& nb_rows, int& nb_cols);



	template<typename T>
	void add_row(IntervalMatrix& A, Array<const ExprNode> &b, const T* bi, map<const ExprNode*, int>& xmap,
		int& nb_rows, int& nb_cols);

	template<typename T>
	Array<Ctc> getEmbeddedLinearSystems(CtcDag& dag_ctc, bool is_mult, int ctc_type, bool extended, int nb_nodes=-1);


}


#endif /* SRC_CONTRACTOR_IBEX_CTCSYSTEMS_H_ */
