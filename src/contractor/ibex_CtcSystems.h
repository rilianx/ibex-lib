/*
 * ibex_CtcSystems.h
 *
 *  Created on: 13-10-2017
 *      Author: victor
 */

#ifndef SRC_CONTRACTOR_IBEX_CTCSYSTEMS_H_
#define SRC_CONTRACTOR_IBEX_CTCSYSTEMS_H_


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


}


#endif /* SRC_CONTRACTOR_IBEX_CTCSYSTEMS_H_ */
