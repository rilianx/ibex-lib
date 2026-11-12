//============================================================================
//                                  I B E X
// File        : ibex_CtcSC.h
// Author      : Ignacio Araya
// License     : See the LICENSE file
// Created     : Nov 08, 2026
//============================================================================

#ifndef __IBEX_CTC_SC_H__
#define __IBEX_CTC_SC_H__

#include "ibex_Ctc.h"
#include "ibex_Linearizer.h"
#include "ibex_LPSolver.h"
#include "ibex_IntervalMatrix.h"
#include "ibex_SCSimplex.h"

#include <vector>

namespace ibex {

/**
 * \ingroup contractor
 *
 * \brief Self-certifying simplex (SC): a polytope-hull contractor.
 *
 * Same role as CtcPolytopeHull: for each variable x_k, the bounds of x_k over
 * the linear relaxation given by a Linearizer are computed and intersected
 * with the box. Instead of an external LP solver, the 2n bound LPs
 *
 *     min / max z_k  s.t.  Abar z = 0,  z in [z]      (Abar = [A | -I])
 *
 * are solved by an internal dual simplex (SCSimplex) that shares one basis
 * among them, and each bound is certified from the multipliers lambda of the
 * final basis by the Neumaier-Shcherbina argument: with gamma = Abar^T lambda
 * computed in interval arithmetic, gamma . z = 0 for every feasible z, so
 *     z_k in -(sum_{j!=k} gamma_j [z_j]) / gamma_k
 * and 0 notin gamma . [z] proves that the relaxation is empty. The bound is
 * valid for any lambda, so the simplex never needs to verify its solution.
 *
 * One call = one linearization and one pass over the bounds:
 *  - twin rows (identical, or proportional within a relative tolerance of
 *    1e-9) are merged; the polytope does not change;
 *  - the first LP is min x_0; then, among the bounds not yet attained by the
 *    last primal point (relative distance >= 1e-8), the nearest one to that
 *    point is solved next; each LP is warm-started from the previous basis;
 *  - an LP stopped by the iteration limit (100 pivots by default) still gives a
 *    certified bound; an infeasibility whose ray does not certify emptiness
 *    is solved once more from a cold start;
 *  - the box is contracted after each LP, and the pass stops as soon as the
 *    relaxation is proved empty.
 *
 * Reference: I. Araya, Certificate-driven polytope-hull contraction for
 * interval branch and bound (submitted).
 */
class CtcSC : public Ctc {
public:
	/**
	 * \brief Warm start between two consecutive LPs of a pass.
	 *
	 * DUAL:     relocate the nonbasic variables, then dual steps only
	 *           (best on systems of equations);
	 * PRIMAL:   keep the nonbasic variables, primal steps when needed
	 *           (best on optimization problems);
	 * ADAPTIVE: alternate both on the first calls and keep the one with the
	 *           lower CPU time per LP, trying the other again every 25 calls
	 *           (primal steps with devex pricing). Not deterministic.
	 */
	typedef enum { DUAL, PRIMAL, ADAPTIVE } WarmStart;

	/**
	 * \param lr        the linearizer (linear relaxation)
	 * \param ws        warm-start strategy
	 * \param tol       feasibility tolerance of the simplex (1e-9; 1e-6 is
	 *                  better for relaxations with nearly parallel rows)
	 * \param max_iter  pivots per LP; -1 means 100, the limit Ibex gives SoPlex
	 */
	CtcSC(Linearizer& lr, WarmStart ws=DUAL, double tol=1e-9, int max_iter=-1);

	virtual ~CtcSC();

	virtual void contract(IntervalVector& box) override;

	virtual void contract(IntervalVector& box, ContractContext& context) override;

	virtual void add_property(const IntervalVector& init_box, BoxProperties& prop) override;

	/** \brief Number of calls that ran a pass. */
	long nb_passes() const { return n_passes; }
	/** \brief Number of LPs solved. */
	long nb_lps() const { return n_lps; }
	/** \brief Number of simplex pivots. */
	long nb_pivots() const { return n_pivots; }

protected:
	/** Linearize and build Abar, [z]. Returns m, 0 or -1 (empty). */
	int linearize(const IntervalVector& box, ContractContext& context);

	/** The pass over the 2n bounds. */
	void pass(IntervalVector& box);

	/** Choose the strategy of the next pass (ADAPTIVE). */
	void choose_strategy();

	Linearizer& lr;
	const int nb_var;
	const WarmStart ws;
	const int max_iter;

	LPSolver mylinearsolver;     // receives the rows of the linearizer
	SCSimplex sx;

	IntervalMatrix A;            // Abar = [A | -I]
	Matrix Af;                   // Abar in floating point
	IntervalVector z;            // [x ; b]
	IntervalVector gamma;        // certificate Abar^T lambda
	std::vector<char> done;      // bounds already solved or skipped
	std::vector<double> xstar;   // last primal point

	// ADAPTIVE: moving average of the CPU time per LP of each warm start
	double ema[2];
	long n_adapt;
	int current;                 // strategy of the current pass (0 dual, 1 primal)

	long n_passes, n_lps, n_pivots;
};

} // namespace ibex

#endif // __IBEX_CTC_SC_H__
