//============================================================================
//                                  I B E X
// File        : ibex_SCSimplex.h
// Author      : Ignacio Araya
// License     : See the LICENSE file
// Created     : Nov 08, 2026
//============================================================================

#ifndef __IBEX_SC_SIMPLEX_H__
#define __IBEX_SC_SIMPLEX_H__

#include "ibex_Matrix.h"
#include "ibex_IntervalVector.h"
#include "ibex_SCBasis.h"

#include <vector>
#include <utility>

namespace ibex {

/**
 * \ingroup numeric
 *
 * \brief Bounded-variable simplex for the 2n bound LPs of CtcSC.
 *
 * Solves  min (or max) z_k  s.t.  Abar z = 0,  z in [z]
 * where Abar = [A | -I] is m x na (na = nx + m: nx structural columns and the
 * m slack columns). The 2n LPs of a node share one basis: each solve can be
 * warm-started from the final basis of the previous one.
 *
 * Revised method with a sparse LU of the basis (SCBasis). The first LP starts
 * from the slack basis, which is dual feasible for any objective +-e_k, and
 * is solved by the dual method (no phase 1).
 *
 * Warm starts:
 *  - dual: every nonbasic variable is relocated to the bound given by the
 *    sign of its new reduced cost; then dual steps only;
 *  - primal: nonbasic variables stay where they are; dual steps repair primal
 *    infeasibility, primal steps (Dantzig, or primal devex) repair dual
 *    infeasibility.
 *
 * Leaving variable: largest relative violation. Entering variable: ratio test
 * with Harris' relative band (1e-3), the largest pivot within the band.
 * x_B is recomputed every 32 pivots.
 *
 * The multipliers returned by #lambda() are always valid input for the
 * Neumaier-Shcherbina certificate computed by the caller, whatever the
 * status: on OPTIMAL and ITER_LIMIT they are those of the final basis, on
 * INFEASIBLE they are the row of B^-1 of the infeasible basic variable (a
 * Farkas ray).
 */
class SCSimplex {
public:
	enum Status { OPTIMAL, INFEASIBLE, ITER_LIMIT, UNBOUNDED, SINGULAR };

	/**
	 * \param feas_tol  primal feasibility tolerance (relative to 1+|bound|)
	 * \param dual_tol  dual feasibility tolerance
	 */
	SCSimplex(double feas_tol=1e-9, double dual_tol=1e-9);

	/**
	 * \brief Set the warm-start strategy used by the next #load.
	 *
	 * \param dual          dual warm start (relocation)
	 * \param primal_steps  allow primal steps (required by the primal warm start)
	 * \param primal_devex  devex pricing in primal steps (Dantzig otherwise)
	 */
	void set_strategy(bool dual, bool primal_steps, bool primal_devex);

	/** \brief Load Abar (m x na) whose first nx columns are structural. */
	void load(const Matrix& Abar, int nx);

	/**
	 * \brief Solve LP_k (min z_k, or max z_k if \a maximize).
	 *
	 * \param warm  start from the current basis if there is one (cold
	 *              start from the slack basis otherwise)
	 */
	Status solve(int k, bool maximize, const IntervalVector& z, int max_iter, bool warm=true);

	/** \brief Multipliers of the last solve (m entries). */
	const std::vector<double>& lambda() const { return y; }

	/** \brief Primal point of the current basis (na entries). */
	void primal_solution(const IntervalVector& z, std::vector<double>& out);

	/** \brief Pivots of the last solve. */
	int iterations() const { return iters; }

	/** \brief Number of calls to solve since construction. */
	long nb_solves() const { return n_solves; }

private:
	SCSimplex(const SCSimplex&);
	SCSimplex& operator=(const SCSimplex&);

	double tol_primal(double lo, double hi, double esc) const;

	bool ensure_fact();
	bool y_from_factorization(int k, bool maximize);
	bool refresh(int k, bool maximize, const IntervalVector& z, bool relocate);
	int  relocate_nonbasic(int k);
	void set_objective(int k, bool maximize, bool relocate);
	void repair_infinite_bounds(const IntervalVector& z);
	void compute_basics(const IntervalVector& z);
	void compute_basics_revised(const IntervalVector& z);
	void pivot_row(int r);
	bool pivot(int r, int col);
	bool dual_step(const IntervalVector& z, bool& infeasible);
	bool primal_step(const IntervalVector& z, bool& unbounded);

	const double feas_tol, dual_tol;

	// strategy of the current node (fixed by load)
	bool next_dual, next_primal_steps, next_primal_devex;
	bool dual_node, pure_node, primal_steps, primal_devex;

	int m, nx, na, nc;
	std::vector<std::vector<std::pair<int,double> > > arow;   // rows of Abar

	std::vector<double> d;      // reduced costs (na)
	std::vector<double> ylam;   // -y (m)
	/* Rows of the tableau. Only the pivot row is kept up to date (by BTRAN);
	 * the others are used by the fallback paths when the LU fails. */
	std::vector<std::vector<int> >    tidx;
	std::vector<std::vector<double> > tval;
	std::vector<int>  basic;     // basic column of each row 1..m
	std::vector<int>  rowof;     // row of each basic column, or -1
	std::vector<char> atupper;   // nonbasic at its upper bound

	std::vector<double> xN, xB, vbuf, xhat, y;
	std::vector<double> rho_buf, rowA_buf, colq_buf, lam_buf, alpha_piv, wp;
	std::vector<int>    cols_buf;

	int  iters;
	bool basis_ready;
	SCBasis fact;
	bool fact_ok;
	bool xb_dirty;
	int  piv_since_refresh;
	int  col_primal_ready;
	bool pivot_failed;
	int  r_inf;
	long n_solves;
};

} // namespace ibex

#endif // __IBEX_SC_SIMPLEX_H__
