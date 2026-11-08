//============================================================================
//                                  I B E X
// File        : ibex_SCBasis.h
// Author      : Ignacio Araya
// License     : See the LICENSE file
// Created     : Nov 08, 2026
//============================================================================

#ifndef __IBEX_SC_BASIS_H__
#define __IBEX_SC_BASIS_H__

#include "ibex_Matrix.h"

#include <vector>
#include <utility>

namespace ibex {

/**
 * \ingroup numeric
 *
 * \brief Sparse LU factorization of a simplex basis (used by CtcSC).
 *
 * The constraint matrix is [Abar | I] (m rows, na+m columns): columns
 * 0..na-1 come from Abar, columns na..na+m-1 are the unit columns.
 * The factorization is the one of SoPlex (SLUFactor), updated after each
 * basis change by Forrest-Tomlin and recomputed from scratch every
 * #refactor_period updates. SoPlex is hidden behind a private implementation,
 * so that this header does not depend on SoPlex.
 *
 * Every SoPlex exception is caught: a failed operation only returns false
 * (or zeros), which the simplex treats as a degraded factorization.
 */
class SCBasis {
public:
	/** \brief Number of Forrest-Tomlin updates between two refactorizations. */
	static constexpr int refactor_period = 50;

	SCBasis();
	~SCBasis();

	/** \brief Load the matrix Abar (m x na). */
	void load_matrix(const Matrix& Abar);

	/** \brief Set and factorize the basis (cols[j]: basic column of row j). */
	bool set_basis(const std::vector<int>& cols);

	/** \brief Column \a col replaces the basic column of row \a j. */
	bool change_basis(int j, int col);

	/**
	 * \brief FTRAN of column \a col (out = B^-1 a_col), kept for the update
	 * that follows (#apply_change).
	 */
	bool prepare_change(int col, double* out);

	/** \brief Basis change that reuses the FTRAN of #prepare_change if possible. */
	bool apply_change(int j, int col);

	/** \brief Column prepared by #prepare_change, or -1. */
	int prepared_column() const { return col_prepared; }

	/** \brief out = e_j^T B^-1 (BTRAN). */
	void btran_row(int j, double* out);

	/** \brief out = B^-1 v (FTRAN). */
	void ftran_vec(const double* v, double* out);

	/** \brief out = lambda^T Abar (na entries). */
	void row_of_A(const double* lambda, double* out) const;

	/** \brief out = Abar u (m entries). */
	void A_times(const double* u, double* out) const;

private:
	SCBasis(const SCBasis&);            // not copyable
	SCBasis& operator=(const SCBasis&);

	bool factorize();

	int m, na;
	std::vector<std::vector<std::pair<int,double> > > row_nz;
	std::vector<int> basic;        // basic column of each row
	std::vector<int> where_basic;  // row of each basic column, or -1
	int updates_since_factor;
	int col_prepared;

	struct Impl;                   // SoPlex objects
	Impl* impl;
};

} // namespace ibex

#endif // __IBEX_SC_BASIS_H__
