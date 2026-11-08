//============================================================================
//                                  I B E X
// File        : ibex_SCBasis.cpp
// Author      : Ignacio Araya
// License     : See the LICENSE file
// Created     : Nov 08, 2026
//============================================================================

#include "ibex_SCBasis.h"

#include "soplex.h"
#include "soplex/slufactor.h"

using namespace soplex;

namespace ibex {

constexpr int SCBasis::refactor_period;

namespace {

/* SLUFactor with its internal timer switched off (the timer costs a system
 * call per solve and is never read). */
struct SLUFactorNoTimer : public SLUFactor {
	SLUFactorNoTimer() {
		solveTime = TimerFactory::switchTimer(solveTime, Timer::OFF);
		timerType = Timer::OFF;
	}
};

} // anonymous namespace

struct SCBasis::Impl {
	std::vector<DSVector> col_nz;   // columns of Abar
	std::vector<DSVector> unit_nz;  // columns of the identity
	SLUFactorNoTimer lu;
	DVector rhs, sol;
	SSVector* eta;

	Impl() : rhs(1), sol(1), eta(NULL) { }
	~Impl() { delete eta; }

	const SVector& column(int c, int na) const {
		return (c < na) ? (const SVector&) col_nz[c] : (const SVector&) unit_nz[c - na];
	}
};

SCBasis::SCBasis() : m(0), na(0), updates_since_factor(0), col_prepared(-1), impl(new Impl()) {
}

SCBasis::~SCBasis() {
	delete impl;
}

void SCBasis::load_matrix(const Matrix& Abar) {
	m  = Abar.nb_rows();
	na = Abar.nb_cols();

	row_nz.assign(m, std::vector<std::pair<int,double> >());
	impl->col_nz.clear();  impl->col_nz.reserve(na);
	impl->unit_nz.clear(); impl->unit_nz.reserve(m);

	for (int c = 0; c < na; ++c) impl->col_nz.push_back(DSVector(m));
	for (int j = 0; j < m; ++j) {
		impl->unit_nz.push_back(DSVector(1));
		impl->unit_nz[j].add(j, 1.0);
	}

	for (int j = 0; j < m; ++j) {
		for (int c = 0; c < na; ++c) {
			const double v = Abar[j][c];
			if (v == 0.0) continue;
			row_nz[j].push_back(std::make_pair(c, v));
			impl->col_nz[c].add(j, v);
		}
	}

	impl->rhs.reDim(m); impl->sol.reDim(m);
	delete impl->eta; impl->eta = new SSVector(m);
}

bool SCBasis::set_basis(const std::vector<int>& cols) {
	if ((int) cols.size() != m) return false;
	basic.assign(m, 0);
	where_basic.assign(na + m, -1);
	for (int j = 0; j < m; ++j) {
		const int c = cols[j];
		if (c < 0 || c >= na + m) return false;
		basic[j] = c;
		where_basic[c] = j;
	}
	updates_since_factor = 0;
	return factorize();
}

/* SoPlex exceptions do not derive from std::exception: catch everything and
 * report a factorization that cannot be used. */
bool SCBasis::factorize() {
	std::vector<const SVector*> cols(m);
	for (int j = 0; j < m; ++j) cols[j] = &impl->column(basic[j], na);
	updates_since_factor = 0;
	col_prepared = -1;   // a column prepared with the former LU is no longer valid
	try {
		return impl->lu.load(&cols[0], m) == SLinSolver::OK;
	} catch (...) {
		return false;
	}
}

bool SCBasis::change_basis(int j, int col) {
	if (j < 0 || j >= m || col < 0 || col >= na + m) return false;

	const int out = basic[j];
	if (out == col) return true;
	where_basic[out] = -1;
	basic[j] = col;
	where_basic[col] = j;

	if (updates_since_factor >= refactor_period)
		return factorize();

	const SVector& newcol = impl->column(col, na);
	try {
		impl->lu.solveRight4update(*impl->eta, newcol);
		SLinSolver::Status st = impl->lu.change(j, newcol, impl->eta);
		++updates_since_factor;
		if (st != SLinSolver::OK) return factorize();
	} catch (...) {
		return factorize();
	}
	return true;
}

bool SCBasis::prepare_change(int col, double* out) {
	col_prepared = -1;
	if (col < 0 || col >= na + m) return false;
	try {
		const SVector& newcol = impl->column(col, na);
		impl->lu.solveRight4update(*impl->eta, newcol);
		for (int r = 0; r < m; ++r) out[r] = (*impl->eta)[r];
		col_prepared = col;
		return true;
	} catch (...) {
		for (int r = 0; r < m; ++r) out[r] = 0.0;
		return false;
	}
}

bool SCBasis::apply_change(int j, int col) {
	if (col_prepared != col) { col_prepared = -1; return change_basis(j, col); }
	col_prepared = -1;
	if (j < 0 || j >= m || col < 0 || col >= na + m) return false;
	const int out = basic[j];
	if (out == col) return true;
	where_basic[out] = -1;
	basic[j] = col;
	where_basic[col] = j;
	if (updates_since_factor >= refactor_period)
		return factorize();
	try {
		SLinSolver::Status st = impl->lu.change(j, impl->column(col, na), impl->eta);
		++updates_since_factor;
		if (st != SLinSolver::OK) return factorize();
	} catch (...) {
		return factorize();
	}
	return true;
}

void SCBasis::btran_row(int j, double* out) {
	try {
		for (int r = 0; r < m; ++r) impl->rhs[r] = 0.0;
		impl->rhs[j] = 1.0;
		impl->lu.solveLeft(impl->sol, impl->rhs);
		for (int r = 0; r < m; ++r) out[r] = impl->sol[r];
	} catch (...) {
		for (int r = 0; r < m; ++r) out[r] = 0.0;
	}
}

void SCBasis::ftran_vec(const double* v, double* out) {
	try {
		for (int r = 0; r < m; ++r) impl->rhs[r] = v[r];
		impl->lu.solveRight(impl->sol, impl->rhs);
		for (int r = 0; r < m; ++r) out[r] = impl->sol[r];
	} catch (...) {
		for (int r = 0; r < m; ++r) out[r] = 0.0;
	}
}

void SCBasis::row_of_A(const double* lambda, double* out) const {
	for (int c = 0; c < na; ++c) out[c] = 0.0;
	for (int j = 0; j < m; ++j) {
		const double l = lambda[j];
		if (l == 0.0) continue;
		const std::vector<std::pair<int,double> >& r = row_nz[j];
		for (size_t t = 0; t < r.size(); ++t) out[r[t].first] += l * r[t].second;
	}
}

void SCBasis::A_times(const double* u, double* out) const {
	for (int j = 0; j < m; ++j) {
		double acc = 0.0;
		const std::vector<std::pair<int,double> >& r = row_nz[j];
		for (size_t t = 0; t < r.size(); ++t) acc += r[t].second * u[r[t].first];
		out[j] = acc;
	}
}

} // namespace ibex
