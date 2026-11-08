//============================================================================
//                                  I B E X
// File        : ibex_SCSimplex.cpp
// Author      : Ignacio Araya
// License     : See the LICENSE file
// Created     : Nov 08, 2026
//============================================================================

#include "ibex_SCSimplex.h"

#include <cmath>
#include <algorithm>

namespace ibex {

namespace {

const double PIV_TOL     = 1e-9;   // smallest acceptable pivot
const double HARRIS_BAND = 1e-3;   // relative band of Harris' ratio test
const int    XB_REFRESH  = 32;     // pivots between two recomputations of x_B

} // anonymous namespace

SCSimplex::SCSimplex(double feas_tol, double dual_tol) :
		feas_tol(feas_tol), dual_tol(dual_tol),
		next_dual(true), next_primal_steps(false), next_primal_devex(false),
		dual_node(true), pure_node(true), primal_steps(false), primal_devex(false),
		m(0), nx(0), na(0), nc(0), iters(0), basis_ready(false), fact_ok(false),
		xb_dirty(true), piv_since_refresh(0), col_primal_ready(-1), pivot_failed(false),
		r_inf(-1), n_solves(0) {
}

void SCSimplex::set_strategy(bool dual, bool primal_steps_, bool primal_devex_) {
	next_dual = dual;
	next_primal_steps = primal_steps_;
	next_primal_devex = primal_devex_;
}

inline double SCSimplex::tol_primal(double, double, double esc) const {
	return feas_tol * esc;
}

void SCSimplex::load(const Matrix& Abar, int nx_) {
	wp.clear();   // primal devex weights: new reference frame per relaxation
	m  = Abar.nb_rows();
	nx = nx_;
	na = Abar.nb_cols();
	nc = na + m;
	dual_node = pure_node = next_dual;
	primal_steps = next_primal_steps;
	primal_devex = next_primal_devex;

	arow.assign(m, std::vector<std::pair<int,double> >());
	for (int j = 0; j < m; ++j)
		for (int c = 0; c < na; ++c) {
			const double v = Abar[j][c];
			if (v != 0.0) arow[j].push_back(std::make_pair(c, v));
		}
	d.assign(na, 0.0);
	ylam.assign(m, 0.0);
	y.assign(m, 0.0);
	tidx.assign(m + 1, std::vector<int>());
	tval.assign(m + 1, std::vector<double>());
	basic.assign(m + 1, -1);
	rowof.assign(nc, -1);
	atupper.assign(nc, 0);
	xN.assign(nc, 0.0);
	xB.assign(m + 1, 0.0);
	lam_buf.assign(m, 0.0);
	vbuf.assign(m, 0.0);
	xhat.assign(na, 0.0);
	cols_buf.assign(m, 0);
	fact.load_matrix(Abar);
	fact_ok = false;
	basis_ready = false;
}

void SCSimplex::repair_infinite_bounds(const IntervalVector& z) {
	for (int c = 0; c < na; ++c) {
		if (rowof[c] >= 0) continue;
		const double v = atupper[c] ? z[c].ub() : z[c].lb();
		if (v != POS_INFINITY && v != NEG_INFINITY) continue;
		const double o = atupper[c] ? z[c].lb() : z[c].ub();
		if (o != POS_INFINITY && o != NEG_INFINITY) atupper[c] = atupper[c] ? 0 : 1;
	}
}

/* Fallback when the LU is not usable: x_B from the stored tableau rows. */
void SCSimplex::compute_basics(const IntervalVector& z) {
	if (fact_ok) { compute_basics_revised(z); return; }
	xb_dirty = false;
	for (int c = 0; c < nc; ++c) {
		if (rowof[c] >= 0) { xN[c] = 0.0; continue; }
		if (c >= na) { xN[c] = 0.0; continue; }
		double v = atupper[c] ? z[c].ub() : z[c].lb();
		if (v == POS_INFINITY || v == NEG_INFINITY) {
			const double o = atupper[c] ? z[c].lb() : z[c].ub();
			if (o == POS_INFINITY || o == NEG_INFINITY) v = 0.0;
			else { v = o; atupper[c] = atupper[c] ? 0 : 1; }
		}
		xN[c] = v;
	}
	for (int r = 1; r <= m; ++r) {
		double acc = 0.0;
		const std::vector<int>&    ix = tidx[r];
		const std::vector<double>& vl = tval[r];
		for (size_t t = 0; t < ix.size(); ++t) {
			const int c = ix[t];
			if (c >= na) break;
			if (rowof[c] >= 0) continue;
			acc += vl[t] * xN[c];
		}
		xB[r] = -acc;
	}
}

/* x_B = B^-1 (-Abar_N x_N), nonbasic variables at their bounds. */
void SCSimplex::compute_basics_revised(const IntervalVector& z) {
	xb_dirty = false;
	piv_since_refresh = 0;
	for (int c = 0; c < na; ++c) {
		if (rowof[c] >= 0) { xN[c] = 0.0; xhat[c] = 0.0; continue; }
		double v = atupper[c] ? z[c].ub() : z[c].lb();
		if (v == POS_INFINITY || v == NEG_INFINITY) {
			const double o = atupper[c] ? z[c].lb() : z[c].ub();
			if (o == POS_INFINITY || o == NEG_INFINITY) v = 0.0;
			else { v = o; atupper[c] = atupper[c] ? 0 : 1; }
		}
		xN[c] = v; xhat[c] = v;
	}
	fact.A_times(&xhat[0], &vbuf[0]);
	for (int j = 0; j < m; ++j) vbuf[j] = -vbuf[j];
	fact.ftran_vec(&vbuf[0], &lam_buf[0]);
	for (int r = 1; r <= m; ++r) xB[r] = lam_buf[r-1];
}

/* Objective row from the stored tableau rows (used at a cold start, where the
 * basis is the slack basis, and as a fallback when the LU is not usable). */
void SCSimplex::set_objective(int k, bool maximize, bool relocate) {
	xb_dirty = true;
	std::fill(d.begin(), d.end(), 0.0);
	std::fill(ylam.begin(), ylam.end(), 0.0);
	d[k] = maximize ? -1.0 : 1.0;
	for (int r = 1; r <= m; ++r) {
		const int b = basic[r];
		const double f = (b < na) ? d[b] : ylam[b - na];
		if (f == 0.0) continue;
		const std::vector<int>&    ix = tidx[r];
		const std::vector<double>& vl = tval[r];
		for (size_t t = 0; t < ix.size(); ++t) {
			const int c = ix[t];
			if (c < na) d[c]       -= f * vl[t];
			else        ylam[c-na] -= f * vl[t];
		}
		if (b < na) d[b] = 0.0; else ylam[b - na] = 0.0;
	}
	if (relocate)
		for (int c = 0; c < na; ++c) {
			if (rowof[c] >= 0) continue;
			atupper[c] = (d[c] < -dual_tol) ? 1 : 0;
		}
}

/* Row r of B^-1 Abar by BTRAN, stored in tidx[r]/tval[r] (entries below
 * 1e-14 times the largest are dropped). */
void SCSimplex::pivot_row(int r) {
	if ((int) rho_buf.size() != m) rho_buf.assign(m, 0.0);
	if ((int) rowA_buf.size() != na) rowA_buf.assign(na, 0.0);
	fact.btran_row(r - 1, &rho_buf[0]);
	fact.row_of_A(&rho_buf[0], &rowA_buf[0]);
	double amax = 0.0;
	for (int c = 0; c < na; ++c) if (std::fabs(rowA_buf[c]) > amax) amax = std::fabs(rowA_buf[c]);
	const double floor = amax * 1e-14;
	tidx[r].clear(); tval[r].clear();
	for (int c = 0; c < na; ++c)
		if (std::fabs(rowA_buf[c]) > floor) { tidx[r].push_back(c); tval[r].push_back(rowA_buf[c]); }
}

/* Column col enters at row r (pivot_row(r) must have been called). Updates
 * d and y with the pivot row and the LU with the entering column. */
bool SCSimplex::pivot(int r, int col) {
	xb_dirty = true;
	const bool reuse = (col_primal_ready == col) && (fact.prepared_column() == col);
	col_primal_ready = -1;
	if (col >= na) return false;
	const double p = rowA_buf[col];
	if (std::fabs(p) < PIV_TOL) return false;
	const double theta = d[col] / p;
	if (theta != 0.0) {
		const std::vector<int>& ix = tidx[r]; const std::vector<double>& vl = tval[r];
		for (size_t t = 0; t < ix.size(); ++t) d[ix[t]] -= theta * vl[t];
		for (int j = 0; j < m; ++j) ylam[j] -= theta * rho_buf[j];   // y += theta*rho
	}
	if ((int) colq_buf.size() != m) colq_buf.assign(m, 0.0);
	const bool col_ok = reuse ? true : fact.prepare_change(col, &colq_buf[0]);
	if ((int) alpha_piv.size() != m+1) alpha_piv.assign(m+1, 0.0);
	alpha_piv[0] = 0.0;
	for (int i = 1; i <= m; ++i) alpha_piv[i] = col_ok ? colq_buf[i-1] : 0.0;
	if (!col_ok) alpha_piv[r] = 0.0;
	const int out = basic[r];
	d[col] = 0.0;
	rowof[out] = -1;
	basic[r] = col;
	rowof[col] = r;
	for (int rr = 1; rr <= m; ++rr) { const int b = basic[rr]; if (b != out) d[b] = 0.0; }
	if (!fact.apply_change(r - 1, col)) {
		for (int rr = 1; rr <= m; ++rr) cols_buf[rr-1] = basic[rr];
		fact_ok = fact.set_basis(cols_buf);
	}
	++piv_since_refresh;
	return fact_ok;
}

bool SCSimplex::dual_step(const IntervalVector& z, bool& infeasible) {
	if (piv_since_refresh >= XB_REFRESH) xb_dirty = true;
	if (xb_dirty) compute_basics(z);

	infeasible = false;
	// leaving variable: largest violation relative to 1+|bound|
	int r_out = -1; double worst = -1.0; bool below = false;
	for (int r = 1; r <= m; ++r) {
		const int c = basic[r];
		const double lo = z[c].lb(), hi = z[c].ub();
		const double esc = 1.0 + ((std::fabs(lo) > std::fabs(hi)) ? std::fabs(lo) : std::fabs(hi));
		const double tol = tol_primal(lo, hi, esc);
		if (xB[r] < lo - tol) {
			const double v = (lo - xB[r]) / esc;
			if (v > worst) { worst = v; r_out = r; below = true; }
		} else if (xB[r] > hi + tol) {
			const double v = (xB[r] - hi) / esc;
			if (v > worst) { worst = v; r_out = r; below = false; }
		}
	}
	if (r_out < 0) return false;   // primal feasible

	// entering variable: Harris' ratio test with a relative band
	int j_in = -1; double best = 0.0, best_piv = 0.0;
	pivot_row(r_out);
	{
		const std::vector<int>&    ix = tidx[r_out];
		const std::vector<double>& vl = tval[r_out];
		double amax = 0.0;
		for (size_t t = 0; t < ix.size(); ++t) {
			if (ix[t] >= na) break;
			const double av = std::fabs(vl[t]);
			if (av > amax) amax = av;
		}
		const double piv_min = (amax > 0.0) ? amax * 1e-11 : PIV_TOL;
		for (size_t t = 0; t < ix.size(); ++t) {
			const int c = ix[t];
			if (c >= na) break;
			if (rowof[c] >= 0) continue;
			const double a = vl[t];
			const double av = std::fabs(a);
			if (av < piv_min || av < PIV_TOL) continue;
			const bool up = atupper[c] ? (a > 0.0) : (a < 0.0);
			if (below != up) continue;
			const double ratio = std::fabs(d[c] / a);
			if (j_in < 0 || ratio < best * (1.0 - HARRIS_BAND)) {
				best = ratio; best_piv = av; j_in = c;
			} else if (ratio <= best * (1.0 + HARRIS_BAND) && av > best_piv) {
				best_piv = av; j_in = c;
				if (ratio < best) best = ratio;
			}
		}
	}
	if (j_in < 0) { r_inf = r_out; infeasible = true; return false; }

	const int c_out = basic[r_out];
	const double v_out = below ? z[c_out].lb() : z[c_out].ub();
	const double xb_r  = xB[r_out];
	const double x_in  = xN[j_in];
	if (!pivot(r_out, j_in)) {
		pivot_failed = true;
		infeasible = false; return false;
	}
	atupper[c_out] = below ? 0 : 1;
	// incremental update of x_B with the entering column
	if ((int) alpha_piv.size() == m+1) {
		const double ar = alpha_piv[r_out];
		if (ar != 0.0 && v_out == v_out && x_in == x_in &&
			v_out != POS_INFINITY && v_out != NEG_INFINITY &&
			x_in  != POS_INFINITY && x_in  != NEG_INFINITY) {
			const double delta = (xb_r - v_out) / ar;
			for (int i = 1; i <= m; ++i)
				if (i != r_out && alpha_piv[i] != 0.0)
					xB[i] -= alpha_piv[i] * delta;
			xB[r_out] = x_in + delta;
			xN[c_out] = v_out;
			xN[j_in]  = 0.0;
			repair_infinite_bounds(z);
			xb_dirty = false;
		} else {
			xb_dirty = true;
		}
	}
	return true;
}

bool SCSimplex::primal_step(const IntervalVector& z, bool& unbounded) {
	unbounded = false;
	// pricing: Dantzig, or primal devex (d_j^2 / w_j)
	if (primal_devex && (int) wp.size() != na) wp.assign(na, 1.0);
	int q = -1; double bestd = 0.0;
	for (int c = 0; c < na; ++c) {
		if (rowof[c] >= 0) continue;
		const double dc = d[c];
		double v;
		if (!atupper[c] && dc < -dual_tol) v = -dc;
		else if (atupper[c] && dc > dual_tol) v = dc;
		else continue;
		if (primal_devex) v = v * v / wp[c];
		if (v > bestd) { bestd = v; q = c; }
	}
	if (q < 0) return false;   // dual feasible
	if (!ensure_fact()) { pivot_failed = true; return false; }
	if (piv_since_refresh >= XB_REFRESH) xb_dirty = true;
	if (xb_dirty) compute_basics(z);

	if ((int) colq_buf.size() != m) colq_buf.assign(m, 0.0);
	if (!fact.prepare_change(q, &colq_buf[0])) { pivot_failed = true; return false; }
	const double sgn = atupper[q] ? -1.0 : 1.0;
	double amax = 0.0;
	for (int r = 1; r <= m; ++r) { const double av = std::fabs(colq_buf[r-1]); if (av > amax) amax = av; }
	const double piv_min = (amax > 0.0) ? amax*1e-11 : PIV_TOL;
	// ratio test (pivots below the threshold are filtered out first)
	double lim = POS_INFINITY; int r_out = -1; bool leaves_up = false; double piv_out = 0.0;
	for (int r = 1; r <= m; ++r) {
		const double a = colq_buf[r-1];
		if (a == 0.0) continue;
		const double av = std::fabs(a);
		if (av < piv_min || av < PIV_TOL) continue;
		const int c = basic[r];
		const double lo = z[c].lb(), hi = z[c].ub();
		const double der = -a * sgn;
		double bound; bool up;
		if (der > 0.0) { bound = (hi - xB[r]) / der; up = true; }
		else           { bound = (lo - xB[r]) / der; up = false; }
		if (bound < 0.0) bound = 0.0;
		if (r_out < 0 || bound < lim*(1.0 - 1e-9)) { lim = bound; r_out = r; leaves_up = up; piv_out = av; }
		else if (bound <= lim*(1.0 + 1e-9) && av > piv_out) { r_out = r; leaves_up = up; piv_out = av; if (bound < lim) lim = bound; }
	}
	const double x_q = xN[q];
	const bool fin_q = (x_q == x_q) && x_q != POS_INFINITY && x_q != NEG_INFINITY;
	const double range = z[q].diam();
	if (range < lim) {   // bound flip of the entering variable
		atupper[q] = atupper[q] ? 0 : 1;
		if (fin_q && range == range && range < POS_INFINITY) {
			for (int r = 1; r <= m; ++r) { const double a = colq_buf[r-1]; if (a != 0.0) xB[r] -= a * sgn * range; }
			xN[q] = x_q + sgn * range;
			xb_dirty = false;
		} else xb_dirty = true;
		return true;
	}
	if (r_out < 0) { unbounded = true; return false; }
	pivot_row(r_out);
	const int c_out = basic[r_out];
	const double v_out = leaves_up ? z[c_out].ub() : z[c_out].lb();
	if (primal_devex && (int) wp.size() == na) {
		const double arq = rowA_buf[q];
		if (arq != 0.0) {
			const double wq = wp[q];
			const std::vector<int>& ix = tidx[r_out]; const std::vector<double>& vl = tval[r_out];
			double wmax = 0.0;
			for (size_t t = 0; t < ix.size(); ++t) {
				const int c = ix[t];
				if (c == q || c >= na || rowof[c] >= 0) continue;
				const double f = vl[t] / arq;
				const double cand = f * f * wq;
				if (cand > wp[c]) wp[c] = cand;
				if (wp[c] > wmax) wmax = wp[c];
			}
			if (c_out < na) { const double nw = wq / (arq * arq); wp[c_out] = nw > 1.0 ? nw : 1.0; }
			if (!(wmax <= 1e8)) std::fill(wp.begin(), wp.end(), 1.0);   // reset the reference frame
		}
	}
	col_primal_ready = q;
	std::vector<double>& alpha = colq_buf;
	if (!pivot(r_out, q)) { pivot_failed = true; return false; }
	atupper[c_out] = leaves_up ? 1 : 0;
	if (fin_q && lim == lim && lim < POS_INFINITY &&
		v_out == v_out && v_out != POS_INFINITY && v_out != NEG_INFINITY) {
		for (int r = 1; r <= m; ++r) {
			if (r == r_out) continue;
			const double a = alpha[r-1];
			if (a != 0.0) xB[r] -= a * sgn * lim;
		}
		xB[r_out] = x_q + sgn * lim;
		xN[c_out] = v_out;
		xN[q] = 0.0;
		repair_infinite_bounds(z);
		xb_dirty = false;
	} else xb_dirty = true;
	return true;
}

bool SCSimplex::ensure_fact() {
	return fact_ok;
}

/* y = c_B^T B^-1 = +-(row of B^-1 where z_k is basic), or 0 if z_k is nonbasic. */
bool SCSimplex::y_from_factorization(int k, bool maximize) {
	if (m <= 0) return false;
	const int r = rowof[k];
	if (r < 0) {
		std::fill(y.begin(), y.end(), 0.0);
		return true;
	}
	if (!ensure_fact()) return false;
	fact.btran_row(r - 1, &lam_buf[0]);
	const double ck = maximize ? -1.0 : 1.0;
	for (int j = 0; j < m; ++j) y[j] = ck * lam_buf[j];
	return true;
}

/* New objective on the current basis: y, d and x_B recomputed from the LU.
 * With relocate (dual warm start) nonbasic variables move to the bound that
 * makes the basis dual feasible. */
bool SCSimplex::refresh(int k, bool maximize, const IntervalVector& z, bool relocate) {
	xb_dirty = true;
	if (m <= 0 || !ensure_fact()) return false;
	if (!y_from_factorization(k, maximize)) return false;
	const double ck = maximize ? -1.0 : 1.0;
	std::fill(d.begin(), d.end(), 0.0);
	d[k] = ck;
	for (int j = 0; j < m; ++j) {
		const double yj = y[j];
		if (yj == 0.0) continue;
		const std::vector<std::pair<int,double> >& a = arow[j];
		for (size_t t = 0; t < a.size(); ++t) d[a[t].first] -= yj * a[t].second;
	}
	for (int j = 0; j < m; ++j) ylam[j] = -y[j];
	for (int r = 1; r <= m; ++r) { const int b = basic[r]; if (b < na) d[b] = 0.0; }
	if (relocate) relocate_nonbasic(k);
	for (int c = 0; c < na; ++c) {
		if (rowof[c] >= 0) { xhat[c] = 0.0; continue; }
		double v = atupper[c] ? z[c].ub() : z[c].lb();
		if (v == POS_INFINITY || v == NEG_INFINITY) {
			const double o = atupper[c] ? z[c].lb() : z[c].ub();
			if (o == POS_INFINITY || o == NEG_INFINITY) v = 0.0;
			else { v = o; atupper[c] = atupper[c] ? 0 : 1; }
		}
		xhat[c] = v;
	}
	fact.A_times(&xhat[0], &vbuf[0]);
	for (int j = 0; j < m; ++j) vbuf[j] = -vbuf[j];
	fact.ftran_vec(&vbuf[0], &lam_buf[0]);
	for (int r = 1; r <= m; ++r) xB[r] = lam_buf[r-1];
	return true;
}

int SCSimplex::relocate_nonbasic(int k) {
	const bool keep_k = !pure_node;
	int n_mov = 0;
	for (int c = 0; c < na; ++c) {
		if (rowof[c] >= 0) continue;
		if (keep_k && c == k) continue;
		const char before = atupper[c];
		if (d[c] < -dual_tol)      atupper[c] = 1;
		else if (d[c] > dual_tol)  atupper[c] = 0;
		if (atupper[c] != before) ++n_mov;
	}
	return n_mov;
}

SCSimplex::Status SCSimplex::solve(int k, bool maximize, const IntervalVector& z, int max_iter, bool warm) {
	++n_solves;
	r_inf = -1;
	xb_dirty = true;
	pivot_failed = false;
	iters = 0;
	if (m <= 0 || k < 0 || k >= na) return SINGULAR;

	const bool basis_ready_before = basis_ready;
	if (!warm || !basis_ready) {
		// slack basis B = -I: dual feasible for any objective +-e_k
		std::fill(rowof.begin(), rowof.end(), -1);
		std::fill(atupper.begin(), atupper.end(), 0);
		for (int r = 1; r <= m; ++r) {
			tidx[r].clear(); tval[r].clear();
			const std::vector<std::pair<int,double> >& a = arow[r-1];
			const int lamcol = na + (r-1);
			for (size_t t = 0; t < a.size(); ++t) {
				tidx[r].push_back(a[t].first);
				tval[r].push_back(-a[t].second);
			}
			tidx[r].push_back(lamcol); tval[r].push_back(-1.0);
			basic[r] = nx + (r-1);
			rowof[nx + (r-1)] = r;
		}
		for (int r = 1; r <= m; ++r) cols_buf[r-1] = basic[r];
		fact_ok = fact.set_basis(cols_buf);
		basis_ready = true;
	}
	const bool warm_start = warm && basis_ready_before;
	const bool relocate = dual_node;
	if (!warm_start || !refresh(k, maximize, z, relocate)) {
		set_objective(k, maximize, !warm_start);
		if (warm_start && relocate) relocate_nonbasic(k);
	}
	repair_infinite_bounds(z);

	Status st = ITER_LIMIT;
	while (iters < max_iter) {
		bool infeasible = false, unbounded = false;
		if (dual_step(z, infeasible)) { ++iters; continue; }
		if (infeasible) { st = INFEASIBLE; break; }
		if (primal_steps && primal_step(z, unbounded)) { ++iters; continue; }
		if (unbounded) { st = UNBOUNDED; break; }
		if (pivot_failed) { st = SINGULAR; break; }
		st = OPTIMAL; break;
	}
	if (st != OPTIMAL && st != ITER_LIMIT) basis_ready = false;

	// multipliers: the Farkas ray of the infeasible row, or c_B^T B^-1
	if (st == INFEASIBLE && r_inf >= 1 && r_inf <= m && ensure_fact()) {
		fact.btran_row(r_inf - 1, &lam_buf[0]);
		for (int j = 0; j < m; ++j) y[j] = lam_buf[j];
	} else if (!y_from_factorization(k, maximize)) {
		for (int r = 0; r < m; ++r) y[r] = -ylam[r];
	}
	repair_infinite_bounds(z);
	return st;
}

void SCSimplex::primal_solution(const IntervalVector& z, std::vector<double>& out) {
	out.assign(na, 0.0);
	if (m <= 0) return;
	compute_basics(z);
	for (int c = 0; c < na; ++c) {
		if (rowof[c] >= 0) out[c] = xB[rowof[c]];
		else               out[c] = atupper[c] ? z[c].ub() : z[c].lb();
	}
}

} // namespace ibex
