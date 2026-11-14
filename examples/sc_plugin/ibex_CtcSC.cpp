//============================================================================
//                                  I B E X
// File        : ibex_CtcSC.cpp
// Author      : Ignacio Araya
// License     : See the LICENSE file
// Created     : Nov 08, 2026
//============================================================================

#include "ibex_CtcSC.h"

#include <cmath>
#include <ctime>
#include <map>

namespace ibex {

namespace {

const double PREC        = 1e-8;   // a bound closer than that to the last primal point is skipped
const double TWIN_TOL    = 1e-9;   // relative tolerance of proportional rows
const double COEF_MAX    = 1e100;  // rows with a larger coefficient are dropped
const long   ADAPT_EXPLORE = 5;    // ADAPTIVE: calls per strategy before choosing
const long   ADAPT_RETRY   = 25;   // ADAPTIVE: the other strategy is tried every 25 calls

double cpu_time() {
#if defined(CLOCK_PROCESS_CPUTIME_ID)
	struct timespec ts;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
	return ts.tv_sec + 1e-9*ts.tv_nsec;
#else
	return ((double) std::clock()) / CLOCKS_PER_SEC;
#endif
}

/* Bound of z_k from gamma, in interval arithmetic: gamma . z = 0 for every
 * feasible z, so z_k in -(sum_{c!=k} gamma_c [z_c]) / gamma_k. */
Interval bound_from_gamma(const IntervalVector& gamma, const IntervalVector& z, int k) {
	if (gamma[k].contains(0.0)) return Interval::all_reals();
	const int na = gamma.size();
	const Interval gk = gamma[k];
	Interval acc(0.0);
	for (int c = 0; c < na; ++c) {
		if (c == k) continue;
		acc += gamma[c] * z[c];
	}
	return (-acc) / gk;
}

/* Measures the CPU time per LP of a pass (ADAPTIVE). */
struct AdaptGuard {
	bool on; double t0; long lp0; const SCSimplex& sx; double* ema; int strat;
	AdaptGuard(bool on, const SCSimplex& sx, double* ema, int strat) :
		on(on), t0(on ? cpu_time() : 0.0), lp0(on ? sx.nb_solves() : 0), sx(sx), ema(ema), strat(strat) { }
	~AdaptGuard() {
		if (!on || strat < 0) return;
		const long lps = sx.nb_solves() - lp0;
		if (lps <= 0) return;
		const double c = (cpu_time() - t0) / lps;
		double& e = ema[strat];
		e = (e < 0.0) ? c : 0.9 * e + 0.1 * c;
	}
};

} // anonymous namespace

CtcSC::CtcSC(Linearizer& lr, WarmStart ws, double tol, int max_iter) :
		Ctc(lr.nb_var()), lr(lr), nb_var(lr.nb_var()), ws(ws), max_iter(max_iter),
		mylinearsolver(lr.nb_var(), LPSolver::Mode::NotCertified),
		sx(tol, tol), A(1,1), Af(1,1), z(1), gamma(1),
		n_adapt(0), current(-1), n_passes(0), n_lps(0), n_pivots(0), n_nan(0), n_dropped(0) {
	ema[0] = ema[1] = -1.0;
	switch (ws) {
	case DUAL:     sx.set_strategy(true, false, false); break;
	case PRIMAL:   sx.set_strategy(false, true, false); break;
	case ADAPTIVE: sx.set_strategy(true, false, true); break;
	}
}

CtcSC::~CtcSC() {
}

void CtcSC::add_property(const IntervalVector& init_box, BoxProperties& prop) {
	lr.add_property(init_box, prop);
}

void CtcSC::contract(IntervalVector& box) {
	ContractContext context(box);
	contract(box, context);
}

void CtcSC::contract(IntervalVector& box, ContractContext& context) {
	if (box.is_unbounded()) return;

	const int m = linearize(box, context);
	if (m == -1) {                 // the relaxation is infeasible
		box.set_empty();
		mylinearsolver.clear_constraints();
		return;
	}
	if (m == 0) return;

	++n_passes;
	pass(box);
	mylinearsolver.clear_constraints();
	context.prop.update(BoxEvent(box, BoxEvent::CONTRACT));
}

int CtcSC::linearize(const IntervalVector& box, ContractContext& context) {
	mylinearsolver.clear_constraints();
	if (box.is_unbounded()) return 0;

	const int m0 = lr.linearize(box, mylinearsolver, context.prop);
	if (m0 <= 0) return m0;

	// rows of the linearizer (the first nb_var rows of the LPSolver are the bounds)
	const Matrix rows0 = mylinearsolver.rows();
	const IntervalVector lhs_rhs0 = mylinearsolver.lhs_rhs();

	/* 1. Consecutive identical rows are merged (the affine linearizer writes
	 * each equation as two one-sided rows). */
	int m = 0;
	std::vector<int> src;   // row of the LPSolver giving the coefficients of each row
	IntervalVector sides(m0);
	for (int i = 0; i < m0; ++i) sides[i] = lhs_rhs0[nb_var+i];
	for (int i = 0; i < m0; ++i) {
		if (m > 0) {
			const int p = src[m-1];
			bool same = true;
			for (int j = 0; j < nb_var && same; ++j) same = (rows0[nb_var+i][j] == rows0[nb_var+p][j]);
			const Interval& a = sides[m-1]; const Interval& b = lhs_rhs0[nb_var+i];
			const double lo = std::max(a.lb(), b.lb()), hi = std::min(a.ub(), b.ub());
			if (same && lo <= hi) {
				sides[m-1] = Interval(lo, hi);
				continue;
			}
		}
		src.push_back(i); sides[m] = lhs_rhs0[nb_var+i]; ++m;
	}

	/* 2. Proportional rows (same support, ratios equal within TWIN_TOL) are
	 * merged into one row whose side is the intersection of the scaled sides. */
	if (m > 1) {
		std::vector<char> alive(m, 1);
		std::map<std::vector<int>, std::vector<int> > groups;
		std::vector<std::vector<int> > supp(m);
		for (int i = 0; i < m; ++i) {
			const Vector& r = rows0[nb_var+src[i]];
			for (int j = 0; j < nb_var; ++j)
				if (r[j] != 0.0) supp[i].push_back(j);
		}
		for (int i = 0; i < m; ++i) {
			if (supp[i].empty()) continue;
			std::vector<int>& g = groups[supp[i]];
			for (size_t t = 0; t < g.size() && alive[i]; ++t) {
				const int q = g[t];
				if (!alive[q]) continue;
				const Vector& a = rows0[nb_var+src[q]];
				const Vector& b = rows0[nb_var+src[i]];
				const int j0 = supp[i][0];
				const double c = b[j0] / a[j0];
				bool prop = true;
				for (size_t u = 0; u < supp[i].size(); ++u) {
					const int j = supp[i][u];
					if (std::fabs(b[j] - c*a[j]) > TWIN_TOL * std::max(std::fabs(b[j]), std::fabs(c*a[j]))) prop = false;
				}
				if (prop) {
					/* b = c*a + d with d = b - c*a (enclosed): b.x in S_i gives
					 * a.x in (S_i - d.x)/c over the box. The residual d.x is not
					 * negligible when the box is large. */
					Interval res(0.0);
					for (size_t u = 0; u < supp[i].size(); ++u) {
						const int j = supp[i][u];
						res += (Interval(b[j]) - Interval(c) * a[j]) * box[j];
					}
					const Interval inter = sides[q] & ((sides[i] - res) / Interval(c));
					if (!inter.is_empty()) { sides[q] = inter; alive[i] = 0; }
				}
			}
			if (alive[i]) g.push_back(i);
		}
		int k = 0;
		for (int i = 0; i < m; ++i)
			if (alive[i]) { src[k] = src[i]; sides[k] = sides[i]; ++k; }
		m = k;
	}

	/* 3. Each row of Abar has the coefficient -1 of its slack. A row with a
	 * coefficient above COEF_MAX next to it (e.g. the slope 1/(lb*ub) of 1/w,
	 * with lb near 1e-308) cannot be pivoted in double precision: B^-1 and x_B
	 * overflow and the simplex returns NaN multipliers. Such a row is dropped,
	 * which only enlarges the polytope. Tiny coefficients are harmless. */
	{
		int k = 0;
		for (int i = 0; i < m; ++i) {
			const Vector& r = rows0[nb_var+src[i]];
			double amax = 0.0;
			for (int j = 0; j < nb_var; ++j) {
				const double v = std::fabs(r[j]);
				if (!(v <= amax)) amax = v;   // also catches NaN
			}
			if (!(amax <= COEF_MAX)) { ++n_dropped; continue; }
			src[k] = src[i]; sides[k] = sides[i]; ++k;
		}
		m = k;
	}

	// 4. Abar = [A | -I] and [z] = [box ; slack intervals]
	A.resize(m, nb_var + m);
	A.clear();
	z.resize(nb_var + m);
	for (int i = 0; i < nb_var; ++i) z[i] = box[i];

	for (int i = 0; i < m; ++i) {
		Interval enclosure(0.0);
		for (int j = 0; j < nb_var; ++j)
			enclosure += Interval(rows0[nb_var+src[i]][j]) * box[j];
		Interval side = sides[i];
		if (side.lb() <= -1e100) side = Interval(NEG_INFINITY, side.ub());
		if (side.ub() >=  1e100) side = Interval(side.lb(), POS_INFINITY);
		const Interval bi = side & enclosure;
		z[nb_var+i] = bi.is_empty() ? enclosure : bi;
	}
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < nb_var; ++j) A[i][j] = rows0[nb_var+src[i]][j];
		A[i][nb_var+i] = Interval(-1.0);
	}
	return m;
}

void CtcSC::choose_strategy() {
	const long nn = ++n_adapt;
	int est;
	if (nn <= 2*ADAPT_EXPLORE || ema[0] < 0.0 || ema[1] < 0.0) est = (int) (nn % 2);
	else {
		const int best = (ema[0] <= ema[1]) ? 0 : 1;
		est = (nn % ADAPT_RETRY == 0) ? 1 - best : best;
	}
	current = est;
	if (est == 0) sx.set_strategy(true, false, true);
	else          sx.set_strategy(false, true, true);
}

void CtcSC::pass(IntervalVector& box) {
	AdaptGuard guard(ws == ADAPTIVE, sx, ema, ws == ADAPTIVE ? -2 : -1);

	const int m  = A.nb_rows();
	const int na = nb_var + m;
	const int nc = 2*nb_var;
	const int tope = max_iter > 0 ? max_iter : 100;

	Af.resize(m, na);
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < na; ++j) {
			const double l = A[i][j].lb(), u = A[i][j].ub();
			Af[i][j] = (l == u) ? l : A[i][j].mid();
		}
	if (ws == ADAPTIVE) { choose_strategy(); guard.strat = current; }
	sx.load(Af, nb_var);

	done.assign(nc, 0);
	xstar.assign(na, 0.0);
	gamma.resize(na);

	for (int step = 0; step < nc; ++step) {
		// choice of the next bound: 2j = lower bound of x_j, 2j+1 = upper bound
		int idx = -1;
		if (step == 0) idx = 0;
		else {
			double best = POS_INFINITY;
			for (int j = 0; j < nb_var; ++j)
				for (int side = 0; side < 2; ++side) {
					const int q = 2*j + side;
					if (done[q]) continue;
					const double bd = side ? z[j].ub() : z[j].lb();
					const double dj = std::fabs(xstar[j] - bd);
					const double rel = (std::fabs(bd) < 1.0) ? dj : std::fabs(dj/bd);
					if (rel < PREC) { done[q] = 1; continue; }   // already attained
					if (dj < best) { best = dj; idx = q; }       // nearest first
				}
		}
		if (idx < 0) break;
		done[idx] = 1;
		const int k = idx/2, side = idx%2;

		++n_lps;
		SCSimplex::Status st = sx.solve(k, side == 1, z, tope, true);

		/* An infeasibility whose ray does not certify emptiness: solve once
		 * more from the slack basis. */
		if (st == SCSimplex::INFEASIBLE) {
			const std::vector<double>& l0 = sx.lambda();
			IntervalVector g0(na, Interval(0.0));
			for (int j = 0; j < m; ++j) {
				if (l0[j] == 0.0) continue;
				const Interval lj(l0[j]);
				for (int c = 0; c < na; ++c)
					if (A[j][c] != Interval(0.0)) g0[c] += lj * A[j][c];
			}
			Interval s0(0.0);
			for (int c = 0; c < na; ++c) s0 += g0[c] * z[c];
			if (s0.is_empty() || s0.contains(0.0)) {
				n_pivots += sx.iterations();
				st = sx.solve(k, side == 1, z, tope, false);
			}
		}
		n_pivots += sx.iterations();

		// certificate gamma = Abar^T lambda, in interval arithmetic
		for (int c = 0; c < na; ++c) gamma[c] = Interval(0.0);
		const std::vector<double>& lam = sx.lambda();
		for (int j = 0; j < m; ++j) {
			if (lam[j] == 0.0) continue;
			const Interval lj(lam[j]);
			for (int c = 0; c < na; ++c)
				if (A[j][c] != Interval(0.0)) gamma[c] += lj * A[j][c];
		}
		Interval sum(0.0);
		for (int c = 0; c < na; ++c) sum += gamma[c] * z[c];
		/* An empty sum comes from non-finite multipliers (NaN): it is not a
		 * certificate. Only a nonempty sum without 0 proves emptiness. */
		if (sum.is_empty()) { ++n_nan; continue; }
		if (!sum.contains(0.0)) { box.set_empty(); return; }   // emptiness proved

		const bool optimal = (st == SCSimplex::OPTIMAL);
		if (!optimal && st != SCSimplex::ITER_LIMIT) continue;
		if (optimal) sx.primal_solution(z, xstar);

		const Interval nb = bound_from_gamma(gamma, z, k);
		if (nb.is_empty()) continue;
		if (!z[k].intersects(nb)) { box.set_empty(); return; }
		z[k] &= nb;
	}
	for (int i = 0; i < nb_var; ++i) box[i] = z[i];
}

} // namespace ibex
