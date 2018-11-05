//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderDefault.h
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Jul 09, 2017
//============================================================================

#ifndef __IBEX_LOUP_FINDER_DEFAULT_H__
#define __IBEX_LOUP_FINDER_DEFAULT_H__

#include "ibex_LoupFinder.h"
#include "ibex_System.h"
#include "ibex_LoupFinderXTaylor.h"
#include "ibex_LoupFinderTrustRegion.h"

namespace ibex {

/**
 * \ingroup optim
 *
 * \brief Default upper-bounding algorithm (for inequalities only).
 *
 * The algorithm uses two upper-bounding techniques:
 * - one based on the search of an inner box:
 *      simple sampling/line probing or in-HC4
 * - one based on the search of an inner polytope:
 *      XTaylor restriction.
 *
 * Note: currently, line probing is disabled unless there is
 * a constraint-free NLP problem (a simple sampling is done
 * otherwise).
 *
 * \note Only works with inequality constraints.
 */
class LoupFinderDefault : public LoupFinder {
public:

	enum mode{xt, abst, both, trustregion} m;

	/**
	 * \brief Create the algorithm for a given system.
	 *
	 * \param sys   - The NLP problem.
	 * \param inHC4 - Flag for building inner boxes. If true, apply inHC4 (inner arithmetic).
	 *                Otherwise, use forward/backward contractor on reversed inequalities.
	 *                Drawbacks of the current implement of inHC4:
	 *                1/ does not work with vector/matrix constraints
	 *                2/ generates symbolically components of the main function (heavy)
	 *
	 */
	LoupFinderDefault(const System& sys,const IntervalVector& initial_box, double alpha, bool inHC4=true, mode m=xt);

	/**
	 * \brief Delete this.
	 */
	virtual ~LoupFinderDefault();

	/**
	 * \brief Find a new loup in a given box.
	 *
	 * \see comments in LoupFinder.
	 */
	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup);

	/*
	 * Loup finder using inner boxes.
	 *
	 * Either HC4 or CtcUnion (of CtcFwdBwd).
	 */
	LoupFinder& finder_probing;

	/**
	 * Loup finder using inner polytopes.
	 */
	LoupFinderXTaylor finder_x_taylor;

	/**
	 * Loup finder using inner polytopes.
	 */
	LoupFinderXTaylor finder_abs_taylor;
	LoupFinderTrustRegion finder_trustregion;
	static string foundby;
	const IntervalVector& initial_box;
	double alpha;

  void print_results(){
		std::cout << (double) finder_x_taylor.nb_opts / finder_x_taylor.nb_calls << " " <<
		(double) finder_abs_taylor.nb_opts / finder_abs_taylor.nb_calls << std::endl;
	}

};

} /* namespace ibex */

#endif /* __IBEX_LOUP_FINDER_DEFAULT_H__ */
