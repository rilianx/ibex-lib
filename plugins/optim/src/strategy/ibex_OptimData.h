//============================================================================
//                                  I B E X                                   
// File        : ibex_OptimData.h
// Author      : Jordan Ninin, Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Oct 18, 2014
//============================================================================

#ifndef __IBEX_OPTIM_DATA_H__
#define __IBEX_OPTIM_DATA_H__

#include "ibex_BoxProp.h"
#include "ibex_Interval.h"
#include "ibex_System.h"

namespace ibex {

/**
 * \ingroup strategy
 *
 * \brief Data required for the Optimizer
 */
class OptimData : public BoxProp {
public:
	/**
	 * \brief Constructor for the root node (followed by a call to init_root).
	 */
	OptimData();

	/**
	 * \brief Delete *this.
	 */
	~OptimData();

	/**
	 * \brief Create a copy
	 */
	virtual BoxProp* copy() const;

	/**
	 * \brief Initialize the value of "pf"
	 *
	 * \see #pf
	 */
	void compute_pf(const Function& goal, const IntervalVector& box);

	/**
	 * \brief Initialize the value of "pu"
	 *
	 * \see #pu
	 */
	void compute_pu(const System& sys, const IntervalVector& box);

	/**
	 * \brief Casado criterion
	 *
	 * Image of the objective on the current box
	 */
	Interval pf;

	/**
	 * \brief Casado criterion
	 *
	 * Constraint factor of the current box : between 0 infeasible and 1 for all constraints satisfied.
	 */
	double pu;

	static const long prop_id;

protected:

	/**
	 * \brief Constructor by copy.
	 */
	explicit OptimData(const OptimData& e);


};

/*================================== inline implementations ========================================*/

inline BoxProp* OptimData::copy() const {
	return new OptimData(*this);
}

} // end namespace ibex

#endif // __IBEX_OPTIM_DATA_H__
