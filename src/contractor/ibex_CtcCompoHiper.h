/*
 * ibex_CtcCompoHiper.h
 *
 *  Created on: 06-11-2017
 *      Author: victor
 */

#ifndef SRC_CONTRACTOR_IBEX_CTC_COMPO_HIPER_H_
#define SRC_CONTRACTOR_IBEX_CTC_COMPO_HIPER_H_


#include "ibex_Ctc.h"
#include "ibex_Array.h"
#include <vector>

namespace ibex {

/** \ingroup contractor
 * \brief Composition of contractors
 *
 * For a box [x] the composition of {c_0,...c_n} performs
 * c_n(...(c_1(c_0([x])))).
 */
class CtcCompoHiper : public Ctc {
public:
	/**
	 * \brief build a composition
	 *
	 * If incremental is true, manages the impact.
	 */
	CtcCompoHiper(const Array<Ctc>& list);

	/**
	 * \brief build the composition of c1 and c2.
	 *
	 * If incremental is true, manages the impact.
	 */
	CtcCompoHiper(Ctc& c1, Ctc& c2);

	/**
	 * \brief Create the composition of 3 contractors
	 */
	CtcCompoHiper(Ctc& c1, Ctc& c2, Ctc& c3);

	/**
	 * \brief Create the composition of 4 contractors
	 */
	CtcCompoHiper(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4);

	/**
	 * \brief Create the composition of 5 contractors
	 */
	CtcCompoHiper(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5);

	/**
	 * \brief Create the composition of 6 contractors
	 */
	CtcCompoHiper(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6);
	CtcCompoHiper(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7);
	/**
	 * \brief Delete *this.
	 */
	~CtcCompoHiper();

	/**
	 * \brief Contract a box.
	 */
	virtual void contract(IntervalVector& box);

	/** The list of sub-contractors */
	Array<Ctc> list;
	/*outputs*/
	std::vector<int> ctc_contractions;
	int llamadas;
	std::vector<int> ctc_calls;
	int all_ctc_contract;

protected:
	void init_impacts();

};

} // end namespace ibex




#endif /* SRC_CONTRACTOR_IBEX_CTCCOMPOHIPER_H_ */
