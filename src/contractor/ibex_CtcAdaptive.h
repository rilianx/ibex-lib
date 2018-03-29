//============================================================================
//                                  I B E X
// File        : Composition of contractors
// Author      : Ignacio Araya
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Mar 23, 2018
// Last Update :
//============================================================================

#ifndef __IBEX_CTC_ADAPTIVE_H__
#define __IBEX_CTC_ADAPTIVE_H__

#include "ibex_Ctc.h"
#include "ibex_Array.h"
#include <map>

using namespace std;

namespace ibex {



class AdaptCell : public Backtrackable {
	public:
		/**
		 * \brief Constructor for the root node (followed by a call to init_root).
		 */
		AdaptCell() : p(NULL) {}



		/**
		 * \brief
		 */
		std::pair<Backtrackable*,Backtrackable*> down(const BisectionPoint&){
			AdaptCell* c1 = new AdaptCell();
			AdaptCell* c2 = new AdaptCell();
			c1->p = this;
			c2->p = this;

			return std::pair<Backtrackable*,Backtrackable*>(c1,c2);
		}

		/** parent of the node **/
		AdaptCell* p;

		/** map from ctc,c to period of application**/
		map<pair<int,int>, int > T;

		/** map from ctc,c to number of failures**/
		map<pair<int,int>, int > F;

	protected:

	/**
	 * \brief Create a copy
	 */
	Backtrackable* copy() const { return new AdaptCell(*this);};
};


/** \ingroup contractor
 * \brief Adaptive contractor
 */
class CtcAdaptive : public Ctc {
public:
	/**
	 * \brief build a composition
	 *
	 * If incremental is true, manages the impact.
	 */
	CtcAdaptive(const Array<Ctc>& list, int m, int L=8, bool bf=true);


	/**
	 * \brief Delete *this.
	 */
	~CtcAdaptive() { }

	/**
	 * \brief Contract a box.
	 */
	virtual void contract(IntervalVector& box){
		std::cout << "error: CtcAdaptive::contract(IntervalVector&) is not implemented, you should use CtcAdaptive::contract(Cell&)" << std::endl;
		exit(0);
	}

	/**
	 * \brief Contract a box.
	 */
	virtual void contract(Cell& c);

	AdaptCell* closest_informed_ancestor(Cell& cell, int ctc, int c);

	/** The list of sub-contractors */
	Array<Ctc> list;

	int* gcalls;
	int* calls;
	int* effective_calls;
	int* nb_input_ctr;
	int* nb_act_ctr;

protected:
	int iter;

	/** Number of constraints */
	int nb_ctr;


	int L;

  bool bf;

};

} // end namespace ibex
#endif // __IBEX_CTC_COMPO_H__
