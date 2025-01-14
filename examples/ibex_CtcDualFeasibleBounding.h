#ifndef __IBEX_CTC_DFB_H__
#define __IBEX_CTC_DFB_H__

#include "ibex_Ctc.h"

namespace ibex {

/**
 * \ingroup contractor
 *
 * \brief Propagation contractor.
 *
 * This class is an implementation of the classical interval variant of the AC3 constraint propagation
 * algorithm.
 *
 */
class CtcDFB : public Ctc {
public:

	/**
	 * \brief Create a DFB contractor for contracting bound(x) in A.x=0
	 */
	CtcDFB(int nb_var, int k, bool upper_contract=false, bool contract_all=false): A(1,1),
    upper_contract(upper_contract), contract_all(contract_all), k(k), Ctc(nb_var) { };    

	/**
	 * \brief Contract a box.
	 */
	virtual void contract(IntervalVector& x_new){

    }

    
    void init(IntervalMatrix A){
        this->A = A;
    }


    IntervalMatrix A;
    bool upper_contract;
    bool contract_all;
    int k;

};

} // namespace ibex
#endif // __IBEX_CTC_DFB_H__
