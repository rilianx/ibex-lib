#ifndef __IBEX_CTC_DFB_H__
#define __IBEX_CTC_DFB_H__

#include "ibex_Ctc.h"
#include <map>
#include <list>

using namespace std;

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

    private:
    int max_iters;
    map<int, int> identity_rows;
    
public:

	/**
	 * \brief Create a DFB contractor for contracting bound(x) in A.x=0
	 */
	CtcDFB(int nb_var, int k, bool upper_contract=false, bool contract_all=false, int max_iters=2): A(1,1),
    x_ref(1), upper_contract(upper_contract), contract_all(contract_all), k(k), 
    Ctc(nb_var), max_iters(max_iters) {  };   

	/**
	 * \brief Contract a box.
	 */
	virtual void contract(IntervalVector& x_new) ;

    /**
     * \brief Virtual destructor to clean up resources.
     */
    virtual ~CtcDFB() {
        //cout << "[CtcDFB] Destroying instance with k=" << k << endl;
        identity_rows.clear();
        A.clear();
        //cout << "[CtcDFB] Instance destroyed successfully." << endl;
    }
    
    
    void init(IntervalMatrix& A, IntervalVector& x_ref);

    std::pair<IntervalVector, IntervalVector> calculateImpacts(
        const IntervalMatrix& A, const IntervalVector& x_new, const IntervalVector& gamma);
    
    std::tuple<int, Interval, Interval> largestImpact(
        const IntervalMatrix& A, const IntervalVector& x_new, const IntervalVector& gamma);

    std::pair<Interval, int> calculateAlpha(
        const IntervalVector& Aj, const IntervalVector& gamma, const Interval& direction);

    void changeSigns(IntervalMatrix& A, IntervalVector& x_new);

    int makeColumnIdentity(IntervalMatrix& A, const int k, bool interchange, int j);

    void makeColumnsIdentity(std::list<int> columns);

    Interval gaussSeidel(IntervalVector& x, int k, IntervalVector& gamma);

    std::pair<Interval, int> getMaxValue(const IntervalVector& vector);

    void regenerateA(IntervalMatrix& Aref);

    double get_virtual_bound();

    double get_perc_impr(int iter);

    double get_Aerror();

    double real_impact(Interval& x_k, int var, double eps=0.01);

    IntervalMatrix A;
    IntervalVector x_ref;
    bool upper_contract;
    bool contract_all;
    int k;

    Interval virtual_x;
    list<double> perc_imprs; //porcentajes de mejora

    //enum State
    enum State {
        INITIAL,
        CONTRACTING,
        FINAL
    };
    State state;

};

} // namespace ibex
#endif // __IBEX_CTC_DFB_H__
