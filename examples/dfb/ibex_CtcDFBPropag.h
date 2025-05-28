#ifndef __IBEX_CTC_ALL_PROPAG_H__
#define __IBEX_CTC_ALL_PROPAG_H__

#include <queue>
#include <map>
#include <set>
#include <vector> // Para std::vector
#include "ibex_Ctc.h"
#include "ibex_Linearizer.h"
#include "ibex_CtcDualFeasibleBounding.h"
#include "ibex_DirectedHyperGraph.h"
#include "ibex_CtcFwdBwd.h"
#include "ibex_ExtendedSystem.h"
#include "ibex_LinearizerXTaylor.h"


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
class CtcDFBPropag : public Ctc {
public:

    CtcDFBPropag(ExtendedSystem& sys, Linearizer& lr, double ratio=0.1, bool stand_alone=true, bool only_hc4=false, bool only_dfb=false);

    virtual ~CtcDFBPropag();

    void update_ref(const IntervalVector& box);

    void linearize(const IntervalVector& box, IntervalMatrix& A, IntervalVector& x);

    virtual void contract(IntervalVector& box);

    void init_dfb_contractors(IntervalMatrix& A, IntervalVector& x_ref);

    double compute_rhs_ub(b_constraint& b_ctr, const IntervalVector& box);

    

    /**
     * \brief Contract a box.
     */
    virtual void contract(IntervalVector& box, ContractContext& context);


    /**
     * \brief The linearization technique
     */
    Linearizer& lr;

    /**
     * \brief  The linear solver that will be used
     */
    LPSolver mylineardummysolver;

    std::vector<CtcDFB*> dfb_ctc; // Vector de punteros a CtcDFB
    std::vector<Ctc*> hc4_ctc;    // Vector de punteros a Ctc

    double ratio;

    IntervalMatrix refA;
    IntervalVector refbox;

    DirectedHyperGraph g; // constraint network (hypergraph)

    //ctc2id
    std::map<Ctc*, int> ctc2id;

    bool stand_alone;

    int count_dfb;
    int count_hc4;

    list< pair<int, double> > history; //dfb contractors

    ExtendedSystem& sys;

    map<int, list <pair<int, b_constraint*> > > adj_b; // x -> list(b_id,b_constraint)

    static bool b_contraction;


};

} // namespace ibex
#endif // __IBEX_CTC_ALL_PROPAG_H__
