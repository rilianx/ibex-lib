#ifndef __IBEX_CTC_DFB_MANAGER_H__
#define __IBEX_CTC_DFB_MANAGER_H__

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
#include "ibex_CtcDFBPropag.h"


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
class CtcDFBManager : public Ctc {
    CtcDFBPropag& dfb_propag;
    Ctc& cid;
    Linearizer& lr;
    LPSolver mylineardummysolver;

public:

    CtcDFBManager(CtcDFBPropag& dfb_propag, Ctc& cid, Linearizer& lr) :
        Ctc(lr.nb_var()), dfb_propag(dfb_propag), cid(cid), lr(lr), mylineardummysolver(nb_var, LPSolver::Mode::Certified) {
        cout << "[CtcDFBManager] Initializing with LPSolver in Certified mode" << endl;


    }

    virtual ~CtcDFBManager() {
        cout << "[CtcDFBManager] Destroying instance" << endl;
        cout << "[CtcDFBManager] Instance destroyed successfully." << endl;
    }
    
    void linearize(const IntervalVector& box, IntervalMatrix& A, IntervalVector& x){
        cout << "[CtcDFBManager] Linearizing box: " << box << endl;
    
        ContractContext context(box);
        int m = lr.linearize(box, mylineardummysolver, context.prop);
    
        cout << "[CtcDFBManager] Linearizer returned m=" << m << endl;
    
        Matrix rows = mylineardummysolver.rows();
        IntervalVector lhs_rhs = mylineardummysolver.lhs_rhs();
    
        //IntervalMatrix A: m*n (n=box.size()+m)
        A.resize(m, nb_var+m);
        //initialize A to 0
        A.clear();
    
        x.resize(nb_var+m); // x U b
    
        //b
        for (int i=0; i<m; i++){
            x[nb_var+i] = lhs_rhs[nb_var+i];
            if (x[nb_var+i].lb() < -1e50)
                x[nb_var+i] = Interval(-1e50, x[nb_var+i].ub());
            if (x[nb_var+i].ub() > 1e50)
                x[nb_var+i] = Interval(x[nb_var+i].lb(), 1e50);
        }
    
    
        for (int i=0; i<m; i++){
            for (int j=0; j<nb_var; j++)
                A[i][j] = rows[nb_var+i][j];
            A[i][nb_var+i] = Interval(-1.0);
        }
    
    
        cout << "[CtcDFBManager] Linearization complete. A dimensions: " << A.nb_rows() << "x" << A.nb_cols() << endl;
    }

    void contract(IntervalVector& box ) {
        ContractContext context(box);
        contract(box,context);
    }

    void contract(IntervalVector& box, ContractContext& context) {
        IntervalMatrix A(1,1);
        IntervalVector x(box);
        linearize(box, A, x);

        dfb_propag.init_dfb_contractors(A, x);

        cout << "[CtcDFBManager] dimension of A: " << dfb_propag.refA.nb_rows() << "x" << dfb_propag.refA.nb_cols() << endl;

        cout << "[CtcDFBManager] Contracting box with CID: " << box << endl;
        // Call the contract method of cid
        cid.contract(box, context);

        cout << "[CtcDFBManager] Contracting complete, final box: " << box << endl;
    }
    

};
}

#endif // __IBEX_CTC_DFB_MANAGER_H__
