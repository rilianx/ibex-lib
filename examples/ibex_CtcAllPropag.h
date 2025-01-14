#ifndef __IBEX_CTC_ALL_PROPAG_H__
#define __IBEX_CTC_ALL_PROPAG_H__

#include <queue>
#include <map>
#include <set>
#include "ibex_Ctc.h"
#include "ibex_Linearizer.h"
#include "ibex_CtcDualFeasibleBounding.h"

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
class CtcAllPropag : public Ctc {
public:

	/**
	 * \brief 
	 */
	CtcAllPropag(ExtendedSystem& sys, Linearizer& lr, double ratio=0.1): Ctc(lr.nb_var()), lr(lr), 
    mylineardummysolver(nb_var, LPSolver::Mode::Certified), refbox(1), refA(1,1), ratio(0.01),  g(sys.nb_ctr, sys.nb_var) { 

        cout << sys << endl;
        for (int i=0; i<lr.nb_var(); i++){
            dfb_ctc.add(*(new CtcDFB(nb_var, i, false, false)));
            dfb_ctc.add(*(new CtcDFB(nb_var, i, true, false)));
        }

	    for (int i=0; i<sys.nb_ctr; i++) {
		    hc4_ctc.add(*(new CtcFwdBwd(sys,i)));
            ctc2id[&hc4_ctc[i]]=i;
            for (int j=0; j<nb_var; j++) {
                if (hc4_ctc[i].input && (*hc4_ctc[i].input)[j]) g.add_arc(i,j,true);
                if (hc4_ctc[i].input && (*hc4_ctc[i].output)[j]) g.add_arc(i,j,false);
            }
        }


        //sizes
        cout << "dfb_ctc.size()=" << dfb_ctc.size() << endl;
        cout << "hc4_ctc.size()=" << hc4_ctc.size() << endl;
	
    } 

    void update_ref(const IntervalVector& box){
        linearize(box, refA, refbox);
    }

    void linearize(const IntervalVector& box, IntervalMatrix& A, IntervalVector& x){
        ContractContext context(box);
		int m = lr.linearize(box, mylineardummysolver, context.prop);

        Matrix rows = mylineardummysolver.rows();
        IntervalVector lhs_rhs = mylineardummysolver.lhs_rhs();

        //IntervalMatrix A: m*n (n=box.size()+m)
        A.resize(m, nb_var+m);
        //initialize A to 0
        for (int i=0; i<m; i++)
            for (int j=0; j<nb_var+m; j++)
                A[i][j] = Interval(0.0);

        x.resize(nb_var+m); // x U b

        for (int i=0; i<nb_var; i++) x=box[i];

        //b
        for (int i=0; i<m; i++){
            x[nb_var+i] = lhs_rhs[nb_var+i];
            if (x[nb_var+i].lb() < -1e20)
                x[nb_var+i] = Interval(NEG_INFINITY, x[nb_var+i].ub());
            if (x[nb_var+i].ub() > 1e20)
                x[nb_var+i] = Interval(x[nb_var+i].lb(), POS_INFINITY);
        }

        for (int i=0; i<m; i++){
            for (int j=0; j<nb_var; j++)
                A[i][j] = rows[nb_var+i][j];
            A[i][nb_var+i] = Interval(-1.0);
        }
    }

    virtual void contract(IntervalVector& box) {
        ContractContext context(box);
        contract(box,context);
    }
    

	/**
	 * \brief Contract a box.
	 */
	virtual void contract(IntervalVector& box, ContractContext& context){
        //matrix initialization (dfb contractors)
        if(refbox.size()!=box.size() || !refbox.is_superset(box)) {
            update_ref(box);
            for (int i=0; i<dfb_ctc.size(); i++)
                dfb_ctc[i].init(refA);   
        }


        // Comparador para ordenar por el primer valor del par
        auto cmp = [](const std::pair<double, Ctc*>& left, const std::pair<double, Ctc*>& right) {
            return left.first < right.first;  // Orden ascendente (mayor a menor)
        };

        std::priority_queue<std::pair<double, Ctc*>, std::vector<std::pair<double, Ctc*>>, decltype(cmp)> pq(cmp);
        set<int> pq_ctrs;

        //initialize the priority queue
        for (int i=0; i<dfb_ctc.size(); i++){
            //se debería priorizar dfb con bounds activos, es decir donde el bound del intervalo sea igual al objetivo
            //pq.push(std::make_pair(1.5, &dfb_ctc[i]));
            //sino
            pq.push(std::make_pair(0.0, &dfb_ctc[i]));
        }

        for (int i=0; i<hc4_ctc.size(); i++){
            pq.push(std::make_pair(1.0, &hc4_ctc[i]));
            pq_ctrs.insert(i);
        }

        IntervalVector old_box(box);

        while (!pq.empty()) {            

            if(CtcDFB* ctc=dynamic_cast<CtcDFB*>(pq.top().second)){
                pq.pop();
                ctc->contract(box);
                if (old_box[ctc->k].ratiodelta(box[ctc->k])>=ratio){
                    set<int> ctrs=g.output_ctrs(ctc->k);
                    for (set<int>::iterator c=ctrs.begin(); c!=ctrs.end(); c++) {
                        //si c no está en pq
                        if(pq_ctrs.find(*c)==pq_ctrs.end()){
                            pq.push(std::make_pair(1.0, &hc4_ctc[*c]));
                            pq_ctrs.insert(*c);
                        }
                    }
                    old_box[ctc->k] = box[ctc->k];
                 }

            }else{
                Ctc* ctc_ = pq.top().second;  pq.pop();
                ctc_->contract(box, context);
                pq_ctrs.erase(ctc2id[ctc_]);

                for (int v=0; v<nb_var; v++){
                    if (old_box[v].ratiodelta(box[v])>=ratio) {
                        //cout << "v=" << v << endl;
                        set<int> ctrs=g.output_ctrs(v);
                        for (set<int>::iterator c=ctrs.begin(); c!=ctrs.end(); c++) {
                            //si c no está en pq
                            //cout << "c=" << *c << endl;
                            if(pq_ctrs.find(*c)==pq_ctrs.end()){
                                pq.push(std::make_pair(1.0, &hc4_ctc[*c]));
                                pq_ctrs.insert(*c);
                            }
                        
                            old_box[v] = box[v];
                        }
                        
                    }
		        }

            }

           

        }
        


    }



	/**
	 * \brief The linearization technique
	 */
	Linearizer& lr;

    /**
	 * \brief  The linear solver that will be used
	 */
	LPSolver mylineardummysolver;

    Array< CtcDFB > dfb_ctc;
    Array< Ctc > hc4_ctc;

    double ratio;

    IntervalMatrix refA;
    IntervalVector refbox;

    DirectedHyperGraph g; // constraint network (hypergraph)

    //ctc2id
    std::map<Ctc*, int> ctc2id;

};

} // namespace ibex
#endif // __IBEX_CTC_DFB_H__
