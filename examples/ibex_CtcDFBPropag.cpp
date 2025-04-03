#include "ibex_CtcDFBPropag.h"
#include <cmath>
#include <tuple>

using namespace std;
using namespace ibex;


CtcDFBPropag::CtcDFBPropag(ExtendedSystem& sys, Linearizer& lr, double ratio, bool stand_alone): Ctc(lr.nb_var()), lr(lr), 
mylineardummysolver(nb_var, LPSolver::Mode::Certified), refbox(1), refA(1,1),
 ratio(ratio),  g(sys.nb_ctr, sys.nb_var), stand_alone(stand_alone) { 

    cout << "[CtcDFBPropag] Initializing with LPSolver in Certified mode" << endl;

    cout << sys << endl;
    for (int i=0; i<lr.nb_var(); i++){
        // Crear punteros dinámicos y almacenarlos en el vector
        dfb_ctc.push_back(new CtcDFB(nb_var, i, false, false));
        dfb_ctc.push_back(new CtcDFB(nb_var, i, true, false));
    }

    for (int i=0; i<sys.nb_ctr; i++) {
        // Crear punteros dinámicos y almacenarlos en el vector
        auto* hc4 = new CtcFwdBwd(sys, i);
        hc4_ctc.push_back(hc4);
        ctc2id[hc4] = i;
        for (int j=0; j<nb_var; j++) {
            if (hc4->input && (*hc4->input)[j]) g.add_arc(i, j, true);
            if (hc4->input && (*hc4->output)[j]) g.add_arc(i, j, false);
        }
    }

    cout << "[CtcDFBPropag] Initialization complete. dfb_ctc.size()=" << dfb_ctc.size()
            << ", hc4_ctc.size()=" << hc4_ctc.size() << endl;
}

CtcDFBPropag::~CtcDFBPropag() {
    cout << "[CtcDFBPropag] Destroying instance" << endl;
    for (auto* obj : dfb_ctc) {
        cout << "[CtcDFBPropag] Deleting CtcDFB instance" << endl;
        delete obj;
    }
    for (auto* obj : hc4_ctc) {
        cout << "[CtcDFBPropag] Deleting Ctc instance" << endl;
        delete obj;
    }
    cout << "[CtcDFBPropag] Instance destroyed successfully." << endl;
}

void CtcDFBPropag::update_ref(const IntervalVector& box){
    linearize(box, refA, refbox);
}

void CtcDFBPropag::linearize(const IntervalVector& box, IntervalMatrix& A, IntervalVector& x){
    cout << "[CtcDFBPropag] Linearizing box: " << box << endl;

    ContractContext context(box);
    int m = lr.linearize(box, mylineardummysolver, context.prop);

    cout << "[CtcDFBPropag] Linearizer returned m=" << m << endl;

    Matrix rows = mylineardummysolver.rows();
    IntervalVector lhs_rhs = mylineardummysolver.lhs_rhs();

    //IntervalMatrix A: m*n (n=box.size()+m)
    A.resize(m, nb_var+m);
    //initialize A to 0
    A.clear();

    x.resize(nb_var+m); // x U b
    for (int i=0; i<nb_var; i++) x[i]=box[i];

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


    cout << "[CtcDFBPropag] Linearization complete. A dimensions: " << A.nb_rows() << "x" << A.nb_cols() << endl;
}

void CtcDFBPropag::contract(IntervalVector& box) {
    ContractContext context(box);
    contract(box,context);
}

void CtcDFBPropag::init_dfb_contractors(IntervalMatrix& A, IntervalVector& x_ref){
    refA.resize(A.nb_rows(), A.nb_cols());
    refA=A;
    refbox.resize(x_ref.size());
    refbox=x_ref;

    if (A.nb_rows() > 0){
        cout << "[CtcDFBPropag] Initializing DFB contractors with refA" << endl;
        for (auto* dfb : dfb_ctc) dfb->init(A, x_ref); 
    }
}

/**
 * \brief Contract a box.
 */
void CtcDFBPropag::contract(IntervalVector& box, ContractContext& context){
    cout << "[CtcDFBPropag] Contracting box: " << box << endl;
    
    if(stand_alone){ //when the contractor is applied stand_alone
        //matrix initialization (dfb contractors)
        update_ref(box);

        if (refA.nb_rows() > 0){
            cout << "[CtcDFBPropag] Initializing DFB contractors with refA" << endl;
            for (auto* dfb : dfb_ctc) dfb->init(refA, refbox);
            
        }
    }
    //refbox contains variables x and b

    //cout << "dimension of A: " << refA.nb_rows() << "x" << refA.nb_cols() << endl;

    for (auto* dfb : dfb_ctc) dfb->state = CtcDFB::INITIAL;
    
    // Comparador para ordenar por el primer valor del par
    auto cmp = [](const std::tuple<double, size_t, Ctc*>& a, 
                    const std::tuple<double, size_t, Ctc*>& b) {
        // Min-heap: el que tenga menor prioridad sale primero
        if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) > std::get<0>(b);
        return std::get<1>(a) > std::get<1>(b); // desempate por orden de llegada
    };

    std::priority_queue<std::tuple<double, size_t, Ctc*>, 
                        std::vector<std::tuple<double, size_t, Ctc*>>, decltype(cmp)> pq(cmp);
    set<int> pq_ctrs;

    size_t pq_order = 0;

    //initialize the priority queue
    if (refA.nb_rows() > 0){
        for (auto* dfb : dfb_ctc) {
            //cout << "[CtcDFBPropag] Adding DFB contractor to priority queue " << dfb->k << endl;
            pq.push({1.0, pq_order++, dfb});
        }
    }

    for (size_t i = 0; i < hc4_ctc.size(); i++) {
        pq.push({0.0, pq_order++, hc4_ctc[i]});
        pq_ctrs.insert(i);
    }

    IntervalVector old_box(box);

    while (!pq.empty()) {            

        if(CtcDFB* ctc=dynamic_cast<CtcDFB*>(std::get<2>(pq.top()))){
            pq.pop();
            cout << "ctc->k=" << ctc->k << endl;
            ctc->contract(box);
            if (box.is_empty()) return;

            cout << "new box=" << box << endl;
            if (old_box[ctc->k].ratiodelta(box[ctc->k])>=ratio){
                set<int> ctrs=g.output_ctrs(ctc->k);
                for (set<int>::iterator c=ctrs.begin(); c!=ctrs.end(); c++) {
                    //si c no está en pq
                    if(pq_ctrs.find(*c)==pq_ctrs.end()){
                        pq.push({1.0, pq_order++, hc4_ctc[*c]});
                        pq_ctrs.insert(*c);
                    }
                }
                old_box[ctc->k] = box[ctc->k];
                }

                if (ctc->state == CtcDFB::CONTRACTING)
                    pq.push({0.0, pq_order++, ctc});
                

        }else{
            Ctc* ctc_ = std::get<2>(pq.top());  pq.pop();
            ctc_->contract(box, context);
            if (box.is_empty()) return;
            
            pq_ctrs.erase(ctc2id[ctc_]);

            for (int v=0; v<nb_var; v++){
                if (old_box[v].ratiodelta(box[v])>=ratio) {
                    //cout << "v=" << v << endl;
                    set<int> ctrs=g.output_ctrs(v);
                    for (set<int>::iterator c=ctrs.begin(); c!=ctrs.end(); c++) {
                        //si c no está en pq
                        //cout << "c=" << *c << endl;
                        if(pq_ctrs.find(*c)==pq_ctrs.end()){
                            pq.push({1.0, pq_order++, hc4_ctc[*c]});
                            pq_ctrs.insert(*c);
                        }
                    
                        old_box[v] = box[v];
                    }
                    
                }
            }

        }

    }
    //cout << "[CtcDFBPropag] Contracting complete, final box: " << box << endl;
    
}


