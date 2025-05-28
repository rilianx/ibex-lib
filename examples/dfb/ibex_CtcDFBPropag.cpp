#include "ibex_CtcDFBPropag.h"
#include "ibex_LinearizerXTaylor.h"
#include <cmath>
#include <tuple>

using namespace std;
using namespace ibex;


CtcDFBPropag::CtcDFBPropag(ExtendedSystem& sys, Linearizer& lr, double ratio, bool stand_alone, bool only_hc4, bool only_dfb): Ctc(lr.nb_var()), lr(lr), 
mylineardummysolver(nb_var, LPSolver::Mode::Certified), refbox(1), refA(1,1),
 ratio(ratio),  g(sys.nb_ctr, sys.nb_var), stand_alone(stand_alone), sys(sys) { 

    cout << "[CtcDFBPropag] Initializing with LPSolver in Certified mode" << endl;

    cout << sys << endl;
    if(!only_hc4){
        for (int i=0; i<lr.nb_var(); i++){
            // Crear punteros dinámicos y almacenarlos en el vector
            dfb_ctc.push_back(new CtcDFB(nb_var, i, false, false, 1));
            dfb_ctc.push_back(new CtcDFB(nb_var, i, true, false, 1));
        }
    }

    if(!only_dfb){
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
    }

    cout << "[CtcDFBPropag] Initialization complete. dfb_ctc.size()=" << dfb_ctc.size()
            << ", hc4_ctc.size()=" << hc4_ctc.size() << endl;
}

CtcDFBPropag::~CtcDFBPropag() {
    //cout << "[CtcDFBPropag] Destroying instance" << endl;
    for (auto* obj : dfb_ctc) {
        //cout << "[CtcDFBPropag] Deleting CtcDFB instance" << endl;
        delete obj;
    }
    for (auto* obj : hc4_ctc) {
        //cout << "[CtcDFBPropag] Deleting Ctc instance" << endl;
        delete obj;
    }
    //cout << "[CtcDFBPropag] Instance destroyed successfully." << endl;
}

void CtcDFBPropag::update_ref(const IntervalVector& box){
    linearize(box, refA, refbox);
}

double CtcDFBPropag::compute_rhs_ub(b_constraint& b_ctr, const IntervalVector& box){
    //cout << "[CtcDFBPropag] Evaluating corner" << endl;
    IntervalVector corner = box;
    for (int i=0; i<nb_var; i++){   
        if (b_ctr.inf[i]==true)
            corner[i] = box[i].lb();
        else
            corner[i] = box[i].ub();
    }
    
    if(b_ctr.c<0){
        Interval g_corner = -sys.f_ctrs[-b_ctr.c-1].eval(corner);
        return (-g_corner + b_ctr.a*corner).ub();
    }else{
        Interval g_corner = sys.f_ctrs[b_ctr.c-1].eval(corner);
        return (-g_corner + b_ctr.a*corner).ub();
    }
    
}

void CtcDFBPropag::linearize(const IntervalVector& box, IntervalMatrix& A, IntervalVector& x){
    cout << "[CtcDFBPropag] Linearizing box: " << box << endl;

    ContractContext context(box);
    int m = lr.linearize(box, mylineardummysolver, context.prop);


    //show list <dynamic_cast>(LinearizerXTaylor) lr.b_ctrs;
    int k=0;
    for (b_constraint& b_ctr : dynamic_cast<LinearizerXTaylor*> (&lr)->b_ctrs) {
        int c = b_ctr.c;
        if (c<0) c = -b_ctr.c;
        c--;

        adj_b[c].push_back(make_pair(nb_var+k, &b_ctr));
        k++;
    }

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
        cout << "b[" << i << "]=" << x[nb_var+i] << endl;
    }

    //coefficients
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
    if (A.nb_rows() > 1){

        refA.resize(A.nb_rows(), A.nb_cols());
        refA=A;
        refbox.resize(x_ref.size());
        refbox=x_ref;

        cout << "[CtcDFBPropag] Initializing DFB contractors with refA" << endl;
        for (auto* dfb : dfb_ctc) dfb->init(A, x_ref);
    }else
        refA.resize(1,1);
    
}

/**
 * \brief Contract a box.
 */
bool CtcDFBPropag::b_contraction = false;
void CtcDFBPropag::contract(IntervalVector& box, ContractContext& context){    
    adj_b.clear();
    if(stand_alone){ //when the contractor is applied stand_alone
        //matrix initialization (dfb contractors)
        update_ref(box);

        if (refA.nb_rows() > 1){
            cout << "[CtcDFBPropag] Initializing DFB contractors with refA" << endl;
            for (auto* dfb : dfb_ctc) dfb->init(refA, refbox);
            
        }
    }

    if (refA.nb_rows() <=1 ) return;
    cout << "[CtcDFBPropag] Contracting box: " << box << endl;

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
    set< pair<int,bool> > dfb_ctrs;
    
    

    size_t pq_order = 0;

    //initialize the priority queue
    if (refA.nb_rows() > 0){
        for (auto* dfb : dfb_ctc) {
            //cout << "[CtcDFBPropag] Adding DFB contractor to priority queue " << dfb->k << endl;
            pq.push({1.0, pq_order++, dfb});
            dfb_ctrs.insert({dfb->k, dfb->upper_contract});
        }
    }

    for (size_t i = 0; i < hc4_ctc.size(); i++) {
        pq.push({0.0, pq_order++, hc4_ctc[i]});
        pq_ctrs.insert(i);
    }

    IntervalVector old_box(box);

    count_dfb = 0;
    count_hc4 = 0;
    while (!pq.empty()) {            

        if(CtcDFB* ctc=dynamic_cast<CtcDFB*>(std::get<2>(pq.top()))){
            pq.pop();
            dfb_ctrs.erase({ctc->k, ctc->upper_contract});
            cout << "ctc->k=" << ctc->k << endl;
            double error = ctc->get_Aerror();
            if (error > 1e-4) {
                cout << "Aerror > 1e-4, regenerating A" << endl;
                ctc->regenerateA(refA);
                cout << "A error after regeneration: " << ctc->get_Aerror() << endl;
            }
            old_box = box;
            for (int i = nb_var; i < refbox.size(); ++i){
                ctc->x_ref[i]=refbox[i]; //restore b
            }

            ctc->contract(box);
            count_dfb+=ctc->iters;
            cout << ctc->get_virtual_bound() << endl;
            
           // cout << ctc->get_perc_impr(3) << endl;

            if (box.is_empty()) break;

            cout << "new box=" << box << endl;
            if (old_box[ctc->k].ratiodelta(box[ctc->k])>=ratio){
                history.push_back(make_pair(count_dfb, box.perimeter()));
                set<int> ctrs=g.output_ctrs(ctc->k);
                for (set<int>::iterator c=ctrs.begin(); c!=ctrs.end(); c++) {
                    //si c no está en pq
                    if(pq_ctrs.find(*c)==pq_ctrs.end()){
                        pq.push({0.0, pq_order++, hc4_ctc[*c]});
                        pq_ctrs.insert(*c);
                    }
                }
                old_box[ctc->k] = box[ctc->k];

                if (ctc->state == CtcDFB::CONTRACTING){
                    pq.push({1.0, pq_order++, ctc});
                    dfb_ctrs.insert({ctc->k, ctc->upper_contract});
                }
            }else{
                if (ctc->state == CtcDFB::CONTRACTING){// && ctc->get_perc_impr(10) > 0.01){
                    pq.push({1.0, pq_order++, ctc});
                    dfb_ctrs.insert({ctc->k, ctc->upper_contract});
                }
            }
                

        }else{
            Ctc* ctc_ = std::get<2>(pq.top());  pq.pop();
            old_box=box;
            ctc_->contract(box, context);
            count_hc4++;
            if (box.is_empty()) break;
            
            pq_ctrs.erase(ctc2id[ctc_]);
            
            //b contraction
            if(b_contraction)
                for(auto p : adj_b[ctc2id[ctc_]])
                    refbox[p.first] = Interval(-1e50, compute_rhs_ub(*p.second, box));  
                

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

                    for(int j=0; j<dfb_ctc.size(); j++){
                        double impact = dfb_ctc[j]->real_impact(box[dfb_ctc[j]->k], v, 0.01);
                        if (impact > 0.01 && dfb_ctrs.find({dfb_ctc[j]->k, dfb_ctc[j]->upper_contract})==dfb_ctrs.end()){
                            cout << "impact:"<< impact << ", k=" << dfb_ctc[j]->k << endl;
                            pq.push({1.0, pq_order++, dfb_ctc[j]});
                            dfb_ctrs.insert({dfb_ctc[j]->k, dfb_ctc[j]->upper_contract});
                            //dfb_ctc[k]->state = CtcDFB::CONTRACTING;
                        }
                    }
                    
                }
            }

        }

    }
    cout << "count_dfb=" << count_dfb << " count_hc4=" << count_hc4 << endl;


    //cout << "[CtcDFBPropag] Contracting complete, final box: " << box << endl;
    
}


