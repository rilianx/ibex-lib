#include "ibex_CtcDFBPropag.h"
#include <cstdlib>
#include "ibex_LinearizerXTaylor.h"
#include <cmath>
#include <tuple>

using namespace std;
using namespace ibex;


CtcDFBPropag::CtcDFBPropag(ExtendedSystem& sys, Linearizer& lr, double ratio, bool stand_alone, bool only_hc4, bool only_dfb): Ctc(lr.nb_var()), lr(lr), 
mylineardummysolver(nb_var, LPSolver::Mode::Certified, LPSolver::default_tolerance, LPSolver::default_timeout, LPSolver::default_max_iter), refbox(1), refA(1,1), refAf(1,1),
 ratio(ratio),  g(sys.nb_ctr, sys.nb_var), stand_alone(stand_alone), sys(sys),
 linearization_infeasible(false), lin_box(1), has_lin(false), call_count(0),
 sweeper(NULL) { 

   // cout << "[CtcDFBPropag] Initializing with LPSolver in Certified mode" << endl;

    //cout << sys << endl;
    sweeper = NULL;
    if (batch_sweep && !only_hc4)
        sweeper = new CtcDFB(nb_var, 0, false, false, sweep_pivots);

    if(!only_hc4){
        for (int i=0; i<lr.nb_var(); i++){
            // Crear punteros dinámicos y almacenarlos en el vector
            /* lazy_pivots: pivotes por visita; -1 = hasta el punto fijo de
             * esa cota en una sola visita. */
            dfb_ctc.push_back(new CtcDFB(nb_var, i, false, false, lazy_pivots));
            dfb_ctc.push_back(new CtcDFB(nb_var, i, true, false, lazy_pivots));
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

  //  cout << "[CtcDFBPropag] Initialization complete. dfb_ctc.size()=" << dfb_ctc.size()
   //         << ", hc4_ctc.size()=" << hc4_ctc.size() << endl;
}

CtcDFBPropag::~CtcDFBPropag() {
    for (size_t c = 0; c < simplex_por_ctc.size(); ++c) delete simplex_por_ctc[c];
    if (sweeper) delete sweeper;
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
//    cout << "[CtcDFBPropag] Linearizing box: " << box << endl;

    linearization_infeasible = false;

    /* Mismo protocolo que CtcPolytopeHull::contract: limpiar el solver antes de
     * linealizar y no linealizar cajas no acotadas (en el sistema extendido de
     * un optimizador la variable objetivo arranca sin cotas). */
    mylineardummysolver.clear_constraints();
    if (box.is_unbounded()) { A.resize(1,1); return; }

    ContractContext context(box);   
    int m = lr.linearize(box, mylineardummysolver, context.prop);

    if (m == -1) {
        /* La relajacion lineal es infactible: la caja no contiene soluciones. */
        linearization_infeasible = true;
        A.resize(1,1);
        return;
    }
    if (m == 0) { A.resize(1,1); return; }   /* nada que linealizar */


    //show list <dynamic_cast>(LinearizerXTaylor) lr.b_ctrs;
    int k=0;
    for (b_constraint& b_ctr : dynamic_cast<LinearizerXTaylor*> (&lr)->b_ctrs) {
        int c = b_ctr.c;
        if (c<0) c = -b_ctr.c;
        c--;

        adj_b[c].push_back(make_pair(nb_var+k, &b_ctr));
        k++;
    }
    

  //  cout << "[CtcDFBPropag] Linearizer returned m=" << m << endl;

    Matrix rows = mylineardummysolver.rows();
//    cout << rows << endl;
    IntervalVector lhs_rhs = mylineardummysolver.lhs_rhs();
//    cout << lhs_rhs << endl;

    //IntervalMatrix A: m*n (n=box.size()+m)
    A.resize(m, nb_var+m);
    //initialize A to 0
    A.clear();

    x.resize(nb_var+m); // x U b
    for (int i=0; i<nb_var; i++) x[i]=box[i];

    //b
    for (int i=0; i<m; i++){
        /* b_i esta definida por b_i = fila_i . x, asi que su rango valido es la
         * evaluacion por intervalos de la fila sobre la caja, intersectada con
         * la cota que reporta el solver LP. Ambas son validas para toda
         * solucion, de modo que su interseccion tambien lo es (y si fuera
         * vacia, eso probaria que la caja no contiene soluciones; aqui se cae
         * al enclosure, que siempre es valido, para no propagar vacuidad desde
         * la linealizacion).
         *
         * Antes se recortaba la cota del solver a +-1e50 cuando era infinita.
         * Eso APRIETA una cota no acotada hasta un valor finito arbitrario y
         * por lo tanto puede excluir soluciones: con coeficientes del orden de
         * 1e72 (linealizaciones de Brown-* sobre cajas de radio 1e9), el valor
         * verdadero de b_i cae muy por fuera de 1e50 y la solucion se perdia. */
        Interval enclosure(0.0);
        for (int j=0; j<nb_var; j++)
            enclosure += Interval(rows[nb_var+i][j]) * box[j];

        Interval bi = lhs_rhs[nb_var+i] & enclosure;
        x[nb_var+i] = bi.is_empty() ? enclosure : bi;
  //      cout << "b[" << i << "]=" << x[nb_var+i] << endl;
    }

    //coefficients
    for (int i=0; i<m; i++){
        for (int j=0; j<nb_var; j++)
            A[i][j] = rows[nb_var+i][j];
        A[i][nb_var+i] = Interval(-1.0);
    }


//    cout << "[CtcDFBPropag] Linearization complete. A dimensions: " << A.nb_rows() << "x" << A.nb_cols() << endl;
   // cout << "[CtcDFBPropag] A Linearized =" << A << endl;
}


void CtcDFBPropag::contract(IntervalVector& box) {
    ContractContext context(box);
    contract(box,context);
}

/* M5(a): los contractores se crean en pares (cota inferior, cota superior) de
 * la misma variable, y el segundo se deduce del primero sin repetir la
 * eliminacion. Ver CtcDFB::init_from_lower. */
void CtcDFBPropag::init_all_dfb(IntervalMatrix& A, IntervalVector& x){
    if (CtcDFB::light_mode) {
        /* Una sola copia en doubles para los 2n contractores. */
        refAf.resize(A.nb_rows(), A.nb_cols());
        for (int i = 0; i < A.nb_rows(); ++i)
            for (int j = 0; j < A.nb_cols(); ++j) refAf[i][j] = A[i][j].mid();
        for (size_t c = 0; c < dfb_ctc.size(); ++c)
            dfb_ctc[c]->init_light(refAf, A, x);
        return;
    }
    if (CtcDFB::simplex_mode) {
        /* Una sola copia flotante de Abar para los 2n contractores, y el
         * simplex se carga una vez por linealizacion. Cada contractor resuelve
         * su propia cota desde la base de las b, asi que no hay estado por
         * contractor que mantener. */
        refAf.resize(A.nb_rows(), A.nb_cols());
        for (int i = 0; i < A.nb_rows(); ++i)
            for (int j = 0; j < A.nb_cols(); ++j) refAf[i][j] = A[i][j].mid();
        static const bool propio = (getenv("DFB_SX_PROPIO") != NULL);
        if (propio) {
            /* Uno por contractor: cada cota conserva SU base entre llamadas. */
            if (simplex_por_ctc.size() != dfb_ctc.size()) {
                for (size_t c = 0; c < simplex_por_ctc.size(); ++c)
                    delete simplex_por_ctc[c];
                simplex_por_ctc.assign(dfb_ctc.size(), (DFBSimplex*)NULL);
                for (size_t c = 0; c < dfb_ctc.size(); ++c)
                    simplex_por_ctc[c] = new DFBSimplex();
            }
            for (size_t c = 0; c < dfb_ctc.size(); ++c) {
                simplex_por_ctc[c]->load(refAf, nb_var);
                dfb_ctc[c]->sx = simplex_por_ctc[c];
            }
        } else {
            simplex.load(refAf, nb_var);
            for (size_t c = 0; c < dfb_ctc.size(); ++c) dfb_ctc[c]->sx = &simplex;
        }
        /* Y ademas la inicializacion normal, porque contract_simplex se cae a
         * la regla vieja cuando el simplex no llega al optimo. */
    }

    for (size_t c = 0; c < dfb_ctc.size(); ++c) {
        bool pareado = CtcDFB::pair_init && c > 0 &&
                       dfb_ctc[c]->upper_contract && !dfb_ctc[c-1]->upper_contract &&
                       dfb_ctc[c]->k == dfb_ctc[c-1]->k;
        if (CtcDFB::float_pivoting) {
            if (pareado) dfb_ctc[c]->init_from_lower_float(*dfb_ctc[c-1], x);
            else         dfb_ctc[c]->init_float(A, x);
        }
        else if (pareado)             dfb_ctc[c]->init_from_lower(*dfb_ctc[c-1], x);
        else                          dfb_ctc[c]->init(A, x);
    }
}

void CtcDFBPropag::init_dfb_contractors(IntervalMatrix& A, IntervalVector& x_ref){
    if (A.nb_rows() > 1){

        refA.resize(A.nb_rows(), A.nb_cols());
        refA=A;
        refbox.resize(x_ref.size());
        refbox=x_ref;

    //    cout << "[CtcDFBPropag] Initializing DFB contractors with refA" << endl;
        init_all_dfb(A, x_ref);
    }else
        refA.resize(1,1);
    
}

/**
 * \brief Contract a box.
 */
bool CtcDFBPropag::b_contraction = false;
bool CtcDFBPropag::trust_linearizer_infeasible = (getenv("DFB_NO_LIN_EMPTY") == NULL);
int  CtcDFBPropag::max_propag_steps = 0;
bool   CtcDFBPropag::batch_sweep   = (getenv("DFB_SWEEP") != NULL);
int    CtcDFBPropag::sweep_pivots  = (getenv("DFB_SWEEP_PIVOTS") ?
                                      atoi(getenv("DFB_SWEEP_PIVOTS")) : 8);
int    CtcDFBPropag::sweep_rounds  = (getenv("DFB_SWEEP_ROUNDS") ?
                                      atoi(getenv("DFB_SWEEP_ROUNDS")) : 3);
long   CtcDFBPropag::n_sweeps      = 0;
int    CtcDFBPropag::dfb_period    = (getenv("DFB_PERIOD") ?
                                      atoi(getenv("DFB_PERIOD")) : 1);
bool   CtcDFBPropag::goal_priority = (getenv("DFB_GOAL_PRIO") != NULL);
int    CtcDFBPropag::goal_iters    = (getenv("DFB_GOAL_ITERS") ?
                                      atoi(getenv("DFB_GOAL_ITERS")) : 4);
bool   CtcDFBPropag::reuse_linearization = (getenv("DFB_REUSE_LIN") != NULL);
double CtcDFBPropag::relin_shrink        = (getenv("DFB_RELIN_SHRINK") ?
                                            atof(getenv("DFB_RELIN_SHRINK")) : 0.5);
long   CtcDFBPropag::n_relin = 0;
long   CtcDFBPropag::n_reuse = 0;
long   CtcDFBPropag::n_reuse_denied = 0;
long   CtcDFBPropag::n_revive = 0;
long   CtcDFBPropag::n_impact_calls = 0;
/* 5 y no 1: la curva en ibexopt es unimodal con maximo en 5 (-7 % celdas y
 * -8 % tiempo en media geometrica sobre 107 instancias, MEDICIONES_PEREZOSO.md
 * §2.3). Y no -1: sin tope por visita el pivoteo flotante puede ciclar y
 * contract_float no retorna, con lo que el limite de tiempo de ibexopt --que se
 * comprueba entre nodos-- no puede actuar; eso costo 72 procesos matados. */
int    CtcDFBPropag::lazy_pivots   = (getenv("DFB_MAXITERS") ?
                                      atoi(getenv("DFB_MAXITERS")) : 5);
int    CtcDFBPropag::pivot_budget  = (getenv("DFB_PIVOT_BUDGET") ?
                                      atoi(getenv("DFB_PIVOT_BUDGET")) : -1);
double CtcDFBPropag::pivot_budget_n = (getenv("DFB_PIVOT_BUDGET_N") ?
                                      atof(getenv("DFB_PIVOT_BUDGET_N")) : 0.0);


/* ================== Barrido por lotes: una base en cadena ================= *
 * Ver batch_sweep en el header. El barrido usa UN tableau para todas las
 * cotas: renormaliza la columna k en la fila 0 y pivotea desde donde quedo la
 * cota anterior. La certificacion garantiza que cualquier lambda de una cota
 * valida, asi que una heuristica imperfecta puede dar cotas debiles pero nunca
 * falsas.
 * ========================================================================= */

void CtcDFBPropag::dfb_sweep(IntervalVector& box) {
    if (!sweeper) return;
    n_sweeps++;

    /* Se arma el tableau una sola vez por barrido (una copia, no 2n). */
    sweeper->init(refA, refbox);
    if (!sweeper->init_ok) return;

    for (int k = 0; k < nb_var; ++k) {
        if (k > 0 && !sweeper->renormalize_to(k)) continue;   /* cambia el objetivo */
        sweeper->k = k;
        sweeper->state = CtcDFB::INITIAL;
        for (int i = nb_var; i < refbox.size(); ++i) sweeper->x_ref[i] = refbox[i];

        /* Los dos lados de x_k desde la misma base: primero se pivotea hacia
         * la cota inferior, luego se niega la fila 0 (flip_side) y se pivotea
         * hacia la superior. Es lo que en la cola hacian dos contractores. */
        for (int lado = 0; lado < 2; ++lado) {
            if (lado == 1) { sweeper->flip_side(); sweeper->state = CtcDFB::INITIAL; }
            for (int r = 0; r < 4 && sweeper->state != CtcDFB::FINAL; ++r) {
                sweeper->contract(box);
                count_dfb += sweeper->iters;
                if (box.is_empty()) return;
                if (sweeper->iters == 0) break;
            }
        }
    }
}

void CtcDFBPropag::hc4_fixpoint(IntervalVector& box, ContractContext& context) {
    if (hc4_ctc.empty()) return;

    std::queue<int> q;
    std::vector<char> en_cola(hc4_ctc.size(), 0);
    for (size_t i = 0; i < hc4_ctc.size(); ++i) { q.push((int)i); en_cola[i] = 1; }

    const long cap = (max_propag_steps > 0) ? max_propag_steps
                     : 50L * (long)(hc4_ctc.size() + 1);
    long pasos = 0;
    IntervalVector old_box(box);

    while (!q.empty()) {
        if (++pasos > cap) break;
        const int c = q.front(); q.pop(); en_cola[c] = 0;

        old_box = box;
        hc4_ctc[c]->contract(box, context);
        count_hc4++;
        if (box.is_empty()) return;

        for (int v = 0; v < nb_var; ++v) {
            if (old_box[v].ratiodelta(box[v]) >= ratio) {
                set<int> ctrs = g.output_ctrs(v);
                for (set<int>::iterator it = ctrs.begin(); it != ctrs.end(); ++it) {
                    if (!en_cola[*it]) { q.push(*it); en_cola[*it] = 1; }
                }
            }
        }
    }
}

void CtcDFBPropag::contract_sweep(IntervalVector& box, ContractContext& context) {
    /* Alternancia gruesa: barrido DFB, punto fijo de HC4, y se repite mientras
     * la caja siga moviendose. Reemplaza el intercalado fino de la cola. */
    for (int ronda = 0; ronda < sweep_rounds; ++ronda) {
        const double antes = box.perimeter();

        dfb_sweep(box);
        if (box.is_empty()) return;

        hc4_fixpoint(box, context);
        if (box.is_empty()) return;

        const double despues = box.perimeter();
        if (!(antes > 0) || !(despues < antes) ||
            (antes - despues) <= ratio * antes) break;
    }
}
void CtcDFBPropag::contract(IntervalVector& box, ContractContext& context){
    adj_b.clear();

    /* Periodo de aplicacion de la parte DFB (ver dfb_period). En las llamadas
     * intermedias no se linealiza ni se inicializan los 2n contractores: solo
     * corre la propagacion HC4, que es lo barato. */
    const bool turno_dfb = (dfb_period <= 1) || (call_count % dfb_period == 0);
    call_count++;

    if(stand_alone && turno_dfb){ //when the contractor is applied stand_alone
        /* Warm start: se reusa la linealizacion vigente solo si la caja actual
         * esta CONTENIDA en la que se linealizo (entonces la relajacion sigue
         * siendo valida) y no encogio tanto como para que quede demasiado
         * floja. Si no, se relineariza y se re-inicializan los contractores. */
        bool reuse = false;
        if (reuse_linearization && has_lin && refA.nb_rows() > 1
            && lin_box.size() == box.size()) {
            if (lin_box.is_superset(box)) {
                const double p0 = lin_box.perimeter(), p1 = box.perimeter();
                reuse = !(p0 > 0 && p1 < relin_shrink * p0);
            } else {
                n_reuse_denied++;   /* caja de otra rama: no es subcaja */
            }
        }

        if (reuse) {
            n_reuse++;
        } else {
            n_relin++;
            //matrix initialization (dfb contractors)
            update_ref(box);

            if (linearization_infeasible && trust_linearizer_infeasible) {
                /* lr.linearize probo que la relajacion es infactible. */
                box.set_empty();
                return;
            }

            if (refA.nb_rows() > 1){
        //        cout << "[CtcDFBPropag] Initializing DFB contractors with refA" << endl;
                init_all_dfb(refA, refbox);
            }
            lin_box.resize(box.size());
            lin_box = box;
            has_lin = (refA.nb_rows() > 1);
        }
    }

    /* Los contractores DFB necesitan la matriz linealizada; los HC4 no.
     * Antes se devolvia la caja intacta cuando la linealizacion no servia
     * (caja no acotada, sin restricciones, o infactible), con lo que la
     * propagacion HC4 no se ejecutaba en absoluto. En un optimizador eso pasa
     * en cada nodo donde la variable objetivo aun no tiene cotas, y el efecto
     * es que --filtering=dfb no filtraba nada: ibexopt terminaba con
     * "possibly unbounded objective" en instancias que hc4 resuelve. */
    const bool dfb_usable = turno_dfb && (refA.nb_rows() > 1);
    if (!dfb_usable && hc4_ctc.empty()) return;

    if (batch_sweep && sweeper) {
        count_dfb = 0; count_hc4 = 0;
        if (dfb_usable) { contract_sweep(box, context); return; }
        hc4_fixpoint(box, context);
        return;
    }
  //  cout << "[CtcDFBPropag] Contracting box: " << box << endl;

    //refbox contains variables x and b

    //cout << "dimension of A: " << refA.nb_rows() << "x" << refA.nb_cols() << endl;

    if (dfb_usable)
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
    if (dfb_usable){
        /* La cola es un min-heap: menor prioridad sale primero. Los HC4 entran
         * con 0.0 y los DFB con 1.0. */
        const int gv = goal_priority ? sys.goal_var() : -1;
        for (auto* dfb : dfb_ctc) {
            //cout << "[CtcDFBPropag] Adding DFB contractor to priority queue " << dfb->k << endl;
            const bool es_objetivo = (gv >= 0 && dfb->k == gv && !dfb->upper_contract);
            pq.push({es_objetivo ? -1.0 : 1.0, pq_order++, dfb});
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

    const long step_cap = (max_propag_steps > 0) ? max_propag_steps
                          : 50L * (long)(hc4_ctc.size() + dfb_ctc.size() + 1);
    long steps = 0;

    /* Tope total de pivotes DFB en esta caja (metodo perezoso). Negativo o cero
     * sin multiplicador = sin tope. */
    budget_spent = 0;
    long budget = -1;
    if (pivot_budget >= 0)              budget = pivot_budget;
    else if (pivot_budget_n > 0.0)      budget = (long)(pivot_budget_n*nb_var + 0.5);

    while (!pq.empty()) {
        if (++steps > step_cap) break;   /* ver max_propag_steps */

        if(CtcDFB* ctc=dynamic_cast<CtcDFB*>(std::get<2>(pq.top()))){
            pq.pop();
            dfb_ctrs.erase({ctc->k, ctc->upper_contract});
            /* Presupuesto agotado: se descarta la cota y se sigue solo con
             * HC4. Lo que ya se contrajo queda: cada pivote dio una cota
             * certificada, no hace falta llegar al punto fijo. */
            if (budget >= 0 && budget_spent >= budget) continue;
    //        cout << "ctc->k=" << ctc->k << endl;
            double error = ctc->get_Aerror();
            if (error > 1e-4) {
     //           cout << "Aerror > 1e-4, regenerating A" << endl;
                ctc->regenerateA(refA);
    //            cout << "A error after regeneration: " << ctc->get_Aerror() << endl;
            }
            old_box = box;
            for (int i = nb_var; i < refbox.size(); ++i){
                ctc->x_ref[i]=refbox[i]; //restore b
            }

            ctc->contract(box);
            count_dfb+=ctc->iters;
            budget_spent += ctc->iters;
        //    cout << ctc->get_virtual_bound() << endl;
            
           // cout << ctc->get_perc_impr(3) << endl;

            if (box.is_empty()) break;

            //cout << "new box=" << box << endl;
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
                    }
                    /* La actualizacion estaba DENTRO del bucle anterior, asi
                     * que no ocurria cuando la variable v no es salida de
                     * ninguna restriccion (por ejemplo la variable objetivo del
                     * sistema extendido). En ese caso la condicion
                     * old_box[v].ratiodelta(box[v])>=ratio quedaba verdadera
                     * para siempre y la cola se realimentaba sin fin: contract
                     * no retornaba y ni el limite de tiempo de ibexopt podia
                     * actuar, porque se comprueba entre nodos. */
                    old_box[v] = box[v];

                    /* Solo tiene sentido reencolar contractores DFB si la
                     * linealizacion sirve: con refA 1x1 sus matrices tienen una
                     * sola columna y real_impact indexa A[0][v] fuera de rango. */
                    for(int j=0; dfb_usable && j<dfb_ctc.size(); j++){
                        double impact = dfb_ctc[j]->real_impact(box[dfb_ctc[j]->k], v, 0.01);
                        n_impact_calls++;
                        if (impact > 0.01 && dfb_ctrs.find({dfb_ctc[j]->k, dfb_ctc[j]->upper_contract})==dfb_ctrs.end()){
                            n_revive++;
                        //    cout << "impact:"<< impact << ", k=" << dfb_ctc[j]->k << endl;
                            pq.push({1.0, pq_order++, dfb_ctc[j]});
                            dfb_ctrs.insert({dfb_ctc[j]->k, dfb_ctc[j]->upper_contract});
                            //dfb_ctc[k]->state = CtcDFB::CONTRACTING;
                        }
                    }
                    
                }
            }

        }

    }
   // cout << "count_dfb=" << count_dfb << " count_hc4=" << count_hc4 << endl;


    //cout << "[CtcDFBPropag] Contracting complete, final box: " << box << endl;
    
}


