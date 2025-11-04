#include "ibex_CtcDualFeasibleBounding.h"
#include <cmath>
#include <tuple>
#include <set>

using namespace std;
using namespace ibex;

void CtcDFB::init(IntervalMatrix& A, IntervalVector& x_ref){
    //cout << "[CtcDFB] Initializing with k=" << k << ", A dimensions: " 
    //     << A.nb_rows() << "x" << A.nb_cols() << endl;
    this->A.resize(A.nb_rows(), A.nb_cols());
    this->A = A;
    this->x_ref.resize(x_ref.size());
    this->x_ref = IntervalVector(x_ref);
    state = INITIAL;

    if (upper_contract) 
        for (int i = 0; i < A.nb_rows(); ++i) this->A[i][k] = -this->A[i][k];

    makeColumnIdentity(this->A, k, true, 0);
    identity_rows.clear();
 

    //cout << "[CtcDFB] Initialization complete for k=" << k << endl;
}

double CtcDFB::get_virtual_bound(){
    return (upper_contract)? -virtual_x.lb():virtual_x.lb();
}

// impacto de var en la contracción de k (solo si bound de var k es activo)
double CtcDFB::real_impact(Interval& x_k, int var, double eps){
    if (var==k) return 0.0;
    if ((upper_contract && (virtual_x.lb()+eps >= -x_k.ub())) || (!upper_contract && (virtual_x.lb()+eps >= x_k.lb()))){
        return std::abs(A[0][var].mid());
    }

    return 0.0;
}

double CtcDFB::get_Aerror(){
    double error = 0;
    //A[0][0] + A[0][1] + A[0][2] + A[0][3]...
    for (int i=0; i<A[0].size(); i++)
        error += A[0][i].diam();
    
    return error;
}

void CtcDFB::regenerateA(IntervalMatrix& Aref){
    //cout << "[CtcDFB] Regenerating A" << endl;

    list<int> columns;
    for (int i = 0; i < A[0].size(); ++i)
        if (A[0][i]==Interval(0)) 
            columns.push_back(i);

    A=Aref;
    if (upper_contract) 
        for (int i = 0; i < A.nb_rows(); ++i) A[i][k] = -A[i][k];

    makeColumnIdentity(A, k, true, 0);
    identity_rows.clear();

    makeColumnsIdentity(columns);

}

void CtcDFB::makeColumnsIdentity(std::list<int> columns) {
    std::set<int> identity_rows_local;

    columns.push_front(k);
    for (int col : columns) {
        // Buscar la mejor fila para hacer la columna col identidad
        int best_row = -1;
        double max_abs_val = 0.0;

        for (int i = 0; i < A.nb_rows(); ++i) {
            if (identity_rows_local.count(i)) continue;
            double val = A[i][col].mid();  // Puedes usar diam() si prefieres
            if (std::abs(val) > max_abs_val) {
                max_abs_val = std::abs(val);
                best_row = i;
            }
        }

        if (best_row == -1) {
            cout << "[CtcDFB::makeColumnsIdentity] No available row to pivot column " << col << endl;
            continue;
        }

        // Hacer la columna col una identidad en la fila best_row
        makeColumnIdentity(A, col, false, best_row);

        // Guardar fila usada
        identity_rows_local.insert(best_row);
        identity_rows[best_row] = col;
        
        //cout << "[CtcDFB::makeColumnsIdentity] Made column " << col << " identity at row " << best_row << endl;
    }
}


//get delta_impr
double CtcDFB::get_perc_impr(int k){
    //sum last k perc_imprs
    double sum = 0;
    int count = 0;
    for (auto it = perc_imprs.rbegin(); it != perc_imprs.rend() && count < k; ++it, ++count) {
        sum += *it;
    }
    if (count < k) 
        return 1;
    else 
        return sum;
}

void CtcDFB::contract(IntervalVector& x_new) {
    //cout << "[CtcDFB] Contracting box with k=" << k << ", box: " << x_new << endl;
    //if(upper_contract) cout << "[CtcDFB] Contracting upper bound" << endl;
    //else cout << "[CtcDFB] Contracting lower bound" << endl;

    x_new.resize(x_ref.size()); //extended dimension for including vector b
    for (int i = nb_var; i < x_ref.size(); ++i) 
        x_new[i] = x_ref[i]; 

    if(state==INITIAL){
        perc_imprs.clear();
        if (upper_contract)  x_ref[k] = -x_ref[k]; // changeSigns(A, x_new);
        virtual_x = gaussSeidel(x_ref, k, A[0]);
        virtual_x = Interval(virtual_x.lb(), x_ref[k].ub());
        if (upper_contract) x_ref[k] = -x_ref[k]; // changeSigns(A, x_new);
        if (x_ref.is_empty()) return;
    }
    
    int i;
    int j;
    Interval alpha;
    Interval delta, direction;
    state = CONTRACTING;

    if (upper_contract) x_new[k] = -x_new[k]; // changeSigns(A, x_new);        

    iters = 0;
   // cout << "x_new before DFB: " << x_new << endl;

    while (max_iters == -1 || iters < max_iters) {
        tie(j, delta, direction) = largestImpact(A, x_new, A[0]);
       
        if (j == -1 || delta == Interval(0)) {
            if (upper_contract) x_new[k] = -x_new[k]; // changeSigns(A, x_new);  
            state = FINAL;

            x_new.resize(nb_var); //original dimension
            cout << iters << endl;
            return;
        } 

        tie(alpha, i) = calculateAlpha(A[j], A[0], direction);

        if (i == -1) {
            if (upper_contract) x_new[k] = -x_new[k]; // changeSigns(A, x_new);  
            
            x_new.resize(nb_var); //original dimension
            x_new.set_empty(); 
            return;
            //continue;
        }

        A[0] = A[0] + alpha * A[j];
        
        A[0][i] = Interval(0);

        makeColumnIdentity(A, i, false, j);
        A[0][k] = Interval(1);

        double x_lb = virtual_x.lb();

        double old_size = virtual_x.diam();
        virtual_x = gaussSeidel(x_new, k, A[0]);
        virtual_x = Interval(virtual_x.lb(), x_new[k].ub());
        perc_imprs.push_back((virtual_x.lb() - x_lb)/old_size);

        //cout << "x_lb[" << k << "] = " << x_lb << endl;

        identity_rows[j] = i;
        if (contract_all) {
            for (const auto& r : identity_rows) {
                gaussSeidel(x_new, r.second, A[r.first]);
                
            }
        }
        //x_lb += std::abs(alpha.mid() * delta.mid());

        //x_new[k] = x_new[k] & Interval(x_lb, x_new[k].ub());
        ++iters;
    }

    if (upper_contract && !x_new.is_empty()) x_new[k] = -x_new[k]; // changeSigns(A, x_new); 
    
    x_new.resize(nb_var); //original dimension
    
}

std::pair<IntervalVector, IntervalVector> CtcDFB::calculateImpacts(
    const IntervalMatrix& A, const IntervalVector& x_new, const IntervalVector& gamma) {
   
            //cout << "[CtcDFB] Calculating impacts. A dimensions: " << A.nb_rows() << "x" << A.nb_cols()
            //     << ", x_new size: " << x_new.size() << ", gamma size: " << gamma.size() << endl;
    
  //  cout << "x_new" << x_new << endl;
    assert(A.nb_cols() == x_new.size() && "Matrix column count must match x_new size");
            
    int m = A.nb_rows();
    int n = A.nb_cols();

    IntervalVector grad_incr = IntervalVector(m, Interval(0));
    IntervalVector grad_decr = IntervalVector(m, Interval(0));

    for (int j = 1; j < m; ++j) {
        IntervalVector x_prime_incr = IntervalVector(n, Interval(0));
        IntervalVector x_prime_decr = IntervalVector(n, Interval(0));

        for (int i = 0; i < n; ++i) {
            //Interval signOfGammai = sign(gamma[i]);
       //     cout << "gamma[" << i << "] = " << gamma[i] << endl;
            if (gamma[i].lb() >= 0.00001) {
                x_prime_incr[i] = Interval(x_new[i].ub());
                x_prime_decr[i] = Interval(x_new[i].ub());
            } else if (gamma[i].ub() <= -0.00001) {
                x_prime_incr[i] = Interval(x_new[i].lb());
                x_prime_decr[i] = Interval(x_new[i].lb());
            }
            else{
                //Interval signOfAji = sign(A[j][i]);
       //         cout << "A[" << j << "][" << i << "] = " << A[j][i] << endl;
                if (A[j][i].lb() >= 0.00001) {
                    x_prime_incr[i] = Interval(x_new[i].ub());
                    x_prime_decr[i] = Interval(x_new[i].lb());
     //               cout << "case 1" << endl;
                } else if (A[j][i].ub() <= -0.00001) {
                    x_prime_incr[i] = Interval(x_new[i].lb());
                    x_prime_decr[i] = Interval(x_new[i].ub());
       //             cout << "case 2" << endl;
                }
         //       cout << "x_prime_incr[" << i << "] = " << x_prime_incr[i] << endl;
           //     cout << "x_prime_decr[" << i << "] = " << x_prime_decr[i] << endl;
            }
        }

        grad_incr[j] = Interval(A[j] * x_prime_incr);
        grad_incr[j] = - grad_incr[j];
        grad_decr[j] = Interval(A[j] * x_prime_decr);
    }

    return {grad_incr, grad_decr};
}

std::tuple<int, Interval, Interval> CtcDFB::largestImpact(
    const IntervalMatrix& A, const IntervalVector& x_new, const IntervalVector& gamma) {
    IntervalVector grad_incr, grad_decr;
    tie(grad_incr, grad_decr) = calculateImpacts(A, x_new, gamma);
    //cout << "grad_incr: " << grad_incr << endl;
    //cout << "grad_decr: " << grad_decr << endl;

    Interval delta_incr, delta_decr;
    int j_incr, j_decr;
    tie(delta_incr, j_incr) = getMaxValue(grad_incr);
    tie(delta_decr, j_decr) = getMaxValue(grad_decr);

    Interval delta = Interval(0);
    int j;
    Interval direction;

  //  cout << "delta_incr: " << delta_incr << " (j=" << j_incr << "), delta_decr: " 
    if (delta_incr.lb() > delta_decr.lb()) {
        delta = delta_incr;
        j = j_incr;
        direction = Interval(1);
    } else {
        delta = delta_decr;
        j = j_decr;
        direction = Interval(-1);
    }

    //Interval signOfDelta = sign(delta);
    
    
    if (delta.ub() <= 0.000001) {
        j = -1;
    }

    return {j, delta, direction};
}

std::pair<Interval, int> CtcDFB::calculateAlpha(
    const IntervalVector& Aj, const IntervalVector& gamma, const Interval& direction) {
    Interval min_alpha(1e20);
    int min_index = -1;
    int n = Aj.size();

    for (int i = 0; i < n; ++i) {
        if (Aj[i].lb() != 0 && Aj[i].ub() != 0) {
            Interval alpha = (gamma[i] / Aj[i]) * direction;
    //        cout << "Alpha candidate for i=" << i << ": " << alpha << endl;
            //Interval signOfAlpha = sign(alpha);
            if (alpha.ub() <= -0.00001) {
                if (std::abs(alpha.lb()) < std::abs(min_alpha.lb())) {
                    min_alpha = alpha;
                    min_index = i;
                }
            }
        }
    }
  
    min_alpha = min_alpha * direction * Interval(-1.0);
    return {min_alpha, min_index};
}

void CtcDFB::changeSigns(IntervalMatrix& A, IntervalVector& x_new) {
    
    //A[0] = -A[0];
    //A[0][this->k] = -A[0][this->k];

    for (int i = 0; i < A.nb_rows(); ++i) 
        A[i][k] = -A[i][k];
    

    x_new[k] = -x_new[k];
}

                                // AUX FUNCTIONS
int CtcDFB::makeColumnIdentity(IntervalMatrix& A, const int k, bool interchange, 
                               int j){
    int m = A.nb_rows();
    int n = A.nb_cols();

    // Paso 1: Encontrar la fila adecuada para el pivoteo
    if (interchange && (A[j][k].lb() == 0 || A[j][k].ub() == 0)) {
        bool found = false;
        for (int jAux = 0; jAux < m; ++jAux) {
            if (jAux == j) continue; // Ya la estamos evaluando
            if (!(A[jAux][k].lb() == 0 || A[jAux][k].ub() == 0)) {
                // Intercambiar las filas j y jAux
                for (int i = 0; i < n; ++i) {
                    std::swap(A[j][i], A[jAux][i]);
                }
                //j = jAux; // Actualizar j con la fila válida
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("No hay ninguna fila con valor no nulo en la columna k");
        }
    }
    

    // Paso 2: Normalizar la fila j respecto del valor en posición k
    if (A[j][k].lb() == 0 || A[j][k].ub() == 0) {
        throw std::invalid_argument("The divider has a bound equal to 0. Which is not allowed.");
    }
    
    //if (interchange){
        A[j] = (Interval(1) / A[j][k]) * A[j];

        A[j][k] = Interval(1);
    //}

    // Paso 3: Hacer ceros los demás elementos en la columna k
    for (int jj = 0; jj < m; ++jj) {
        if (jj == j) continue;
        Interval factor = A[jj][k];///A[j][k];
        for (int i = 0; i < n; ++i) {
            A[jj][i] -= factor * A[j][i];
        }
        A[jj][k] = Interval(0);
    }

    return j;
}

// Aplica la fórmula de Gauss-Seidel sobre el vector x
Interval CtcDFB::gaussSeidel(IntervalVector& x, int k, IntervalVector& gamma){
    double epsilon = 1e-6;
    int n = gamma.size();
    

    if (k == -1 || k >= n){
        throw std::invalid_argument("Invalid k.");
    }

//    if (gamma[k] != Interval(1)){
//        cout << "gamma[k] = " << gamma[k] << endl;
//        throw std::invalid_argument("Gamma[k] is not 1.");
//    }

    Interval tmp = Interval(gamma[k]);
    gamma[k] = Interval(0);

    Interval xContract = -(IntervalVector(gamma) * IntervalVector(x));
    gamma[k] = Interval(tmp);
    //cout << -gamma << "*" << x << " = " << xContract << endl;
    if (gamma[k] != Interval(1))
        xContract = xContract / gamma[k];

    if (x[k].intersects(xContract))
    {
        if (contract_all){
    //        std::cout << "x[k] = " << x[k] << std::endl;
    //        std::cout << "xContract = " << xContract << std::endl;
        }
        x[k] = x[k] & xContract;
        if (contract_all){
    //        std::cout << "x[k] = " << x[k] << std::endl;
        }
    }
    if (x[k].is_empty()) x.set_empty();

    return xContract;
}

// Encuentra el índice del valor máximo en un vector
std::pair<Interval, int> CtcDFB::getMaxValue(const IntervalVector& vector) {
    Interval max_value = Interval(vector[0]);
    int max_index = 0;
    int n = vector.size();

    for (int i = 1; i < n; ++i) {
        if (vector[i].lb() > max_value.lb()) {
            max_value = Interval(vector[i]);
            max_index = i;
        }
    }

    return {max_value, max_index};
}
