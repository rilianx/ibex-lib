/*
 * ibex_Conditioners.h
 *
 *  Created on: 13-10-2017
 *      Author: victor
 */

#ifndef SRC_NUMERIC_IBEX_CONDITIONERS_H_
#define SRC_NUMERIC_IBEX_CONDITIONERS_H_

#include <list>
#include <set>
#include "ibex_IntervalMatrix.h"
#include "ibex_Matrix.h"
#include "ibex_Linear.h"
#include "ibex_Expr.h"
#include <algorithm>
using namespace std;

namespace ibex {

    /**
     * TODO: add comments, cambiar prototipo de funcion
     * \brief This function performs the Gauss-Jordan elimination to the matrix A. All the row operations are stored
     * on the matrix P .
     *
     */

	pair<int,int> find_next_pivot(Matrix & PA, IntervalVector x,set<int> & ban_rows, set<int> & ban_cols);
	void best_gauss_jordan_coll (IntervalMatrix A, IntervalVector x, vector<IntervalMatrix> & perm_list,
							vector <vector <pair <int,int> > > & proj_vars, double prec);
	Matrix gauss_jordan (IntervalMatrix& A, double prec);

	void gauss_jordan_all (IntervalMatrix& A, vector<Matrix>& permutations,vector < vector < pair <int, int> > > &pair_contr_all , double prec);

	Matrix gauss_jordan_test (IntervalMatrix& A, vector<pair <double, int > > order_cols,vector<pair <int, int> >& pair_contr,double prec);

	void new_pseudoinverse(Array <const ExprNode>& xn, Array <const ExprNode>& bn, Matrix PA_aux,
			IntervalMatrix & A, double prec);
    /**
     * \brief This function computes the pseudoinverse of the Matrix A. It also obtains the conditioner matrix P.
     */
	bool pseudoinverse(Matrix A, IntervalMatrix& P);

	void combinatorial(IntervalMatrix A,int cols,int rows,std::vector< std::vector <int> > & comb_piv);

	Matrix best_gauss_jordan (IntervalMatrix& A, IntervalVector box, double prec);

	bool max_piv(Matrix B, set<int> & ban_rows, set<int> & ban_cols, pair<int,int> & max_values,double prec);

}



#endif /* SRC_NUMERIC_IBEX_CONDITIONERS_H_ */
