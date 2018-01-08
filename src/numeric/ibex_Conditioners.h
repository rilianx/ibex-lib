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

using namespace std;

namespace ibex {

    /**
     * TODO: add comments, cambiar prototipo de funcion
     * \brief This function performs the Gauss-Jordan elimination to the matrix A. All the row operations are stored
     * on the matrix P.
     *
     */
	Matrix gauss_jordan (IntervalMatrix& A, double prec=1e-7);

	void gauss_jordan_all (IntervalMatrix& A, double prec);

	void new_pseudoinverse(Array <const ExprNode> & xn,Array <const ExprNode> & bn, Matrix PA_aux, Matrix PA,
			IntervalMatrix & P, IntervalMatrix & A,IntervalVector & b,double prec);
    /**
     * \brief This function computes the pseudoinverse of the Matrix A. It also obtains the conditioner matrix P.
     */
	bool pseudoinverse(Matrix A, IntervalMatrix& P);

	void combinatorial(IntervalMatrix A,int cols,int rows,std::vector< std::vector <int> > & comb_piv);

}



#endif /* SRC_NUMERIC_IBEX_CONDITIONERS_H_ */
