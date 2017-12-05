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
using namespace std;

namespace ibex {

    /**
     * TODO: add comments, cambiar prototipo de funcion
     * \brief This function performs the Gauss-Jordan elimination to the matrix A. All the row operations are stored
     * on the matrix P.
     *
     */
	Matrix gauss_jordan (IntervalMatrix& A, double prec=1e-7);

	Matrix gauss_jordan_test (IntervalMatrix& A, double prec=1e-7);
    /**
     * \brief This function computes the pseudoinverse of the Matrix A. It also obtains the conditioner matrix P.
     */
	bool pseudoinverse(Matrix A, IntervalMatrix& P);

}



#endif /* SRC_NUMERIC_IBEX_CONDITIONERS_H_ */
