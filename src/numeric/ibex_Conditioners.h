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

class Conditioners{
	public:
	Matrix gauss_jordan (IntervalMatrix& A, double prec=1e-7);
	void gauss_jordan_collection (IntervalMatrix &A, list<Matrix> &list_P, list<IntervalMatrix> &list_PA, double prec=1e-7);
	bool pseudoinverse(Matrix A, IntervalMatrix& P);
};



}



#endif /* SRC_NUMERIC_IBEX_CONDITIONERS_H_ */
