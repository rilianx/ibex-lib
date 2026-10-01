//============================================================================
//                                  I B E X                                   
// File        : ibex_LargestFirst.cpp
// Author      : Bertrand Neveu, Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : May 19, 2012
// Last Update : Dec 25, 2017
//============================================================================

#include "ibex_LargestFirst.h"
#include "ibex_NoBisectableVariableException.h"

#include <cstdlib>
using namespace std;

namespace ibex {

LargestFirst::LargestFirst(double prec, double ratio1) : Bsc(prec), ratio(ratio1) {

}

LargestFirst::LargestFirst(const Vector& prec, double ratio1) : Bsc(prec), ratio(ratio1) {

}

  bool LargestFirst::nobisectable(const IntervalVector & box, int i) const {
    return too_small (box, i);
  }


/* Variables con precision 0 primero, comparadas por su diametro; si ninguna
 * es bisecable, las demas por diam/prec. Con un vector de precisiones que
 * contiene ceros, `diam/prec` daba +inf para todas esas variables: empataban,
 * y se bisecaba siempre la primera. Es lo que pasa en ibexopt con eps_x=0. */
BisectionPoint LargestFirst::choose_var(const Cell& cell) {

	const IntervalVector& box=cell.box;

	int var =-1;
	double l=0.0;
	static const bool viejo = (getenv("IBEX_LF_VIEJO") != NULL);  /* restaura el criterio anterior */
	if (!viejo && !uniform_prec()) {
		for (int i=0; i< box.size(); i++)
			if (prec(i)==0 && !nobisectable(box,i) && (var==-1 || box[i].diam()>l)) {
				var=i;
				l=box[i].diam();
			}
		if (var!=-1) return BisectionPoint(var,ratio,true);
	}
	for (int i=0; i< box.size(); i++)	{
	  if (!(nobisectable (box,i))){
			if (var==-1) {
				var=i;
				l = uniform_prec()? box[i].diam() : (box[i].diam()/prec(i));
			}
			else {
				double l_tmp = uniform_prec()? box[i].diam() : (box[i].diam()/prec(i));
				if (l_tmp>l) {
					var = i;
					l = l_tmp;
				}
			}
		}
	}

	if (var !=-1){
		return BisectionPoint(var,ratio,true);
	}
	else {
		throw NoBisectableVariableException();
	}

}

} // end namespace ibex
