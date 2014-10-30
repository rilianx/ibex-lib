//============================================================================
//                                  I B E X                                   
// File        : ibex_CellSet.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 12, 2012
// Last Update : May 12, 2012
//============================================================================

#include "ibex_CellSet.h"

namespace ibex {

template <class T>
void CellSet<T>::flush() {
	while (!cset.empty()) {
		delete *cset.begin();
		cset.erase(cset.begin());
	}
}

template <class T>
int CellSet<T>::size() const {
	return cset.size();
}

template <class T>
bool CellSet<T>::empty() const {
	return cset.empty();
}

template <class T>
void CellSet<T>::push(Cell* cell) {
	if (capacity>0 && size()==capacity) throw CellBufferOverflow();
	cset.insert(cell);
}

template <class T>
Cell* CellSet<T>::pop() {
	Cell* c = *cset.begin();
	cset.erase(cset.begin());
	return c;
}

template <class T>
Cell* CellSet<T>::top() const {
	return *cset.begin();
}

  // E.g.: called in Optimizer in case of a new upper bound
// on the objective ("loup"). This function then removes (and deletes) from
// the heap all the cells with a cost greater than loup.
template <class T>
void CellSet<T>::contract_heap(double loup, int y){
  std::set<Cell*>::iterator it0=cset.begin();

  while(it0!=cset.end()) {
	  if((*it0)->box[y].lb() > loup){
	     delete *it0;
         cset.erase(it0++);
      }else
         it0++;
  }
}

template class CellSet<minLB>;

} // end namespace ibex
