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

CellBS::CellBS() : lb(0.0), depth(0), id(0) {

}

CellBS::CellBS(const CellBS& e) : lb(e.lb), depth(e.depth) {

}

std::pair<Backtrackable*,Backtrackable*> CellBS::down() {
	CellBS* c1=new CellBS(*this);
	CellBS* c2=new CellBS(*this);
	c1->depth=depth+1;
	c2->depth=depth+1;
	return std::pair<Backtrackable*,Backtrackable*>(c1,c2);
}

CellBS::~CellBS(){ }

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

template <class T>
Cell* CellSet<T>::random_pop(double lb_max) {
	std::set<Cell*>::iterator it=cset.begin();
	double roulette=0.0;
	int first_depth=(*it)->get<CellBS>().depth;
	while(it!=cset.end() && (*it)->get<CellBS>().lb <=lb_max){
		roulette += std::pow(2.0, (double)(first_depth-(*it)->get<CellBS>().depth));
		it++;
	}
	roulette=(double) (rand()/RAND_MAX) * roulette;
	
	it=cset.begin();
	while(it!=cset.end() && roulette>=0.0){
	  roulette -= std::pow(2.0, (double)(first_depth-(*it)->get<CellBS>().depth));
	  it++;
    }
	
	it--;
	Cell* c = *it;
	cset.erase(it);
	return c;
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
