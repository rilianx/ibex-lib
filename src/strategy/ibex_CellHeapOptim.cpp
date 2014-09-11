//============================================================================
//                                  I B E X                                   
// File        : ibex_CellHeapOptim.cpp
// Author      : Gilles Chabert, Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 14, 2012
// Last Update : Apr 7, 2014
//============================================================================

// Implementation with 2 heaps : the first one with Comparatorlb, the second with another comparator 

#include "ibex_CellHeapOptim.h"
#include "ibex_Optimizer.h"
#include <algorithm>
using namespace std;

namespace ibex {


namespace {
/* we need these dummy classes just because
 * the xxx_heap functions take the comparator
 * argument by copy (hence, we cannot give
 * "*this" since the class is abstract) */



// the classical best first search comparator, based on minimizing the lower bound of the cost estimate of the cell with the upper bound of the cost for breaking the ties.
// this comparator is used in the first heap (buffer of Optimizer   crit==LB)
struct CellComparatorFirstCrit {
	bool operator()(const pair<OptimCell*,pair<double,double> >& c1, const pair<OptimCell*,pair<double,double> >& c2) {
		if( c1.second.first !=  c2.second.first)
			return c1.second.first >= c2.second.first;
		else
			return c1.second.second >= c2.second.second;
	}
};


// the other comparators  used in the second heap  (buffer2  of Optimizer)
struct CellComparatorSecondCrit {
	bool operator()(const pair<OptimCell*,pair<double,double> >& c1, const pair<OptimCell*,pair<double,double> >& c2) {
		if( c1.second.second !=  c2.second.second)
			return c1.second.second >= c2.second.second;
		else
			return c1.second.first >= c2.second.first;
	}
};



CellHeapOptim::CellHeapOptim(const int y, criterion crit) : y(y) , crit(crit){;}


double CellHeapOptim::cost_first(const OptimCell& c) const {
	switch (crit)	{
	case LB  : return  c.box[y].lb();
	case UB  : return  c.box[y].ub();
	case C3  : return -((c.loup - c.pf.lb()) / c.pf.diam() );
	case C5  : return -(c.pu * (c.loup - c.pf.lb()) / c.pf.diam());
	case C7  : return  c.box[y].lb()/(c.pu*(c.loup-c.pf.lb())/c.pf.diam());
	case PU  : return -c.pu	;
	case PF  : return  c.pf.lb();
	}

}

/** The second criterion */
double CellHeapOptim::cost_second(const OptimCell& c) const {
	switch (crit)	{
	case LB  : 	return c.box[y].ub();
	case PF  : 	return c.pf.ub();
	case UB  : case PU  : case C7  : case C5  : case C3  :
				return c.box[y].lb();
	}

}


pair<double,double> CellHeapOptim::cost(const OptimCell& c) const {
	return pair<double,double>( cost_first(c),cost_second(c));
}


/* removes from the top of the heap the cells that have already been removed from the other heap (heap_present <2) */
void CellHeapOptim::cleantop()  {
	while (!empty() && (lopt.front().first)->heap_present < 2)
	{ OptimCell *c = pop();
	if (c->heap_present == 0) delete c;
	}


}

/* "heap destruction" made by another comparator and reconstruction of the heap with its comparator : useful for diversification by breaking the ties another way*/
void CellHeapOptim:: makeheap()
{
	make_heap (lopt.begin(),lopt.end(),CellComparatorSecondCrit());
	make_heap (lopt.begin(),lopt.end(),CellComparatorFirstCrit());
}

void CellHeapOptim::flush() {
	for (vector<pair<OptimCell*,pair<double,double> > >::iterator it=lopt.begin(); it!=lopt.end(); it++)
	{ OptimCell* cell=it->first;
	cell->heap_present--;
	if (cell->heap_present==0)	  delete cell;
	// delete it->second;  // not need now?
	}
	lopt.clear();
}


// E.g.: called in Optimizer in case of a new upper bound
// on the objective ("loup"). This function then removes (and deletes) from
// the heap all the cells with a cost greater than loup.
void CellHeapOptim::contract_heap(double loup)  {

	//	cout << " before contract heap  "  << lopt.size() <<  endl ;
// sort the heap with the criterion LB
	switch (crit)	{
	case LB  : 	case PF  :
		sort_heap(lopt.begin(),lopt.end(),CellComparatorFirstCrit()); break; //TODO JN: is it really necessary? The heap should be already sort.
	case UB  :  case C3  :  case C5  :  case C7  :  case PU  :
		sort_heap(lopt.begin(),lopt.end(),CellComparatorSecondCrit()); break;
	}

// detect the number of OptimCell to remove
	vector<pair<OptimCell*,pair<double,double> > >::iterator it0=lopt.begin();
	int k=0;
	switch (crit)	{
	case LB  : 	case PF  :
		while (it0!=lopt.end() && it0->second.first > loup) { it0++; k++; }  break;
	case UB  :  case C3  :  case C5  :  case C7  :  case PU  :
		while (it0!=lopt.end() && it0->second.second > loup) { it0++; k++; } break;
	}

// remove this OptimCell
	for (int i=0;i<k;i++) {
		//	cout << " i " << " heap_present " << (lopt[i].first)->heap_present << endl;
		(lopt[i].first)->heap_present--;
		if ((lopt[i].first)->heap_present==0) {
			delete lopt[i].first;
		}
		//delete lopt[i].second;  //not need now?
	}
	if (k>0) lopt.erase(lopt.begin(),it0);

// update the value of the criterion if necessary
	switch (crit)	{

	case C3  : case C5  : case C7  : {
		for (unsigned int i=0;i<lopt.size();i++) {
			(lopt[i].first)->loup=loup;
			lopt[i].second.first = cost_first(*(lopt[i].first));
		}
		break;
	}
	case LB  : case UB  :  	case PF  : case PU  : default : break;
	}

// rebuild the heap with the first criterion
	make_heap(lopt.begin(),lopt.end(),CellComparatorFirstCrit());

}


// remove the cell from the buffer and decrements its heap_present counter
OptimCell* CellHeapOptim::pop() {
	OptimCell* c = lopt.front().first;
	c->heap_present--;
	pop_heap(lopt.begin(), lopt.end(), CellComparatorFirstCrit());
	// delete lopt.back().second; // not need now?
	lopt.pop_back(); // removes the "best"
	return c;     // and returns it
}

void CellHeapOptim::push(OptimCell* cell) {
	if (capacity>0 && size()==capacity) throw CellBufferOverflow();

	lopt.push_back(pair<OptimCell*,pair<double,double> >(cell,&cost(*cell)));

	push_heap(lopt.begin(), lopt.end(), CellComparatorFirstCrit());

	cell->heap_present++;

}


// returns the cell on the top of the heap without modifying the heap
OptimCell* CellHeapOptim::top() const {
	return lopt.front().first;
}

double CellHeapOptim::minimum()  {
	return lopt.front().second.first;
}

int CellHeapOptim::size() const {
	return lopt.size();
}
bool CellHeapOptim::empty() const {
	return lopt.empty();
}

ostream& operator<<(ostream& os, const CellHeapOptim& heap) {
	os << "[ ";
	for (vector<pair<OptimCell*,pair<double,double> > >::const_iterator it=heap.lopt.begin(); it!=heap.lopt.end(); it++)
		os << it->first->box << " ";
	return os << "]";
}


} // end namespace ibex
