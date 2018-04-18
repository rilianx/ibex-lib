/*
 * ibex_TreeCell.h
 *
 *  Created on: 05-02-2018
 *      Author: cduartes
 */

#ifndef PSO_SRC_STRATEGY_IBEX_TREECELLOPT_H_
#define PSO_SRC_STRATEGY_IBEX_TREECELLOPT_H_

#include "ibex_CellBufferOptim.h"
#include "ibex_Backtrackable.h"
#include "ibex_ExtendedSystem.h"
#include "ibex_CellCostFunc.h"
#include "ibex_Interval.h"
#include "ibex_Random.h"
#include "ibex_System.h"

#include <map>

using namespace std;

namespace ibex {
class CellPSO : public Backtrackable {
	public:
		/**
		 * \brief Constructor for the root node (followed by a call to init_root).
		 */
		CellPSO() : p(NULL), left(NULL), right(NULL) {}

		/**
		 * \brief Duplicate the structure into the left/right nodes
		 */
		std::pair<Backtrackable*,Backtrackable*> down(){
			CellPSO* c1 = new CellPSO();
			CellPSO* c2 = new CellPSO();

			return std::pair<Backtrackable*,Backtrackable*>(c1,c2);
		}

		/** parent of the node **/
		Cell* p;

		/** childrens of the node **/
		Cell* left;
		Cell* right;
};
//TODO create base class TreeCell and extend TreeCellOpt
class TreeCellOpt {
	public:
	TreeCellOpt(System* orig_sys) : orig_sys(orig_sys), root(NULL), min(NEG_INFINITY){}

	/*
	 * Is tree empty?
	 */
	bool is_empty();

	/*
	 * Insert a new node into the tree
	 * If root exists, give it as children of last_node.
	 */
	void insert(Cell* cell, Cell* last_node=NULL);

	/*
	 * Trim the tree
	 * Delete every node that has no children left starting from last_node.
	 * Return true if the minlb_node was removed
	 */
	bool trim(Cell* last_node);

	/*
	 * Random selection of leaf-node.
	 */
	Cell* random_node(Cell* node=NULL);

  /*
  * Return the node with the min lb
		*/
	//Cell* minlb_node();

	Cell* diving_node(Cell* minlb);

	/*
	 * Return true if Vector is contained by any node
	 */
	bool isContained(const Vector& x){
		Cell* aux = search(x);
		std::cout << "found in (" << aux << ")" << endl;
		if(aux)
			return true;
		else
			return false;
	}

	/*
	 * Return the node who contains Vector
	 */
	Cell* search(const Vector& x) const;

	/*
	virtual bool hasGB(Cell* cell){
		bool have = true;
		for(int i=0; i<orig_sys->box.size(); i++){
			if(!cell->box[i].contains(gb[i])){
				have = false;
				break;
			}
		}
		cout << "(" << cell << ") have::" << have << ":" << gb << endl;
		return have;
	}*/

	/*
	 * Deletes every node from tree.
	 */
	void cleanTree(){
		stack<Cell*> S;
		if(root) S.push(root);
		else return;

		min = POS_INFINITY;
		while(!S.empty()){
			Cell* c = S.top(); S.pop();
			if(c->get<CellPSO>().left || c->get<CellPSO>().right){
				if(c->get<CellPSO>().left) S.push(c->get<CellPSO>().left);
				if(c->get<CellPSO>().right) S.push(c->get<CellPSO>().right);
			}else{
				trim(c);
			}
		}
	}

	/*
	 * Contract or removes leave nodes with ub<loup
	 */
	void contract(double loup);

	double minimum(){return min;}

	~TreeCellOpt();

	Cell* root;

protected:
	//Cell* last_node;

	System* orig_sys;

	double min;
};

} /* namespace ibex */

#endif /* PSO_SRC_STRATEGY_IBEX_TREECELLOPT_H_ */