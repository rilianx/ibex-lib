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
	TreeCellOpt(System* orig_sys) : orig_sys(orig_sys), root(NULL){}

	void add_backtrackable(Cell& root) {
		  root.add<CellPSO>();
	}
	/*
	 * Is tree empty?
	 */
	bool empty(){
		//std::cout << "root: " << root;
		if(root == NULL){
			cout << " [NULL?]" << endl;
			return true;
		}else{
			//std::cout << " [NOTNULL]" << endl;
			return false;
		}
	}
	/*
	 * Introduces root into tree;
	 * If root exist, give it as children to the last cell top-ed.
	 */
	void add(Cell* cell, Cell* last_node){
		cout << "push(" << cell << ")" << endl;
		if(!root){
			root = cell;
		}else if (last_node){
			if(!last_node->get<CellPSO>().left){
				//std::cout << "left of " << last_node << endl;
				last_node->get<CellPSO>().left = cell;
				cell->get<CellPSO>().p = last_node;
			}else if(!last_node->get<CellPSO>().right){
				//std::cout << "right of " << last_node << endl;
				last_node->get<CellPSO>().right = cell;
				cell->get<CellPSO>().p = last_node;
			}else{
				cout << "no-where to assign cell" << endl;
			}
		}
	}

	/*
	 * Trim this tree
	 * Delete every node that has no children left starting from last_node.
	 */
	bool trim(Cell* last_node){
		int count = 0;
		if(last_node){
			Cell* father = last_node->get<CellPSO>().p;
			Cell* aux = last_node;
			while(aux){
				if(!aux->get<CellPSO>().left && !aux->get<CellPSO>().right){
					father = aux->get<CellPSO>().p;
					if(father){
						//std::cout << "(" << aux << ") has father (" << father  << ")" << endl;
						//delete last_node from his father
						if(father->get<CellPSO>().left == aux){
							father->get<CellPSO>().left = NULL;
						}else if(father->get<CellPSO>().right == aux){
							father->get<CellPSO>().right = NULL;
						}
						std::cout << "trim(" << aux << ") : " << aux->box << endl;
						count++;
						delete aux;
						aux = father;
						last_node = NULL;
					}else{
						std::cout << "trim(" << aux << ") : " << aux->box << endl;
						count++;
						delete root;
						root = NULL;
						last_node = NULL;
					}
				}else{
					std::cout << "(" << aux << ") has childrens";
					if(aux->get<CellPSO>().left)
						std::cout << " (" << aux->get<CellPSO>().left << ")";
					if(aux->get<CellPSO>().right)
						std::cout << " (" << aux->get<CellPSO>().right << ")";
					std::cout << endl;
					break;
				}
			}
			std::cout << "trimed(" << count << ")" << endl;
		}else{
			std::cout << "Deleting last_node more than once / it's first time" << endl;
		}
		if(count > 0)
			return true;
		else
			return false;
	}

	/*
	 * Random selection of leaf-node.
	 */
	Cell* nodeSelection(){
		Cell* aux = root;
		if(aux == NULL) { std::cout << "root is NULL" << endl; return NULL; }
		while(aux){
			float sel = RNG::rand(0,1);
			if(sel > 0.5){
				if(aux->get<CellPSO>().right)
					aux = aux->get<CellPSO>().right;
				else if (aux->get<CellPSO>().left)
					aux = aux->get<CellPSO>().left;
				else
					break;
			}else{
				if(aux->get<CellPSO>().left)
					aux = aux->get<CellPSO>().left;
				else if(aux->get<CellPSO>().left)
					aux = aux->get<CellPSO>().left;
				else
					break;
			}
		}
		return aux;
	}

	/*
	 * Return true if Vector is contained by any node
	 */
	bool isContained(const Vector& x){
		Cell* aux = nodeContainer(x);
		std::cout << "found in (" << aux << ")" << endl;
		if(aux)
			return true;
		else
			return false;
	}

	/*
	 * Return the node who contains Vector
	 */
	Cell* nodeContainer(const Vector& x) const{
		Cell* aux = root;
		bool bool_gb_found;
		bool bool_have;
		if(!aux){
			std::cout << "root_null" << endl;
			return NULL;
		}
		while(aux){
			bool_gb_found = false;
			if(!aux->get<CellPSO>().left && !aux->get<CellPSO>().right){
				// is a leaf
				//std::cout << "(" << aux << ") is leaf" << endl;
				return aux;
			}else{
				if(aux->get<CellPSO>().left){
					//std::cout << "(" << aux << ") has left (" << aux->get<CellPSO>().left << ")" << endl;
					bool_have = true;
					for(int i=0; i<orig_sys->box.size(); i++){
						if(!aux->get<CellPSO>().left->box[i].contains(x[i])){
							bool_have = false;
							break;
						}
					}
					if(bool_have){
						aux = aux->get<CellPSO>().left;
						//std::cout << "left has it" << endl;
						bool_gb_found = true;
					}
				}
				if(!bool_gb_found && aux->get<CellPSO>().right){
					//std::cout << "(" << aux << ") has right (" << aux->get<CellPSO>().right << ")" << endl;
					bool_have = true;
					for(int i=0; i<orig_sys->box.size(); i++){
						if(!aux->get<CellPSO>().right->box[i].contains(x[i])){
							bool_have = false;
							break;
						}
					}
					if(bool_have){
						aux = aux->get<CellPSO>().right;
						//std::cout << "right has it" << endl;
						bool_gb_found = true;
					}
				}
				if(!bool_gb_found){
					aux = NULL;
					std::cout << "no one has it" << endl;
					break;
				}
			}
		}
		return aux;
	}

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
		//TODO better approach
		//root = NULL;
	}

	void contract(double loup) {
		//TODO contraction method
	}

	~TreeCellOpt();

protected:
	//Cell* last_node;
	Cell* root;
	System* orig_sys;
};

} /* namespace ibex */

#endif /* PSO_SRC_STRATEGY_IBEX_TREECELLOPT_H_ */
