/*
 * ibex_TreeCellOpt.cpp
 *
 *  Created on: 26 feb. 2018
 *      Author: iaraya
 */


#include "ibex_TreeCellOpt.h"

using namespace std;

namespace ibex {

/*
 * Is tree empty?
 */
bool TreeCellOpt::is_empty(){
	//std::cout << "root: " << root;
	if(root == NULL){
		cout << " [NULL?]" << endl;
		return true;
	}else{
		//std::cout << " [NOTNULL]" << endl;
		return false;
	}
}

void TreeCellOpt::insert(Cell* cell, Cell* last_node){
	cout << "push(" << cell << ")" << endl;
	if(!root){
		root = cell;
	}else if (last_node){
		if(!last_node->get<CellPSO>().left){
			std::cout << "left of " << last_node << endl;
			last_node->get<CellPSO>().left = cell;
			cell->get<CellPSO>().p = last_node;
		}else if(!last_node->get<CellPSO>().right){
			std::cout << "right of " << last_node << endl;
			last_node->get<CellPSO>().right = cell;
			cell->get<CellPSO>().p = last_node;
		}else{
			cout << "no-where to assign cell" << endl;
		}
	}
}

bool TreeCellOpt::trim(Cell* last_node){
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

Cell* TreeCellOpt::random_node(){
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
			else if(aux->get<CellPSO>().right)
				aux = aux->get<CellPSO>().right;
			else
				break;
		}
	}
	return aux;
}

Cell* TreeCellOpt::search(const Vector& x) const{
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
			std::cout << "(" << aux << ") is leaf" << endl;
			return aux;
		}else{
			if(aux->get<CellPSO>().left){
				std::cout << "(" << aux << ") has left (" << aux->get<CellPSO>().left << ")" << endl;
				bool_have = true;
				for(int i=0; i<orig_sys->box.size(); i++){
					if(!aux->get<CellPSO>().left->box[i].contains(x[i])){
						bool_have = false;
						break;
					}
				}
				if(bool_have){
					aux = aux->get<CellPSO>().left;
					std::cout << "left has it" << endl;
					bool_gb_found = true;
				}
			}
			if(!bool_gb_found && aux->get<CellPSO>().right){
				std::cout << "(" << aux << ") has right (" << aux->get<CellPSO>().right << ")" << endl;
				bool_have = true;
				for(int i=0; i<orig_sys->box.size(); i++){
					if(!aux->get<CellPSO>().right->box[i].contains(x[i])){
						bool_have = false;
						break;
					}
				}
				if(bool_have){
					aux = aux->get<CellPSO>().right;
					std::cout << "right has it" << endl;
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

void TreeCellOpt::contract(double loup) {
	stack<Cell*> S;
	S.push(root);

	while(!S.empty()){
		Cell* c = S.top(); S.pop();
		if(c->get<CellPSO>().left || c->get<CellPSO>().right){
			if(c->get<CellPSO>().left) S.push(c->get<CellPSO>().left);
			if(c->get<CellPSO>().right) S.push(c->get<CellPSO>().right);
		}else{
			if(c->box[c->box.size()-1].lb() > loup) trim(c);
			else c->box[c->box.size()-1] &=Interval(NEG_INFINITY, loup);
		}
	}
}

}
