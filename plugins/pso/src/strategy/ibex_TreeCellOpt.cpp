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
		//cout << " [NULL?]" << endl;
		return true;
	}else{
		//std::cout << " [NOTNULL]" << endl;
		return false;
	}
}

void TreeCellOpt::insert(Cell* cell, Cell* last_node){
	//cout << "push(" << cell << ")" << endl;
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
	}else
		cout << "no last_node" << endl;
}

bool TreeCellOpt::trim(Cell* last_node, Cell* minlb){
	bool ret = false;
	if(last_node){
		Cell* father = last_node->get<CellPSO>().p;
		Cell* aux = last_node;
		while(aux){
			if(!aux->get<CellPSO>().left && !aux->get<CellPSO>().right){
				father = aux->get<CellPSO>().p;
				if(father){
					if(father->get<CellPSO>().left == aux) father->get<CellPSO>().left = NULL;
					else if(father->get<CellPSO>().right == aux) father->get<CellPSO>().right = NULL;

					if(aux==minlb) ret=true;
					delete aux;
					aux = father;
					last_node = NULL;
				}else{
					if(root==minlb) ret=true;
					delete root;
					root = NULL;
					last_node = NULL;
					aux=NULL;
				}
			}else
				break;
		}
	}else{
		std::cout << "Deleting last_node more than once / it's first time" << endl;
	}

	return ret;
}

Cell* TreeCellOpt::diving_node(Cell* minlb){
	Cell* aux = minlb;

	if(aux == NULL) {
		//std::cout << "root is NULL" << endl;
		return NULL;
	}
  int n = aux->box.size()-1;

	while(aux){

		if(!aux->get<CellPSO>().left && !aux->get<CellPSO>().right){
			return aux;
		}else if(!aux->get<CellPSO>().left)
			aux=aux->get<CellPSO>().right;
		else if(!aux->get<CellPSO>().right)
		  aux=aux->get<CellPSO>().left;
		else{
			if(aux->get<CellPSO>().left->box[n].lb() < aux->get<CellPSO>().right->box[n].lb())
			   aux=aux->get<CellPSO>().left;
			else
		  	 aux=aux->get<CellPSO>().right;
		}
  }
	return aux;
}

Cell* TreeCellOpt::minlb_node(){
	Cell* aux = root;
	Cell* minlb = NULL;
	if(aux == NULL) {
		//std::cout << "root is NULL" << endl;
		return NULL;
	}
  int n = aux->box.size()-1;

	stack<Cell*> S;
	S.push(aux);

	while(!S.empty()){
		aux=S.top(); S.pop();

		if(!aux->get<CellPSO>().left && !aux->get<CellPSO>().right){
			if(!minlb || aux->box[n].lb() < minlb->box[n].lb())
			   minlb=aux;
		}else{
			if(aux->get<CellPSO>().left && (!minlb || aux->get<CellPSO>().left->box[n].lb() < minlb->box[n].lb()))
				S.push(aux->get<CellPSO>().left);

			if(aux->get<CellPSO>().right && (!minlb || aux->get<CellPSO>().right->box[n].lb() < minlb->box[n].lb()))
				S.push(aux->get<CellPSO>().right);
	  }

  }
	return minlb;
}

Cell* TreeCellOpt::random_node(Cell* node){

	Cell* aux = root;
  if(node) aux = node;
	if(aux == NULL) {
		//std::cout << "root is NULL" << endl;
		return NULL;
	}
	while(aux){
		float sel = RNG::rand(0,1);
		if(sel > 0.5){
			if(aux->get<CellPSO>().right){
				aux = aux->get<CellPSO>().right;
			}else if (aux->get<CellPSO>().left){
				aux = aux->get<CellPSO>().left;
			}else
				break;
		}else{
			if(aux->get<CellPSO>().left){
				aux = aux->get<CellPSO>().left;
			}else if(aux->get<CellPSO>().right){
				aux = aux->get<CellPSO>().right;
			}else
				break;
		}
	}
	return aux;
}

Cell* TreeCellOpt::search(const Vector& x) const{
	//cout << "search" << endl;
	Cell* aux = root;
	bool bool_have;
	if(!aux){
		//std::cout << "root_null" << endl;
		//cout << "NULL1" << endl;
		return NULL;
	}

	stack<Cell*> S;
	S.push(aux);

	while(!S.empty()){
		aux=S.top(); S.pop();

		if(!aux->get<CellPSO>().left && !aux->get<CellPSO>().right){
			// is a leaf
			//std::cout << "(" << aux << ") is leaf" << endl;

			bool_have = true;
			for(int i=0; i<orig_sys->box.size(); i++){
				if(!aux->box[i].contains(x[i])){
					bool_have = false;
					break;
				}
			}
			if(bool_have){
				S.push(aux->get<CellPSO>().left);
				//std::cout << "this leaf has it" << endl;
				return aux;
			}else{
				return NULL;
			}
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
					S.push(aux->get<CellPSO>().left);
					//aux = aux->get<CellPSO>().left;
					//std::cout << "left has it" << endl;
					//continue;
				}
			}

			if(aux->get<CellPSO>().right){
				//std::cout << "(" << aux << ") has right (" << aux->get<CellPSO>().right << ")" << endl;
				bool_have = true;
				for(int i=0; i<orig_sys->box.size(); i++){
					if(!aux->get<CellPSO>().right->box[i].contains(x[i])){
						bool_have = false;
						break;
					}
				}
				if(bool_have){
					S.push(aux->get<CellPSO>().right);
					//aux = aux->get<CellPSO>().right;
					//std::cout << "right has it" << endl;
					//continue;
				}
			}



			//cout << "NULL2" << endl;
			//return NULL;
		}
	}

	//cout << "NULL" << endl;
	return NULL;
}

void TreeCellOpt::contract(double loup) {

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
			if(c->box[c->box.size()-1].lb() > loup) trim(c);
			else{
				c->box[c->box.size()-1] &=Interval(NEG_INFINITY, loup);
				if(c->box[c->box.size()-1].lb() < min) min=c->box[c->box.size()-1].lb();
			}
		}
	}

}

}
