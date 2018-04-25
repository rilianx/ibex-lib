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

bool TreeCellOpt::trim(Cell* last_node){
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

					delete aux;
					aux = father;
					last_node = NULL;
				}else{

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

Cell* TreeCellOpt::minlb(Cell* minlb){
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

/*
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
}*/

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

bool TreeCellOpt::contains(const IntervalVector& box, const Vector& x) const{

	for(int i=0; i<x.size(); i++){
		if(!box[i].contains(x[i]))
			return false;
	}
	return true;
}


Cell* TreeCellOpt::search(const Vector& x) const{
	//cout << "search" << endl;
	Cell* aux = root;
	if(!aux) return NULL;

	stack<Cell*> S;
	S.push(aux);

	while(!S.empty()){
		aux=S.top(); S.pop();

		if(!aux->get<CellPSO>().left && !aux->get<CellPSO>().right){

			if( contains(aux->box, x) )
				return aux;


		}else{
			if(aux->get<CellPSO>().left){

				if(contains(aux->get<CellPSO>().left->box, x))
					S.push(aux->get<CellPSO>().left);

			}

			if(aux->get<CellPSO>().right){

				if(contains(aux->get<CellPSO>().right->box, x))
					S.push(aux->get<CellPSO>().right);

			}

		}
	}

	return NULL;
}

double TreeCellOpt::manhattan_distance(const Vector& x, const IntervalVector& box) const{
	double distance = 0.0;
	for(int i=0; i<x.size(); i++){
		if(!box[i].contains(x[i]))
			distance+= (x[i]<box[i].lb())? box[i].lb()-x[i] : x[i] - box[i].ub();

	}
	return distance;

}

Cell* TreeCellOpt::closest_node(const Vector& x) const{
	//cout << "search" << endl;
	Cell* aux = root;
	Cell* closest=NULL;
	if(!aux) return NULL;


	map<double, Cell*> M; // lb_distance --> node

	double min_distance= POS_INFINITY;

	M.insert(make_pair(0, root));

	while(!M.empty()){
		pair<double, Cell*> pair =*M.begin();
		if(pair.first >= min_distance) return closest;

		aux=pair.second;
		M.erase(M.begin());

		if(!aux->get<CellPSO>().left && !aux->get<CellPSO>().right){

			double d=manhattan_distance(x, aux->box);

			if( d==0 )
				return aux;
			else if(d<min_distance){
				closest=aux;
				min_distance=d;
			}

		}else{
			if(aux->get<CellPSO>().left){

				double lb_d=manhattan_distance(x, aux->get<CellPSO>().left->box);
				M.insert(make_pair(lb_d,aux->get<CellPSO>().left));

			}

			if(aux->get<CellPSO>().right){

				double lb_d=manhattan_distance(x, aux->get<CellPSO>().right->box);
				M.insert(make_pair(lb_d,aux->get<CellPSO>().right));

			}

		}
	}

	return closest;
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
