//============================================================================
//                                  I B E X                                   
// File        : ibex_CellHeap_2.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 12, 2012
// Last Update : May 12, 2012
//============================================================================

#include "ibex_CellHeap_2.h"
#include <algorithm>

using namespace std;

namespace ibex {

CellHeap_2::~CellHeap_2() {
	delete root;
}

void CellHeap_2::flush() {
	delete root;
}

int CellHeap_2::size() const {
	return nb_cells;
}

// E.g.: called by manage_buffer in Optimizer in case of a new upper bound
// on the objective ("loup"). This function then removes (and deletes) from
// the heap all the cells with a cost greater than loup.
void CellHeap_2::contract_heap(double loup)
{
	//TODO
	//cout << " before contract heap  " << l.size() << endl;

	sort_heap(l.begin(),l.end(),CellComparator());
	vector<pair<Cell*,double> >::iterator it0=l.begin();

	int k=0;
	while (it0!=l.end() && it0->second > loup) { it0++; k++; }

	for (int i=0;i<k;i++) {
		delete l[i].first;
	}

	if (k>0) l.erase(l.begin(),it0);

	make_heap(l.begin(), l.end() ,CellComparator());

	//cout << " after contract heap " << l.size() << endl;

}


bool CellHeap_2::empty() const {
	return (nb_cells==0);
}

void CellHeap_2::push(Cell* cell) {
	if (capacity>0 && size()==capacity) throw CellBufferOverflow();
	//TODO
}

Cell* CellHeap_2::pop() {

	//TODO
	Cell* cell;
	double crit_cell;
	Cell* c_return = root->elt;

	if (nb_cells==1) {
		root->elt=NULL;
		delete root;
		nb_cells--;
		return c_return;

	} else {
		// RECUPERATION DU DERNIER ELEMENT et suppression du noeud associe
		// calcul de la hauteur
		int hauteur = 0;
		int aux = nb_cells;
		while(aux>1) {
			aux = aux / 2;
			hauteur++;
		}

		// pt = pointeur sur l element courant
		HeapNode * pt = *root;
		HeapNode * pt2;
		for (int pos=hauteur-1; 0; pos--) {
			if ( nb_cells & (1 << pos)) {  // test the "pos"th bit of the integer nb_cells
				//	   print *,  "DROITE"
				// on doit arriver sur une feuille ssi la pile est vide
				if (pos!=0) {
					if (pt->right) {
						ibex_error("MOD_GEST_MINIMIER: erreur extraction D 1");
					} else {
						pt = pt->right;
					}
				} else {
					if (pt->right) {
						ibex_error("MOD_GEST_MINIMIER: erreur extraction D 2");
					} else {
						pt2= pt;
						pt = pt->right;
						pt2->right=NULL;

						cell = pt->elt;
						crit_cell= pt->crit;

						pt->father->right=NULL;
						pt->right=NULL;
						pt->left=NULL;
						pt->father=NULL;
						pt->elt=NULL;
						delete pt;
						pt=NULL;
					}


				}
			} else {
				//  print *, "GAUCHE"
				// on doit arriver sur une feuille ssi la pile est vide
				if (pos!=0) {
					if (pt->left) {
						ibex_error("MOD_GEST_MINIMIER: erreur extraction G 1");
					} else {
						pt = pt->left;
					}

				} else {
					if (pt->left) {
						ibex_error("MOD_GEST_MINIMIER: erreur extraction G 2");
					} else {
						pt2= pt;
						pt = pt->left;
						pt2->left = NULL;

						cell = pt->elt;
						crit_cell= pt->crit;

						pt->father->left=NULL;
						pt->right=NULL;
						pt->left=NULL;
						pt->father=NULL;
						pt->elt=NULL;
						delete pt;
						pt=NULL;
					}
				}

			}


		}
		// CEL pointe vers l element qu on remet en tete
		// PT2 pointe vers le noeud supprime
		root->elt = cell;
		root->crit = crit_cell;

		// PRECOLATION
		bool b=true;
		pt = root;

		while (b&&(pt->left)) {
			if (pt->right) {
				if (pt->is_sup(pt->left)) {
					if (pt->right->is_sup(pt->left)) {
						// le gauche est le plus petit
						// on echange l'element de pt avec son fils gauche
						cell = pt->left->elt;
						crit_cell = pt->left->crit;

						pt->left->elt = pt->elt;
						pt->left->crit= pt->crit;

						pt->elt = cell;
						pt->crit = crit_cell;

						// au suivant
						pt = pt->left;
					} else {
						// le droit est le plus petit
						// on echange l'element de pt avec son fils droit
						cell = pt->right->elt;
						crit_cell = pt->right->crit;

						pt->right->elt = pt->elt;
						pt->right->crit= pt->crit;

						pt->elt = cell;
						pt->crit = crit_cell;

						// au suivant
						pt = pt->right;
					}

				} else {
					if (pt->is_sup(pt->right)) {
						// le droit est le plus petit
						// on echange l'element de pt avec son fils droit
						cell = pt->right->elt;
						crit_cell = pt->right->crit;

						pt->right->elt = pt->elt;
						pt->right->crit= pt->crit;

						pt->elt = cell;
						pt->crit = crit_cell;

						// au suivant
						pt = pt->right;
					} else {
						// le noeud courant est le plus petit des 3 donc on arrete
						b=false;

					}
				}
			} else { // il n'y a plus de fils droit mais il y a un fils gauche
				if (pt->is_sup(pt->left)) {
					// le gauche est le plus petit
					// on echange l'element de pt avec son fils gauche
					cell = pt->left->elt;
					crit_cell = pt->left->crit;

					pt->left->elt = pt->elt;
					pt->left->crit= pt->crit;

					pt->elt = cell;
					pt->crit = crit_cell;

					// au suivant
					pt = pt->left;
				} else {
					// le noeud courant est le plus petit des 3 donc on arrete , on doit meme avoir touche le fond ;-)
					b=false;
				}
			}
		}

		nb_cells--;
		return c_return ;
	}

}

Cell* CellHeap_2::top() const {
	return root->elt;
}


ostream& operator<<(ostream& os, const CellHeap_2& heap) {
	os << "[ ";
	os << heap.root << " ";
	return os << "]";
}




HeapNode::HeapNode(): elt(NULL), crit(NEG_INFINITY), right(NULL), left(NULL), father(NULL) { }

HeapNode::HeapNode(Cell* elt, double crit): elt(elt), crit(crit), father(NULL), right(NULL), left(NULL) { }

HeapNode::HeapNode(Cell* elt, double crit, HeapNode * father): elt(elt), crit(crit), father(father), right(NULL), left(NULL) { }

HeapNode::~HeapNode() {
	delete elt;
	if (right) delete right;
	if (left) delete left;
}

bool HeapNode::is_sup(HeapNode *n2) const {
	return crit >= n2->crit;
}

std::ostream& operator<<(std::ostream& os, const HeapNode& node) {
	os << *(node.elt) << " ";
	os <<  *(node.right) << *(node.left);
	return os;
}


}  // end namespace ibex
