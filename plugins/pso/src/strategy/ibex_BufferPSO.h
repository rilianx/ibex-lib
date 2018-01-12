/*
 * ibex_BufferPSO.h
 *
 *  Created on: 9 ene. 2018
 *      Author: iaraya
 */

#ifndef PSO_SRC_STRATEGY_IBEX_BUFFERPSO_H_
#define PSO_SRC_STRATEGY_IBEX_BUFFERPSO_H_


#include "ibex_CellBufferOptim.h"
#include "ibex_ExtendedSystem.h"
#include "ibex_CellCostFunc.h"
#include "ibex_Interval.h"
#include "ibex_Random.h"

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
			const Cell* p;

			/** childrens of the node **/
			Cell* left;
			Cell* right;
	};

	class BufferPSO : public CellBufferOptim {
		public:
			BufferPSO();

			~BufferPSO();

			virtual void add_backtrackable(Cell& root) {
				  root.add<CellPSO>();
			}

			/*
			 * Deletes every node from tree.
			 */
			virtual void flush(){
				//es importante, pero para el final
				//TODO Make a better solution
				delete root;
			}

			/*
			 *TODO Count every node in tree?
			 */
			virtual unsigned int size() const{
				//no es importante (retornar 1)
				return 1;
			}

			/*
			 * Is tree empty?
			 */
			virtual bool empty(){
				if(!root) return true;
				else return false;
			}

			/*
			 * Allows to add the root of the tree
			 */
			virtual void push(Cell* cell){
				if(!root) root = cell;
			}

			/*
			 * Deletes last node accessed.
			 * Additionally deletes father knowledge of this child.
			 */
			virtual Cell* pop(){
				if(last_node == NULL) { cout << "Deleting last_node more than once or last_node NULL." << endl;}
				const Cell* father = last_node->get<CellPSO>().p;
				Cell* aux;
				if(father == NULL) {
					root = NULL;
				}else{
					if(father->get<CellPSO>().right == last_node)
						father->get<CellPSO>().right == NULL;
					if(father->get<CellPSO>().left == last_node)
						father->get<CellPSO>().left == NULL;
				}
				aux = last_node;
				last_node = NULL;
				return aux;
			}

			/*
			 * Returns node who contains gBest.
			 */
			virtual Cell* top() {
				return nodeContainer(gb);
			}

			virtual std::ostream& print(std::ostream& os) const;

			virtual double minimum() const {return NEG_INFINITY;}

			//TODO: no es importante!! pero hay que hacer algo al respecto
			virtual void contract(double loup) { ;}

			virtual void set_gbest(Vector& gbest) { gb = gbest; }

			/*
			 * Return the node who contains Vector
			 */
			virtual Cell* nodeContainer(Vector x) {
				Cell* aux = root;
				if(!aux->box.contains(x)) { return NULL; }
				while(aux){
					if(aux->get<CellPSO>().right && aux->get<CellPSO>().right->box.contains(x)){
						aux = aux->get<CellPSO>().right;
						last_node = aux;
					}else if(aux->get<CellPSO>().left && aux->get<CellPSO>().left->box.contains(x)){
						aux = aux->get<CellPSO>().left;
						last_node = aux;
					}else{
						last_node = aux;
						break;
					}
				}
				return aux;
			}


			/*
			 * Return true if Vector is contained by any node
			 */
			static bool isContained(BufferPSO* buffer, Vector x){
				Cell* aux = buffer->nodeContainer(x);
				if(aux != NULL)
					return false;
				else
					return true;
			}

			virtual Cell* nodeSelection(){
				double sel = RNG::rand(0,1);

				Cell* aux = root;
				if(!aux) { return NULL; }
				while(aux){
					if(sel > 0.5)
						if(aux->get<CellPSO>().right)
							aux = aux->get<CellPSO>().right;
						else if (aux->get<CellPSO>().left)
							aux = aux->get<CellPSO>().left;
					else
						if(aux->get<CellPSO>().left)
							aux = aux->get<CellPSO>().left;
						else if(aux->get<CellPSO>().left)
							aux = aux->get<CellPSO>().left;
				}
				return aux;
			}

		protected:
			Cell *root;
			Vector gb;
			Cell* last_node;
	};

} // namespace ibex

#endif
