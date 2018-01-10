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
#include <map>

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

			//es importante, pero para el final
			virtual void flush();

			//no es importante (retornar 1)
			virtual unsigned int size() const{
				return 1;
			}

			//root = NULL
			virtual bool empty() const;

			virtual void push(Cell* cell){
				//nada mas que hacer
				if(!root) root = cell;
			}

			virtual Cell* pop(){
				//Eliminar last_nodo y actualizar padre (recursiva)
				return last_node;
			}

			virtual Cell* top() {
				//retornar nodo que contiene gb
				//Cell* aux = root->get<CellPSO>().left;
				Cell* aux = root;
				if(!aux->box.contains(gb)) { return NULL; }
				while(aux){
					if(aux->get<CellPSO>().right && aux->get<CellPSO>().right->box.contains(gb)){
						aux = aux->get<CellPSO>().right;
						last_node = aux;
					}else if(aux->get<CellPSO>().left && aux->get<CellPSO>().left->box.contains(gb)){
						aux = aux->get<CellPSO>().left;
						last_node = aux;
					}else{
						break;
					}
				}

				return aux;
			}

			virtual std::ostream& print(std::ostream& os) const;

			virtual double minimum() const {return NEG_INFINITY;}

			//TODO: no es importante!! pero hay que hacer algo al respecto
			virtual void contract(double loup) { ;}

			void set_gbest(Vector& gbest) { gb=gbest; }

		protected:
			Cell *root;
			Vector gb;

			Cell* last_node;
	};

} // namespace ibex

#endif
