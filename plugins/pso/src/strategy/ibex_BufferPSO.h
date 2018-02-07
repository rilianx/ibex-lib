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
#include "ibex_TreeCellOpt.h"
#include "ibex_PSOSwarm.h"
#include "ibex_Interval.h"
#include "ibex_Random.h"
#include <map>

using namespace std;

namespace ibex {
	class BufferPSO : public CellBufferOptim {
		public:
			BufferPSO(System* orig_sys, double c1, double c2, int particles, int iterations, double p) :
				last_node(NULL), swarm(new PSOSwarm(tree,orig_sys,c1,c2,particles,iterations,p)), loup_point(NULL),
				loup_value(NULL), tree(new TreeCellOpt(orig_sys)), loup_cell(NULL){
			}
			virtual ~BufferPSO(){}

			/*
			 * Deletes everything from buffer.
			 */
			virtual void flush(){
				tree->cleanTree();
			}

			/*
			 * Returns amount of nodes in buffer.
			 */
			virtual unsigned int size() const{
				//TODO better approach
				if(!empty())
					return 1;
				else
					return 0;
			}

			virtual bool empty() const{
				return tree->empty();
			}

			/*
			 * Pushes a cell into buffer.
			 */
			virtual void push(Cell* cell){
				tree->add(cell);
 			}

			virtual Cell* pop() { return NULL; }

			/*
			 * Clean tree.
			 * PSO movement
			 * Returns node who contains gBest and update last_node.
			 */
			virtual Cell* top() const{
				Cell* aux;
				if(tree->trim()){
					loup_cell = tree->nodeContainer(loup_point.mid());
					if(!tree->nodeContainer(swarm->getGBestPosition()))
						swarm->initializePSO();
				}
				std::cout << "get gb: " << loup_point.mid() << endl;
				if(loup_cell && loup_value < swarm->getGBestValue())
					aux = loup_cell;
				else
					aux = tree->nodeContainer(swarm->getGBestPosition());
				last_node = aux;
				return aux;
			}

			virtual std::ostream& print(std::ostream& os) const { return os;}

			virtual double minimum() const {return NEG_INFINITY;}

			/*
			 * Update cell who contains loup_point, then contract tree.
			 */
			virtual void contract(double loup) {
				loup_cell = tree->nodeContainer(loup_point.mid());
				tree->contract(loup);
			}

			virtual void set_loup_point(IntervalVector &loup_point){
				this->loup_point = loup_point;
				Cell* aux = tree->nodeContainer(swarm->getGBestPosition());
				if(aux)
					loup_cell = aux;
			}

			virtual void set_loup_value(double &loup_value){
				this->loup_value = loup_value;
			}


		protected:
			/* loup info */
			IntervalVector &loup_point;
			double &loup_value;
			mutable Cell* loup_cell;

			/* structures */
			mutable Cell* last_node;
			TreeCellOpt* tree;
			PSOSwarm* swarm;
	};

} // namespace ibex

#endif
