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
#include "ibex_System.h"
#include <map>

using namespace std;

namespace ibex {
	class BufferPSO : public CellBufferOptim {
		public:
			BufferPSO(System* orig_sys, double c1, double c2, int particles, int iterations, double p) :
				last_node(NULL), loup_point(new IntervalVector(0)), loup_cell(NULL), tree(new TreeCellOpt(orig_sys)){
				double aux = 0.0;
				loup_value = &aux;
				swarm = new PSOSwarm(tree,orig_sys,c1,c2,particles,iterations,p);
			}
			virtual ~BufferPSO(){}

			virtual void add_backtrackable(Cell& root) {
				  root.add<CellPSO>();
			}

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
				tree->add(cell,last_node);
 			}

			virtual Cell* pop() { return NULL; }

			/*
			 * Clean tree.
			 * PSO movement
			 * Returns node who contains gBest and update last_node.
			 */
			virtual Cell* top() const{
				Cell* aux;
				if(!swarm->isInitialized()){
					swarm->resetPSO();
				}else if(tree->trim(last_node)){
					if(!tree->nodeContainer(swarm->getGBestPosition())){
						swarm->resetPSO();
					}
				}
				//swarm->executePSO();
				loup_cell = tree->nodeContainer(loup_point->mid());
				std::cout << "get gb: " << loup_point->mid() << endl;
				if(loup_cell && *loup_value < swarm->getGBestValue()){
					std::cout << "loup_cell" << endl;
					aux = loup_cell;
				}else{
					std::cout << "get gbest: " << swarm->getGBestPosition() << endl;
					aux = tree->nodeContainer(swarm->getGBestPosition());
				}
				std::cout << "ln (" << last_node << ") ret (" << aux  << ")" << endl;
				last_node = aux;
				return aux;
			}

			virtual std::ostream& print(std::ostream& os) const { return os;}

			virtual double minimum() const {return NEG_INFINITY;}

			/*
			 * Update cell who contains loup_point, then contract tree.
			 */
			virtual void contract(double loup) {
				loup_cell = tree->nodeContainer(loup_point->mid());
				tree->contract(loup);
			}

			virtual void set_loup_point(IntervalVector* loup_point){
				this->loup_point = loup_point;
				Cell* aux = tree->nodeContainer(swarm->getGBestPosition());
				if(aux)
					loup_cell = aux;
			}

			virtual void set_loup_value(double* loup_value){
				this->loup_value = loup_value;
			}

			void set_values(IntervalVector* lp_point, double* lp_value){
				loup_point = lp_point;
				loup_value = lp_value;
			}


		protected:
			/* loup info */
			IntervalVector* loup_point;
			double* loup_value;
			mutable Cell* loup_cell;

			/* structures */
			mutable Cell* last_node;
			TreeCellOpt* tree;
			PSOSwarm* swarm;
	};

} // namespace ibex

#endif
