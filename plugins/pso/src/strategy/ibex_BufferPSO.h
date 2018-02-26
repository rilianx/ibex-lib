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
				last_node(NULL), tree(new TreeCellOpt(orig_sys)){
				double aux = 0.0;
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
				return tree->is_empty();
			}

			/*
			 * Pushes a cell into buffer.
			 */
			virtual void push(Cell* cell){
				tree->insert(cell,last_node);
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
					if(!tree->search(swarm->getGBestPosition())){
						swarm->resetPSO();
					}
				}

				//swarm->executePSO();


				std::cout << "get gbest: " << swarm->getGBestPosition() << endl;
				aux = tree->search(swarm->getGBestPosition());


				//el gbest debeería estar en un nodo
				assert(aux);

				std::cout << "ln (" << last_node << ") ret (" << aux  << ")" << endl;
				last_node = aux;
				return aux;
			}

			virtual std::ostream& print(std::ostream& os) const { return os;}

			virtual double minimum() const {return NEG_INFINITY;}

			virtual pair<double, Vector> get_gbest()  {
				if(swarm->getGBestCountViolations()==0){
					return make_pair(swarm->getGBestValue(), swarm->getGBestPosition());
				}else
					return make_pair(POS_INFINITY, Vector(1));
			}

			virtual void update_gbest(Vector loup_point, double loup){
				if(loup<swarm->getGBestValue()){
					if (tree->search(loup_point))
						swarm->setGBest(loup_point);
				}
			}

			/*
			 * Contract the tree.
			 */
			virtual void contract(double loup) {
				tree->contract(loup);
			}




		protected:
			/* structures */
			mutable Cell* last_node;
			TreeCellOpt* tree;
			PSOSwarm* swarm;
	};

} // namespace ibex

#endif
