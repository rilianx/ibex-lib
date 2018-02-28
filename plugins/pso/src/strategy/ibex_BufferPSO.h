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
			BufferPSO(PSOSwarm* swarm) :
				last_node(NULL), tree(swarm->getTree()), swarm(swarm){
				double aux = 0.0;
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

			virtual Cell* pop() { tree->trim(last_node); last_node=NULL; return NULL; }


			virtual Cell* top() const{ cout << "deprecated, you should use Cell*  top(double loup) instead" << endl;}
			/*
			 * Clean tree.
			 * PSO movement
			 * Returns node who contains gBest and update last_node.
			 */
			virtual Cell* top(double loup) const{
				Cell* aux;
				if(!swarm->isInitialized()){
					swarm->resetPSO(loup);
				}

				if(!tree->search(swarm->getGBestPosition())){
				//	cout << "gbest removed!: resetPSO" << endl;
					swarm->resetPSO(loup);
				}

				//cout << "box_gb:" <<  tree->search(swarm->getGBestPosition()) << endl;
				swarm->executePSO(loup);

				//std::cout << "gbest: " << swarm->getGBestPosition() << endl;
				//std::cout << "gbest: " << swarm->getGBestValue() << "+" << swarm->getGBestPenalty() << endl;
				aux = tree->search(swarm->getGBestPosition());


				//el gbest debeería estar en un nodo
				assert(aux);

				//std::cout << "ln (" << last_node << ") ret (" << aux  << ")" << endl;
				last_node = aux;
				return aux;
			}

			virtual std::ostream& print(std::ostream& os) const { return os;}

			virtual double minimum() const {return tree->minimum();}

			virtual pair<double, Vector> get_gbest()  {
				if(swarm->getGBestPenalty()==0){
					return make_pair(swarm->getGBestValue(), swarm->getGBestPosition());
				}else
					return make_pair(POS_INFINITY, swarm->getGBestPosition());
			}

			virtual void update_gbest(Vector loup_point, double loup){
				//cout << "changing the gbest " <<  swarm->getGBestValue() << "-->" << loup << endl;
				swarm->update_gbest(loup_point, loup);
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
