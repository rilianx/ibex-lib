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
				char file_name[] = "output_innode.txt";
				swarm->startPlot(file_name);
				double aux = 0.0;
			}
			virtual ~BufferPSO(){}

			virtual void add_backtrackable(Cell& root) {
				 // root.add<CellPSO>();
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

			virtual Cell* pop() { cout << "BufferPSO::pop() function is not implemented, use trim instead" << endl; return NULL; }

			virtual void trim() {
				if(tree->trim(last_node, swarm->minlb_node )) swarm->minlb_node=NULL;
				last_node=NULL;
			}


			virtual Cell* top() const{ cout << "deprecated, you should use Cell*  top(double loup) instead" << endl; return NULL;}
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

				//gbest should be inside a node.
				if(!tree->search(swarm->getGBestPosition())){ //  ||
			  //!swarm->minlb_node) {
					//cout << "gbest removed!: resetPSO" << endl;
				  //swarm->resetGBest(loup);
					swarm->resetPSO(loup);

				}


				swarm->executePSO(loup);

				//std::cout << "gbest: " << swarm->getGBestPosition() << endl;
				//std::cout << "gbest: " << swarm->getGBestValue() << "+" << swarm->getGBestPenalty() << endl;
				//if(swarm->getGBestPenalty() < 0.1)
				aux = tree->search(swarm->getGBestPosition());
				//else
				  //aux = tree->diving_node(swarm->minlb_node);
				//cout << aux << endl;


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


			mutable Cell* last_node;

		protected:
			/* structures */

			TreeCellOpt* tree;
			PSOSwarm* swarm;
	};

} // namespace ibex

#endif
