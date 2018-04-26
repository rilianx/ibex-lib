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
#include "ibex_Interval.h"
#include "ibex_Random.h"
#include "ibex_System.h"
#include <map>

using namespace std;

namespace ibex {

	class Part{

		public:

		Part(const Vector& pos) : position(pos){
			position.resize(orig_sys->box.size());
			cost=orig_sys->goal->eval(position).ub();
			viols=computePenalty();
		}

		void reinitialize(const Vector& pos) {
		  position=pos;
		  position.resize(orig_sys->box.size());
			cost=orig_sys->goal->eval(position).ub();
			viols=computePenalty();
		}


		double computePenalty(){
			double sum = 0;
			for(int i=0; i<orig_sys->ctrs.size();i++){
				Interval eval = orig_sys->ctrs[i].f.eval(position);
				if((orig_sys->ctrs[i].op == LEQ || orig_sys->ctrs[i].op == LT) && eval.ub()>0.0){
					sum += eval.ub();
				}
				else if((orig_sys->ctrs[i].op == GEQ || orig_sys->ctrs[i].op == GT) && eval.lb()<0.0){
					sum += -eval.ub();
				}
				else if(orig_sys->ctrs[i].op == EQ ){
					sum += std::max(0.0, abs(eval).ub() /*- eqeps*/);
				}
			}
			return sum;
		}

		double fitness(double current_ub){
			if(current_ub==POS_INFINITY) return -viols;
			else if(viols>0.0) return -(viols + std::max(0.0, cost - current_ub));
			else	return -(cost - current_ub);
		}

		static System* orig_sys;
		Vector position;
		double cost;
		double viols;

	};



	class BufferPSO : public CellBufferOptim {
		public:
			BufferPSO(TreeCellOpt* tree, int popsize=10) :
				last_node(NULL), tree(tree), particles(NULL), popsize(popsize){

					particles = new Part*[popsize];
					for(int i=0;i<popsize; i++) particles[i]=NULL;
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

      /*
			* Removes every node with no children in the last_node branch
			*/
			virtual void trim();


			virtual Cell* top() const{ cout << "deprecated, you should use Cell*  top(double loup) instead" << endl; return NULL;}

			virtual Cell* top(double loup) const;

			virtual std::ostream& print(std::ostream& os) const { return os;}

			virtual double minimum() const {return tree->minimum();}

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

			Part** particles;

			int popsize;
	};

} // namespace ibex

#endif
