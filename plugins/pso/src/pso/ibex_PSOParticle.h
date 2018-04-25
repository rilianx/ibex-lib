/*
 * ibex_PSOParticle.h
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#ifndef PSO_SRC_IBEX_PSOPARTICLE_H_
#define PSO_SRC_IBEX_PSOPARTICLE_H_
#include "ibex_System.h"
#include "ibex_TreeCellOpt.h"

namespace ibex{
	class PSOParticle {
		public:
			PSOParticle(TreeCellOpt* tree, System* orig_sys, double eqeps=1e-8);
			virtual ~PSOParticle();

			void updateVelocityAndPosition(TreeCellOpt* tree, System* orig_sys, Vector gBest, double c1, double c2, double x);

			static double computePenalty(System* orig_sys, Vector& position);
			void selectBestInternal(TreeCellOpt* tree);
			bool isFeasible();
			bool isBestFeasible();
			Vector getPosition();
			Vector getBestPosition();
			double getValue();
			double getBestValue();
			double getPenalty();
			double getBestPenalty();


			void initialize(TreeCellOpt* tree, System* orig_sys){
				// Randomize particle position into valid place
				if(tree->is_empty()) return;

				Cell* c = tree->random_node(NULL);
				//Cell* c = minlb_node;


				position = c->box.random();
				velocity = Vector::zeros(orig_sys->box.size());
				position.resize(orig_sys->box.size());

				value = orig_sys->goal->eval(position).ub();
				//cout << value << endl;
				penalty = computePenalty(orig_sys, position);
				//cout << penalty << endl;

				// store best values
				pBest = position;
				vBest = value;
				peBest = penalty;

				//cout << "pbest: " <<  vBest << "+" << peBest  << endl;

				if(!tree->search(position)){
				   cout << "ERROR: The new particle "<< position <<" is not in the tree" << endl;
				}
			}

			void update_pBest(TreeCellOpt* tree, System* orig_sys, double loup){
				//cout << position << endl;
				if(tree->search(position)){
					value = orig_sys->goal->eval(position).ub();
					penalty = computePenalty(orig_sys, position);


					if(compute_fitness(value, penalty, loup) < compute_fitness(vBest, peBest, loup)){
						//cout << "is inside -> update pBest" << endl;
						pBest = position;
						vBest = value;
						peBest = penalty;
						//cout << "new pbest: " <<  vBest << "+" << peBest << " --> " << compute_fitness(vBest, peBest, loup) << endl;
					}
				}
				//cout << "is NOT inside" << endl;

			}

			static double compute_fitness(double val, double viol, double loup){
				if(loup==POS_INFINITY) return viol;

				if(viol>0.0){
					return viol + std::max(0.0,val-loup);
				}else
					return val - loup;
			}


		protected:
			double value;
			double penalty;
			Vector position;

			double vBest;
			double peBest;
			Vector pBest;

			double eqeps;

			Vector velocity;

	};
}
#endif /* PSO_SRC_IBEX_PSOPARTICLE_H_ */
