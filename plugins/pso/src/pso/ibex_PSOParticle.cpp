/*
 * ibex_PSOParticle.cpp
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#include "ibex_PSOParticle.h"
#include "ibex_Interval.h"
#include "ibex_Random.h"
#include <climits>
#include "stdlib.h"

using namespace std; //delete
namespace ibex {
	PSOParticle::PSOParticle(TreeCellOpt* tree, System* orig_sys, double c1, double c2) : position(orig_sys->box.random()), pBest(position), velocity(Vector::zeros(orig_sys->box.size())){
		// Randomize particle position into valid place
		if(!tree->is_empty()){
			std::cout << "select node"<< endl;
			Cell* c = tree->random_node();
			std::cout << "created in (" << c << ")" << endl;
			position = c->box.random();
			position.resize(orig_sys->box.size());
		}
		pBest = position;
		calculateFitness(orig_sys); //needs position and pBest
		// store best values
		vBest = value;
		peBest = penalty;
	}

	/*
	 * \param function(Vector pos) external function to check if position is inside a node from a tree
	 * updates current velocity and position of a particle
	 */
	void PSOParticle::updateVelocityAndPosition(TreeCellOpt* tree, System* orig_sys, Vector gBest, double c1, double c2, double p) {
		double maxVel;	//max velocity value

		// ** Update Velocity **
		double rand1 = RNG::rand(0,1);
		double rand2 = RNG::rand(0,1);
		Vector newV = velocity + (c1*rand1*(pBest - position) + (c2*rand2*(gBest-position)));
		velocity = newV;

		// check for max velocity exceeded
		for(int i=0; i<orig_sys->box.size(); i++){
			maxVel = p*orig_sys->box[i].diam(); 		// set maxVel for each variable
			if(velocity[i] > maxVel)
				velocity[i] = maxVel;
			else if(velocity[i] < -maxVel)
				velocity[i] = -maxVel;
		}

		// ** Update Position based on new Velocity **
		position = position + velocity;


		// check if buffer's tree exist (B&B-PSO strategy)
		if(!tree){
			// check for position exceed box size, stop particle and change velocity's direction
			for(int i=0; i<orig_sys->box.size(); i++){
				if(position[i] > orig_sys->box[i].ub()){
					position[i] = orig_sys->box[i].ub();	// stick to edge (ub)
					velocity[i] = -velocity[i];				// change direction
				}
				if(position[i] < orig_sys->box[i].lb()){
					position[i] = orig_sys->box[i].lb();	// stick to edge (lb)
					velocity[i] = -velocity[i]; 			// change direction
				}
			}
		}
	}

	/*
	 * update value (fitness) and violations for current position
	 * fitness = 0 + error; f(x) < f(bestx)
	 * fitness = (f(x)-f(bestx)) + error; f(x) >= f(bestx)
	 */
	void PSOParticle::calculateFitness(System* orig_sys){
		double sum = 0, M = 10;
		for(int i=0; i<orig_sys->ctrs.size();i++){
			Interval eval = orig_sys->ctrs[i].f.eval(position);
			if((orig_sys->ctrs[i].op == LEQ || orig_sys->ctrs[i].op == LT) && eval.ub()>0.0){
				sum += eval.ub();
			}
			else if((orig_sys->ctrs[i].op == GEQ || orig_sys->ctrs[i].op == GT) && eval.lb()<0.0){
				sum += -eval.ub();
			}
			else if(orig_sys->ctrs[i].op == EQ ){
				sum += abs(eval).ub();
			}
		}

		// evaluate
		double fx = orig_sys->goal->eval(position).ub();
		double bx = orig_sys->goal->eval(pBest).ub();

		//calculate fitness
		double fitness;
		if(fx < bx)
			fitness = 0;
		else
			fitness = fx - bx;

		//store results
		value = fitness + sum;
		penalty = sum;
	}

	/*
	 * Select best position based on fitness previous calculation
	 * Considering pBest can't be outside a box
	 */
	void PSOParticle::selectBestInternal(TreeCellOpt* tree){
		//TODO DELETE COUT
		//if(tree){
			if(tree->search(position)){
				//cout << "position inside -> update pBest" << endl;
				if(value < vBest){
					pBest = position;
					vBest = value;
					peBest = penalty;
				} //else: keep
			}
		/*}else{
			if(value < vBest){
				pBest = position;
				vBest = value;
				peBest = penalty;
			} //else: keep
		}*/
	}

	/*
	 * Another strategy to select best particle based on Feasibility rule (K. Deb - 2000)
	 */
	/*void PSOParticle::selectBestInternal(){
		if(peBest == 0 && penalty == 0){		// both feasible: select better fitness.
			if(value < vBest){
				pBest = position;
				vBest = value;
				peBest = penalty;
			}	// else: keep
		}else if(peBest > 0 && penalty == 0){	// only one is feasible: select feasible.
				pBest = position;
				vBest = value;
				peBest = penalty;
		}else if(peBest > 0 && penalty > 0){	// both infeasible: select the one with less penalty.
			if(penalty < peBest){
				pBest = position;
				vBest = value;
				peBest = penalty;
			}	// else: keep
		}
	}*/

	/*
	 * check for violations of best position known, returns boolean
	 */
	bool PSOParticle::isBestFeasible(){
		if(peBest > 0)
			return false;
		else
			return true;
	}

	/*
	 * check for violations of current position, returns boolean
	 */
	bool PSOParticle::isFeasible(){
		if(penalty > 0)
			return false;
		else
			return true;
	}

	/*
	 * GETTER AND SETTER
	 */
	Vector PSOParticle::getPosition(){
		return position;
	}

	Vector PSOParticle::getBestPosition(){
		return pBest;
	}

	double PSOParticle::getBestValue(){
		return vBest;
	}

	double PSOParticle::getValue(){
		return value;
	}

	double PSOParticle:: getPenalty(){
		return penalty;
	}

	double PSOParticle:: getBestPenalty(){
		return peBest;
	}

	PSOParticle::~PSOParticle() {
		delete this;
	}

}

