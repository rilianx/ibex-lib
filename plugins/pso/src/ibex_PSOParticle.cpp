/*
 * ibex_PSOParticle.cpp
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#include "ibex_PSOParticle.h"
#include <climits>
#include "stdlib.h"

using namespace std; //DELETE THIS
namespace ibex {
	PSOParticle::PSOParticle(System* orig_sys, double c1, double c2) : position(orig_sys->box.random()), pBest(position), velocity(Vector::zeros(orig_sys->box.size())){
		calculateFitness(orig_sys);
		// store best position
		pBest = position;
		vBest = value;
		peBest = penalty;
	}

	/*
	 * updates current velocity and position of a particle
	 */
	void PSOParticle::updateVelocityAndPosition(System* orig_sys, PSOParticle* gBest, double c1, double c2, double p) {
		double maxVel;	//max velocity value

		// ** Update Velocity **
		double rand1 = (double)random()/(double)RAND_MAX;
		double rand2 = (double)random()/(double)RAND_MAX;
		Vector newV = velocity + (c1*rand1*(pBest - position) + (c2*rand2*(gBest->getBestPosition()-position)));
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

	/*
	 * update value (fitness) and violations for current position
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

		// store results
		value = orig_sys->goal->eval(position).ub();
		penalty = sum;
	}

	/*
	 * Select best position based on Feasibility rule (K. Deb - 2000)
	 */
	void PSOParticle::selectBestInternal(){
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
	}

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

