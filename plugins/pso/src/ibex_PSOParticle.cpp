/*
 * ibex_PSOParticle.cpp
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#include "ibex_PSOParticle.h"
//#include "ibex_PSOSwarm.h"
#include <iostream> //DELETE THIS
#include <climits>
#include "stdlib.h"

using namespace std; //DELETE THIS
namespace ibex {
	PSOParticle::PSOParticle(System* orig_sys, double c1, double c2) : position(orig_sys->box.random()), pBest(position), velocity(Vector::zeros(orig_sys->box.size())){
		calculateFitness(orig_sys);
		vBest = value;
	}

	void PSOParticle::updateVelocityAndPosition(System* orig_sys, PSOParticle* gBest, double c1, double c2, double p) {
		double maxVel;	//max velocity value
		// ** Update Velocity **
		double rand1 = (double)random()/(double)RAND_MAX;
		double rand2 = (double)random()/(double)RAND_MAX;
		Vector newV = velocity + (c1*rand1*(pBest - position) + (c2*rand2*(gBest->getBestPosition()-position)));
		velocity = newV;

		// check for max velocity exceeded
		for(int i=0; i<orig_sys->box.size(); i++){
			maxVel = p*orig_sys->box[i].diam();
			if(velocity[i] > maxVel)
				velocity[i] = maxVel;
			else if(velocity[i] < -maxVel)
				velocity[i] = -maxVel;
		}

		// ** Update Position based on new Velocity **
		Vector newP = position + velocity;
		position = newP;

		// If position exceed box size, stop particle and change velocity's direction
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

		calculateFitness(orig_sys);	//update value (fitness)

		// Feasibility rule (K. Deb - 2000)
		// both feasible: select better fitness.
		if(vioBest == 0 && violations == 0){
			if(value < vBest){
				pBest = position;
				vBest = value;
				vioBest = violations;
			}
		}
		// only one is feasible: select feasible.
		if(vioBest > 0 && violations == 0){
			pBest = position;
			vBest = value;
			vioBest = violations;
		}
		// both infeasible: select the one with less violations.
		if(vioBest > 0 && violations > 0){
			if(violations < vioBest){
				pBest = position;
				vBest = value;
				vioBest = violations;
			}
		}
	}

	void PSOParticle::calculateFitness(System* orig_sys){
		int sum = 0, M = 100;
		int countViolations = 0;
		double fitness;
		for(int i=0; i<orig_sys->ctrs.size();i++){
			Interval eval = orig_sys->ctrs[i].f.eval(position);
			if((orig_sys->ctrs[i].op == LEQ || orig_sys->ctrs[i].op == LT) && eval.ub()>0.0){
				sum += eval.ub();
				countViolations++;
				//cout << "error("<< i <<"):" << eval.ub() << endl;
			}
			else if((orig_sys->ctrs[i].op == GEQ || orig_sys->ctrs[i].op == GT) && eval.lb()<0.0){
				sum += -eval.ub();
				countViolations++;
				//cout << "error("<< i <<"):" << -eval.lb() << endl;
			}
			else if(orig_sys->ctrs[i].op == EQ ){
				sum += abs(eval).ub();
				countViolations++;
				//cout << "error("<< i <<"):" << abs(eval).ub() << endl;
			}
		}
		if(sum > 0){
			fitness = sum*M;
		}else{
			fitness = orig_sys->goal->eval(position).ub();
		}

		violations = countViolations;
		value = fitness;
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

	bool PSOParticle::isBestFeasible(){
		if(vioBest > 0)
			return false;
		else
			return true;
	}

	int PSOParticle:: getViolations(){
		return violations;
	}

	int PSOParticle:: getBestViolations(){
		return vioBest;
	}

	bool PSOParticle::isFeasible(){
		if(violations > 0)
			return false;
		else
			return true;
	}

	PSOParticle::~PSOParticle() {
		delete this;
	}

}

