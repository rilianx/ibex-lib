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
	PSOParticle::PSOParticle(System* orig_sys, double c1, double c2) : position(orig_sys->box.random()), pBest(position), velocity(orig_sys->box.random()){
		value = calculateFitness(orig_sys);
		constrictor = 1;
	}

	void PSOParticle::updateVelocityAndPosition(System* orig_sys, PSOParticle* gBest, double c1, double c2, bool bounce, double constrictorMod) {
		double oldFitness = calculateFitness(orig_sys);

		// ** Update Velocity **
		double rand1 = (double)random()/(double)RAND_MAX;
		double rand2 = (double)random()/(double)RAND_MAX;
		Vector newV = constrictor*(velocity + (c1*rand1*(pBest - position)) + (c2*rand2*(gBest->getBestPosition()-position)));
		velocity = newV;

		// ** Update Position based on new Velocity
		Vector newP = position + velocity;
		position = newP;

		// If position exceed box size, stop particle and change velocity's direction.
		// Plus, apply factor of constriction to velocity
		for(int i=0; i<orig_sys->box.size(); i++){
			if(position[i] > orig_sys->box[i].ub()){
				if(bounce == true){
					position[i] = orig_sys->box[i].ub() - (position[i] - orig_sys->box[i].ub());	// bounces with lost of 60%
					velocity[i] = -velocity[i];														// change direction
					if(position[i] < orig_sys->box[i].lb()){	// if bounce exceed other edge
						position[i] = orig_sys->box[i].lb();	// stick to edge
						velocity[i] = -velocity[i];				// change to original direction
					}
				}else{
					position[i] = orig_sys->box[i].ub();	// stick to edge
					velocity[i] = -velocity[i];				// change direction
				}
				constrictor = constrictor*constrictorMod;
			}
			if(position[i] < orig_sys->box[i].lb()){
				if(bounce == true){
					position[i] = orig_sys->box[i].lb() + (orig_sys->box[i].lb() - position[i]);	// bounces with lost of 60%
					velocity[i] = -velocity[i];														// change direction
					if(position[i] > orig_sys->box[i].ub()){	// if bounce exceed other edge
						position[i] = orig_sys->box[i].ub();	// stick to edge
						velocity[i] = -velocity[i];				// change to original direction
					}
				}else{
					position[i] = orig_sys->box[i].lb();	// stick to edge
					velocity[i] = -velocity[i]; 			// change direction
				}
				constrictor = constrictor*constrictorMod;
			}
		}
		if(calculateFitness(orig_sys) < oldFitness)
			pBest = position;
	}

	double PSOParticle::calculateFitness(System* orig_sys){
		int countPenalty = 0;
		double fitness;
		for(int i=0; i<orig_sys->ctrs.size();i++){
			Interval eval = orig_sys->ctrs[i].f.eval(position);
			if((orig_sys->ctrs[i].op == LEQ || orig_sys->ctrs[i].op == LT) && eval.ub()>0.0){
				//cout << "error("<< i <<"):" << eval.ub() << endl;
				countPenalty++;
			}
			else if((orig_sys->ctrs[i].op == GEQ || orig_sys->ctrs[i].op == GT) && eval.lb()<0.0){
				//cout << "error("<< i <<"):" << -eval.lb() << endl;
				countPenalty++;
			}
			else if(orig_sys->ctrs[i].op == EQ ){
				//cout << "error("<< i <<"):" << abs(eval).ub() << endl;
				countPenalty++;
			}
		}
		if(countPenalty > 0){
			fitness = INT_MAX;
		}else
			fitness = orig_sys->goal->eval(position).ub();

		return fitness;
	}

	Vector PSOParticle::getBestPosition(){
		return pBest;
	}

	PSOParticle::~PSOParticle() {
		delete this;
	}

}

