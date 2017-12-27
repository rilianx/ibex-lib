/*
 * ibex_PSOParticle.cpp
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#include "ibex_PSOParticle.h"
//#include "ibex_PSOSwarm.h"
#include <iostream> //DELETE THIS
#include "stdlib.h"

using namespace std; //DELETE THIS
namespace ibex {
	PSOParticle::PSOParticle(System* orig_sys, double c1, double c2) : position(orig_sys->box.random()), pBest(position), velocity(Vector::zeros(orig_sys->box.size())){
		value = calculateFitness(orig_sys);
		constrictor = 1;
	}

	void PSOParticle::updateVelocityAndPosition(System* orig_sys, PSOParticle* gBest, double c1, double c2) {
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
				position[i] = orig_sys->box[i].ub() - (position[i] - orig_sys->box[i].ub());
				//position[i] = orig_sys->box[i].ub();
				velocity[i] = -velocity[i];
				constrictor = constrictor*0.85;
			}
			if(position[i] < orig_sys->box[i].lb()){
				position[i] = orig_sys->box[i].lb() + (orig_sys->box[i].lb() - position[i]);
				//position[i] = orig_sys->box[i].lb();
				velocity[i] = -velocity[i];
				constrictor = constrictor*0.85;
			}
		}
		if(calculateFitness(orig_sys) < oldFitness)
			pBest = position;
	}

	double PSOParticle::calculateFitness(System* orig_sys){
		bool constrained = true;
		int countPenalty = 0;
		Interval eval;
		if(constrained == false){
			eval = orig_sys->goal->eval(position);
			return eval.ub();
		}else if(constrained == true){
			for(int i=0; i<orig_sys->ctrs.size();i++){
				eval = orig_sys->ctrs[i].f.eval(position);
				if((orig_sys->ctrs[i].op == LEQ || orig_sys->ctrs[i].op == LT) && eval.ub()>0.0){
					cout << "error("<< i <<"):" << eval.ub() << endl;
					countPenalty++;
				}
				else if((orig_sys->ctrs[i].op == GEQ || orig_sys->ctrs[i].op == GT) && eval.lb()<0.0){
					cout << "error("<< i <<"):" << -eval.lb() << endl;
					countPenalty++;
				}
				else if(orig_sys->ctrs[i].op == EQ ){
					cout << "error("<< i <<"):" << abs(eval).ub() << endl;
					countPenalty++;
				}
			}
		}
		if(countPenalty != 0){
			return RAND_MAX;
		}
	}

	Vector PSOParticle::getBestPosition(){
		return pBest;
	}

	PSOParticle::~PSOParticle() {
		delete this;
	}

}

