/*
 * ibex_PSOParticle.cpp
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#include "ibex_PSOParticle.h"
//#include "stdlib.h"

namespace ibex {

	PSOParticle::PSOParticle(System* orig_sys, double c1, double c2){
		velocity = 0;
		position = orig_sys->box.random();
		value = calculateFitness(orig_sys);
		pBest = position;
	}

	void PSOParticle::updateVelocityAndPosition(System* orig_sys, PSOParticle gBest, double c1, double c2) {
		double oldFitness = calculateFitness(orig_sys);

		// ** Update Velocity **
		double rand1 = (double)random()/(double)RAND_MAX;
		double rand2 = (double)random()/(double)RAND_MAX;
		velocity = velocity + c1*rand1*(pBest - position) + c2*rand2*(gBest - position);

		// ** Update Position based on new Velocity
		position = position + velocity;
		if(calculateFitness(orig_sys) > oldFitness)
			pBest = position;
	}

	double PSOParticle::calculateFitness(System* orig_sys){
		double fitness = orig_sys->goal->eval(position);
		return fitness;
	}

	Vector PSOParticle::getPosition(){
		return this->position;
	}

	PSOParticle::~PSOParticle() {
		delete this;
	}

}

