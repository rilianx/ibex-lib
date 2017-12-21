/*
 * ibex_PSOParticle.cpp
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#include "ibex_PSOParticle.h"
#include "stdlib.h"

namespace ibex {

	PSOParticle::PSOParticle(IntervalVector box, System* orig_sys, double c1, double c2){
		this->velocity = 0;
		this->position = box.random();
		this->value = this->calculateFitness(orig_sys);
		this->pBest = this->position;
		this->vBest = this->value;
		this->c1 = c1;
		this->c2 = c2;
	}

	void PSOParticle::updateVelocity(PSOParticle gBest) {
		double rand1 = (double)random()/(double)RAND_MAX;
		double rand2 = (double)random()/(double)RAND_MAX;
		double newV = velocity + c1*rand1*(pBest - position) + c2*rand2*(gBest - position);
		velocity = newV;
	}
	void PSOParticle::updatePosition(){
		Vector newP = position + velocity;
		position = newP;
	}
	void PSOParticle::setBestPosition(Vector position, double value){
		vBest = value;
		pBest = position;
	}
	double PSOParticle::calculateFitness(System* orig_sys){
		double newValue = orig_sys->goal->eval(position);
		return newValue;
	}

	PSOParticle::~PSOParticle() {
		// TODO Auto-generated destructor stub
	}

}

