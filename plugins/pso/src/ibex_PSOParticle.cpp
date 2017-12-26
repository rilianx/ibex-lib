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
	}

	void PSOParticle::updateVelocityAndPosition(System* orig_sys, PSOParticle* gBest, double c1, double c2) {
		bool trace = true;
		if(trace) cout << "Update Velocity and Position" << endl;
		double oldFitness = calculateFitness(orig_sys);

		// ** Update Velocity **
		double rand1 = (double)random()/(double)RAND_MAX;
		double rand2 = (double)random()/(double)RAND_MAX;

		if(trace) cout << "bestPosition: " << pBest << endl;
		if(trace) cout << "globalBest: " << gBest->getPosition() << endl;
		Vector newV = velocity + (c1*rand1*(pBest - position)) + (c2*rand2*(gBest->getPosition()-position));
		//Vector newV = velocity;
		if(trace) cout << "New velocity: " << newV << endl;
		velocity = newV;

		// ** Update Position based on new Velocity
		Vector newP = position + velocity;
		//Vector newP = position;
		position = newP;

		for(int i=0; i<orig_sys->box.size(); i++){
			if(position[i] > orig_sys->box[i].ub())
				position[i] = orig_sys->box[i].ub();
			if(position[i] < orig_sys->box[i].lb())
				position[i] = orig_sys->box[i].lb();
		}

		if(calculateFitness(orig_sys) > oldFitness)
			pBest = position;
	}

	double PSOParticle::calculateFitness(System* orig_sys){
		Interval eval = orig_sys->goal->eval(position);
		double fitness = eval.ub();
		return fitness;
	}

	Vector PSOParticle::getPosition(){
		return position;
	}

	PSOParticle::~PSOParticle() {
		delete this;
	}

}

