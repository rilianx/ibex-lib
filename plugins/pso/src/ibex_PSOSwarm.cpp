/*
 * ibex_PSOSwarm.cpp
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#include "ibex_PSOSwarm.h"

namespace ibex{
	PSOSwarm::PSOSwarm(IntervalVector box, System* orig_sys, double c1, double c2, int nParticles, int limit){
		this->c1 = c1;
		this->c2 = c2;
		this->nParticles = nParticles;
		this->particles[nParticles];
		this->limit = limit;

		double currentBest;

		//initializate particles
		for(int i=0; i<nParticles; i++){
			particles[i] = new PSOParticle(box, orig_sys, c1, c2);
			currentBest = particles[i].calculateFitness(orig_sys);
			if(currentBest > gBest.calculateFitness(orig_sys))
				updateGBest(particles[i]);//update currentBest particle

		}
	}

	PSOSwarm::updateParticles(System* orig_sys, double c1, double c2){
		double currentBest;
		for(int i=0; i<nParticles; i++){
			particles[i].updateVelocity(gBest,c1,c2);
			//MISSING LIMIT MIN-MAX VEL
			particles[i].updatePosition();
			currentBest = particles[i].calculateFitness(orig_sys);
			if(currentBest > gBest.calculateFitness(orig_sys))
				updateGBest(particles[i]);//update currentBest particle
		}
	}
	PSOSwarm::updateGBest(PSOParticle particle){
		this->gBest = particle;
	}
	Vector PSOSwarm::getGBestPosition(){
		return gBest.getPosition();
	}

	PSOSwarm::~PSOSwarm() {
		// TODO Auto-generated destructor stub
	}
}
