/*
 * ibex_PSOSwarm.cpp
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#include "ibex_PSOSwarm.h"

namespace ibex{
	PSOSwarm::PSOSwarm(double c1, double c2, int nParticles, int limit) : c1(c1), c2(c2), gBest(NULL){
		this->nParticles = nParticles;
		particles = new PSOParticle*[nParticles];
		this->limit = limit;
	}

	Vector PSOSwarm::executePSO(System* orig_sys){
		double currentFitness;
		int iterations = 0;
		bool exit = false;

		// ** Initialize particles on random places**
		for(int i=0; i<nParticles; i++){
			particles[i] = new PSOParticle(orig_sys, c1, c2);
			currentFitness = particles[i].calculateFitness(orig_sys);
			if(currentFitness > gBest.calculateFitness(orig_sys))
				gBest = particles[i];	//update currentBest particle
		}

		// ** Iterations **
		do{
			iterations++;
			// ** Move particles **
			for(int i=0; i<nParticles; i++){
				// ** Update velocity and position of every particle **
				//MISSING LIMIT MIN-MAX VEL
				particles[i].updateVelocityAndPosition(orig_sys,gBest,c1,c2);
				// ** Calculate fitness of every particle **
				currentFitness = particles[i].calculateFitness(orig_sys);
				// ** Choose better fitness
				if(currentFitness > gBest.calculateFitness(orig_sys))
					gBest = particles[i];	//update currentBest particle
			}
		}while(exit == false || iterations >= limit);
		return getGBestPosition();
	}

	Vector PSOSwarm::getGBestPosition(){
		return gBest.getPosition();
	}

	PSOSwarm::~PSOSwarm() {
		for(int i=0; i<nParticles; i++){
			delete particles[i];
		}
		delete this;
	}
}
