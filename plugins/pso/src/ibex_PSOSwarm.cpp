/*
 * ibex_PSOSwarm.cpp
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#include "ibex_PSOSwarm.h"
#include <iostream> //DELETE THIS

using namespace std; //DELETE THIS
namespace ibex{
	bool PSOSwarm::trace = false;

	PSOSwarm::PSOSwarm(double c1, double c2, int nParticles, int limit) : c1(c1), c2(c2), gBest(NULL){
		this->nParticles = nParticles;
		particlesArray = new PSOParticle*[nParticles];
		this->limit = limit;
	}

	Vector PSOSwarm::executePSO(System* orig_sys, bool bounce, double constrictorMod){
		double currentFitness;
		int iterations = 0;
		bool exit = false;
		if (trace) cout << "Particles: " << nParticles << endl;
		if (trace) cout << "Iterations limit: " << limit << endl;
		if (trace) cout << "Constrictor modificator: " << constrictorMod << endl;
		if (trace) cout << "Cognitive parameter (c1): " << c1 << endl;
		if (trace) cout << "Social parameter (c2): " << c2 << endl;

		// ** Initialize particles on random places**
		if (trace) cout << "\033[0;33m# Initialize population of particles with random position and velocity" << endl;
		for(int i=0; i<nParticles; i++){
			particlesArray[i] = new PSOParticle(orig_sys, c1, c2);
			currentFitness = particlesArray[i]->calculateFitness(orig_sys);
			if(gBest){
				if(currentFitness < gBest->calculateFitness(orig_sys,gBest->getBestPosition()))
					gBest = particlesArray[i];	//update globalBest particle
			}else{
				gBest = particlesArray[i];
			}
		}
		if (trace) cout << "\033[0;32mgBest fitness: " << gBest->calculateFitness(orig_sys,gBest->getBestPosition()) << endl;
		if (trace) cout << "\033[0;31mAt: " << gBest->getBestPosition() << endl;

		// ** Iterations **
		if (trace) cout << "\033[0;33m# ITERATIONS" << endl;
		do{
			iterations++;
			for(int i=0; i<nParticles; i++){
				currentFitness = 0;

				// # Update velocity and position of every particle.
				particlesArray[i]->updateVelocityAndPosition(orig_sys,gBest,c1,c2, bounce, constrictorMod);

				// # Evaluate objective (fitness) of every particle.
				currentFitness = particlesArray[i]->calculateFitness(orig_sys);

				// # Select the particle with best fitness (min) and save as gBest.
				if(currentFitness < gBest->vBest){
					gBest = particlesArray[i];	//update globalBest particle
					if (trace) cout << "\033[0;32mgBest fitness:" << gBest->calculateFitness(orig_sys,gBest->getBestPosition()) << endl;
					if (trace) cout << "\033[0;31mgAt: " << gBest->getBestPosition() << endl;
				}
			}
		}while(exit == false && iterations < limit);
		if (trace) cout << "\033[0;33m# END ITERATIONS" << endl;
		return getGBestPosition();
	}

	Vector PSOSwarm::getGBestPosition(){
		return gBest->getBestPosition();
	}

	PSOSwarm::~PSOSwarm() {
		for(int i=0; i<nParticles; i++){
			delete particlesArray[i];
		}
		delete this;
	}
}
