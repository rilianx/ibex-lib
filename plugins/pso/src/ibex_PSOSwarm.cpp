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

	Vector PSOSwarm::executePSO(System* orig_sys, bool bounce){
		double currentFitness = 0;
		int iterations = 0;
		bool exit = false;
		if (trace) cout << "Particles: " << nParticles << endl;
		if (trace) cout << "Iterations limit: " << limit << endl;
		if (trace) cout << "Coef c1: " << c1 << endl;
		if (trace) cout << "Coef c2: " << c2 << endl;

		// ** Initialize particles on random places**
		if (trace) cout << "# Initialize population of particles with random position and velocity" << endl;
		for(int i=0; i<nParticles; i++){
			particlesArray[i] = new PSOParticle(orig_sys, c1, c2);
			currentFitness = particlesArray[i]->calculateFitness(orig_sys);
			if(gBest){
				if(currentFitness > gBest->calculateFitness(orig_sys))
					gBest = particlesArray[i];	//update currentBest particle
			}else{
				gBest = particlesArray[i];
			}
		}
		if (trace) cout << "Current gBest fitness: " << gBest->calculateFitness(orig_sys) << endl;
		if (trace) cout << "At: " << gBest->getBestPosition() << endl;

		// ** Iterations **
		if (trace) cout << "# ITERATIONS" << endl;
		do{
			iterations++;
			// ** Move particles **
			for(int i=0; i<nParticles; i++){
				// ** #Update velocity and position of every particle **
				particlesArray[i]->updateVelocityAndPosition(orig_sys,gBest,c1,c2, bounce);
				// ** #Calculate fitness of every particle **
				currentFitness = particlesArray[i]->calculateFitness(orig_sys);
				// ** #Choose better fitness
				if(currentFitness < gBest->calculateFitness(orig_sys)){
					gBest = particlesArray[i];	//update currentBest particle
					if (trace) cout << "Current gBest fitness: " << gBest->calculateFitness(orig_sys) << endl;
					if (trace) cout << "At: " << gBest->getBestPosition() << endl;
				}
			}
		}while(exit == false && iterations < limit);
		if (trace) cout << "# END ITERATIONS" << endl;
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
