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

	Vector PSOSwarm::executePSO(System* orig_sys){
		double currentFitness = 0;
		int iterations = 0;
		bool exit = false;

		// ** Initialize particles on random places**
		if (trace) cout << "#Initialize particles" << endl;
		for(int i=0; i<nParticles; i++){
			if (trace) cout << "particle: " << i << endl;
			particlesArray[i] = new PSOParticle(orig_sys, c1, c2);
			if (trace) cout << "pos: " << particlesArray[i]->getPosition() << endl;
			currentFitness = particlesArray[i]->calculateFitness(orig_sys);
			if (trace) cout << "fitness: " << particlesArray[i]->calculateFitness(orig_sys) << endl;
			if(gBest){
				if(currentFitness > gBest->calculateFitness(orig_sys))
					gBest = particlesArray[i];	//update currentBest particle
			}else{
				gBest = particlesArray[i];
			}
		}

		// ** Iterations **
		if (trace) cout << "#Iterations" << endl;
		do{
			iterations++;
			// ** Move particles **
			if (trace) cout << "#Move particles" << endl;
			for(int i=0; i<nParticles; i++){
				if (trace) cout << "particle: " << i << endl;
				// ** Update velocity and position of every particle **
				//MISSING LIMIT MIN-MAX VEL
				//if (trace) cout << "oldPos: " << particlesArray[i]->getPosition() << endl;
				particlesArray[i]->updateVelocityAndPosition(orig_sys,gBest,c1,c2);
				//if (trace) cout << "newPos: " << particlesArray[i]->getPosition() << endl;
				// ** Calculate fitness of every particle **
				currentFitness = particlesArray[i]->calculateFitness(orig_sys);
				// ** Choose better fitness
				//if(currentFitness > gBest->calculateFitness(orig_sys))
					//gBest = particlesArray[i];	//update currentBest particle
			}
		}while(exit == false || iterations >= limit);
		return getGBestPosition();
	}

	Vector PSOSwarm::getGBestPosition(){
		return gBest->getPosition();
	}

	PSOSwarm::~PSOSwarm() {
		for(int i=0; i<nParticles; i++){
			delete particlesArray[i];
		}
		delete this;
	}
}
