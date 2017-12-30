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

	Vector PSOSwarm::executePSO(System* orig_sys, double p){
		double fitness;
		int iterations = 0;
		bool exit = false;
		if (trace) cout << "Particles: " << nParticles << endl;
		if (trace) cout << "Iterations limit: " << limit << endl;
		if (trace) cout << "Box's diameter ponderator (p): " << p << endl;
		if (trace) cout << "Cognitive parameter (c1): " << c1 << endl;
		if (trace) cout << "Social parameter (c2): " << c2 << endl;

		// ** Initialize particles on random places**
		if (trace) cout << "\033[0;33m# Initialize population of particles with random position and velocity" << endl;
		for(int i=0; i<nParticles; i++){
			particlesArray[i] = new PSOParticle(orig_sys, c1, c2);
			particlesArray[i]->calculateFitness(orig_sys);
			if(gBest){
				selectParticle(particlesArray[i]);
			}else{
				gBest = particlesArray[i];
				if (trace) cout << "\033[0;32mgBest fitness["<< gBest->getBestViolations() <<"]:" << gBest->getBestValue() << "\t First" << endl;
				if (trace) cout << "\033[0;31mAt: " << gBest->getBestPosition() << endl;
			}
		}

		// ** Iterations **
		if (trace) cout << "\033[0;33m# ITERATIONS" << endl;
		do{
			iterations++;
			for(int i=0; i<nParticles; i++){
				// # Update velocity and position of every particle.
				particlesArray[i]->updateVelocityAndPosition(orig_sys,gBest,c1,c2, p);

				// # Evaluate objective (fitness) of every particle.
				particlesArray[i]->calculateFitness(orig_sys);

				// # Select the particle with best fitness (min) and save as gBest.
				selectParticle(particlesArray[i]);
			}
		}while(exit == false && iterations < limit);
		if (trace) cout << "\033[0;33m# END ITERATIONS" << endl;
		if (trace) cout << "\033[0mvBest["<< gBest->getBestValue() <<"] violBest["<< gBest->getBestViolations() <<"]" << endl;
		return getGBestPosition();
	}

	void PSOSwarm::selectParticle(PSOParticle* particle){
		// Feasibility rule (K. Deb - 2000)
		// both feasible: select better fitness.
		if(gBest->isBestFeasible() == true && particle->isFeasible() == true){
			if(particle->getValue() < gBest->getBestValue()){
				gBest = particle;	//update globalBest particle
				if (trace) cout << "\033[0;32mgBest fitness["<< gBest->getBestViolations() <<"]:" << gBest->getBestValue() << "\t Criteria[1]" << endl;
				if (trace) cout << "\033[0;31mAt: " << gBest->getBestPosition() << endl;
			}
		}

		// only one is feasible: select feasible.
		if(gBest->isBestFeasible() == false && particle->isFeasible() == true){
			gBest = particle;	//update globalBest particle
			if (trace) cout << "\033[0;32mgBest fitness["<< gBest->getBestViolations() <<"]:" << gBest->getBestValue() << "\t Criteria[2]" << endl;
			if (trace) cout << "\033[0;31mAt: " << gBest->getBestPosition() << endl;
		}

		// both infeasible: select better fitness.
		if(gBest->isBestFeasible() == false && particle->isFeasible() == false){
			if(particle->getViolations() < gBest->getBestViolations()){
				gBest = particle;	//update globalBest particle
				if (trace) cout << "\033[0;32mgBest fitness["<< gBest->getBestViolations() <<"]:" << gBest->getBestValue() << "\t Criteria[3]" << endl;
				if (trace) cout << "\033[0;31mAt: " << gBest->getBestPosition() << endl;
			}
		}
	}

	Vector PSOSwarm::getGBestPosition(){
		return gBest->getBestPosition();
	}

	double PSOSwarm::getGBestValue(){
		return gBest->getBestValue();
	}

	PSOSwarm::~PSOSwarm() {
		for(int i=0; i<nParticles; i++){
			delete particlesArray[i];
		}
		delete this;
	}
}
