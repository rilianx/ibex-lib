/*
 * ibex_PSOSwarm.cpp
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#include "ibex_PSOSwarm.h"

using namespace std;
namespace ibex{
	bool PSOSwarm::trace = false;
	bool plot = true;

	/*
	 * Initialize swarm
	 * Instance particles and set parameters
	 * Initialize particles
	 */
	//TODO pair array ExtArray<c_i,&vector_i>
	PSOSwarm::PSOSwarm(TreeCellOpt* tree, System* orig_sys, double c1, double c2, double x, int nParticles, int limit) :
			tree(tree), c1(c1), c2(c2), x(x), nParticles(nParticles), orig_sys(orig_sys), limit(limit),
			gBest(Vector::zeros(orig_sys->box.size())), gpenalty(0.0), gValue(0.0), initialized(false){

		particlesArray = new PSOParticle*[nParticles];
		for(int i=0; i<nParticles; i++) particlesArray[i]=NULL;
	}

	void PSOSwarm::resetPSO(double loup){
	  if(trace) cout << "reset_PSO" << endl;
		for(int i=0; i<nParticles; i++){
			if(!particlesArray[i]) particlesArray[i] = new PSOParticle(tree, orig_sys);
			else particlesArray[i]->initialize(tree, orig_sys);
			//cout << "particle["<<i<<"]:" << particlesArray[i]->getBestValue() << "+" << particlesArray[i]->getBestPenalty()<< endl;

			// Update the gBest if the pBest is better than the gBest
			if(i==0 || PSOParticle::compute_fitness(particlesArray[i]->getBestValue(),particlesArray[i]->getBestPenalty(), loup)
			< PSOParticle::compute_fitness(gValue,gpenalty, loup) ){
				gBest = particlesArray[i]->getBestPosition();
				gValue = particlesArray[i]->getBestValue();
				gpenalty = particlesArray[i]->getBestPenalty();

			}

		}

    initialized=true;
	}

	/*
	 * Reset particles data and information of swarm
	 */
	void PSOSwarm::resetGBest(double loup){
	  if(trace) cout << "resetGBest" << endl;
		for(int i=0; i<nParticles; i++){
			if(!tree->search(particlesArray[i]->getBestPosition()))
				particlesArray[i]->initialize(tree,orig_sys);

			// Update the gBest if the pBest is better than the gBest
			if(i==0 || PSOParticle::compute_fitness(particlesArray[i]->getBestValue(),particlesArray[i]->getBestPenalty(), loup)
			< PSOParticle::compute_fitness(gValue,gpenalty, loup) ){
				gBest = particlesArray[i]->getBestPosition();
				gValue = particlesArray[i]->getBestValue();
				gpenalty = particlesArray[i]->getBestPenalty();

			}

		}
		//cout << gValue << endl;
		initialized = true;
		//cout << "initialized" << endl;
	}

	/*
	 * Execute PSO algorithm
	 */
	void PSOSwarm::executePSO(double loup){
		double fitness;
		int iterations = 0;

		if (trace) cout << "Particles: " << nParticles << endl;
		if (trace) cout << "Iterations limit: " << limit << endl;
		if (trace) cout << "Cognitive parameter (c1): " << c1 << endl;
		if (trace) cout << "Social parameter (c2): " << c2 << endl;

		//cout << gValue << endl;
		//vericar que todos los pbest se encuentren dentro de una caja
		for(int i=0; i<nParticles; i++){
			if(!tree->search(particlesArray[i]->getBestPosition()))
				particlesArray[i]->initialize(tree,orig_sys);
		}



		// ** Iterations **
		if (trace) cout << "\033[0;33m# ITERATIONS" << endl;

		while(iterations < limit){


			if (plot) iterationPlot();		//dump iteration to file
			iterations++;
			for(int i=0; i<nParticles; i++){
				// Update velocity and position of every particle.
				particlesArray[i]->updateVelocityAndPosition(tree,orig_sys,gBest,c1,c2,x);

				// Update the pBest if the new position is better than the previous one
				// AND the pBest is located into some box
				particlesArray[i]->update_pBest(tree,orig_sys, loup);

				// Update the gBest if the pBest is better than the gBest
				if(PSOParticle::compute_fitness(particlesArray[i]->getBestValue(),particlesArray[i]->getBestPenalty(), loup)
							< PSOParticle::compute_fitness(gValue,gpenalty, loup) ){
					gBest = particlesArray[i]->getBestPosition();
					gValue = particlesArray[i]->getBestValue();
					gpenalty = particlesArray[i]->getBestPenalty();
					if(trace) cout << "new gbest: " <<  gValue << "+" << gpenalty << endl;
				}


			}

		}

		if (trace) cout << "\033[0;33m# END ITERATIONS" << endl;
		if (trace) cout << "\033[0mvBest["<< getGBestValue() <<"] violBest["<< getGBestPenalty() <<"]" << endl;
	}



	void PSOSwarm::startPlot(char* file_name){
		cout << "to file " << file_name << endl;
		output.open(file_name);
		//output << "[";
	}

	void PSOSwarm::iterationPlot(){
		for(int i=0; i < nParticles ; i++){
			for(int a=0; a < orig_sys->box.size(); a++){
				output << particlesArray[i]->getPosition()[a] << ",";
			}
			output << particlesArray[i]->getPosition()[orig_sys->nb_var] << ";";
		}
		output << getGBestPosition()[0] << "," << getGBestPosition()[1] << endl;
	}

	void PSOSwarm::closePlot(){
		cout << "close file output.txt" << endl;
		//output << "]" << endl;
		output.close();
		//system("python3 plot.py");
	}

	/*
	 * GETTER AND SETTER
	 */


	Vector PSOSwarm::getGBestPosition(){
		return gBest;
	}

	double PSOSwarm::getGBestValue(){
		 return gValue;
	}

	double PSOSwarm::getGBestPenalty(){
		return gpenalty;
	}


	bool PSOSwarm::isGBestFeasible(){
		if (gpenalty > 0)
			return false;
		else
			return true;
	}

	bool PSOSwarm::isInitialized(){
		return initialized;
	}

	PSOSwarm::~PSOSwarm() {
		for(int i=0; i<nParticles; i++){
			delete particlesArray[i];
		}
		delete this;
	}
}
