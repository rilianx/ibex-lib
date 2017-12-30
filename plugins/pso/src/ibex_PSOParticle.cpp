/*
 * ibex_PSOParticle.cpp
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#include "ibex_PSOParticle.h"
//#include "ibex_PSOSwarm.h"
#include <iostream> //DELETE THIS
#include <climits>
#include "stdlib.h"

using namespace std; //DELETE THIS
namespace ibex {
	PSOParticle::PSOParticle(System* orig_sys, double c1, double c2) : position(orig_sys->box.random()), pBest(position), velocity(Vector::zeros(orig_sys->box.size())){
		value = calculateFitness(orig_sys);
		vBest = value;
	}

	void PSOParticle::updateVelocityAndPosition(System* orig_sys, PSOParticle* gBest, double c1, double c2, double p) {
		double maxVel;	//max velocity value
		// ** Update Velocity **
		double rand1 = (double)random()/(double)RAND_MAX;
		double rand2 = (double)random()/(double)RAND_MAX;
		Vector newV = velocity + (c1*rand1*(pBest - position) + (c2*rand2*(gBest->getBestPosition()-position)));
		velocity = newV;

		// check for max velocity exceeded
		for(int i=0; i<orig_sys->box.size(); i++){
			maxVel = p*orig_sys->box[i].diam();
			if(velocity[i] > maxVel)
				velocity[i] = maxVel;
			else if(velocity[i] < -maxVel)
				velocity[i] = -maxVel;
		}

		// ** Update Position based on new Velocity **
		Vector newP = position + velocity;
		position = newP;

		// If position exceed box size, stop particle and change velocity's direction
		for(int i=0; i<orig_sys->box.size(); i++){
			if(position[i] > orig_sys->box[i].ub()){
				position[i] = orig_sys->box[i].ub();	// stick to edge (ub)
				velocity[i] = -velocity[i];				// change direction
			}
			if(position[i] < orig_sys->box[i].lb()){
				position[i] = orig_sys->box[i].lb();	// stick to edge (lb)
				velocity[i] = -velocity[i]; 			// change direction
			}
		}
		value = calculateFitness(orig_sys);	//update value (fitness)
		if(value < vBest){					//update best position and best fitness
			pBest = position;
			vBest = value;
		}
	}

	double PSOParticle::calculateFitness(System* orig_sys){
		int sum = 0;
		double fitness;
		for(int i=0; i<orig_sys->ctrs.size();i++){
			Interval eval = orig_sys->ctrs[i].f.eval(position);
			if((orig_sys->ctrs[i].op == LEQ || orig_sys->ctrs[i].op == LT) && eval.ub()>0.0){
				sum += eval.ub();
				//cout << "error("<< i <<"):" << eval.ub() << endl;
			}
			else if((orig_sys->ctrs[i].op == GEQ || orig_sys->ctrs[i].op == GT) && eval.lb()<0.0){
				sum += -eval.ub();
				//cout << "error("<< i <<"):" << -eval.lb() << endl;
			}
			else if(orig_sys->ctrs[i].op == EQ ){
				sum += abs(eval).ub();
				//cout << "error("<< i <<"):" << abs(eval).ub() << endl;
			}
		}
		if(sum > 0){
			fitness = sum;
			feasible = false;
		}else{
			fitness = orig_sys->goal->eval(position).ub();
			feasible = true;
		}

		return fitness;
	}

	Vector PSOParticle::getBestPosition(){
		return pBest;
	}

	double PSOParticle::getBestValue(){
		return vBest;
	}

	PSOParticle::~PSOParticle() {
		delete this;
	}

}

