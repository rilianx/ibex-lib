/*
 * ibex_PSOParticle.cpp
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#include "ibex_PSOParticle.h"
#include "ibex_Interval.h"
#include "ibex_Random.h"
#include <climits>
#include "stdlib.h"

using namespace std; //delete
namespace ibex {
	PSOParticle::PSOParticle(TreeCellOpt* tree, System* orig_sys, double eqeps) :
			position(orig_sys->box.random()), pBest(position),
			velocity(Vector::zeros(orig_sys->box.size())), eqeps(eqeps){
		initialize(tree, orig_sys);
	}

	/*
	 * updates current velocity and position of a particle
	 */
	void PSOParticle::updateVelocityAndPosition(TreeCellOpt* tree, System* orig_sys, Vector gBest, double c1, double c2, double x) {
		double maxVel;	//max velocity value

		// ** Update Velocity **

		double rand1 = RNG::rand(0,1);
		double rand2 = RNG::rand(0,1);


		//variante distinto random por componente del vector
		//Vector newV (position.size());
		//for(int i=0;i<newV.size();i++)
			//newV[i] = x*(velocity[i] + (c1*RNG::rand(0,1)*(pBest[i] - position[i]) +
		//			(c2*RNG::rand(0,1)*(gBest[i]-position[i]))));

		Vector newV = x*(velocity + (c1*rand1*(pBest - position) + (c2*rand2*(gBest-position))));
		//if (velocity == newV) cout << "same speed" << endl;

		velocity = newV;
		//if (velocity == Vector::zeros(orig_sys->box.size())) cout << "speed zero! checkout!" << endl;
		//cout << velocity << endl;
		//velocity = x * orig_sys->box.random();

		//cout << position << endl;

		// ** Update Position based on new Velocity **
		position = position + velocity;
		//cout << position[0] << "," << position[1] << endl;

	}


	double PSOParticle::computePenalty(System* orig_sys, Vector& position){
		double sum = 0;
		for(int i=0; i<orig_sys->ctrs.size();i++){
			Interval eval = orig_sys->ctrs[i].f.eval(position);
			if((orig_sys->ctrs[i].op == LEQ || orig_sys->ctrs[i].op == LT) && eval.ub()>0.0){
				sum += eval.ub();
			}
			else if((orig_sys->ctrs[i].op == GEQ || orig_sys->ctrs[i].op == GT) && eval.lb()<0.0){
				sum += -eval.ub();
			}
			else if(orig_sys->ctrs[i].op == EQ ){
				sum += std::max(0.0, abs(eval).ub() /*- eqeps*/);
			}
		}

    
		return sum;
	}



	/*
	 * check for violations of best position known, returns boolean
	 */
	bool PSOParticle::isBestFeasible(){
		if(peBest > 0)
			return false;
		else
			return true;
	}

	/*
	 * check for violations of current position, returns boolean
	 */
	bool PSOParticle::isFeasible(){
		if(penalty > 0)
			return false;
		else
			return true;
	}

	/*
	 * GETTER AND SETTER
	 */
	Vector PSOParticle::getPosition(){
		return position;
	}

	Vector PSOParticle::getBestPosition(){
		return pBest;
	}

	double PSOParticle::getBestValue(){
		return vBest;
	}

	double PSOParticle::getValue(){
		return value;
	}

	double PSOParticle:: getPenalty(){
		return penalty;
	}

	double PSOParticle:: getBestPenalty(){
		return peBest;
	}

	PSOParticle::~PSOParticle() {
		delete this;
	}

}
