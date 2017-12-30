/*
 * ibex_PSOParticle.h
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#ifndef PSO_SRC_IBEX_PSOPARTICLE_H_
#define PSO_SRC_IBEX_PSOPARTICLE_H_
#include "ibex_System.h"

namespace ibex{
	class PSOParticle {
	public:
		PSOParticle(System* orig_sys, double c1, double c2);
		virtual ~PSOParticle();

		void updateVelocityAndPosition(System* orig_sys, PSOParticle* gBest,double c1, double c2, double p);
		double calculateFitness(System* orig_sys);
		Vector getBestPosition();
		double getBestValue();

	protected:
		double value;
		double vBest;
		Vector position;
		Vector pBest;
		Vector velocity;
		bool feasible;
	};
}
#endif /* PSO_SRC_IBEX_PSOPARTICLE_H_ */
