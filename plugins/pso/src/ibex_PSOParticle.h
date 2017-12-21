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
		PSOParticle(IntervalVector box, System* orig_sys, double c1, double c2);
		virtual ~PSOParticle();

		void updateVelocity(PSOParticle gBest);
		void updatePosition();
		void setBestPosition(Vector position, double value);
		double calculateFitness(System* orig_sys);

	protected:
		double value;
		double vBest;
		Vector position;	// must be box to solve
		double pBest;		// must be box to solve
		double velocity;
		double c1;
		double c2;
	};
}

#endif /* PSO_SRC_IBEX_PSOPARTICLE_H_ */
