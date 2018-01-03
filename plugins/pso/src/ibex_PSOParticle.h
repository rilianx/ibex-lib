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
		void calculateFitness(System* orig_sys);
		void selectBestInternal();
		bool isFeasible();
		bool isBestFeasible();

		Vector getPosition();
		Vector getBestPosition();
		double getValue();
		double getBestValue();
		int getViolations();
		int getBestViolations();

	protected:
		double value;
		double vBest;
		Vector position;
		Vector pBest;
		Vector velocity;
		int violations;
		int vioBest;
	};
}
#endif /* PSO_SRC_IBEX_PSOPARTICLE_H_ */
