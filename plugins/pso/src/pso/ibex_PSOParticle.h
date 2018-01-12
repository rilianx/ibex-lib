/*
 * ibex_PSOParticle.h
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#ifndef PSO_SRC_IBEX_PSOPARTICLE_H_
#define PSO_SRC_IBEX_PSOPARTICLE_H_
#include "ibex_System.h"
#include "ibex_BufferPSO.h"

namespace ibex{
	class PSOParticle {
	public:
		PSOParticle(BufferPSO* buffer, System* orig_sys, double c1, double c2);
		virtual ~PSOParticle();

		void updateVelocityAndPosition(BufferPSO* buffer, System* orig_sys, Vector gBest,double c1, double c2, double p);
		void calculateFitness(System* orig_sys);
		void selectBestInternal(BufferPSO* buffer);
		bool isFeasible();
		bool isBestFeasible();
		Vector getPosition();
		Vector getBestPosition();
		double getValue();
		double getBestValue();
		double getPenalty();
		double getBestPenalty();
		bool isValid(BufferPSO* buffer);

	protected:
		double value;
		double vBest;
		Vector position;
		Vector pBest;
		Vector velocity;
		double penalty;
		double peBest;
	};
}
#endif /* PSO_SRC_IBEX_PSOPARTICLE_H_ */
