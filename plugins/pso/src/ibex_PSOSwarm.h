/*
 * ibex_PSOSwarm.h
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

#ifndef PSO_SRC_IBEX_PSOSWARM_H_
#define PSO_SRC_IBEX_PSOSWARM_H_
#include "ibex_System.h"
#include "ibex_PSOParticle.h"

namespace ibex{
	class PSOSwarm {
	public:
		PSOSwarm(IntervalVector box, System* orig_sys, double c1, double c2, int nParticles, int limit);
		virtual ~PSOSwarm();

		void updateParticles(System* orig_sys, double c1, double c2);
		void updateGBest(PSOParticle particle);
		Vector getGBestPosition();

	protected:
		int nParticles;
		int limit;
		PSOParticle particles[];
		PSOParticle gBest;
		double c1;
		double c2;
	};
}
#endif /* PSO_SRC_IBEX_PSOSWARM_H_ */
