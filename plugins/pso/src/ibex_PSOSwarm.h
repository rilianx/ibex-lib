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
		PSOSwarm(System* orig_sys);
		virtual ~PSOSwarm();


	protected:
		int numParticles;
		PSOParticle particles;
	};
}
#endif /* PSO_SRC_IBEX_PSOSWARM_H_ */
