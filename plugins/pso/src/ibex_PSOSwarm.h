/*
 * ibex_PSOSwarm.h
 *
 *  Created on: 21-12-2017
 *      Author: cduartes
 */

/*
 * PSO Algorithm
 * # Initialize population of particles with random position and velocity.
 * # Evaluate objective (fitness) of every particle.
 * 		- Select this position as best position for particle.
 * # Select the particle with best fitness (min) and save as gBest.
 * ITERATIONS
 * # Update velocity and position of every particle.
 * # Evaluate objective (fitness) of every particle.
 * 		- Select this position as best position for particle.
 * # Select the particle with best fitness (min) and save as gBest.
 * END ITERATIONS
 *
 * Notes:
 * 		This version limit particles to box boundaries simulating stick to edges.
 * 		Velocity of particles is limited by "p", a ponderation factor for box's diameter, setting as max velocity p*diam.
 *		Iterations ends when limit are exceeded.
 */

#ifndef PSO_SRC_IBEX_PSOSWARM_H_
#define PSO_SRC_IBEX_PSOSWARM_H_
#include "ibex_System.h"
#include "ibex_PSOParticle.h"

namespace ibex{
	class PSOSwarm {
	public:
		PSOSwarm(double c1, double c2, int nParticles, int limit);
		virtual ~PSOSwarm();

		Vector executePSO(System* orig_sys, double p);
		Vector getGBestPosition();

		static bool trace;

	protected:
		int nParticles;
		int limit;
		//Stock** stockArray = new Stock*[n];
		PSOParticle** particlesArray; //array of particles
		PSOParticle* gBest;
		double c1;
		double c2;
	};
}
#endif /* PSO_SRC_IBEX_PSOSWARM_H_ */
