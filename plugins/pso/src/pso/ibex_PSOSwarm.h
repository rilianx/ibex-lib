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
 * 		- Select best position of particle and save it.
 * # Select the particle with best fitness (min) and save as gBest.
 * END ITERATIONS
 *
 * Notes:
 * 		This version limit particles to box boundaries simulating stick to edges.
 * 		Velocity of particles is limited by "p", a ponderation factor for box's diameter, setting as max velocity p*diam.
 *		Iterations ends when limit are exceeded.
 *		with --trace PSO will show in screen when gBest is updated and exporting a file to plot with python
 */

#ifndef PSO_SRC_IBEX_PSOSWARM_H_
#define PSO_SRC_IBEX_PSOSWARM_H_
#include "ibex_System.h"
#include <fstream> // plot file library
#include <iostream> // screen text library
#include <stdlib.h>
#include <iomanip>
#include "ibex_PSOParticle.h"

namespace ibex{
	class PSOSwarm {
	public:
		PSOSwarm(System* orig_sys, double c1, double c2, int nParticles, int limit);
		virtual ~PSOSwarm();

		Vector executePSO(double p);
		void startPlot();
		void iterationPlot();
		void closePlot();
		Vector getGBestPosition();
		double getGBestValue();
		double getGBestPenalty();
		bool isGBestFeasible();
		int getGBestCountViolations();
		static bool trace;

	protected:
		int nParticles;
		int limit;
		PSOParticle** particlesArray; //array of particles
		Vector gBest;
		double c1;
		double c2;
		std::ofstream output;
		System* orig_sys;

		void selectParticle(PSOParticle* particle);
	};
}
#endif /* PSO_SRC_IBEX_PSOSWARM_H_ */
