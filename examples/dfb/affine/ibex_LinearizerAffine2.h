/* ============================================================================
 * I B E X - Definition of the Linearizer of Affine2 forms
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENCE.
 *
 * Author(s)   : Jordan Ninin
 * Created     : June 6, 2020
 * ---------------------------------------------------------------------------- */


#ifndef __IBEX_LINEARIZER_AFFINE2_H__
#define __IBEX_LINEARIZER_AFFINE2_H__

#include "ibex_Affine.h"
#include "ibex_System.h"
#include "ibex_Linearizer.h"

//#include <vector>

namespace ibex {

/**
 * \ingroup numeric
 * \brief Affine-based linearization
 *
 * This class is an implementation of the ART algorithm
 * \author Jordan Ninin
 * \date May 2013
 */

void lin_perf_volcar();

class LinearizerAffine2 : public Linearizer {

public:

	LinearizerAffine2 (const System& sys);

	~LinearizerAffine2 ();

	/**
	 * \biref  ART iteration.
	 *
	 *  Linearize the system and performs 2n calls to Simplex in order to reduce the 2 bounds of each variable
	 */
	virtual int linearize(const IntervalVector& box, LPSolver& lp_solver);

	/**
	 * \brief Linealizacion CON IDENTIDAD DE FILAS.
	 *
	 * Igual que `linearize` pero en vez de empujar las restricciones al
	 * solver devuelve, para CADA restriccion y cada lado, la fila en espacio
	 * de `x`, el rango admisible de `fila.x`, y una identidad estable
	 * `2*ctr + lado` (lado 0 = cota superior, 1 = cota inferior).
	 *
	 * Dos diferencias con `linearize`, y las dos hacen falta para poder
	 * CONGELAR la matriz entre llamadas:
	 *
	 *  - emite tambien las restricciones **redundantes** en esta caja, que
	 *    `linearize` descarta; si no, la fila `i` de una llamada no es la
	 *    misma restriccion que la fila `i` de la siguiente;
	 *  - devuelve la fila y la cota por separado, para poder quedarse con una
	 *    fila vieja y corregir solo la cota.
	 *
	 * Devuelve el numero de filas emitidas, o -1 si detecta infactibilidad.
	 */
	int linearize_id(const IntervalVector& box, std::vector<Vector>& rows,
	                 std::vector<Interval>& rango, std::vector<int>& ids,
	                 const std::vector<char>* evaluar = NULL);

	/** \brief Cuantas restricciones tiene el sistema. */
	int nb_ctr() const;

	/** \brief Restriccion de origen de cada fila emitida por la ultima
	 *  llamada a `linearize`, en el orden en que se agregaron al solver.
	 *  Lo usa la sonda del lagrangiano (`DFBH_LAGR`). */
	std::vector<int> fila_ctr;

	const System& sistema() const { return sys; }

	/** \brief Forma afin de la restriccion `ctr` sobre `box`; vacia si la
	 *  evaluacion no la produce. Solo lectura, para sondas. */
	Affine2 forma_afin(int ctr, const IntervalVector& box);

	/** \brief Variables de las que depende la restriccion `ctr`. Permite
	 *  re-linealizar solo las restricciones cuyas variables se movieron. */
	const std::vector<int>& used_vars(int ctr) const;

	/**
	 * \brief Generation of a linear approximation of the inner region
	 *
	 */
	int inlinearization(const IntervalVector& box, LPSolver& lp_solver);

	/**
	 * \brief Generation of a linear approximation of the linear objective function
	 *
	 */
	bool goal_linearization(const IntervalVector& box, LPSolver& lp_solver);

private:
	/**
	 * \brief The system
	 */
	const System& sys;

	/**
	 * \brief Affine evaluator for the goal function (if any)
	 */
	AffineEval<AF_Default>* goal_af_evl;

	/**
	 * \brief Affine evaluators for the constraints functions
	 */
	AffineEval<AF_Default>** ctr_af_evl;
};

} // end namespace ibex

#endif /* __IBEX_LINEAR_RELAX_AFFINE2_H__ */

