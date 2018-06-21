//============================================================================
//                                  I B E X                                   
// File        : Round-robin bisector
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : May 8, 2012
// Last Update : Dec 25, 2017
//============================================================================

#ifndef __IBEX_ROUND_ROBIN_H__
#define __IBEX_ROUND_ROBIN_H__

#include "ibex_Bsc.h"
#include "ibex_Backtrackable.h"

namespace ibex {

/**
 * \ingroup bisector
 *
 * \brief Round-robin bisector.
 *
 * See #choose_var(const Cell& cell) for details.
 *
 */
class RoundRobin : public Bsc {
public:
	/**
	 * \brief Create a bisector with round-robin heuristic.
	 *
	 * \param prec             - allow to skip variables whose domain is too small
	 * \param ratio (optional) - the ratio between the diameters of the left and the right parts of the
	 *                           bisected interval. Default value is 0.45.
	 */
  RoundRobin(double prec, double ratio=Bsc::default_ratio());

	/**
	 * \brief Create a bisector with round-robin heuristic.
	 *
	 * \param prec             - allow to skip variables whose domain is too small
	 * \param ratio (optional) - the ratio between the diameters of the left and the right parts of the
	 *                           bisected interval. Default value is 0.45.
	 */
  RoundRobin(const Vector& prec, double ratio=Bsc::default_ratio());

	/**
	 * \brief Return next variable to be bisected.
	 *
	 * called by Bsc::bisect(...)
	 *
	 * The next variable is the first one after the last bisected (using variable ordering based
	 * based on indices, and in a circular way) with a domain sufficiently large. Hence, if the
	 * domain of the var immediatly following the last one is wider than \a w, it will be bisected.
	 * If the cell is the root cell (no previous bisection), the variable inspected first is the
	 * first one. In case all the variables have domain diameter less than \a w, a
	 * NoBisectableVariableException is raised.
	 */
	virtual BisectionPoint choose_var(const Cell& cell);

	/**
	 * \brief Add an instance of #ibex::BisectedVar to the backtrackable data of the root cell.
	 */
	void add_backtrackable(Cell& root);

	/**
	 * \brief Ratio to choose the split point.
	 *
	 * Ratio between the diameters of the left and right parts of a bisected
	 * interval.
	 */
	const double ratio;
	/*
  protected :
  int var_select(const IntervalVector& box, int n, int &last_var);
	*/
};



class RoundRobinNvar : public RoundRobin {
  public :
  /**
   * \brief Create a bisector with round-robin heuristic.on the first nbvars variables
   *
   * \param nbvars : number of variables to bisect ; the bisection will be done only on the first nbvars variables
   * \param w - the minimum width (diameter) an interval must have to be bisected (see #w).
   *            Default value is 1e-07.
   * \param ratio (optional) - the ratio between the diameters of the left and the right parts of the
   *            bisected interval. Default value is 0.45.
   */


  RoundRobinNvar (int nbvars, double prec, double ratio=Bsc::default_ratio())
      : RoundRobin (prec,ratio), nbvars(nbvars) {;};

  RoundRobinNvar (int nbvars,const Vector& prec, double ratio=Bsc::default_ratio())
      : RoundRobin (prec,ratio), nbvars(nbvars) {;};

  BisectionPoint choose_var(const Cell& cell);
    //    pair<IntervalVector,IntervalVector> bisect(Cell& cell);
  /*
    virtual  std::pair<IntervalVector,IntervalVector> bisect(const IntervalVector& box, int& last_var);
 protected :
  */
    int nbvars;

  };

/*

   
class RoundRobinQInter : public RoundRobin {
  public :
 
    
  RoundRobinQInter (double prec, double ratio=Bsc::default_ratio())
      : RoundRobin ( prec,ratio) {;};

  RoundRobinQInter (const Vector& prec, double ratio=Bsc::default_ratio())
      : RoundRobin (prec,ratio) {;};
  //  std::pair<IntervalVector,IntervalVector> bisect(Cell& cell);
};
*/
} // end namespace ibex

#endif // __IBEX_ROUND_ROBIN_H__


