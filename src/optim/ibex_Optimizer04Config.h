//============================================================================
//                                  I B E X                                   
// File        : ibex_Optimizer04Config.h
// Author      : Bertrand Neveu, Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Dec 11, 2014
// Last Update : Dec 04, 2020
//============================================================================

#ifndef __IBEX_OPTIMIZER_04_CONFIG_
#define __IBEX_OPTIMIZER_04_CONFIG_

#include "ibex_OptimizerConfig.h"
#include "ibex_OptimMemory.h"
#include "ibex_Linearizer.h"
#include "ibex_ExtendedSystem.h"

namespace ibex {


class Optimizer04Config : public OptimizerConfig, protected Memory {
public:

	//propags
	static double hc4_ratio_propag;
	static double hc44cid_propag;
	static double hc44xn_propag;
	static double hc43bcidhc4_propag;
	static double poly_hc4_propag;
	static double all_propag;

	static double ratio_bisect;

	//3bcid
	static double s3b;
	static double scid;
	static double var_min_width;
	//only acid
	static bool optim;
	static double ct_ratio;

	Optimizer04Config(int argc, char** argv);

	//lr=xn
	static double relax_ratio;
	static int poly_maxiter;
	static double poly_eps;

	//finder
	static bool in_hc4;

	/** The system */
	System *sys;

protected:

	// ============================================================================
	virtual void load_sys(const char* filename);

	virtual unsigned int nb_var();

	virtual Ctc& get_ctc();

	virtual Bsc& get_bsc();

	virtual Linearizer* get_linear_relax();

	virtual LoupFinder& get_loup_finder();

	virtual CellBufferOptim& get_cell_buffer();

	virtual int goal_var();
	// ============================================================================


	NormalizedSystem *norm_sys;
	ExtendedSystem *ext_sys;
	std::string filename;
	std::string filtering;
	std::string linearrelaxation;
	std::string bisection;
	std::string strategy;

	int beamsize;
};

} // namespace ibex

#endif // __IBEX_OPTIMIZER_04_CONFIG_
