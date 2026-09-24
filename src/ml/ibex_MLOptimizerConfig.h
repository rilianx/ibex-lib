//============================================================================
//                                  I B E X
// File        : ibex_MLOptimizerConfig.h
// Author      : Ignacio Araya
// License     : See the LICENSE file
// Created     : 2026
//============================================================================

#ifndef __IBEX_ML_OPTIMIZER_CONFIG_H__
#define __IBEX_ML_OPTIMIZER_CONFIG_H__

#include "ibex_DefaultOptimizerConfig.h"

#include <string>

namespace ibex {

/**
 * \ingroup ml
 *
 * \brief IbexOpt's default configuration, with a selectable bisector.
 *
 * Everything else -- the contractor, the loup finder, the cell buffer -- stays
 * exactly what `ibexopt` uses, so that a comparison between bisectors varies
 * one thing only.
 *
 * \remark This class exists as a separate level of inheritance on purpose.
 * #DefaultOptimizerConfig::get_bsc() is called from Optimizer's constructor;
 * a class deriving from *both* the config and Optimizer cannot override it,
 * because at that point its own vptr is not installed yet. Sitting between the
 * two, this class is a fully constructed base by the time Optimizer's
 * constructor runs, so the override is dispatched.
 */
class MLOptimizerConfig : public DefaultOptimizerConfig {
public:

	/**
	 * \brief The bisectors that can be compared.
	 *
	 * All of them work on the extended system and all but #BSC_ROUNDROBIN
	 * delegate to OptimLargestFirst when their own criterion does not apply --
	 * the same fallback chain IbexOpt uses. #BSC_LSMEAR_GUARD runs IbexOpt's
	 * own LSmear and hands over to RoundRobin, for good, once LSmear keeps
	 * bisecting its parent's variable (see BscHijackGuard).
	 */
	typedef enum {
		BSC_LSMEAR_MG,     //!< LSmear, variant LSMEAR_MG: IbexOpt's default
		BSC_LSMEAR,        //!< LSmear, variant LSMEAR (Jacobian over the box)
		BSC_SMEARSUMREL,   //!< SmearSumRelative: LSmear's own fallback, alone
		BSC_SMEARSUM,      //!< SmearSum (Hansen)
		BSC_SMEARMAX,      //!< SmearMax (Kearfott)
		BSC_SMEARMAXREL,   //!< SmearMaxRelative
		BSC_LARGESTFIRST,  //!< OptimLargestFirst: widest domain, objective aware
		BSC_ROUNDROBIN,    //!< RoundRobin: the naive baseline
		BSC_LSMEAR_GUARD   //!< LSmear until hijacked, then RoundRobin (BscHijackGuard)
	} Bisector;

	/**
	 * \brief The linear relaxation the X-Newton step is built on.
	 *
	 * IbexOpt relaxes each constraint by a Taylor form taken at a corner of the
	 * box (#RELAX_XTAYLOR). Affine arithmetic instead carries a linear form
	 * plus an accumulated error term through the whole evaluation, which keeps
	 * the dependency between occurrences of the same variable that an interval
	 * Taylor form throws away.
	 */
	typedef enum {
		RELAX_XTAYLOR,   //!< LinearizerXTaylor: what ibexopt uses
		RELAX_AFFINE,    //!< LinearizerAffine2
		RELAX_BOTH       //!< both, composed: more constraints, tighter hull
	} Relaxation;

	MLOptimizerConfig(const System& sys, double rel_eps_f, double abs_eps_f,
			double eps_h, bool rigor, bool inHC4, bool kkt,
			double random_seed, const Vector& eps_x,
			Bisector bisector=BSC_LSMEAR_MG,
			Relaxation relaxation=RELAX_XTAYLOR);

	/** \brief Which bisector this configuration builds. */
	Bisector get_bisector() const;

	/** \brief Its name, as accepted on the command line. */
	static const char* bisector_name(Bisector b);

	/**
	 * \brief Parse a bisector name.
	 *
	 * \return false if the name is unknown.
	 */
	static bool parse_bisector(const std::string& name, Bisector& out);

	/** \brief All the accepted names, comma separated. */
	static std::string bisector_names();

	/** \brief Which relaxation this configuration builds. */
	Relaxation get_relaxation() const;

	/** \brief Its name, as accepted on the command line. */
	static const char* relaxation_name(Relaxation r);

	/** \brief Parse a relaxation name. \return false if unknown. */
	static bool parse_relaxation(const std::string& name, Relaxation& out);

	/** \brief All the accepted names, comma separated. */
	static std::string relaxation_names();

protected:

	virtual Bsc& get_bsc() override;
	virtual Ctc& get_ctc() override;

	Bisector bisector;
	Relaxation relaxation;
	Bsc* bsc_cache;
	Ctc* ctc_cache;
};

inline MLOptimizerConfig::Bisector MLOptimizerConfig::get_bisector() const {
	return bisector;
}

inline MLOptimizerConfig::Relaxation MLOptimizerConfig::get_relaxation() const {
	return relaxation;
}

} // end namespace ibex

#endif // __IBEX_ML_OPTIMIZER_CONFIG_H__
