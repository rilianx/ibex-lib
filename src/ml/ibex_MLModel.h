//============================================================================
//                                  I B E X
// File        : ibex_MLModel.h
// Author      : Ignacio Araya
// License     : See the LICENSE file
// Created     : 2026
//============================================================================

#ifndef __IBEX_ML_MODEL_H__
#define __IBEX_ML_MODEL_H__

#include "ibex_Vector.h"

#include <string>
#include <vector>

namespace ibex {

/**
 * \ingroup ml
 *
 * \brief A branching rule evaluated inside the solver.
 *
 * Scores each variable from a fixed per-variable feature vector (see
 * #MLModel::NB_FEATURES) and the bisector picks the argmax among the
 * candidates. Two families are supported, both stored as plain text so that no
 * inference library has to be linked in:
 *
 * - **linear**: score = bias + w . phi
 * - **gbdt**:   score = base + sum of regression trees over phi
 *
 * Without this, a learned rule can only be compared to LSmear on node counts,
 * which says nothing about wall-clock time: the point of evaluating the model
 * in C++ is to pay its real cost inside the search.
 *
 * ### File format
 *
 * \verbatim
 * ibexml-model 1
 * type linear
 * nfeat 14
 * bias 0.1
 * w -0.3 1.2 0 ... (nfeat values)
 * \endverbatim
 *
 * \verbatim
 * ibexml-model 1
 * type gbdt
 * nfeat 14
 * base 0.0
 * tree 5                      # number of nodes in this tree
 * S 3 1.5 1 4                 # node 0: if phi[3] <= 1.5 go to node 1, else 4
 * S 0 -2.0 2 3                # node 1
 * L 0.25                      # node 2: leaf
 * L -0.75                     # node 3
 * L 1.10                      # node 4
 * tree 1                      # next tree
 * L 0.5
 * \endverbatim
 *
 * Blank lines and lines starting with '#' are ignored. Node indices are
 * relative to the tree they belong to. See `export_linear` / `export_gbdt` in
 * `python/ibexml.py`.
 */
class MLModel {
public:

	/**
	 * \brief Number of per-variable features the model is fed with.
	 *
	 * The order is fixed and must match `encode()` in `python/ibexml.py`
	 * (VAR_FEATURES followed by SCORE_FEATURES):
	 *
	 *  0 lb            5 rel_diam      10 smear_sum
	 *  1 ub            6 is_goal       11 smear_max
	 *  2 diam          7 bisectable    12 smear_sum_rel
	 *  3 mid           8 too_small     13 lsmear
	 *  4 mag           9 dual          14 dual_ok
	 *
	 * Feature 14 is the same for every variable of a node. It is there because
	 * 9 and 13 are set to 0 when the LP relaxation did not solve, which a model
	 * would otherwise be unable to tell from a genuine zero -- and that is
	 * exactly where LSmear hands over to its fallback.
	 *
	 * Features 15..29 are *comparative*: choosing a branching variable is a
	 * comparison between the candidates of the node, and a per-variable model
	 * that only ever sees one variable at a time cannot express that. For each
	 * of the five base quantities
	 *
	 *     lsmear (13), smear_sum (10), smear_sum_rel (12), dual (9), rel_diam (5)
	 *
	 * three values are derived over the admissible candidates C of the node:
	 *
	 *   15,18,21,24,27  rank    |{i in C : v_i < v_j}| / (|C|-1), or 0.5 if |C|==1
	 *   16,19,22,25,28  relmax  v_j / max_{i in C} |v_i|   (0 if that max is 0)
	 *   17,20,23,26,29  relsum  v_j / sum_{i in C} |v_i|   (0 if that sum is 0)
	 *
	 * in that order (rank, relmax, relsum for lsmear, then for smear_sum, and
	 * so on). Variables outside C keep 0: they are never scored anyway.
	 *
	 * Monotone per-feature transforms are deliberately *not* provided: trees are
	 * invariant to them and a linear model folds them into its weights. The
	 * comparative features above are the ones that cannot be recovered from a
	 * single variable's vector.
	 *
	 * Non-finite values are sanitized exactly as `encode()` does before the
	 * model sees them: NaN becomes 0 and infinities are clipped to
	 * +/-#CLIP.
	 */
	static const int NB_FEATURES = 30;

	/** \brief Clipping bound applied to infinite features. */
	static const double CLIP;

	/** \brief Index of the first comparative feature. */
	static const int FIRST_DERIVED = 15;

	/** \brief Feature indices the comparative features are derived from. */
	static const int DERIVED_FROM[5];

	typedef enum { LINEAR, GBDT } Type;

	/**
	 * \brief Load a model from a file.
	 *
	 * \throw std::runtime_error if the file is missing or malformed, or if its
	 *        "nfeat" does not match #NB_FEATURES.
	 */
	explicit MLModel(const std::string& filename);

	/**
	 * \brief Score one variable from its feature vector.
	 *
	 * \param phi - #NB_FEATURES values, already sanitized.
	 */
	double score(const double* phi) const;

	/** \brief The model family. */
	Type type() const;

	/** \brief Number of trees (0 for a linear model). */
	int nb_trees() const;

	/** \brief A one-line description, for reporting. */
	std::string description() const;

	/** \brief Replace NaN by 0 and infinities by +/-#CLIP. */
	static double sanitize(double x);

protected:
	/* A regression tree, stored as a flat array of nodes. */
	struct Node {
		bool leaf;
		int feat;      // internal nodes only
		double thr;    // internal nodes: go left iff phi[feat] <= thr
		int left;      // internal nodes: index of the left child
		int right;     // internal nodes: index of the right child
		double value;  // leaves only
	};

	double eval_tree(const std::vector<Node>& tree, const double* phi) const;

	Type _type;
	int _nfeat;
	double _bias;                        // linear: bias; gbdt: base score
	std::vector<double> _w;              // linear only
	std::vector<std::vector<Node> > _trees;  // gbdt only
	std::string _filename;
};

inline MLModel::Type MLModel::type() const { return _type; }

inline int MLModel::nb_trees() const { return (int) _trees.size(); }

} // end namespace ibex

#endif // __IBEX_ML_MODEL_H__
