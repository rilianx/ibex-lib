//============================================================================
//                                  I B E X
// File        : ibex_MLModel.cpp
// Author      : Ignacio Araya
// License     : See the LICENSE file
// Created     : 2026
//============================================================================

#include "ibex_MLModel.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace ibex {

const double MLModel::CLIP = 1e12;
const int MLModel::NB_FEATURES;
const int MLModel::FIRST_DERIVED;
// lsmear, smear_sum, smear_sum_rel, dual, rel_diam
const int MLModel::DERIVED_FROM[5] = { 13, 10, 12, 9, 5 };

double MLModel::sanitize(double x) {
	if (std::isnan(x)) return 0.0;
	if (x > CLIP) return CLIP;
	if (x < -CLIP) return -CLIP;
	return x;
}

namespace {

void fail(const string& file, int line, const string& msg) {
	ostringstream o;
	o << file << ":" << line << ": " << msg;
	throw runtime_error(o.str());
}

} // anonymous namespace

MLModel::MLModel(const string& filename) :
		_type(LINEAR), _nfeat(NB_FEATURES), _bias(0.0), _filename(filename) {

	ifstream in(filename.c_str());
	if (!in.is_open())
		throw runtime_error("cannot open the model file '" + filename + "'");

	string line;
	int lineno = 0;
	bool header = false, type_seen = false;

	while (getline(in, line)) {
		lineno++;

		// strip comments and surrounding blanks
		size_t h = line.find('#');
		if (h != string::npos) line = line.substr(0, h);
		size_t a = line.find_first_not_of(" \t\r\n");
		if (a == string::npos) continue;
		size_t b = line.find_last_not_of(" \t\r\n");
		line = line.substr(a, b - a + 1);

		istringstream ls(line);
		string key;
		ls >> key;

		if (!header) {
			int version = 0;
			if (key != "ibexml-model" || !(ls >> version) || version != 1)
				fail(filename, lineno, "expected the header 'ibexml-model 1'");
			header = true;
			continue;
		}

		if (key == "type") {
			string t;
			ls >> t;
			if (t == "linear") _type = LINEAR;
			else if (t == "gbdt") _type = GBDT;
			else fail(filename, lineno, "unknown model type '" + t + "'");
			type_seen = true;

		} else if (key == "nfeat") {
			if (!(ls >> _nfeat)) fail(filename, lineno, "expected an integer after 'nfeat'");
			if (_nfeat != NB_FEATURES) {
				ostringstream o;
				o << "the model expects " << _nfeat << " features but this build provides "
				  << NB_FEATURES;
				fail(filename, lineno, o.str());
			}

		} else if (key == "bias" || key == "base") {
			if (!(ls >> _bias)) fail(filename, lineno, "expected a real after '" + key + "'");

		} else if (key == "w") {
			_w.clear();
			double v;
			while (ls >> v) _w.push_back(v);
			if ((int) _w.size() != _nfeat) {
				ostringstream o;
				o << "'w' has " << _w.size() << " values, expected " << _nfeat;
				fail(filename, lineno, o.str());
			}

		} else if (key == "tree") {
			int n = 0;
			if (!(ls >> n) || n <= 0)
				fail(filename, lineno, "expected a positive node count after 'tree'");

			vector<Node> tree;
			tree.reserve(n);
			for (int k = 0; k < n; k++) {
				if (!getline(in, line)) fail(filename, lineno, "truncated tree");
				lineno++;
				size_t hh = line.find('#');
				if (hh != string::npos) line = line.substr(0, hh);
				istringstream ns(line);
				string kind;
				if (!(ns >> kind)) { k--; continue; }   // blank line inside a tree

				Node node;
				node.leaf = false; node.feat = 0; node.thr = 0;
				node.left = -1; node.right = -1; node.value = 0;

				if (kind == "L") {
					node.leaf = true;
					if (!(ns >> node.value)) fail(filename, lineno, "expected a leaf value");
				} else if (kind == "S") {
					if (!(ns >> node.feat >> node.thr >> node.left >> node.right))
						fail(filename, lineno, "expected 'S <feat> <thr> <left> <right>'");
					if (node.feat < 0 || node.feat >= _nfeat)
						fail(filename, lineno, "feature index out of range");
				} else
					fail(filename, lineno, "expected 'S' or 'L', got '" + kind + "'");

				tree.push_back(node);
			}

			// every child index must be inside this tree
			for (size_t k = 0; k < tree.size(); k++) {
				if (tree[k].leaf) continue;
				if (tree[k].left < 0 || tree[k].left >= (int) tree.size() ||
					tree[k].right < 0 || tree[k].right >= (int) tree.size())
					fail(filename, lineno, "child index out of range in a tree");
			}
			_trees.push_back(tree);

		} else
			fail(filename, lineno, "unknown keyword '" + key + "'");
	}

	if (!header) throw runtime_error("'" + filename + "' is empty");
	if (!type_seen) throw runtime_error("'" + filename + "' has no 'type' line");
	if (_type == LINEAR && (int) _w.size() != _nfeat)
		throw runtime_error("'" + filename + "' is linear but has no 'w' line");
	if (_type == GBDT && _trees.empty())
		throw runtime_error("'" + filename + "' is a gbdt but has no tree");
}

double MLModel::eval_tree(const vector<Node>& tree, const double* phi) const {
	int i = 0;
	// the depth is bounded by the number of nodes, so a malformed cycle cannot hang
	for (size_t step = 0; step <= tree.size(); step++) {
		const Node& n = tree[i];
		if (n.leaf) return n.value;
		i = (phi[n.feat] <= n.thr) ? n.left : n.right;
	}
	return 0.0;
}

double MLModel::score(const double* phi) const {
	double s = _bias;
	if (_type == LINEAR) {
		for (int j = 0; j < _nfeat; j++) s += _w[j] * phi[j];
	} else {
		for (size_t t = 0; t < _trees.size(); t++) s += eval_tree(_trees[t], phi);
	}
	return s;
}

string MLModel::description() const {
	ostringstream o;
	o << _filename << " (" << (_type == LINEAR ? "linear" : "gbdt");
	if (_type == GBDT) o << ", " << _trees.size() << " trees";
	o << ", " << _nfeat << " features)";
	return o.str();
}

} // end namespace ibex
