//============================================================================
//                                  I B E X
//
//                             ***************
//                                IbexOpt-ML
//                             ***************
//
// Node-level instrumentation of IbexOpt, meant to be driven from Python in
// order to learn a bisection heuristic.
//
// Four modes:
//
//   --collect : run the real (LSmear- or model-driven) search and write one
//               JSON training sample per line.
//   --run     : run the real search but ask the client, over stdin/stdout, for
//               every branching decision. The node selection stays inside, so
//               node counts are comparable with ibexopt.
//   --solve   : run the search to completion under the current rule and report
//               nodes and wall-clock time.
//   (default) : a server reading one JSON command per line on stdin.
//
// Author      : Ignacio Araya
// License     : See the LICENSE file
//============================================================================

#include "ibex.h"
#include "parse_args.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

#ifdef __IBEX_WITH_IPOPT__
namespace ibex { extern double ipopt_diam; extern bool ipopt_defer; extern double ipopt_expansion; extern double ipopt_expansion_first; extern double ipopt_start_max; extern double ipopt_max_time; extern double ipopt_call_limit; extern long ipopt_adaptive; extern bool expr_sqr_rewrite; extern double ipopt_corr_time; extern int ipopt_multistart; extern bool cvx_stats; extern bool cvx_close_on; extern long cvx_backoff_max; extern int cvx_psd_test; extern long cvx_center_indef; extern long optimizer_extra_lb_raised, optimizer_extra_lb_pruned; extern long cvx_nodes, cvx_with_dual, cvx_certified, cvx_close, cvx_better_lb; extern double ipopt_guided; extern bool ipopt_guided_nodes; extern double ipopt_total_time, ipopt_corr_total_time, cvx_total_time; extern int ipopt_guided_patience; extern bool ipopt_guided_local; extern int ipopt_guided_calls; extern int ipopt_guided_improved; extern double ipopt_guided_time; extern int ipopt_ms_improved; extern int ipopt_ms_seed; }   // ibex_LoupFinderIpopt.cpp: Ipopt only on boxes this narrow
#endif
using namespace ibex;

namespace {

/* Answer "{"ok":false,"error":...}" on stdout. */
void fail(const string& msg) {
	JsonOut out(cout);
	out.obj();
	out.kv("ok", false);
	out.kv("error", msg);
	out.end_obj();
	cout << "\n" << flush;
}

/* Read the mandatory "box" argument of a command. */
IntervalVector read_box(MLNodeServer& server, const JsonValue& cmd) {
	const JsonValue* b = cmd.get("box");
	if (b==NULL) throw JsonError("missing 'box'");
	return b->as_box(server.nb_ext_var());
}

/* Read the sampling/diving parameters a command may override. */
MLNodeServer::SampleParams read_sample_params(const JsonValue& cmd,
		const MLNodeServer::SampleParams& def) {
	MLNodeServer::SampleParams p = def;
	p.budget            = cmd.get_long  ("budget",        def.budget);
	p.budget_start      = cmd.get_long  ("budget_start",  def.budget_start);
	p.prune             = cmd.get_bool  ("prune",         def.prune);
	p.max_depth         = cmd.get_int   ("max_depth",     def.max_depth);
	p.include_goal      = cmd.get_bool  ("include_goal",  def.include_goal);
	p.topk              = cmd.get_int   ("topk",          def.topk);
	p.depth             = cmd.get_int   ("depth",         -1);
	p.last_bisected_var = cmd.get_int   ("last_bisected_var", -1);
	return p;
}

/*
 * Run one command. Returns false when the client asked to stop.
 */
bool run_command(MLNodeServer& server, const JsonValue& cmd,
		const MLNodeServer::SampleParams& defaults) {

	string what = cmd.get_string("cmd", "");

	if (what=="quit" || what=="exit") return false;

	JsonOut out(cout);

	if (what=="info") {
		out.obj();
		out.kv("ok", true);
		out.key("info");
		server.write_info(out);
		out.end_obj();

	} else if (what=="reset") {
		const JsonValue* b = cmd.get("box");
		// no box given: restart from the very beginning, not from the already
		// contracted root
		IntervalVector box = (b==NULL || b->is_null()) ? server.initial_box() : b->as_box();
		server.reset(box, cmd.get_double("loup", POS_INFINITY));
		out.obj();
		out.kv("ok", true);
		out.kv("root", server.root_box());
		out.kv("loup", server.current_loup());
		out.end_obj();

	} else if (what=="contract") {
		MLNodeServer::ContractResult r =
				server.contract(read_box(server, cmd), cmd.get_bool("keep_loup", true));
		out.obj();
		out.kv("ok", true);
		out.kv("status", r.status==MLNodeServer::OPEN ? "open" :
				(r.status==MLNodeServer::EPS ? "eps" : "pruned"));
		out.kv("box", r.box);
		out.kv("eps_lb", r.eps_lb);
		out.kv("loup_before", r.loup_before);
		out.kv("loup_after", r.loup_after);
		out.kv("time", r.time);
		out.end_obj();

	} else if (what=="features") {
		out.obj();
		out.kv("ok", true);
		out.key("node");
		server.write_features(out, read_box(server, cmd),
				cmd.get_int("depth", -1), cmd.get_int("last_bisected_var", -1));
		out.end_obj();

	} else if (what=="lsmear" || what=="choose_var") {
		IntervalVector box = read_box(server, cmd);
		out.obj();
		out.kv("ok", true);
		try {
			BisectionPoint bp = server.choose_var(box);
			out.kv("var", (int) bp.var);
			out.kv("pos", bp.pos);
			out.kv("rel_pos", bp.rel_pos);
		} catch (NoBisectableVariableException&) {
			out.knull("var");
			out.knull("pos");
			out.knull("rel_pos");
		}
		out.end_obj();

	} else if (what=="model_var") {
		IntervalVector box = read_box(server, cmd);
		int v = server.model_var(box, cmd.get_bool("include_goal", defaults.include_goal),
				cmd.get_int("topk", defaults.topk));
		out.obj();
		out.kv("ok", true);
		if (v>=0) out.kv("var", v); else out.knull("var");
		out.end_obj();

	} else if (what=="bisect") {
		IntervalVector box = read_box(server, cmd);
		int var = cmd.get_int("var", -1);
		double pos = cmd.get_double("pos", std::numeric_limits<double>::quiet_NaN());
		bool rel = cmd.get_bool("rel_pos", true);

		if (var<0) { // let the bisector decide
			BisectionPoint bp = server.choose_var(box);   // may throw
			var = (int) bp.var;
			pos = bp.pos;
			rel = bp.rel_pos;
		}
		pair<IntervalVector,IntervalVector> p = (rel || std::isnan(pos)) ?
				server.split(box, var) : server.split(box, var, pos);

		out.obj();
		out.kv("ok", true);
		out.kv("var", var);
		out.kv("left", p.first);
		out.kv("right", p.second);
		out.end_obj();

	} else if (what=="dive") {
		const JsonValue* v = cmd.get("var");
		if (v==NULL || v->is_null()) throw JsonError("missing 'var'");
		MLNodeServer::DiveResult d = server.dive(read_box(server, cmd), v->as_int(),
				cmd.get_long("budget", defaults.budget),
				cmd.get_int("max_depth", defaults.max_depth));
		out.obj();
		out.kv("ok", true);
		out.key("dive");
		server.write_dive(out, d);
		out.end_obj();

	} else if (what=="probe") {
		vector<pair<int,MLNodeServer::ProbeResult> > pr = server.hc4_probes(read_box(server, cmd),
				cmd.get_int("dims", 1), cmd.get_int("parts", 4),
				cmd.get_bool("include_goal", defaults.include_goal), cmd.get_int("topk", defaults.topk),
				cmd.get_string("ctc", "hc4")=="lp" ? MLNodeServer::PROBE_LP :
				(cmd.get_string("ctc", "hc4")=="full" ? MLNodeServer::PROBE_FULL :
				(cmd.get_string("ctc", "hc4")=="proc" ? MLNodeServer::PROBE_PROC : MLNodeServer::PROBE_HC4)));
		out.obj();
		out.kv("ok", true);
		out.key("probes").arr();
		for (size_t k=0; k<pr.size(); k++) {
			out.obj();
			out.kv("var", pr[k].first);
			out.kv("empties", pr[k].second.empties);
			out.kv("logvol", pr[k].second.logvol);
			out.kv("time", pr[k].second.time);
			out.end_obj();
		}
		out.end_arr();
		out.end_obj();

	} else if (what=="candidates") {
		vector<int> c = server.candidates(read_box(server, cmd),
				cmd.get_bool("include_goal", defaults.include_goal),
				cmd.get_int("topk", defaults.topk));
		out.obj();
		out.kv("ok", true);
		out.key("vars").arr();
		for (size_t k=0; k<c.size(); k++) out.val(c[k]);
		out.end_arr();
		out.end_obj();

	} else if (what=="sample") {
		MLNodeServer::SampleParams p = read_sample_params(cmd, defaults);
		out.obj();
		out.kv("ok", true);
		out.key("sample");
		server.write_sample(out, read_box(server, cmd), p);
		out.end_obj();

	} else if (what=="features_vector") {
		vector<double> phi;
		server.var_features(read_box(server, cmd), phi);
		int F = MLModel::NB_FEATURES;
		out.obj();
		out.kv("ok", true);
		out.kv("nb_features", F);
		out.key("phi").arr();
		for (size_t j=0; j*F<phi.size(); j++) {
			out.arr();
			for (int k=0; k<F; k++) out.val(phi[j*F+k]);
			out.end_arr();
		}
		out.end_arr();
		out.end_obj();

	} else if (what=="set_loup") {
		server.set_loup(cmd.get_double("loup", POS_INFINITY));
		out.obj();
		out.kv("ok", true);
		out.kv("loup", server.current_loup());
		out.end_obj();

	} else if (what=="state") {
		out.obj();
		out.kv("ok", true);
		out.kv("loup", server.current_loup());
		out.kv("uplo", server.get_uplo());
		out.kv("nb_cells", server.get_nb_cells());
		out.end_obj();

	} else {
		throw JsonError("unknown command '" + what + "'");
	}

	cout << "\n" << flush;
	return true;
}

} // anonymous namespace

int main(int argc, char** argv) {

#ifdef __IBEX_NO_LP_SOLVER__
	ibex_error("ibexopt-ml requires a LP Solver (use -DLP_LIB with cmake)");
	exit(1);
#endif

	args::ArgumentParser parser("********* IbexOpt-ML (branching-data extraction) *********.",
			"Extract node features and dive labels from a Minibex file, run the search "
			"under a learned rule, or serve one node at a time over stdin/stdout.");

	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	args::Flag version(parser, "version", "Display the version of Ibex.", {'v',"version"});

	// --- optimizer settings (same meaning as in ibexopt) ---
	args::ValueFlag<double> rel_eps_f(parser, "float", "Relative precision on the objective.", {'r', "rel-eps-f"});
	args::ValueFlag<double> abs_eps_f(parser, "float", "Absolute precision on the objective.", {'a', "abs-eps-f"});
	args::ValueFlag<double> eps_h(parser, "float", "Equality relaxation value.", {"eps-h"});
	args::ValueFlag<double> eps_x_arg(parser, "float", "Precision on the variables.", {"eps-x"});
	args::ValueFlag<double> random_seed(parser, "float", "Random seed (reproducibility).", {"random-seed"});
	args::ValueFlag<int>    simpl_level(parser, "int", "Expression simplification level (0..3).", {"simpl"});
	args::ValueFlag<double> initial_loup(parser, "float", "Initial \"loup\" (a priori known upper bound).", {"initial-loup"});
	args::Flag rigor(parser, "rigor", "Activate rigor mode.", {"rigor"});
	args::Flag kkt(parser, "kkt", "Activate the Kuhn-Tucker contractor.", {"kkt"});

	// --- modes ---
	args::Flag collect(parser, "collect", "Run the search and write samples (one JSON object per line).", {"collect"});
	args::Flag run(parser, "run", "Run the search, asking the client for every branching variable.", {"run"});
	args::Flag solve(parser, "solve", "Run the search to completion and report nodes and time.", {"solve"});
	args::ValueFlag<string> output_file(parser, "filename", "Where to write the samples in --collect mode (default: stdout).", {'o',"output"});

	// --- the learned rule ---
	args::ValueFlag<string> model_file(parser, "filename", "Branch with this model instead of the bisector (see MLModel for the format).", {"model"});
	args::ValueFlag<double> grasp_alpha(parser, "float", "With --bisector lsmear-grasp: choose among the candidates >= alpha * the best LSmear impact. Default: 0.8.", {"grasp-alpha"});
	args::ValueFlag<int> tabu_tenure(parser, "int", "With --bisector lsmear-tabu: a variable that captures LSmear (it wants the parent's variable again) is tabu for N levels along the branch. With lsmear-recent: the variables of the last N ancestors are tabu. Default: 2.", {"tabu-tenure"});
	args::ValueFlag<double> sb_ratio_arg(parser, "float", "With --oracle-score sb: leave the bisector's choice only if strong branching prunes more children, or as many with at most this fraction of the open volume. Default: 0 (always strong branching).", {"sb-ratio"});
	args::Flag sb_vol_only(parser, "sb-vol-only", "With --sb-ratio: judge by the open volume alone, also when strong branching prunes more children.", {"sb-vol-only"});
	args::ValueFlag<int> probe_dims(parser, "int", "With --oracle-score hc4: dimensions reduced per probe (the candidate plus N-1 partners from LSmear's ranking). Default: 1.", {"probe-dims"});
	args::ValueFlag<int> probe_parts(parser, "int", "With --oracle-score hc4: slices per reduced dimension. Default: 4.", {"probe-parts"});
	args::ValueFlag<string> probe_ctc(parser, "hc4|lp|full|proc", "With --oracle-score hc4: the probe's contractor. hc4 (default); lp: one HC4 pass then one polytope hull (XTaylor+affine); full: the search's contractor; proc: contraction and upper bounding, as a step of a dive.", {"probe-ctc"});
	args::ValueFlag<int> ipopt_first(parser, "int", "With --loup ipoptxn: also call Ipopt at each of the first N loup-finder calls (a multi-start near the root). Default: 0.", {"ipopt-first"});
	args::ValueFlag<double> dive_loup(parser, "float", "Dives (sampling, oracle) prune with this upper bound when it is lower than the search's: an optimum known in advance, so that the labels measure the branching alone. The search itself keeps its own loup.", {"dive-loup"});
	args::Flag sb_need_pruned(parser, "sb-need-pruned", "With --sb-ratio: deviate from the bisector only when the probe prunes (or empties) strictly more than its choice; a smaller volume alone is not enough.", {"sb-need-pruned"});
	args::ValueFlag<int> pc_rel(parser, "int", "With --oracle-score pc: probe a variable until it has N observations (reliability). Default: 4.", {"pc-rel"});
	args::Flag pc_light(parser, "pc-light", "With --oracle-score pc: probe with the search's contractor alone (no upper bounding) instead of fully processing both children.", {"pc-light"});
	args::ValueFlag<int> acid_fixed(parser, "int", "Non-adaptive ACID: 3BCID on N variables at every node (capped at all of them: a large N means all, as in ACID's first tuning phase), no tuning. For comparisons between searches that would tune ACID differently.", {"acid-fixed"});
	args::ValueFlag<double> late_loup(parser, "float", "Give the search this upper bound right after the root is processed (unlike --initial-loup, the root's contraction is unchanged).", {"late-loup"});
	args::ValueFlag<double> ipopt_diam_arg(parser, "float", "With --loup ipopt*: Ipopt is only called on boxes whose widest variable is at most this wide (negative: no threshold, unbounded boxes included). Default: 1e8.", {"ipopt-diam"});
	args::ValueFlag<double> ipopt_expansion_arg(parser, "float", "With --loup ipopt*: half-width (relative to |x| when |x|>1) of the box around Ipopt's point in which a short nested search looks for a certified feasible point. Default: 1e-6.", {"ipopt-expansion"});
	args::ValueFlag<double> ipopt_expansion_first_arg(parser, "float", "With --ipopt-expansion: first try the default 1e-6 box with this CPU timeout (seconds), and widen only if it certifies nothing.", {"ipopt-expansion-first"});
	args::ValueFlag<double> ipopt_start_max_arg(parser, "float", "With --ipopt-diam negative: Ipopt is not called on a box whose starting point has a coordinate larger than this in absolute value. Default: 1e8.", {"ipopt-start-max"});
	args::ValueFlag<double> ipopt_max_time_arg(parser, "float", "With --loup ipopt*: CPU limit of one Ipopt call in seconds (Ipopt's max_cpu_time). Default: none.", {"ipopt-max-time"});
	args::Flag bound_trace_arg(parser, "bound-trace", "With --solve: add \"bounds\": [[time, nodes, loup, uplo], ...] (start, each new loup, uplo moves at geometrically spaced times, end).", {"bound-trace"});
	args::Flag sqr_rewrite_arg(parser, "sqr-rewrite", "Simplify e*e into sqr(e) when both factors are the same expression (exact interval square).", {"sqr-rewrite"});
	args::ValueFlag<long> ipopt_adaptive_arg(parser, "int", "With --loup ipopt*: adaptive Ipopt schedule; the gap between scheduled calls resets to 1 after a call that improves the loup and doubles (up to this value) after one that does not. 0: fixed schedule (default).", {"ipopt-adaptive"});
	args::ValueFlag<double> ipopt_corr_time_arg(parser, "float", "With --ipopt-expansion: CPU timeout (s) of the correction search in the wide box. Default: 1.", {"ipopt-corr-time"});
	args::ValueFlag<double> ipopt_call_limit_arg(parser, "float", "With --loup ipopt*: CPU seconds after which one Ipopt call is made to give up (its evaluations report an error). <= 0: none. Default: 5.", {"ipopt-call-limit"});
	args::ValueFlag<double> ipopt_guided_arg(parser, "float", "With --loup ipopt*: guided restarts, an extra Ipopt run from the current node's box (of 10 random points, the farthest from previous Ipopt solutions), while these runs use less than this fraction of the CPU time.", {"ipopt-guided"});
	args::ValueFlag<int> ipopt_guided_patience_arg(parser, "int", "With --ipopt-guided: after this many guided runs in a row without improving the loup, stop them until the loup improves. Default: 0 (never stop).", {"ipopt-guided-patience"});
	args::Flag ipopt_guided_nodes_arg(parser, "ipopt-guided-nodes", "With --ipopt-guided: budget in loup-finder calls (guided runs <= fraction x nodes) instead of CPU time, so the search is deterministic.", {"ipopt-guided-nodes"});
	args::Flag ipopt_guided_local_arg(parser, "ipopt-guided-local", "With --ipopt-guided: Ipopt bounded by the node's box instead of the root box.", {"ipopt-guided-local"});
	args::Flag cvx_stats_arg(parser, "cvx-stats", "Measure, at each node, whether the Lagrangian of the last Ipopt run is certifiably convex on the box and whether its tangent bound would close it (counters in the JSON).", {"cvx-stats"});
	args::Flag cvx_close_arg(parser, "cvx-close", "Raise each node's lower bound with the tangent bound of the Lagrangian of the last Ipopt run when it is certified convex on the box (interval Hessian, Gershgorin); boxes above the upper bound are pruned.", {"cvx-close"});
	args::ValueFlag<int> cvx_backoff_arg(parser, "int", "With --cvx-close: after a node not certified convex, check 2x later, up to every N nodes; a certified node resets it. Default: 0 (every node).", {"cvx-backoff"});
	args::ValueFlag<int> cvx_test_arg(parser, "int", "PSD test of the interval Hessian for --cvx-close/--cvx-stats: 0 Gershgorin (default), 1 Gershgorin + comparison matrix + midpoint-radius (interval Cholesky).", {"cvx-test"});
	args::ValueFlag<int> ipopt_multistart_arg(parser, "int", "With --loup ipopt*: at the first Ipopt call, this many extra Ipopt runs from random points of the root box.", {"ipopt-multistart"});
	args::ValueFlag<int> ipopt_ms_seed_arg(parser, "int", "Seed of --ipopt-multistart's starting points (default 1).", {"ipopt-ms-seed"});
	args::Flag ipopt_defer_arg(parser, "ipopt-defer", "With --loup ipopt*: a scheduled Ipopt call that falls on a box wider than --ipopt-diam waits for the next box narrow enough instead of being skipped.", {"ipopt-defer"});
	args::ValueFlag<long> guard_horizon(parser, "int", "With --bisector lsmear-guard: only switch to round-robin within the first N decisions. Default: 0 (no horizon).", {"guard-horizon"});
	args::ValueFlag<string> bisector_arg(parser, "name", "Bisector to use. One of: "
			+ MLOptimizerConfig::bisector_names() + ". Default: lsmear (what ibexopt uses).", {"bisector"});
	args::ValueFlag<string> relax_arg(parser, "name", "Linear relaxation behind the X-Newton step. One of: "
			+ MLOptimizerConfig::relaxation_names() + ". Default: xtaylor (what ibexopt uses).", {"relax"});
	args::ValueFlag<double> ratio_arg(parser, "float", "Where a bisected domain is cut. "
			"Default: 0.5 (the middle, what ibexopt uses). Ibex's own Bsc default is "
			"0.45, which is what the hand-assembled strategies in the examples get: "
			"they never pass a ratio. It reaches only the rules that do not compute "
			"their own bisection point.", {"bisect-ratio"});
	args::ValueFlag<string> loup_arg(parser, "name", "Upper bounding. One of: "
			+ MLOptimizerConfig::loup_names() + ". Default: default (what ibexopt uses). "
			"The ipopt* variants change the incumbent, so a run with one is not comparable "
			"with a run without.", {"loup"});
	args::ValueFlag<int> ipopt_freq(parser, "int", "Ipopt is called every N calls to the loup finder, "
			"plus the 10th, 20th and 50th, plus whenever another finder improves the incumbent. Default: 100.", {"ipopt-freq"});
	args::Flag ipopt_qp(parser, "ipopt-qp", "The problem is a QP (quadratic objective, linear constraints): "
			"the Lagrangian hessian is then computed once instead of at every call.", {"ipopt-qp"});

	// --- diving / sampling ---
	args::ValueFlag<long>   budget(parser, "int", "Ceiling on the dive node budget. Default: 200.", {"budget"});
	args::ValueFlag<long>   budget_start(parser, "int", "First dive budget tried; doubled while every candidate stays censored. 0 disables. Default: 25.", {"budget-start"});
	args::Flag oracle_arg(parser, "oracle", "Branch with the one-step oracle: dive on every candidate at every node (--budget, --budget-start, --no-prune apply) and bisect on the one that closes in the fewest nodes. Expensive.", {"oracle"});
	args::ValueFlag<string> oracle_score(parser, "nodes|depth", "With --oracle: rank the dives by their size (nodes, the default), lazily by how deep they went (depth: one fixed --budget for all, no pruning), or by strong branching (sb: contract both children only; more pruned children, then smaller total volume).", {"oracle-score"});
	args::ValueFlag<int> oracle_depth_arg(parser, "int", "With --oracle: use it only at nodes up to this depth; deeper, the bisector decides (0: everywhere).", {"oracle-max-depth"});
	args::Flag no_prune(parser, "no-prune", "Give every candidate the full budget instead of the best result so far.", {"no-prune"});
	args::ValueFlag<int>    max_depth(parser, "int", "Dive depth limit (0: none). Default: 0.", {"max-depth"});
	args::ValueFlag<int>    topk(parser, "int", "Evaluate only the k most promising candidates (0: all). Default: 0.", {"topk"});
	args::Flag no_goal(parser, "no-goal", "Do not evaluate the objective variable as a bisection candidate.", {"no-goal"});

	// --- search control ---
	args::ValueFlag<double> sample_prob(parser, "float", "Probability of sampling a visited node in --collect mode. Default: 0.1.", {"sample-prob"});
	args::Flag stop_at_max_samples(parser, "stop-at-max-samples", "In --collect mode, end the search once --max-samples are written (nothing after the last sample is used). Size --sample-prob so the cap is reached near the end of the tree, or the samples crowd its top.", {"stop-at-max-samples"});
	args::ValueFlag<long>   max_samples(parser, "int", "Stop sampling after that many (the search keeps going; 0: no limit). Default: 1000.", {"max-samples"});
	args::ValueFlag<long>   max_nodes(parser, "int", "Stop the search after that many nodes (0: no limit).", {"max-nodes"});
	args::ValueFlag<double> timeout(parser, "float", "Time limit of the search, in seconds.", {'t',"timeout"});
	args::Flag no_features(parser, "no-features", "In --run mode, do not send the node features (decisions only).", {"no-features"});
	args::Flag progress(parser, "progress", "Report progress on stderr.", {"progress"});

	args::Positional<string> filename(parser, "filename", "The name of the MINIBEX file.");

	try {
		parser.ParseCLI(argc, argv);
	} catch (args::Help&) {
		cout << parser;
		return 0;
	} catch (args::ParseError& e) {
		cerr << e.what() << endl << parser;
		return 1;
	} catch (args::ValidationError& e) {
		cerr << e.what() << endl << parser;
		return 1;
	}

	if (version) {
		cout << "IbexOpt-ML Release " << _IBEX_RELEASE_ << endl;
		return 0;
	}

	if (filename.Get()=="") {
		cerr << "no input file (try ibexopt-ml --help)" << endl;
		return 1;
	}

	System* sys = NULL;
	MLNodeServer* server = NULL;
	MLModel* model = NULL;
	int rc = 0;

	try {
		if (sqr_rewrite_arg) ibex::expr_sqr_rewrite = true;   // before parsing (simplification)
		sys = new System(filename.Get().c_str(),
				simpl_level ? simpl_level.Get() : ExprNode::default_simpl_level);

		if (!sys->goal)
			ibex_error("input file has no goal (it is not an optimization problem).");

		double seed = random_seed ? random_seed.Get() : DefaultOptimizerConfig::default_random_seed;

		Vector eps_x(sys->nb_var, eps_x_arg ? eps_x_arg.Get() : OptimizerConfig::default_eps_x);

		MLOptimizerConfig::Bisector bisector = MLOptimizerConfig::BSC_LSMEAR_MG;
		if (bisector_arg && !MLOptimizerConfig::parse_bisector(bisector_arg.Get(), bisector)) {
			cerr << "unknown bisector '" << bisector_arg.Get() << "'. One of: "
			     << MLOptimizerConfig::bisector_names() << endl;
			delete sys;
			return 1;
		}

		MLOptimizerConfig::Relaxation relaxation = MLOptimizerConfig::RELAX_XTAYLOR;
		if (relax_arg && !MLOptimizerConfig::parse_relaxation(relax_arg.Get(), relaxation)) {
			cerr << "unknown relaxation '" << relax_arg.Get() << "'. One of: "
			     << MLOptimizerConfig::relaxation_names() << endl;
			delete sys;
			return 1;
		}

		double bisect_ratio = ratio_arg ? ratio_arg.Get()
				: DefaultOptimizerConfig::default_bisect_ratio;

		MLOptimizerConfig::LoupFinderKind loup = MLOptimizerConfig::LOUP_DEFAULT;
		if (loup_arg && !MLOptimizerConfig::parse_loup(loup_arg.Get(), loup)) {
			cerr << "unknown upper bounding '" << loup_arg.Get() << "'. One of: "
			     << MLOptimizerConfig::loup_names() << endl;
			delete sys;
			return 1;
		}

		if (grasp_alpha) MLOptimizerConfig::grasp_alpha = grasp_alpha.Get();
		if (tabu_tenure) MLOptimizerConfig::tabu_tenure = tabu_tenure.Get();
#ifdef __IBEX_WITH_IPOPT__
		if (ipopt_diam_arg) ibex::ipopt_diam = ipopt_diam_arg.Get();
		if (ipopt_expansion_arg) ibex::ipopt_expansion = ipopt_expansion_arg.Get();
		if (ipopt_call_limit_arg) ibex::ipopt_call_limit = ipopt_call_limit_arg.Get();
		if (ipopt_adaptive_arg) ibex::ipopt_adaptive = ipopt_adaptive_arg.Get();
		if (ipopt_corr_time_arg) ibex::ipopt_corr_time = ipopt_corr_time_arg.Get();
		if (ipopt_guided_arg) ibex::ipopt_guided = ipopt_guided_arg.Get();
		if (ipopt_guided_patience_arg) ibex::ipopt_guided_patience = ipopt_guided_patience_arg.Get();
		if (ipopt_guided_nodes_arg) ibex::ipopt_guided_nodes = true;
		if (ipopt_guided_local_arg) ibex::ipopt_guided_local = true;
		if (cvx_stats_arg) ibex::cvx_stats = true;
		if (cvx_close_arg) ibex::cvx_close_on = true;
		if (cvx_test_arg) ibex::cvx_psd_test = cvx_test_arg.Get();
		if (cvx_backoff_arg) ibex::cvx_backoff_max = cvx_backoff_arg.Get();
		if (ipopt_multistart_arg) ibex::ipopt_multistart = ipopt_multistart_arg.Get();
		if (ipopt_ms_seed_arg) ibex::ipopt_ms_seed = ipopt_ms_seed_arg.Get();
		if (ipopt_max_time_arg) ibex::ipopt_max_time = ipopt_max_time_arg.Get();
		if (ipopt_start_max_arg) ibex::ipopt_start_max = ipopt_start_max_arg.Get();
		if (ipopt_expansion_first_arg) ibex::ipopt_expansion_first = ipopt_expansion_first_arg.Get();
		if (ipopt_defer_arg) ibex::ipopt_defer = true;
#endif
		if (ipopt_first) MLOptimizerConfig::ipopt_first = ipopt_first.Get();
		MLOptimizerConfig::grasp_seed = (unsigned int) seed;

		server = new MLNodeServer(*sys,
				rel_eps_f ? rel_eps_f.Get() : OptimizerConfig::default_rel_eps_f,
				abs_eps_f ? abs_eps_f.Get() : OptimizerConfig::default_abs_eps_f,
				eps_h ? eps_h.Get() : NormalizedSystem::default_eps_h,
				rigor.Get(),
				DefaultOptimizerConfig::default_inHC4,
				kkt.Get(),
				seed,
				eps_x,
				bisector,
				relaxation,
				loup,
				ipopt_freq ? ipopt_freq.Get() : 100,
				ipopt_qp.Get(),
				bisect_ratio);

		if (guard_horizon && !server->set_guard_horizon(guard_horizon.Get())) {
			cerr << "--guard-horizon needs --bisector lsmear-guard" << endl;
			delete server;
			delete sys;
			return 1;
		}

		if (oracle_arg) server->set_oracle(true);
		if (oracle_depth_arg) server->set_oracle_max_depth(oracle_depth_arg.Get());
		if (sb_ratio_arg) server->set_sb_ratio(sb_ratio_arg.Get());
		if (dive_loup) server->set_dive_loup(dive_loup.Get());
		if (acid_fixed) server->set_acid_fixed(acid_fixed.Get());
		if (late_loup) server->set_late_loup(late_loup.Get());
		if (sb_vol_only) server->set_sb_vol_only(true);
		if (sb_need_pruned) server->set_sb_need_pruned(true);
		if (probe_dims) server->set_probe_dims(probe_dims.Get());
		if (probe_parts) server->set_probe_parts(probe_parts.Get());
		if (probe_ctc) {
			if (probe_ctc.Get()=="hc4") server->set_probe_ctc(MLNodeServer::PROBE_HC4);
			else if (probe_ctc.Get()=="lp") server->set_probe_ctc(MLNodeServer::PROBE_LP);
 			else if (probe_ctc.Get()=="full") server->set_probe_ctc(MLNodeServer::PROBE_FULL);
			else if (probe_ctc.Get()=="proc") server->set_probe_ctc(MLNodeServer::PROBE_PROC);
			else { cerr << "--probe-ctc: hc4, lp, full or proc" << endl; delete server; delete sys; return 1; }
		}
		if (oracle_score) {
			if (oracle_score.Get()=="depth") server->set_oracle_depth(true);
			else if (oracle_score.Get()=="sb") server->set_oracle_sb(true);
			else if (oracle_score.Get()=="hc4") server->set_oracle_hc4(true);
			else if (oracle_score.Get()=="pc") { server->set_oracle_pc(true, pc_rel ? pc_rel.Get() : 4); if (pc_light) server->set_pc_light(true); }
			else if (oracle_score.Get()!="nodes") { cerr << "--oracle-score: nodes, depth, sb, hc4 or pc" << endl; delete server; delete sys; return 1; }
		}

		if (model_file) {
			model = new MLModel(model_file.Get());
			server->set_model(model);
			if (progress) cerr << "[model] " << model->description() << endl;
		}

		double loup0 = initial_loup ? initial_loup.Get() : POS_INFINITY;

		MLNodeServer::SampleParams sp;
		sp.budget       = budget       ? budget.Get()       : 200;
		sp.budget_start = budget_start ? budget_start.Get() : 25;
		sp.prune        = !no_prune;
		sp.max_depth    = max_depth    ? max_depth.Get()    : 0;
		sp.include_goal = !no_goal;
		sp.topk         = topk         ? topk.Get()         : 0;

		MLNodeServer::RunParams rp;
		rp.sample_prob = sample_prob ? sample_prob.Get() : 0.1;
		rp.max_samples = max_samples ? max_samples.Get() : 1000;
		rp.stop_at_max_samples = stop_at_max_samples;
		rp.max_nodes   = max_nodes   ? max_nodes.Get()   : 0;
		rp.timeout     = timeout     ? timeout.Get()     : 0;
		rp.progress    = progress.Get();
		rp.features    = !no_features;

		if (collect) {
			ofstream file;
			ostream* os = &cout;
			if (output_file) {
				file.open(output_file.Get().c_str());
				if (!file.is_open()) ibex_error("cannot open the output file");
				os = &file;
			}

			long nb = server->collect(sys->box, loup0, *os, rp, sp);

			if (output_file) file.close();

			cerr << "[collect] " << nb << " samples, "
			     << server->get_nb_cells() << " nodes, "
			     << server->elapsed() << "s, "
			     << "loup=" << server->get_loup() << " uplo=" << server->get_uplo()
			     << " (" << server->status() << ")" << endl;

		} else if (run) {
			cout.precision(17);
			server->run_interactive(sys->box, loup0, cin, cout, rp, sp);

		} else if (solve) {
			server->solve(sys->box, loup0, rp);
			JsonOut out(cout);
			out.obj();
			out.kv("status", server->status());
			out.kv("loup", server->get_loup());
			out.kv("uplo", server->get_uplo());
			out.kv("nodes", server->get_nb_cells());
			out.kv("time", server->elapsed());
			out.kv("bisector", MLOptimizerConfig::bisector_name(bisector));
			out.kv("relax", MLOptimizerConfig::relaxation_name(relaxation));
			// not "loup": that key already carries the incumbent value above
			out.kv("loup_finder", MLOptimizerConfig::loup_name(loup));
			out.kv("bisect_ratio", bisect_ratio);
			out.kv("rule", model!=NULL ? model->description()
					: string(MLOptimizerConfig::bisector_name(bisector)));
			{
				// when the upper bound reached its final value (within the
				// solver's relative precision, 1e-3, and exactly): what node
				// selection can hope to move earlier
				const std::vector<MLNodeServer::LoupEvent>& H = server->loup_history;
				double fin = server->get_loup();
				long near_n=-1, last_n=-1; double near_t=-1, last_t=-1;
				for (size_t k=0; k<H.size(); k++) {
					if (near_n<0 && H[k].loup <= fin + 1e-3*std::fabs(fin)) { near_n = H[k].nodes; near_t = H[k].time; }
					if (H[k].loup <= fin) { last_n = H[k].nodes; last_t = H[k].time; break; }
				}
				out.kv("loup_updates", (long) H.size()-1);
				out.kv("ipopt_time", ibex::ipopt_total_time);
				out.kv("ipopt_corr_time", ibex::ipopt_corr_total_time);
				out.kv("cvx_time", ibex::cvx_total_time);
				if (ibex::ipopt_multistart>0) out.kv("ms_improved", (long) ibex::ipopt_ms_improved);
				if (ibex::cvx_close_on) { out.kv("cvx_raised", ibex::optimizer_extra_lb_raised); out.kv("cvx_pruned", ibex::optimizer_extra_lb_pruned); }
				if (ibex::cvx_stats) { out.kv("cvx_nodes", ibex::cvx_nodes); out.kv("cvx_with_dual", ibex::cvx_with_dual); out.kv("cvx_certified", ibex::cvx_certified); out.kv("cvx_close", ibex::cvx_close); out.kv("cvx_better_lb", ibex::cvx_better_lb); out.kv("cvx_center_indef", ibex::cvx_center_indef); }
				if (ibex::ipopt_guided>0) { out.kv("guided_calls", (long) ibex::ipopt_guided_calls); out.kv("guided_improved", (long) ibex::ipopt_guided_improved); out.kv("guided_time", ibex::ipopt_guided_time); }
				if (bound_trace_arg) {
					out.key("bounds"); out.arr();
					const std::vector<MLNodeServer::BoundEvent>& B = server->bound_history;
					for (size_t k=0; k<B.size(); k++) { out.arr(); out.val(B[k].time); out.val(B[k].nodes); out.val(B[k].loup); out.val(B[k].uplo); out.end_arr(); }
					out.end_arr();
				}
				if (getenv("IBEX_LOUP_POINT") && server->get_loup()<POS_INFINITY)
					out.kv("loup_point", server->get_loup_point());
				out.kv("loup_near_nodes", near_n);
				out.kv("loup_near_time", near_t);
				out.kv("loup_final_nodes", last_n);
				out.kv("loup_final_time", last_t);
				if (!server->acid_tuning().empty()) out.kv("acid_nbcidvar", (long) server->acid_tuning()[0].nbcidvar);
				if (acid_fixed) out.kv("acid_fixed", (long) acid_fixed.Get());
			}
			if (getenv("IBEX_ACID_TRACE")) for (size_t k=0; k<server->acid_tuning().size(); k++)
				std::cerr << "acid " << k << " nbcidvar=" << server->acid_tuning()[k].nbcidvar
				          << " nbcalls=" << server->acid_tuning()[k].nbcalls
				          << " nbtuning=" << server->acid_tuning()[k].nbtuning << std::endl;
			if (oracle_arg) {
				out.kv("oracle_calls", server->oracle_calls);
				out.kv("oracle_fallbacks", server->oracle_fallbacks);
				out.kv("oracle_deviations", server->oracle_deviations);
				out.kv("pc_probes", server->pc_probes);
			}
			if (server->guard_switched_at()>=-1) {
				out.kv("guard_horizon", server->guard_horizon());
				out.kv("guard_switched_at", server->guard_switched_at());
			}
			out.end_obj();
			cout << endl;

		} else {
			// ---------------- server mode ----------------
			server->reset(sys->box, loup0);
			cout.precision(17);

			string line;
			while (getline(cin, line)) {
				size_t first = line.find_first_not_of(" \t\r\n");
				if (first==string::npos) continue;

				try {
					JsonValue cmd = JsonValue::parse(line);
					if (!run_command(*server, cmd, sp)) break;
				} catch (JsonError& e) {
					fail(e.what());
				} catch (NoBisectableVariableException&) {
					fail("no bisectable variable");
				} catch (LPException&) {
					fail("LP solver failure");
				} catch (std::exception& e) {
					fail(string("error: ") + e.what());
				}
			}
		}

	} catch (ibex::SyntaxError& e) {
		cerr << e << endl;
		rc = 1;
	} catch (std::exception& e) {
		cerr << "error: " << e.what() << endl;
		rc = 1;
	}

	if (server) delete server;
	if (model) delete model;
	if (sys) delete sys;
	return rc;
}
