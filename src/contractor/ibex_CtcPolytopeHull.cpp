#include <string>
#include <cmath>
//============================================================================
//                                  I B E X
// File        : ibex_CtcPolytopeHull.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Oct 31, 2013
// Last update : Aug 01, 2018
//============================================================================

#include "ibex_CtcPolytopeHull.h"
#include <cstdio>
#include <ctime>
#include <cstdlib>

#include "ibex_LinearizerFixed.h"

using namespace std;

namespace ibex {

namespace {

class PolytopeHullEmptyBoxException { };

}

CtcPolytopeHull::CtcPolytopeHull(Linearizer& lr, int max_iter, int time_out, double eps) :
		Ctc(lr.nb_var()), lr(lr),
		mylinearsolver(nb_var, LPSolver::Mode::Certified, eps, time_out, max_iter),
		contracted_vars(BitSet::all(nb_var)), own_lr(false), primal_sols(2*nb_var, nb_var),
		primal_sol_found(2*nb_var) {
	/* NO estaban inicializados: se acumulaba sobre basura. */
	n_soplex_iterations = 0; n_soplex_calls = 0;

}

CtcPolytopeHull::CtcPolytopeHull(const Matrix& A, const Vector& b, int max_iter, int time_out, double eps) :
		Ctc(A.nb_cols()), lr(*new LinearizerFixed(A,b)),
		mylinearsolver(nb_var, LPSolver::Mode::Certified, eps, time_out, max_iter),
		contracted_vars(BitSet::all(nb_var)), own_lr(true), primal_sols(2*nb_var, nb_var),
		primal_sol_found(2*nb_var) {
	n_soplex_iterations = 0; n_soplex_calls = 0;

}

namespace { double ph_t_lin = 0.0, ph_t_total = 0.0;
bool ph_on() { static const bool v = (getenv("IBEX_PH_STATS") != NULL); return v; }
double ph_reloj() { if (!ph_on()) return 0.0;
                    struct timespec ts; clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
                    return ts.tv_sec + 1e-9*ts.tv_nsec; } }
extern double ph_t_solve, ph_t_extrae, ph_t_cert;
extern bool ph_stats();

/* Estados de SoPlex por LP (`IBEX_PH_STATS=1`, resumen al destruir). */
static long ph_st[5] = {0,0,0,0,0};   // OptimalProved, InfeasibleProved, Infeasible, Unknown, otro
/* Holgura de Neumaier-Shcherbina: contraccion que permitiria el optimo
 * flotante de SoPlex contra la que deja la cota certificada, sumadas como
 * fraccion del ancho de la variable. */
static double ph_pot = 0.0, ph_real = 0.0; static long ph_nopt = 0, ph_perdida = 0;
/* Destino de las 2n cotas de cada llamada: resueltas con un LP, saltadas por
 * Achterberg (el punto del LP anterior ya estaba en la cota), o sin resolver
 * porque la pasada se corto (vacio sin prueba, desconocido, tope). */
static long ph_cotas = 0, ph_resueltas = 0, ph_achterberg = 0, ph_sin_resolver = 0, ph_llamadas = 0;
static void ph_holgura(double flot, double cert, const Interval& xi) {
	const double w = xi.diam();
	if (!(w > 0) || !std::isfinite(w) || !std::isfinite(flot)) return;
	const double pot = std::max(0.0, std::min(w, flot - xi.lb())) / w;
	const double re  = std::isfinite(cert) ? std::max(0.0, std::min(w, cert - xi.lb())) / w : 0.0;
	++ph_nopt; ph_pot += pot; ph_real += re;
	static const int vol = getenv("IBEX_PH_VOLCAR") ? atoi(getenv("IBEX_PH_VOLCAR")) : 0;
	if (ph_nopt <= vol) fprintf(stderr, "[ns] x=[%.10g,%.10g] flot=%.10g cert=%.10g pot=%.4f real=%.4f\n", xi.lb(), xi.ub(), flot, cert, pot, re);
	if (pot > 1e-6 && re < pot * 0.9) ++ph_perdida;
}
static void ph_cuenta(LPSolver::Status st) {
	if (st == LPSolver::Status::OptimalProved) ++ph_st[0];
	else if (st == LPSolver::Status::InfeasibleProved) ++ph_st[1];
	else if (st == LPSolver::Status::Infeasible) ++ph_st[2];
	else if (st == LPSolver::Status::Unknown) ++ph_st[3];
	else ++ph_st[4];
}

CtcPolytopeHull::~CtcPolytopeHull() {
	if (own_lr) delete &lr;
	if (getenv("IBEX_PH_STATS"))
		fprintf(stderr, "[phull] estados: optimo=%ld vacio_probado=%ld vacio_sin_prueba=%ld desconocido=%ld otro=%ld\n",
		        ph_st[0], ph_st[1], ph_st[2], ph_st[3], ph_st[4]);
	if (getenv("IBEX_PH_STATS") && ph_cotas)
		fprintf(stderr, "[phull] cotas: %ld en %ld llamadas | resueltas %.1f%%  saltadas por Achterberg %.1f%%  sin resolver (pasada cortada) %.1f%%\n",
		        ph_cotas, ph_llamadas, 100.0*ph_resueltas/ph_cotas, 100.0*ph_achterberg/ph_cotas, 100.0*ph_sin_resolver/ph_cotas);
	if (getenv("IBEX_PH_STATS") && ph_nopt)
		fprintf(stderr, "[phull] holgura NS: %ld LPs optimos | contraccion potencial (flotante) %.5f real (certificada) %.5f -> se pierde %.1f%% | LPs que pierden >10%%: %ld (%.1f%%)\n",
		        ph_nopt, ph_pot/ph_nopt, ph_real/ph_nopt, ph_pot > 0 ? 100.0*(1.0 - ph_real/ph_pot) : 0.0,
		        ph_perdida, 100.0*ph_perdida/ph_nopt);
}

void CtcPolytopeHull::add_property(const IntervalVector& init_box, BoxProperties& map) {
	lr.add_property(init_box, map);
}

void CtcPolytopeHull::contract(IntervalVector& box) {
	ContractContext context(box);
	contract(box,context);
}

void CtcPolytopeHull::contract(IntervalVector& box, ContractContext& context) {
	primal_sol_found.clear();
	mylinearsolver.clear_constraints();

	if (box.is_unbounded()) return;

	try {
		//returns the number of constraints in the linearized system
		const double ph_a = ph_reloj();
		int cont = lr.linearize(box, mylinearsolver, context.prop);
		ph_t_lin += ph_reloj() - ph_a;
		/* OJO: esta traza estaba SIN guardar, en el camino caliente: imprimia
		 * la propagacion en cada nodo del arbol. Ademas de ensuciar la salida,
		 * invalida cualquier comparacion de tiempos contra produccion. */
		//cout << context.prop << endl;

		//cout << "[polytope-hull] end of LR" << endl;

		if (cont==-1) throw PolytopeHullEmptyBoxException();

		if (cont==0) return;

		const double ph_b = ph_reloj();
		optimizer(box);
		ph_t_total += ph_reloj() - ph_b;

		//mylinearsolver.writeFile("LP.lp");
		//system ("cat LP.lp");
		//cout << "[polytope-hull] box after LR: " << box << endl;
	}
	catch(PolytopeHullEmptyBoxException& e) {
		box.set_empty(); // empty the box before exiting
		mylinearsolver.clear_constraints();
	}

	context.prop.update(BoxEvent(box,BoxEvent::CONTRACT));
}

void CtcPolytopeHull::set_contracted_vars(const BitSet& vars) {
	contracted_vars = vars;
}

void CtcPolytopeHull::optimizer(IntervalVector& box) {

	Interval opt(0.0);
	int* inf_bound = new int[nb_var]; // indicator inf_bound = 1 means the inf bound is feasible or already contracted, call to simplex useless (cf Baharev)
	int* sup_bound = new int[nb_var]; // indicator sup_bound = 1 means the sup bound is feasible or already contracted, call to simplex useless
	long ph_lps_esta = 0, ph_fijas = 0;

	for (int i=0; i<nb_var; i++) {

		if (contracted_vars[i]) {
			inf_bound[i]=0;
			sup_bound[i]=0;
		} else {
			inf_bound[i]=1;
			sup_bound[i]=1;
			ph_fijas += 2;
		}
	}

	/* Iteraciones de SoPlex por LP en PRODUCCION, o sea con el arranque tibio
	 * secuencial del unico LPSolver. Es el termino que faltaba para comparar
	 * contra los pivotes de DFB: las sondas resolvian SoPlex en frio, que no es
	 * lo que corre aca. `IBEX_PH_STATS=1` lo imprime cada 2000 LPs. */
	if (([]() { static const bool v_ = (getenv("IBEX_PH_STATS") != NULL); return v_; }()) && n_soplex_calls > 0 && (n_soplex_calls % 2000) < 2*nb_var)
	{
		const double T = ph_t_lin + ph_t_total;
		fprintf(stderr, "[phull] LPs=%d  iter/LP=%.2f | armado %.0f%%  simplex %.0f%%  "
		        "extrae %.0f%%  certif %.0f%%  otro %.0f%% | us/LP=%.1f  us/iter=%.1f "
		        "(total %.2fs)\n",
		        n_soplex_calls, (double)n_soplex_iterations/n_soplex_calls,
		        100*ph_t_lin/T, 100*ph_t_solve/T, 100*ph_t_extrae/T, 100*ph_t_cert/T,
		        100*(T-ph_t_lin-ph_t_solve-ph_t_extrae-ph_t_cert)/T,
		        1e6*T/n_soplex_calls, 1e6*T/n_soplex_iterations, T);
	}

	int nexti=-1;   // the next variable to be contracted
	int infnexti=0; // the bound to be contracted contract  infnexti=0 for the lower bound, infnexti=1 for the upper bound
	LPSolver::Status stat=LPSolver::Status::Unknown;

	// Update the bounds the variables
	mylinearsolver.set_bounds(box);

	for(int ii=0; ii<(2*nb_var); ii++) {  // at most 2*n calls

		int i= ii/2;
		if (nexti != -1) i=nexti;
		//cout << "[polytope-hull]->[optimize] var n°"<< i << " infnexti=" << infnexti << " infbound=" << inf_bound[i] << " supbound=" << sup_bound[i] << endl;
		//cout << "[polytope-hull]->[optimize] box before simplex: " << box << endl;
		if (infnexti==0 && inf_bound[i]==0)  // computing the left bound : minimizing x_i
		{
			inf_bound[i]=1;
			mylinearsolver.set_cost(i, 1.0);
			stat = mylinearsolver.minimize();
			ph_cuenta(stat); ++ph_resueltas; ++ph_lps_esta;
			n_soplex_iterations += mylinearsolver.mysoplex->numIterations();
			n_soplex_calls++;
			

			mylinearsolver.set_cost(i, 0.0);
			//cout << "[polytope-hull]->[optimize] simplex for left bound returns stat:" << stat << endl;
			if (stat == LPSolver::Status::OptimalProved) {
				opt = mylinearsolver.minimum();
				if (ph_stats()) ph_holgura(mylinearsolver.not_proved_primal_sol()[i], opt.lb(), box[i]);
				//std::cout << "opt=" << opt << endl;
				if(opt.lb()>box[i].ub()) {
					delete[] inf_bound;
					delete[] sup_bound;
					throw PolytopeHullEmptyBoxException();
				}
				primal_sols[2*i]=mylinearsolver.not_proved_primal_sol();
				primal_sol_found.add(2*i);

				if(opt.lb() > box[i].lb()) {
					box[i]=Interval(opt.lb(),box[i].ub());
					mylinearsolver.set_bounds(i,box[i]);
					history.push_back(make_pair(n_soplex_iterations, box.perimeter()));
				}

				if (!choose_next_variable(box,nexti,infnexti, inf_bound, sup_bound)) {
					break;
				}
			}
			else if (stat == LPSolver::Status::InfeasibleProved) {
				delete[] inf_bound;
				delete[] sup_bound;
				// the infeasibility is proved, the EmptyBox exception is raised
				throw PolytopeHullEmptyBoxException();
			}

			else if (stat == LPSolver::Status::Infeasible) {
				// the infeasibility is found but not proved, no other call is needed
				break;
			}

			else if (stat == LPSolver::Status::Unknown) {
				int next=-1;
				for (int j=0;j<nb_var;j++) {
					if (inf_bound[j]==0) {
						nexti=j;  next=0;  infnexti=0;
						break;
					}
					else if  (sup_bound[j]==0) {
						nexti=j;  next=0;  infnexti=1;
						break;
					}
				}
				if (next==-1)  break;
			}

		}
		else if (infnexti==1 && sup_bound[i]==0) { // computing the right bound :  maximizing x_i
			sup_bound[i]=1;
			mylinearsolver.set_cost(i, -1.0);
			stat= mylinearsolver.minimize();
			ph_cuenta(stat); ++ph_resueltas; ++ph_lps_esta;
			n_soplex_iterations += mylinearsolver.mysoplex->numIterations();
			n_soplex_calls++;


			mylinearsolver.set_cost(i, 0.0);
			//cout << "[polytope-hull]->[optimize] simplex for right bound returns stat=" << stat << endl;
			if( stat == LPSolver::Status::OptimalProved) {
				opt = -mylinearsolver.minimum();
				if (ph_stats()) ph_holgura(-mylinearsolver.not_proved_primal_sol()[i], -opt.ub(), -box[i]);
				//std::cout << "opt=" << opt << endl;
				if(opt.ub() <box[i].lb()) {
					delete[] inf_bound;
					delete[] sup_bound;
					throw PolytopeHullEmptyBoxException();
				}

				primal_sols[2*i+1]=mylinearsolver.not_proved_primal_sol();
				primal_sol_found.add(2*i+1);

				if (opt.ub() < box[i].ub()) {
					box[i] =Interval( box[i].lb(), opt.ub());
					mylinearsolver.set_bounds(i,box[i]);
					history.push_back(make_pair(n_soplex_iterations, box.perimeter()));
				}

				if (!choose_next_variable(box,nexti,infnexti, inf_bound, sup_bound)) {
					break;
				}
			}
			else if(stat == LPSolver::Status::InfeasibleProved) {
				delete[] inf_bound;
				delete[] sup_bound;
				// the infeasibility is proved,  the EmptyBox exception is raised
				throw PolytopeHullEmptyBoxException();
			}
			else if (stat == LPSolver::Status::Infeasible) {
				// the infeasibility is found but not proved, no other call is needed
				break;
			}

			else if (stat == LPSolver::Status::Unknown) {
				int next=-1;
				for (int j=0;j<nb_var;j++) {
					if (inf_bound[j]==0) {
						nexti=j;  next=0;  infnexti=0;
						break;
					}
					else if (sup_bound[j]==0) {
						nexti=j;  next=0;  infnexti=1;
						break;
					}
				}
				if (next==-1) break;
			}
		}
		else break; // in case of stat==MAX_ITER  we do not recall the simplex on a another variable  (for efficiency reason)
	}
	if (ph_stats()) {
		/* las marcadas sin LP son las de Achterberg; las que quedaron en 0, sin resolver */
		long marcadas = 0, ceros = 0;
		for (int j = 0; j < nb_var; ++j) { marcadas += inf_bound[j] + sup_bound[j]; ceros += (inf_bound[j]==0) + (sup_bound[j]==0); }
		++ph_llamadas; ph_cotas += 2*nb_var - ph_fijas; ph_sin_resolver += ceros;
		ph_achterberg += marcadas - ph_fijas - ph_lps_esta;
	}
	delete[] inf_bound;
	delete[] sup_bound;
}

bool CtcPolytopeHull::choose_next_variable(IntervalVector & box, int & nexti, int & infnexti, int* inf_bound, int* sup_bound) {

	bool found = false;
	// the primal solution : used by choose_next_variable
	Vector primal_solution = mylinearsolver.not_proved_primal_sol();
	//cout << " primal " << primal_solution << endl;

	// The Achterberg heuristic for choosing the next variable (nexti) and its bound (infnexti) to be contracted (cf Baharev paper)
	// and updating the indicators if a bound has been found feasible (with the precision prec_bound)
	// called only when a primal solution is found by the LP solver (use of primal_solution)

	// double prec_bound = mylinearsolver.getEpsilon(); // relative precision for the indicators TODO change with the precision of the optimizer ??
	double prec_bound = 1.e-8; // relative precision for the indicators      :  compatibility for testing  BNE
	/* `IBEX_PH_NOACH=1`: sonda, apaga el salto de Achterberg (prec_bound = 0). */
	static const bool noach = (getenv("IBEX_PH_NOACH") != NULL);
	if (noach) prec_bound = -1.0;
	double delta=1.e100;
	double deltaj=delta;

	for (int j=0; j<nb_var; j++)	{

		/* `IBEX_PH_SALTEO=ancho`: tolerancia relativa al ancho de la caja, no a
		 * |cota|; la de Achterberg saltea cotas contraibles en cajas angostas
		 * lejos del cero. */
		static const bool salteo_ancho = getenv("IBEX_PH_SALTEO") && std::string(getenv("IBEX_PH_SALTEO")) == "ancho";
		const double wj = box[j].diam();
		const bool usa_ancho = salteo_ancho && wj > 0 && wj < POS_INFINITY;
		if (inf_bound[j]==0) {
			deltaj= fabs(primal_solution[j]- box[j].lb());
			if ((usa_ancho && deltaj < prec_bound*wj) ||
			    (!usa_ancho && fabs (box[j].lb()) < 1 && deltaj < prec_bound) ||
					(!usa_ancho && fabs (box[j].lb()) >= 1 && fabs (deltaj /(box[j].lb())) < prec_bound))	{
				inf_bound[j]=1;
			}
			if (inf_bound[j]==0 && deltaj < delta) 	{
				nexti=j; infnexti=0;delta=deltaj; found =true;
			}
		}

		if (sup_bound[j]==0) {
			deltaj = fabs (primal_solution[j]- box[j].ub());


			if ((usa_ancho && deltaj < prec_bound*wj) ||
			    (!usa_ancho && fabs (box[j].ub()) < 1 && deltaj < prec_bound) 	||
					(!usa_ancho && fabs (box[j].ub()) >= 1 && fabs (deltaj/(box[j].ub())) < prec_bound)) {
				sup_bound[j]=1;
			}
			if (sup_bound[j]==0 && deltaj < delta) {
				nexti=j; infnexti=1;delta=deltaj;  found =true;
			}

		}


	}
	return found;
}

} // end namespace ibex
