/* ============================================================================
 * I B E X - Implementation of the Linearizer of Affine2 forms
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENCE.
 *
 * Author(s)   : Jordan Ninin
 * Created     : June 6, 2020
 * ---------------------------------------------------------------------------- */

#include <vector>
#include <math.h>
#include "ibex_LinearizerAffine2.h"
#include "ibex_Exception.h"

#include <time.h>
namespace ibex {

/* Reparto interno de `linearize_id` (`DFBH_LINPERF=1`): cuanto del costo es
 * FIJO por llamada —armar las variables afines, que es O(n^2) sin importar
 * cuantas restricciones se evaluen— y cuanto es proporcional a las
 * restricciones. Decide si linealizar selectivamente puede abaratar algo. */
/* `pow(2,-50)` estaba DENTRO del bucle de variables, o sea `nb_ctr * nb_var`
 * llamadas a `pow` por linealizacion. Es una constante. */
static const double DOS_A_MENOS_50 = 8.8817841970012523e-16;
/* `DFBH_LINSLOW=1` restaura el bucle anterior —`pow(2,-50)` adentro y `rad`/`mid`
 * recalculados por restriccion— para medir el arreglo a igualdad de binario.
 * Los dos caminos dan resultados identicos bit a bit. */
static bool lin_lento() { static const bool v = (getenv("DFBH_LINSLOW") != NULL); return v; }
double lin_t_vars = 0.0, lin_t_eval = 0.0, lin_t_fila = 0.0;
long   lin_n_llam = 0, lin_n_ctr = 0;
static double lin_ahora() {
	struct timespec ts; clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
	return ts.tv_sec + 1e-9*ts.tv_nsec;
}
void lin_perf_volcar() {
	if (!lin_n_llam) return;
	const double T = lin_t_vars + lin_t_eval + lin_t_fila;
	if (T <= 0) return;
	fprintf(stderr, "[linperf] llamadas=%ld  restricciones evaluadas=%ld (%.2f por llamada) | "
	        "total=%.3fs | variables afines (COSTO FIJO) %.1f%%  evaluacion por restriccion %.1f%%  "
	        "extraccion de fila %.1f%% | us por llamada: fijo %.1f  eval %.1f  fila %.1f\n",
	        lin_n_llam, lin_n_ctr, (double)lin_n_ctr/lin_n_llam, T,
	        100*lin_t_vars/T, 100*lin_t_eval/T, 100*lin_t_fila/T,
	        1e6*lin_t_vars/lin_n_llam, 1e6*lin_t_eval/lin_n_llam, 1e6*lin_t_fila/lin_n_llam);
}

// the constructor
LinearizerAffine2::LinearizerAffine2(const System& sys1) :
				Linearizer(sys1.nb_var), sys(sys1),
				goal_af_evl(NULL),
				ctr_af_evl(new Affine2Eval*[sys1.nb_ctr]) {

	if (sys1.goal) {
		goal_af_evl = new Affine2Eval(*sys1.goal);
	}

	for (int i = 0; i < sys.nb_ctr; i++) {
		ctr_af_evl[i] = new Affine2Eval(sys.ctrs[i].f);
	}
}

LinearizerAffine2::~LinearizerAffine2() {
	for (int i = 0; i < sys.nb_ctr; i++) {
		delete ctr_af_evl[i];
	}
	delete[] ctr_af_evl;
}

bool LinearizerAffine2::goal_linearization(const IntervalVector& box, LPSolver& lp_solver) {
	// Linearization of the objective function by AF2

	if (!sys.goal) {
		ibex_error("LinearRelaxAffine2: there is no goal function to linearize.");
	}

	std::pair<Domain*,Affine2Domain*> res = goal_af_evl->eval(box, Affine2Variables(box));
	Affine2 af2 = res.second->i();
	if (af2.is_empty()) {
		return false;
	}
	try {
	if (af2.size() == sys.nb_var) { // if the affine2 form is valid
		// convert the epsilon variables to the original box
		double tmp=0;
		for (int i =0; i <sys.nb_var; i++) {
			tmp = box[i].rad();
			if (tmp==0) { // sensible case to avoid rowconst[i]=NaN
				if (af2.val(i)==0)
					lp_solver.set_cost(i, 0);
				else {
					return false; // sensible case to avoid
				}
			} else {
				lp_solver.set_cost(i, af2.val(i) / tmp);
			}
		}
	} else {
		return false;
	}
	return true;
	} catch (LPException&) {
		return false;
	}
}


int LinearizerAffine2::inlinearization(const IntervalVector& box, LPSolver& lp_solver) {

	Affine2 af2;

	int cont=0;
	Interval ev(0), center(0), err(0);
	Vector rowconst(sys.nb_var);

	// Create the linear relaxation of each constraint
	for (int ctr = 0; ctr < sys.nb_ctr; ctr++) {
		CmpOp op = sys.ctrs[ctr].op;

		std::pair<Domain*,Affine2Domain*> res = ctr_af_evl[ctr]->eval(box, Affine2Variables(box));
		ev  = res.first->i();
		af2 = res.second->i();
		//ev  = ctr_af_evl[ctr]->eval(box).i();
		//af2 = ctr_af_evl[ctr]->af2.top->i();

		//std::cout <<ev<<":::"<< af2<<"  "<<af2.size()<<"  " <<sys.nb_var<< std::endl;

		if (af2.size() == sys.nb_var) { // if the affine2 form is valid
			bool b_abort=false;
			// convert the epsilon variables to the original box
			double tmp=0;
			center =0;
			err =0;
			for (int i =0;(!b_abort) &&(i <sys.nb_var); i++) {
				tmp = box[i].rad();
				if (tmp==0) { // sensible case to avoid rowconst[i]=NaN
					if (af2.val(i)==0)
						rowconst[i]=0;
					else {
						b_abort =true;
					}
				} else {
					rowconst[i] =af2.val(i) / tmp;
					center += rowconst[i]*box[i].mid();
					err += fabs(rowconst[i]) * DOS_A_MENOS_50;
				}
			}
			if (!b_abort) {
				switch (op) {
				case LEQ:
				case LT: {
					if (0.0 < ev.ub()) {
						lp_solver.add_constraint(rowconst, LEQ,	(-(af2.err()+err) - (af2.mid()-center)).lb());
						cont++;
					}
					break;
				}
				case GEQ:
				case GT: {
					if (ev.lb() < 0.0) {
						lp_solver.add_constraint(rowconst, GEQ,	((af2.err()+err) - (af2.mid()-center)).ub());
						cont++;
					}
					break;
				}
				case EQ: {
					not_implemented("LinearRelaxAffine2::inlinearization not implemented for equality constraints");
				}
				default:
					break;
				}
			}
		}

	}

	return -1;
}


int LinearizerAffine2::nb_ctr() const { return sys.nb_ctr; }

Affine2 LinearizerAffine2::forma_afin(int ctr, const IntervalVector& box) {
	Affine2Variables varaf2(box);
	std::pair<Domain*,Affine2Domain*> res = ctr_af_evl[ctr]->eval(box, varaf2);
	Affine2 af2 = res.second->i();
	if (res.first->i().is_empty() || af2.size() != sys.nb_var) af2.set_empty();
	return af2;
}

const std::vector<int>& LinearizerAffine2::used_vars(int ctr) const {
	return sys.ctrs[ctr].f.used_vars;
}

/* Composicion de la relajacion, para explicar por que `linearize` y
 * `linearize_id` no emiten el mismo numero de filas (`DFBH_LINCOMP=1`).
 * `linearize` descarta las redundantes y las LT/GT, y duplica cada EQ;
 * `linearize_id` sin `DFBH_LINMISMO` hace lo contrario en los dos sentidos. */
static long lc_llamadas=0, lc_ctr_vistas=0, lc_ltgt=0, lc_redund=0, lc_utiles=0, lc_eq=0, lc_eq_chica=0;
namespace { struct LcVolcado { ~LcVolcado() {
	if (!getenv("DFBH_LINCOMP") || !lc_llamadas) return;
	const double k = 1.0/ (double)lc_llamadas;
	fprintf(stderr,
	  "[lincomp] por llamada: utiles(LEQ/GEQ) %.2f  redundantes %.2f  LT/GT %.2f"
	  "  EQ %.2f  EQ-chica %.2f\n"
	  "[lincomp] filas linearize  = utiles + 2*EQ                 = %.2f\n"
	  "[lincomp] filas linearize_id = utiles + redund + LTGT + EQ + EQchica = %.2f\n",
	  lc_utiles*k, lc_redund*k, lc_ltgt*k, lc_eq*k, lc_eq_chica*k,
	  (lc_utiles + 2.0*lc_eq)*k,
	  (lc_utiles + lc_redund + lc_ltgt + lc_eq + lc_eq_chica)*k);
} } lc_volcado; }

int LinearizerAffine2::linearize_id(const IntervalVector& box, std::vector<Vector>& rows,
                                    std::vector<Interval>& rango, std::vector<int>& ids,
                                    const std::vector<char>* evaluar) {
	++lc_llamadas;
	rows.clear(); rango.clear(); ids.clear();
	static const bool linperf = (getenv("DFBH_LINPERF") != NULL);
	const double t0 = linperf ? lin_ahora() : 0.0;
	Affine2 af2;
	Affine2Variables varaf2(box);
	/* `rad` y `mid` de la caja no dependen de la restriccion: se calculaban
	 * `nb_ctr` veces cada uno. Se izan fuera del bucle; los valores son los
	 * mismos, asi que el resultado es identico bit a bit. */
	std::vector<double> rad_(sys.nb_var), mid_(sys.nb_var);
	for (int i = 0; i < sys.nb_var; i++) { rad_[i] = box[i].rad(); mid_[i] = box[i].mid(); }
	if (linperf) { lin_t_vars += lin_ahora() - t0; ++lin_n_llam; }
	Vector rowconst(sys.nb_var);
	Interval ev(0.0), center(0.0), err(0.0);

	for (int ctr = 0; ctr < sys.nb_ctr; ctr++) {
		/* Con mascara, solo se evaluan las restricciones marcadas: el resto
		 * conserva la fila y el rango de su ultima evaluacion, que siguen
		 * siendo validos porque la caja solo se encogio. */
		if (evaluar && ctr < (int)evaluar->size() && !(*evaluar)[ctr]) continue;
		const CmpOp op = sys.ctrs[ctr].op;
		const double t1 = linperf ? lin_ahora() : 0.0;
		std::pair<Domain*,Affine2Domain*> res = ctr_af_evl[ctr]->eval(box, varaf2);
		if (linperf) { lin_t_eval += lin_ahora() - t1; ++lin_n_ctr; }
		ev  = res.first->i();
		af2 = res.second->i();
		if (ev.is_empty()) af2.set_empty();
		if (af2.size() != sys.nb_var) continue;   /* sin forma afin: se omite */

		const double t2 = linperf ? lin_ahora() : 0.0;
		bool b_abort = false;
		center = 0; err = 0;
		for (int i = 0; (!b_abort) && (i < sys.nb_var); i++) {
			const double tmp = lin_lento() ? box[i].rad() : rad_[i];
			if (tmp == 0) {
				if (af2.val(i) == 0) rowconst[i] = 0; else b_abort = true;
			} else {
				rowconst[i] = af2.val(i) / tmp;
				center += rowconst[i]*(lin_lento() ? box[i].mid() : mid_[i]);
				err += fabs(rowconst[i]) * (lin_lento() ? pow(2,-50) : DOS_A_MENOS_50);
			}
		}
		if (b_abort) continue;

		/* `g(x)` en espacio de x: fila.x + k + E[-1,1], con k = af2.mid()-center. */
		const Interval E = af2.err() + err;
		const Interval k = af2.mid() - center;
		const double sup = (E - k).ub();    /* fila.x <= sup  si  g <= 0 */
		const double inf = (-E - k).lb();   /* fila.x >= inf  si  g >= 0 */

		/* `DFBH_LINMISMO=1`: emitir EXACTAMENTE las filas que emite
		 * `linearize`, para que la relajacion congelada recien congelada sea
		 * comparable fila a fila con la fresca. Sin esto `linearize_id` emite
		 * de mas (las redundantes sobre la caja, y las LT/GT que `linearize`
		 * solo usa como prueba de infactibilidad) y de menos (una fila
		 * bilatera por EQ en vez de dos unilateras). El poliedro es el mismo,
		 * pero `m` cambia y con el el tableau, el orden de seleccion de cota y
		 * el gasto de pivotes. */
		/* Dos diferencias INDEPENDIENTES contra `linearize`, cada una con su
		 * interruptor, porque se cancelaban entre si en el conteo de filas:
		 *
		 *  `DFBH_LINRED=1`  descarta las filas redundantes sobre la caja y las
		 *     LT/GT, igual que `linearize`. Sin el, `linearize_id` las emite:
		 *     son validas pero no cortan, solo agrandan el tableau.
		 *  `DFBH_LINEQ2=1`  parte cada EQ en dos filas unilateras, igual que
		 *     `linearize`. Sin el, va una sola fila bilatera. El poliedro es el
		 *     mismo en los dos casos, pero dos filas con la MISMA normal
		 *     duplican una cota: el tableau crece y la heuristica gasta una
		 *     optimizacion entera en una cota repetida.
		 *
		 * Los dos juntos reproducen `linearize` fila a fila (`DFBH_LINMISMO=1`).
		 * En el sistema extendido la unica EQ es el objetivo `y = f(x)`, asi que
		 * `LINEQ2` decide si la fila mas importante de la relajacion va una o dos
		 * veces. */
		static const bool mismo  = (getenv("DFBH_LINMISMO") != NULL);
		static const bool sin_red = mismo || (getenv("DFBH_LINRED")  != NULL);
		static const bool eq_dos  = mismo || (getenv("DFBH_LINEQ2") != NULL);
		static const double tol_eq = 2.0 * 1e-9;   /* LPSolver::default_tolerance */

		++lc_ctr_vistas;
		switch (op) {
		case LT: case GT:            ++lc_ltgt;   break;
		case LEQ: if (!(0.0 < ev.ub())) ++lc_redund; else ++lc_utiles; break;
		case GEQ: if (!(ev.lb() < 0.0)) ++lc_redund; else ++lc_utiles; break;
		case EQ:  if (ev.diam() <= tol_eq) ++lc_eq_chica; else ++lc_eq; break;
		default: break;
		}

		switch (op) {
		case LT:
			if (sin_red) { if (ev.lb() == 0.0) return -1; break; }
			/* fall through */
		case LEQ:
			if (0.0 < ev.lb()) return -1;             /* infactible */
			if (sin_red && !(0.0 < ev.ub())) break;   /* redundante sobre la caja */
			rows.push_back(rowconst); rango.push_back(Interval(NEG_INFINITY, sup));
			ids.push_back(2*ctr);
			break;
		case GT:
			if (sin_red) { if (ev.ub() == 0.0) return -1; break; }
			/* fall through */
		case GEQ:
			if (ev.ub() < 0.0) return -1;
			if (sin_red && !(ev.lb() < 0.0)) break;
			rows.push_back(rowconst); rango.push_back(Interval(inf, POS_INFINITY));
			ids.push_back(2*ctr+1);
			break;
		case EQ:
			if (!ev.contains(0.0)) return -1;
			if (eq_dos) {
				if (sin_red && ev.diam() <= tol_eq) break;
				/* dos filas unilateras, en el mismo orden que `linearize` */
				rows.push_back(rowconst); rango.push_back(Interval(inf, POS_INFINITY));
				ids.push_back(2*ctr+1);
				rows.push_back(rowconst); rango.push_back(Interval(NEG_INFINITY, sup));
				ids.push_back(2*ctr);
				break;
			}
			if (sin_red && ev.diam() <= tol_eq) break;
			rows.push_back(rowconst); rango.push_back(Interval(inf, sup));
			ids.push_back(2*ctr);
			break;
		default: break;
		}
		if (linperf) lin_t_fila += lin_ahora() - t2;
	}
	return (int)rows.size();
}

/*********generation of the linearized system*********/
int LinearizerAffine2::linearize(const IntervalVector& box, LPSolver& lp_solver) {

	Affine2 af2;
	Affine2Variables varaf2(box);
	Vector rowconst(sys.nb_var);
	Interval ev(0.0);
	Interval center(0.0);
	Interval err(0.0);
	CmpOp op;
	int cont = 0;
	fila_ctr.clear();
	const bool lento = lin_lento();
	std::vector<double> radL(sys.nb_var), midL(sys.nb_var);
	if (!lento) for (int i = 0; i < sys.nb_var; i++) { radL[i] = box[i].rad(); midL[i] = box[i].mid(); }

	// Create the linear relaxation of each constraint
	for (int ctr = 0; ctr < sys.nb_ctr; ctr++) {

		op  = sys.ctrs[ctr].op;

		std::pair<Domain*,Affine2Domain*> res = ctr_af_evl[ctr]->eval(box, varaf2);
		ev  = res.first->i();
		af2 = res.second->i();
		//ev  = ctr_af_evl[ctr]->eval(box).i();
		//af2 = ctr_af_evl[ctr]->af2.top->i();


		if (ev.is_empty()) {
			af2.set_empty();
		}
		//std::cout <<ev<<":::"<< af2<<"  "<<af2.size()<<"  " <<sys.nb_var<< std::endl;

		if (af2.size() == sys.nb_var) { // if the affine2 form is valid
			bool b_abort=false;
			// convert the epsilon variables to the original box
			double tmp=0;
			center =0;
			err =0;
			for (int i =0;(!b_abort) &&(i <sys.nb_var); i++) {
				tmp = lento ? box[i].rad() : radL[i];
				if (tmp==0) { // sensible case to avoid rowconst[i]=NaN
					if (af2.val(i)==0)
						rowconst[i]=0;
					else {
						b_abort =true;
					}
				} else {
					rowconst[i] =af2.val(i) / tmp;
					center += rowconst[i]*(lento ? box[i].mid() : midL[i]);
					err += fabs(rowconst[i]) * (lento ? pow(2,-50) : DOS_A_MENOS_50);
				}
			}
			if (!b_abort) {
				switch (op) {
				case LT:
					if (ev.lb() == 0.0) return -1;
					break;
				case LEQ:
					if (0.0 < ev.lb()) return -1;
					else if (0.0 < ev.ub()) {
						try {
							lp_solver.add_constraint(rowconst, LEQ,	((af2.err()+err) - (af2.mid()-center)).ub());
							cont++; fila_ctr.push_back(ctr);
						} catch (LPException&) { }
					}
					break;
				case GT:
					if (ev.ub() == 0.0) return -1;
					break;
				case GEQ:
					if (ev.ub() < 0.0) return -1;
					else if (ev.lb() < 0.0) {
						try {
							lp_solver.add_constraint(rowconst, GEQ,	(-(af2.err()+err) - (af2.mid()-center)).lb());
							cont++; fila_ctr.push_back(ctr);
						} catch (LPException&) { }
					}
					break;
				case EQ:
					if (!ev.contains(0.0)) return -1;
					else {
						if (ev.diam()>2*lp_solver.tolerance()) {
							try {
								lp_solver.add_constraint(rowconst, GEQ,	(-(af2.err()+err) - (af2.mid()-center)).lb());
								cont++; fila_ctr.push_back(ctr);
								lp_solver.add_constraint(rowconst, LEQ,	((af2.err()+err) - (af2.mid()-center)).ub());
								cont++; fila_ctr.push_back(ctr);
							} catch (LPException&) { }
						}
					}
					break;
				}
			}
		}

	}
	return cont;

}

//void CtcART::convert_back(IntervalVector & box, IntervalVector & epsilon) {
//
//	for (int i = 0; i < box.size(); i++) {
//		box[i] &= box[i].mid() + (box[i].rad() * epsilon[i]);
//	}
//}

}
