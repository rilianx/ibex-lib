#include "ibex_LPSolver.h"
#include <cstdlib>

namespace ibex {

/** \brief Stream out \a x. */
std::ostream& operator<<(std::ostream& os, const LPSolver::Status x){

	switch(x) {
	case(LPSolver::Status::Optimal) :{
			os << "OPTIMAL";
			break;
	}
	case(LPSolver::Status::Infeasible) :{
			os << "INFEASIBLE";
			break;
	}
	case(LPSolver::Status::OptimalProved) :{
			os << "OPTIMAL_PROVED";
			break;
	}
	case(LPSolver::Status::InfeasibleProved) :{
			os << "INFEASIBLE_PROVED";
			break;
	}
	case(LPSolver::Status::Timeout) :{
			os << "TIME_OUT";
			break;
	}
	case(LPSolver::Status::MaxIter) :{
			os << "MAX_ITER";
			break;
	}
	case(LPSolver::Status::Unbounded) :{
			os << "UNBOUNDED";
			break;
	}
	case(LPSolver::Status::Unknown) :{
		os << "UNKNOWN";
		break;
	}
	}
	return os;
}

bool LPSolver::neumaier_shcherbina_postprocessing() {
	/* Camino original, conservado para comparar: construye rows() (m x n) y su
	 * transpuesta (n x m) en cada resolucion del LP. Con IBEX_NS_LEGACY=1 se
	 * vuelve a el. */
	static const bool legacy = (getenv("IBEX_NS_LEGACY") != NULL);
	if (legacy) {
		Matrix A_trans = rows_transposed();
		IntervalVector b = lhs_rhs();
		IntervalVector rest = A_trans*uncertified_dual_;
		rest -= cost();
		obj_ = uncertified_dual_*b - rest*ivec_bounds_;
		return true;
	}

	/* Misma formula, acumulada por filas y salteando los duales nulos.
	 *
	 *     obj_ = lambda'*b - (A'*lambda - c)'*[x]
	 *
	 * Una fila con lambda_i == 0 no aporta ni a lambda'*b ni a A'*lambda, asi
	 * que se puede omitir: el dual del simplex es disperso (solo las filas
	 * activas en la base son no nulas), de modo que el trabajo pasa de m*n a
	 * nnz(lambda)*n. Y sobre todo desaparecen las dos matrices densas, que se
	 * asignaban y llenaban una vez por cota, o sea 2n veces por caja.
	 *
	 * La acumulacion de A'*lambda queda en punto flotante, igual que en el
	 * camino original (alli Matrix*Vector es un producto de dobles); solo
	 * cambia el orden de asociacion. */
	const int n = nb_vars();
	const int m = nb_rows();

	/* OJO: con IBEX_NS_RIGOROUS=1 el residuo A'*lambda - c se acumula en
	 * aritmetica de intervalos. El camino por omision lo acumula en dobles,
	 * igual que el original, y por eso NO es riguroso: el error de redondeo del
	 * residuo entra multiplicado por el ancho de la caja. El interruptor esta
	 * para medir cuanto costaria cerrar ese agujero. */
	static const bool rigorous = (getenv("IBEX_NS_RIGOROUS") != NULL);

	Interval dual_b(0.0);
	if (rigorous) {
		IntervalVector rest(n, Interval::zero());
		for (int i = 0; i < m; ++i) {
			const double di = uncertified_dual_[i];
			if (di == 0.0) continue;
			dual_b += di*lhs_rhs(i);
			const Vector r = row(i);
			for (int j = 0; j < n; ++j) {
				if (r[j] != 0.0) rest[j] += di*r[j];
			}
		}
		for (int j = 0; j < n; ++j) rest[j] -= cost(j);
		obj_ = dual_b - rest*ivec_bounds_;
		return true;
	}

	Vector rest(n, 0.0);
	for (int i = 0; i < m; ++i) {
		const double di = uncertified_dual_[i];
		if (di == 0.0) continue;
		dual_b += di*lhs_rhs(i);
		const Vector r = row(i);
		for (int j = 0; j < n; ++j) {
			if (r[j] != 0.0) rest[j] += di*r[j];
		}
	}
	rest -= cost();

	obj_ = dual_b - rest*ivec_bounds_;
	return true;
}

bool LPSolver::neumaier_shcherbina_infeasibility_test() {
    ibex::Matrix A_trans = rows_transposed();
    IntervalVector b = lhs_rhs();
    Vector lambda(1);
	// It is possible that the solver does not find an infeasible direction
	// even when the problem is infeasible.
    bool infeasible_dir_found = uncertified_infeasible_dir(lambda);
    if(!infeasible_dir_found) {
        return false;
    }


    IntervalVector rest = A_trans * lambda ;
    Interval d = rest * ivec_bounds_ - lambda * b;

    // if 0 does not belong to d, the infeasibility is proved
    return !d.contains(0.0);
}


void LPSolver::invalidate() {
    status_ = LPSolver::Status::Unknown;
    has_solution_ = false;
    has_infeasible_dir_ = false;
}

bool LPSolver::is_feasible() const {
	return (status_ == LPSolver::Status::Optimal && mode_ == LPSolver::Mode::NotCertified)
	|| (status_ == LPSolver::Status::OptimalProved && mode_ == LPSolver::Mode::Certified);
}

}  // end namespace ibex
