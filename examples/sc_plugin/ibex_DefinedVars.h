//============================================================================
// Variables definidas de un sistema: las que una igualdad v = g(x) expresa
// como funcion de las demas.  Un solver no necesita bisectarlas: cuando las
// variables libres son chicas, HC4 hacia adelante achica tambien las definidas.
//
// Deteccion estructural (no depende de nombres ni de quien escribio el modelo):
//   - restriccion de igualdad f = 0 cuya expresion, aplanada en sus sumandos de
//     primer nivel (+, -, menos unario), tiene un sumando  c*v  con c constante
//     no nula y v una variable escalar que no aparece en los demas sumandos.
//     Ibex guarda "v = g" como v - g = 0, y su simplificacion puede sacar
//     constantes ((v - g) - 2), por eso se aplana en vez de mirar solo v - g;
//   - entre esas candidatas se aceptan solo definiciones sin ciclos: v se acepta
//     cuando todas las variables de g son libres o ya aceptadas.
// Si una variable tiene varias definiciones, basta una aceptada.
//
// Pensado para el Minibex de ibex_cse (auxiliares de CSE), pero sirve para
// cualquier modelo con variables auxiliares escritas como definiciones.
//============================================================================

#ifndef __IBEX_DEFINED_VARS_H__
#define __IBEX_DEFINED_VARS_H__

#include "ibex.h"
#include <vector>

namespace ibex {

namespace defvars_detail {

/* indice de la variable del sistema si e es un simbolo escalar de f, si no -1 */
inline int indice_simbolo(const Function& f, const ExprNode& e) {
	const ExprSymbol* s = dynamic_cast<const ExprSymbol*>(&e);
	if (!s || !s->dim.is_scalar()) return -1;
	for (int k = 0; k < f.nb_arg(); k++)
		if (&f.arg(k) == s) return k;
	return -1;
}

/* variables (indices) que aparecen en e */
inline std::vector<int> variables_de(const Function& f, const ExprNode& e) {
	std::vector<int> vs;
	ExprSubNodes nodos(e);
	for (int i = 0; i < nodos.size(); i++) {
		const int k = indice_simbolo(f, nodos[i]);
		if (k >= 0) vs.push_back(k);
	}
	return vs;
}

/* aplana los sumandos de primer nivel de e: (nodo, signo) */
inline void sumandos(const ExprNode& e, bool neg, std::vector<std::pair<const ExprNode*, bool> >& out) {
	if (const ExprAdd* a = dynamic_cast<const ExprAdd*>(&e)) {
		sumandos(a->left, neg, out); sumandos(a->right, neg, out); return;
	}
	if (const ExprSub* b = dynamic_cast<const ExprSub*>(&e)) {
		sumandos(b->left, neg, out); sumandos(b->right, !neg, out); return;
	}
	if (const ExprMinus* m = dynamic_cast<const ExprMinus*>(&e)) {
		sumandos(m->expr, !neg, out); return;
	}
	out.push_back(std::make_pair(&e, neg));
}

/* si t es  v  o  c*v / v*c  con c constante escalar no nula, devuelve v; si no -1 */
inline int variable_sola(const Function& f, const ExprNode& t) {
	int v = indice_simbolo(f, t);
	if (v >= 0) return v;
	const ExprMul* m = dynamic_cast<const ExprMul*>(&t);
	if (!m) return -1;
	for (int lado = 0; lado < 2; lado++) {
		const ExprConstant* c = dynamic_cast<const ExprConstant*>(lado == 0 ? &m->left : &m->right);
		if (!c || !c->dim.is_scalar() || c->get_value().contains(0.0)) continue;
		v = indice_simbolo(f, lado == 0 ? m->right : m->left);
		if (v >= 0) return v;
	}
	return -1;
}

} // namespace defvars_detail

/**
 * \brief Marca las variables definidas del sistema (true = definida).
 *
 * Solo para sistemas con variables escalares; si hay alguna variable no
 * escalar devuelve todo false (no se detecta nada).
 */
inline std::vector<bool> defined_variables(const System& sys) {
	using namespace defvars_detail;
	const int n = sys.nb_var;
	std::vector<bool> definida(n, false);
	if (sys.args.size() != n) return definida;          // hay variables no escalares

	// candidatas: (v, variables de g)
	std::vector<std::pair<int, std::vector<int> > > cand;
	for (int c = 0; c < sys.nb_ctr; c++) {
		const NumConstraint& ctr = sys.ctrs[c];
		if (ctr.op != EQ || ctr.f.image_dim() != 1 || ctr.f.nb_arg() != n) continue;
		std::vector<std::pair<const ExprNode*, bool> > ts;
		sumandos(ctr.f.expr(), false, ts);
		std::vector<std::vector<int> > vars(ts.size());
		for (size_t t = 0; t < ts.size(); t++) vars[t] = variables_de(ctr.f, *ts[t].first);
		for (size_t t = 0; t < ts.size(); t++) {
			const int v = variable_sola(ctr.f, *ts[t].first);
			if (v < 0) continue;
			std::vector<int> g;
			bool aparece = false;
			for (size_t u = 0; u < ts.size(); u++) {
				if (u == t) continue;
				for (size_t w = 0; w < vars[u].size(); w++) {
					if (vars[u][w] == v) aparece = true;
					g.push_back(vars[u][w]);
				}
			}
			if (!aparece) { cand.push_back(std::make_pair(v, g)); break; }
		}
	}

	// aceptacion sin ciclos: v se acepta si todas las variables de g son libres
	// (no candidatas) o ya aceptadas
	std::vector<bool> es_cand(n, false);
	for (size_t t = 0; t < cand.size(); t++) es_cand[cand[t].first] = true;
	bool cambio = true;
	while (cambio) {
		cambio = false;
		for (size_t t = 0; t < cand.size(); t++) {
			const int v = cand[t].first;
			if (definida[v]) continue;
			bool ok = true;
			for (size_t u = 0; u < cand[t].second.size() && ok; u++) {
				const int w = cand[t].second[u];
				if (es_cand[w] && !definida[w]) ok = false;
			}
			if (ok) { definida[v] = true; cambio = true; }
		}
	}

	// nunca dejar el sistema sin variables libres
	int libres = 0;
	for (int i = 0; i < n; i++) if (!definida[i]) libres++;
	if (libres == 0) definida.assign(n, false);
	return definida;
}

} // namespace ibex

#endif
