/* Stub of the optional external LP engines (SoPlex 7, HiGHS, Gurobi).
 *
 * The contractor CtcDFBHull can solve the bound LPs with these engines
 * (DFBH_CAJA_SOPLEX7, DFBH_CAJA_HIGHS, DFBH_CAJA_GUROBI) for side comparisons
 * mentioned in the paper. None of the ladder rungs R0-R5 uses them. When they
 * are not installed, build_paper.sh links this file instead, and selecting one
 * of these engines aborts with a message. */
#include "ibex_Caja7.h"
#include "ibex_CajaHighs.h"
#include "ibex_CajaGurobi.h"
#include <cstdio>
#include <cstdlib>

namespace {
[[noreturn]] void no_disponible(const char* motor) {
    std::fprintf(stderr, "%s was not linked: rebuild with build_paper.sh --externals\n", motor);
    std::exit(2);
}
}

namespace ibex {
#ifndef CAJA7_ENLAZADA   /* -DCAJA7_ENLAZADA: la caja de SoPlex 7/8 se enlaza aparte */
Caja7::Caja7() : sx(nullptr) { no_disponible("SoPlex 7"); }
Caja7::~Caja7() {}
void Caja7::cargar(int, int, const double*, const double*, const double*, const double*, const double*, double, int) {}
void Caja7::cota(int, double, double) {}
int Caja7::resolver(int, bool, bool, double*, double*, int*) { return -1; }
#endif

CajaHighs::CajaHighs() : h(nullptr), n_(0), m_(0) { no_disponible("HiGHS"); }
CajaHighs::~CajaHighs() {}
void CajaHighs::cargar(int, int, const double*, const double*, const double*, const double*, const double*, double, bool) {}
void CajaHighs::cota(int, double, double) {}
int CajaHighs::resolver(int, bool, bool, double*, double*, int*) { return -1; }

CajaGurobi::CajaGurobi() : env(nullptr), mod(nullptr), n_(0), m_(0) { no_disponible("Gurobi"); }
CajaGurobi::~CajaGurobi() {}
void CajaGurobi::cargar(int, int, const double*, const double*, const double*, const double*, const double*, double, bool) {}
void CajaGurobi::cota(int, double, double) {}
int CajaGurobi::resolver(int, bool, bool, double*, double*, int*) { return -1; }
}
