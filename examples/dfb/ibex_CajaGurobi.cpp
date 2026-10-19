#include "gurobi_c.h"
#include "ibex_CajaGurobi.h"
#include <vector>
#include <cstdio>
#include <cstdlib>
namespace ibex {
#define ENV ((GRBenv*) env)
#define MOD ((GRBmodel*) mod)
CajaGurobi::CajaGurobi() : env(nullptr), mod(nullptr), n_(0), m_(0) {
    GRBenv* e = nullptr;
    if (GRBemptyenv(&e) || GRBsetintparam(e, "OutputFlag", 0) || GRBstartenv(e)) {
        fprintf(stderr, "[gurobi] no se pudo crear el entorno: %s\n", e ? GRBgeterrormsg(e) : "?");
        exit(1);
    }
    GRBsetintparam(e, "Threads", 1);
    GRBsetintparam(e, "Method", getenv("DFBH_GRB_METODO") ? atoi(getenv("DFBH_GRB_METODO")) : 1);
    GRBsetintparam(e, "InfUnbdInfo", 1);
    GRBsetintparam(e, "DualReductions", 0);
    env = e;
}
CajaGurobi::~CajaGurobi() {
    if (mod) GRBfreemodel(MOD);
    if (env) GRBfreeenv(ENV);
}
static double inf_g(double v) { return v >= 1e100 ? GRB_INFINITY : (v <= -1e100 ? -GRB_INFINITY : v); }
void CajaGurobi::cargar(int n, int m, const double* lb, const double* ub, const double* lhs,
                        const double* rhs, const double* A, double tol, bool presolve) {
    if (mod) { GRBfreemodel(MOD); mod = nullptr; }
    std::vector<double> l(n), u(n);
    for (int i = 0; i < n; ++i) { l[i] = inf_g(lb[i]); u[i] = inf_g(ub[i]); }
    GRBmodel* g = nullptr;
    GRBnewmodel(ENV, &g, "caja", n, nullptr, l.data(), u.data(), nullptr, nullptr);
    GRBenv* me = GRBgetenv(g);
    const double t = tol < 1e-9 ? 1e-9 : (tol > 1e-2 ? 1e-2 : tol);
    GRBsetdblparam(me, "FeasibilityTol", t);
    GRBsetdblparam(me, "OptimalityTol", t);
    GRBsetintparam(me, "Presolve", presolve ? -1 : 0);
    GRBsetintparam(me, "ModelSense", GRB_MINIMIZE);
    std::vector<int> ind; std::vector<double> val;
    for (int j = 0; j < m; ++j) {
        const double lo = inf_g(lhs[j]), up = inf_g(rhs[j]);
        ind.clear(); val.clear();
        const bool libre = (lo == -GRB_INFINITY && up == GRB_INFINITY);
        if (!libre)
            for (int i = 0; i < n; ++i) {
                const double a = A[(long) j * n + i];
                if (a != 0.0) { ind.push_back(i); val.push_back(a); }
            }
        const int nz = (int) ind.size();
        int* pi = nz ? ind.data() : nullptr; double* pv = nz ? val.data() : nullptr;
        /* una fila con dos lados finitos distintos va como rango: Gurobi le
         * agrega una holgura propia, que no cambia el numero de filas ni `Pi` */
        if (libre)               GRBaddconstr(g, 0, nullptr, nullptr, GRB_LESS_EQUAL, GRB_INFINITY, nullptr);
        else if (lo == up)       GRBaddconstr(g, nz, pi, pv, GRB_EQUAL, up, nullptr);
        else if (lo == -GRB_INFINITY) GRBaddconstr(g, nz, pi, pv, GRB_LESS_EQUAL, up, nullptr);
        else if (up == GRB_INFINITY)  GRBaddconstr(g, nz, pi, pv, GRB_GREATER_EQUAL, lo, nullptr);
        else                     GRBaddrangeconstr(g, nz, pi, pv, lo, up, nullptr);
    }
    GRBupdatemodel(g);
    mod = g; n_ = n; m_ = m;
}
void CajaGurobi::cota(int i, double lb, double ub) {
    GRBsetdblattrelement(MOD, "LB", i, inf_g(lb));
    GRBsetdblattrelement(MOD, "UB", i, inf_g(ub));
}
int CajaGurobi::resolver(int k, bool maximizar, bool frio, double* lam, double* x, int* iters) {
    GRBmodel* g = MOD;
    GRBsetdblattrelement(g, "Obj", k, maximizar ? -1.0 : 1.0);
    if (frio) GRBreset(g, 0);
    GRBoptimize(g);
    int st = 0; double it = 0;
    GRBgetintattr(g, "Status", &st);
    GRBgetdblattr(g, "IterCount", &it);
    *iters = (int) it;
    int res = 3;
    if (st == GRB_OPTIMAL) {
        if (!GRBgetdblattrarray(g, "Pi", 0, m_, lam) && !GRBgetdblattrarray(g, "X", 0, n_, x)) res = 0;
    } else if (st == GRB_INFEASIBLE) {
        res = GRBgetdblattrarray(g, "FarkasDual", 0, m_, lam) ? 2 : 1;
    }
    GRBsetdblattrelement(g, "Obj", k, 0.0);
    return res;
}
}
