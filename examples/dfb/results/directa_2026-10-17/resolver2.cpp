// Reconstruye por la API un LP guardado por la caja directa (formato .txtlp exacto) y lo resuelve en frio.
// uso: resolver2 archivo.txtlp archivo.set
#include "soplex.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace soplex;
static double leer(FILE* f) { char b[64]; if (fscanf(f, "%63s", b) != 1) exit(2); return strtod(b, 0); }
int main(int argc, char** argv) {
    FILE* f = fopen(argv[1], "r"); int n, m, k, mx;
    if (!f || fscanf(f, "%d %d %d %d", &n, &m, &k, &mx) != 4) return 1;
    SoPlex sx;
    sx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
    sx.loadSettingsFile(argv[2]);
    sx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
    sx.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
    LPColSetReal cols(n); DSVectorReal vacio(0);
    for (int i = 0; i < n; ++i) { double lo = leer(f), up = leer(f); cols.add(i == k ? (mx ? -1.0 : 1.0) : 0.0, lo, vacio, up); }
    sx.addColsReal(cols);
    LPRowSetReal rows(m); double amax = 0, amin = 1e300;
    for (int j = 0; j < m; ++j) {
        double lo = leer(f), up = leer(f); DSVectorReal r(n);
        for (int i = 0; i < n; ++i) { double a = leer(f); if (a != 0) { r.add(i, a); double b = a < 0 ? -a : a; if (b > amax) amax = b; if (b < amin) amin = b; } }
        rows.add(lo, r, up);
    }
    sx.addRowsReal(rows);
    sx.optimize();
    int st = (int) sx.status();
    printf("estado %d iteraciones %d |a| en [%.1e, %.1e]\n", st, sx.numIterations(), amin, amax);
    return 0;
}
