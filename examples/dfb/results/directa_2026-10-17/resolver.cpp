// Resuelve un LP guardado por la caja directa (LP, parametros y base opcional) e informa estado e iteraciones.
// uso: resolver archivo.lp archivo.set [archivo.bas]
#include "soplex.h"
#include <cstdio>
using namespace soplex;
int main(int argc, char** argv) {
    SoPlex sx;
    sx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
    if (!sx.loadSettingsFile(argv[2])) { printf("sin parametros\n"); }
    sx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
    if (!sx.readFile(argv[1])) { printf("no se pudo leer %s\n", argv[1]); return 1; }
    if (argc > 3 && !sx.readBasisFile(argv[3])) printf("base no leida\n");
    sx.optimize();
    int st = (int) sx.status();
    printf("estado %d iteraciones %d objetivo %.12g\n", st, sx.numIterations(), st == 1 ? sx.objValueReal() : 0.0);
    return 0;
}
