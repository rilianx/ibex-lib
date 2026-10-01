# Ejecutables
TARGETS = ibexopt example_dfb

# Compilador y flags
CXX = g++
CXXFLAGS = -O3 -msse3 -frounding-math -std=gnu++11 -g -MMD -MP

# Includes
INCLUDES = \
  -I../../interval_lib_wrapper/gaol \
  -I../../lp_lib_wrapper/soplex \
  -I../../src \
  -I../../src/arithmetic -I../../src/bisector -I../../src/cell -I../../src/combinatorial \
  -I../../src/contractor -I../../src/data -I../../src/function -I../../src/loup \
  -I../../src/numeric -I../../src/optim -I../../src/parser -I../../src/predicate -I../../src/set \
  -I../../src/solver -I../../src/strategy -I../../src/symbolic -I../../src/symbolic/../operators \
  -I../../src/system -I../../src/tools \
  -I../../interval_lib_wrapper/gaol/gaol-4.2.0 \
  -I../../interval_lib_wrapper/gaol/gaol-4.2.0-build \
  -I../../interval_lib_wrapper/gaol/gaol-4.2.0-build/gaol \
  -I../../interval_lib_wrapper/gaol/mathlib-2.1.0/src \
  -I../../interval_lib_wrapper/gaol/mathlib-2.1.0-build/src \
  -I../../lp_lib_wrapper/soplex/soplex-4.0.2/src \
  -Iaffine


LDFLAGS = \
  -L../../src \
  -L../../interval_lib_wrapper/gaol/gaol-4.2.0-build \
  -L../../interval_lib_wrapper/gaol/mathlib-2.1.0-build \
  -L../../lib

LDLIBS = -libex -lsoplex -lgaol -lgdtoa -lultim

# Fuentes por ejecutable
ibexopt_SRCS = ibex_CtcDFBHull.cpp ibexopt.cpp ibex_CtcDualFeasibleBounding.cpp ibex_Optimizer05Config.cpp ibex_CtcDFBPropag.cpp ibex_DFBSparse.cpp ibex_DFBSimplex.cpp ibex_DFBBasis.cpp affine/ibex_Affine2_fAF2.cpp affine/ibex_Affine3_fAFFullI.cpp affine/ibex_AffineVar.cpp affine/ibex_CtcAffineFwdBwd.cpp affine/ibex_LinearizerAffine2.cpp
example_dfb_SRCS = dfb_non_linear_system.cpp ibex_CtcDualFeasibleBounding.cpp ibex_Optimizer05Config.cpp ibex_CtcDFBPropag.cpp ../../src/contractor/ibex_CtcPolytopeHull.cpp ibex_DFBSparse.cpp ibex_DFBSimplex.cpp ibex_DFBBasis.cpp
# Todos los .cpp en un set único
SRCS = $(ibexopt_SRCS) $(example_dfb_SRCS)
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.d)

# Mapeo ejecutables a objetos
ibexopt_OBJS = $(ibexopt_SRCS:.cpp=.o)
example_dfb_OBJS = $(example_dfb_SRCS:.cpp=.o)

all: $(TARGETS)

ibexopt: $(ibexopt_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

example_dfb: $(example_dfb_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

-include $(DEPS)

run-ibexopt: ibexopt
	./ibexopt example.txt --filtering=dfb --lr=no

run-example_dfb: example_dfb
	./example_dfb

clean:
	rm -f $(OBJS) $(TARGETS) $(DEPS)

.PHONY: all clean run-ibexopt run-example_dfb

# --- Banco de pruebas y diagnostico de DFB ---
# Ojo: SRCS/DEPS de arriba no incluyen estos fuentes, asi que sus dependencias
# de cabecera se incluyen aparte (EXTRA_DEPS). Sin esto, make no rehace los .o
# cuando cambia ibex_CtcDualFeasibleBounding.h y quedan objetos con layouts de
# clase distintos, que fallan en tiempo de ejecucion.
bench_dfb_SRCS = bench_dfb.cpp ibex_CtcDualFeasibleBounding.cpp ibex_CtcDFBPropag.cpp ../../src/contractor/ibex_CtcPolytopeHull.cpp ibex_DFBSparse.cpp ibex_DFBSimplex.cpp ibex_DFBBasis.cpp
bench_dfb_OBJS = $(bench_dfb_SRCS:.cpp=.o)

bench_dfb: $(bench_dfb_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

diag_dfb_SRCS = diag_dfb.cpp ibex_CtcDualFeasibleBounding.cpp ibex_CtcDFBPropag.cpp ../../src/contractor/ibex_CtcPolytopeHull.cpp ibex_DFBSparse.cpp ibex_DFBSimplex.cpp ibex_DFBBasis.cpp
diag_dfb_OBJS = $(diag_dfb_SRCS:.cpp=.o)

diag_dfb: $(diag_dfb_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

# Arnes de regresion de los casos de prueba de DFB (los 7 habilitados)
tests_dfb_all_SRCS = tests_dfb_all.cpp ibex_CtcDualFeasibleBounding.cpp ibex_CtcDFBPropag.cpp ../../src/contractor/ibex_CtcPolytopeHull.cpp ibex_DFBSparse.cpp ibex_DFBSimplex.cpp ibex_DFBBasis.cpp
tests_dfb_all_OBJS = $(tests_dfb_all_SRCS:.cpp=.o)

tests_dfb_all: $(tests_dfb_all_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

# Simplex dual de variables acotadas, validado contra SoPlex
test_dfbsimplex_SRCS = test_dfbsimplex.cpp ibex_DFBSimplex.cpp ibex_DFBBasis.cpp
test_dfbsimplex_OBJS = $(test_dfbsimplex_SRCS:.cpp=.o)

test_dfbsimplex: $(test_dfbsimplex_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

# Equivalencia bit a bit entre el tableau denso y el disperso
test_dfbsparse_SRCS = test_dfbsparse.cpp ibex_DFBSparse.cpp
test_dfbsparse_OBJS = $(test_dfbsparse_SRCS:.cpp=.o)

test_dfbsparse: $(test_dfbsparse_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

# Dispersion de los terminos a*z dentro de una fila (limite del escalado)
spread_SRCS = spread.cpp ibex_CtcDualFeasibleBounding.cpp ibex_CtcDFBPropag.cpp ibex_DFBSparse.cpp ibex_DFBSimplex.cpp ibex_DFBBasis.cpp
spread_OBJS = $(spread_SRCS:.cpp=.o)

spread: $(spread_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

# Dispersion del pricing: decide si la factorizacion puede abaratarlo
probe_price_SRCS = probe_price.cpp ibex_CtcDualFeasibleBounding.cpp ibex_CtcDFBPropag.cpp ibex_DFBSparse.cpp ibex_DFBSimplex.cpp ibex_DFBBasis.cpp
probe_price_OBJS = $(probe_price_SRCS:.cpp=.o)

probe_price: $(probe_price_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

# Equivalencia entre el tableau denso y la base factorizada (etapa 1)
test_dfbbasis_SRCS = test_dfbbasis.cpp ibex_DFBBasis.cpp
test_dfbbasis_OBJS = $(test_dfbbasis_SRCS:.cpp=.o)

test_dfbbasis: $(test_dfbbasis_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)


# --- Duelo directo DFB vs SoPlex, sin propagacion ---
# Misma linealizacion, misma caja fija, las mismas 2n cotas en el mismo orden y
# el mismo criterio de parada (hasta el optimo). Es la unica comparacion que
# aisla el METODO del volumen de LPs y del intercalado con HC4.
duelo_dfb_SRCS = duelo_dfb.cpp ibex_CtcDualFeasibleBounding.cpp ibex_CtcDFBPropag.cpp ../../src/contractor/ibex_CtcPolytopeHull.cpp ibex_DFBSparse.cpp ibex_DFBSimplex.cpp ibex_DFBBasis.cpp
duelo_dfb_OBJS = $(duelo_dfb_SRCS:.cpp=.o)

duelo_dfb: $(duelo_dfb_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

probe_lp_SRCS = probe_lp.cpp ibex_DFBSimplex.cpp ibex_DFBBasis.cpp
probe_lp_OBJS = $(probe_lp_SRCS:.cpp=.o)

probe_lp: $(probe_lp_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

EXTRA_SRCS = probe_lp.cpp ibex_CtcDFBHull.cpp duelo_dfb.cpp spread.cpp bench_dfb.cpp diag_dfb.cpp tests_dfb_all.cpp incr_dfb.cpp chain_dfb.cpp \
             test_dfbbasis.cpp ibex_DFBBasis.cpp probe_price.cpp \
             test_dfbsparse.cpp ibex_DFBSparse.cpp \
             test_dfbsimplex.cpp ibex_DFBSimplex.cpp
EXTRA_DEPS = $(EXTRA_SRCS:.cpp=.d)
-include $(EXTRA_DEPS)

.PHONY: clean-dfb-tools
clean-dfb-tools:
	rm -f duelo_dfb bench_dfb diag_dfb tests_dfb_all incr_dfb $(EXTRA_SRCS:.cpp=.o) $(EXTRA_DEPS)

# Metrica incremental: pivotes tibio vs frio tras un cambio de cotas
incr_dfb_SRCS = incr_dfb.cpp ibex_CtcDualFeasibleBounding.cpp ibex_CtcDFBPropag.cpp ../../src/contractor/ibex_CtcPolytopeHull.cpp ibex_DFBSparse.cpp ibex_DFBSimplex.cpp ibex_DFBBasis.cpp
incr_dfb_OBJS = $(incr_dfb_SRCS:.cpp=.o)

incr_dfb: $(incr_dfb_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

# Una base en cadena para las n cotas, contra n contractores independientes
chain_dfb_SRCS = chain_dfb.cpp ibex_CtcDualFeasibleBounding.cpp ibex_CtcDFBPropag.cpp ../../src/contractor/ibex_CtcPolytopeHull.cpp ibex_DFBSparse.cpp ibex_DFBSimplex.cpp ibex_DFBBasis.cpp
chain_dfb_OBJS = $(chain_dfb_SRCS:.cpp=.o)

chain_dfb: $(chain_dfb_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)
