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
  -I../../lp_lib_wrapper/soplex/soplex-4.0.2/src 


LDFLAGS = \
  -L../../src \
  -L../../interval_lib_wrapper/gaol/gaol-4.2.0-build \
  -L../../interval_lib_wrapper/gaol/mathlib-2.1.0-build \
  -L../../lib

LDLIBS = -libex -lsoplex -lgaol -lgdtoa -lultim

# Fuentes por ejecutable
ibexopt_SRCS = ibexopt.cpp ibex_CtcDualFeasibleBounding.cpp ibex_Optimizer05Config.cpp ibex_CtcDFBPropag.cpp
example_dfb_SRCS = dfb_non_linear_system.cpp ibex_CtcDualFeasibleBounding.cpp ibex_Optimizer05Config.cpp ibex_CtcDFBPropag.cpp ../../src/contractor/ibex_CtcPolytopeHull.cpp
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
