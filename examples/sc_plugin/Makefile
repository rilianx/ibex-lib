# Builds CtcSC and two drivers against the in-place Ibex of this repository
# (built with cmake, GAOL and SoPlex 4.0.2). The affine linearizer used by the
# drivers comes from examples/dfb/affine (it is not part of Ibex).

CXX = g++
CXXFLAGS = -O3 -msse3 -frounding-math -std=gnu++11 -MMD -MP

IBEX = ../..
AFF  = ../dfb/affine

INCLUDES = \
  -I. -I$(AFF) \
  -I$(IBEX)/interval_lib_wrapper/gaol \
  -I$(IBEX)/lp_lib_wrapper/soplex \
  -I$(IBEX)/src \
  -I$(IBEX)/src/arithmetic -I$(IBEX)/src/bisector -I$(IBEX)/src/cell -I$(IBEX)/src/combinatorial \
  -I$(IBEX)/src/contractor -I$(IBEX)/src/data -I$(IBEX)/src/function -I$(IBEX)/src/loup \
  -I$(IBEX)/src/numeric -I$(IBEX)/src/optim -I$(IBEX)/src/parser -I$(IBEX)/src/predicate -I$(IBEX)/src/set \
  -I$(IBEX)/src/solver -I$(IBEX)/src/strategy -I$(IBEX)/src/symbolic -I$(IBEX)/src/operators \
  -I$(IBEX)/src/system -I$(IBEX)/src/tools \
  -I$(IBEX)/interval_lib_wrapper/gaol/gaol-4.2.0 \
  -I$(IBEX)/interval_lib_wrapper/gaol/gaol-4.2.0-build \
  -I$(IBEX)/interval_lib_wrapper/gaol/gaol-4.2.0-build/gaol \
  -I$(IBEX)/interval_lib_wrapper/gaol/mathlib-2.1.0/src \
  -I$(IBEX)/interval_lib_wrapper/gaol/mathlib-2.1.0-build/src \
  -I$(IBEX)/lp_lib_wrapper/soplex/soplex-4.0.2/src

LDFLAGS = \
  -L$(IBEX)/src \
  -L$(IBEX)/interval_lib_wrapper/gaol/gaol-4.2.0-build \
  -L$(IBEX)/interval_lib_wrapper/gaol/mathlib-2.1.0-build \
  -L$(IBEX)/lib

LDLIBS = -libex -lsoplex -lgaol -lgdtoa -lultim

SC_SRCS  = ibex_CtcSC.cpp ibex_SCSimplex.cpp ibex_SCBasis.cpp
AFF_SRCS = $(AFF)/ibex_Affine2_fAF2.cpp $(AFF)/ibex_Affine3_fAFFullI.cpp $(AFF)/ibex_AffineVar.cpp \
           $(AFF)/ibex_CtcAffineFwdBwd.cpp $(AFF)/ibex_LinearizerAffine2.cpp

SC_OBJS  = $(SC_SRCS:.cpp=.o)
AFF_OBJS = $(patsubst $(AFF)/%.cpp,aff_%.o,$(AFF_SRCS))

TARGETS = sc_solve sc_opt

all: $(TARGETS)

sc_solve: sc_solve.o $(SC_OBJS) $(AFF_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

sc_opt: sc_opt.o $(SC_OBJS) $(AFF_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

aff_%.o: $(AFF)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f *.o *.d $(TARGETS)

-include $(wildcard *.d)

.PHONY: all clean
