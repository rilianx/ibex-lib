# Builds CtcSC and two drivers (sc_solve, sc_opt) against an installed Ibex.
#
# Ibex must be built with SoPlex as LP library (CtcSC uses its LU):
#   cmake -DLP_LIB=soplex -DCMAKE_INSTALL_PREFIX=<prefix> .. && make && make install
# Then:
#   make                                   # if ibex.pc is found by pkg-config
#   make PKG_CONFIG_PATH=<prefix>/share/pkgconfig
#
# The affine linearizer (affine/, by J. Ninin, LGPL) is used only by the
# drivers: Ibex develop does not ship it.

PKG_CONFIG ?= pkg-config
IBEX_CFLAGS := $(shell PKG_CONFIG_PATH="$(PKG_CONFIG_PATH)" $(PKG_CONFIG) --cflags ibex)
IBEX_LIBS   := $(shell PKG_CONFIG_PATH="$(PKG_CONFIG_PATH)" $(PKG_CONFIG) --libs ibex)
ifeq ($(IBEX_LIBS),)
$(error ibex.pc not found: set PKG_CONFIG_PATH to the folder of ibex.pc)
endif

CXX      ?= g++
CXXFLAGS ?= -O3 -frounding-math
CXXFLAGS += -std=c++14 -MMD -MP
CPPFLAGS += -I. -Iaffine $(IBEX_CFLAGS)
# Ibex installs its libraries in <prefix>/lib (and gaol/soplex in lib/ibex/3rd):
comma := ,
RPATH := $(patsubst -L%,-Wl$(comma)-rpath$(comma)%,$(filter -L%,$(IBEX_LIBS)))

SC_SRCS  = ibex_CtcSC.cpp ibex_SCSimplex.cpp ibex_SCBasis.cpp
AFF_SRCS = $(wildcard affine/*.cpp)
SC_OBJS  = $(SC_SRCS:.cpp=.o)
AFF_OBJS = $(AFF_SRCS:.cpp=.o)

all: sc_solve sc_opt

sc_solve: sc_solve.o $(SC_OBJS) $(AFF_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(IBEX_LIBS) $(RPATH)

sc_opt: sc_opt.o $(SC_OBJS) $(AFF_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(IBEX_LIBS) $(RPATH)

clean:
	rm -f *.o *.d affine/*.o affine/*.d sc_solve sc_opt

-include $(wildcard *.d affine/*.d)

.PHONY: all clean
