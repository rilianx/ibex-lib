#!/bin/bash
# Builds the binaries used in the experiments of the SC paper
# (examples/dfb/paper/main.tex):
#   ibexopt_esc, ibexsolve_dfb_esc   the optimizer and the system solver
#   ibexreplay_esc, ibexreplay_cmp   replay of the banks of passes
#                                    (both are the same program: _cmp is the
#                                    name used by the compo scripts)
#
# Prerequisite: Ibex built in place from the root of this repository, with
# GAOL and SoPlex 4.0.2 (the defaults):
#   cmake -DCMAKE_BUILD_TYPE=Release -DINTERVAL_LIB=gaol -DLP_LIB=soplex . && make
# which leaves src/libibex.a and lib/libsoplex.a.
#
# Usage:
#   ./build_paper.sh                 core build: every rung R0-R5 of the paper
#   ./build_paper.sh --externals     also links SoPlex 7.1.6, HiGHS 1.15.1 and
#                                    Gurobi 12.0.3, used only for the side
#                                    comparisons of Section 5 (not in the ladder)
#     with SX7SRC=<SoPlex 7.1.6 sources>  SX7=<its build dir, configured with
#                                        cmake -DCMAKE_CXX_FLAGS=-Dsoplex=soplex7>
#          HIGHSSRC=<HiGHS 1.15.1 sources> HIGHS=<its build dir>
#          GUROBI=<gurobi1203/linux64>
set -e
cd "$(dirname "$0")"

ext=""; libs=""
if [ "$1" = "--externals" ]; then
  : "${SX7SRC:?}" "${SX7:?}" "${HIGHSSRC:?}" "${HIGHS:?}" "${GUROBI:?}"
  g++ -O3 -std=gnu++14 -Dsoplex=soplex7 -I$SX7SRC/src -I$SX7 -c ibex_Caja7.cpp -o ibex_Caja7.o
  g++ -O3 -std=gnu++17 -I$HIGHSSRC/highs -I$HIGHSSRC -I$HIGHS -I$HIGHSSRC/extern -c ibex_CajaHighs.cpp -o ibex_CajaHighs.o
  g++ -O3 -I$GUROBI/include -c ibex_CajaGurobi.cpp -o ibex_CajaGurobi.o
  ext="ibex_Caja7.o ibex_CajaHighs.o ibex_CajaGurobi.o"
  libs="$SX7/lib/libsoplex.a $HIGHS/lib/libhighs.a -L$GUROBI/lib -lgurobi120 -Wl,-rpath,$GUROBI/lib -lpthread -lz"
else
  make ibex_CajasExternasStub.o
  ext="ibex_CajasExternasStub.o"
fi
# every object of the two programs, from the sources of this commit
objs=$( (make -n -W ibexopt.o ibexopt; make -n -W ibexsolve_dfb.o ibexsolve_dfb) | grep -o '[^ ]*\.o' | sort -u)
make $objs ibex_CajaSPx.o ibexreplay.o
ext="ibex_CajaSPx.o $ext"

# The link lines are those of the Makefile, plus the engines above.
l=$(make -n -W ibexopt.o ibexopt | tail -1)
eval "${l/-o ibexopt /-o ibexopt_esc } $ext $libs"
eval "$(echo "$l" | sed -e 's/ ibexopt\.o / ibexreplay.o /' -e 's/-o ibexopt /-o ibexreplay_esc /') $ext $libs"
cp ibexreplay_esc ibexreplay_cmp
l=$(make -n -W ibexsolve_dfb.o ibexsolve_dfb | tail -1)
eval "${l/-o ibexsolve_dfb /-o ibexsolve_dfb_esc } $ext $libs"
ls -la ibexopt_esc ibexsolve_dfb_esc ibexreplay_esc ibexreplay_cmp
