# CtcSC: self-certifying simplex contractor for Ibex

`CtcSC` is a polytope-hull contractor with the role and the interface of
`ibex::CtcPolytopeHull`. You give it a `Linearizer`. On each call it
linearizes the system over the box and contracts every variable to its bounds
over the relaxation.

Unlike `CtcPolytopeHull`, it does not call an external LP solver. The 2n
bound LPs are solved by its own bounded-variable dual simplex
(`SCSimplex`). All 2n LPs share one basis, and each one is warm-started from
the last basis of the previous LP. Each bound is certified once, from the
multipliers of the final basis, with the Neumaier–Shcherbina argument in
interval arithmetic. Any multipliers give a valid bound, so the simplex never
has to verify its solution:

- an LP stopped at its iteration limit still gives a certified bound;
- an infeasibility whose Farkas ray does not certify emptiness is solved once
  more from a cold start.

What a pass does:

- twin rows of the relaxation are merged: identical rows, and rows
  proportional within a relative tolerance of 1e-9;
- the first LP minimizes x_0;
- a bound already attained by the last primal point (relative distance
  below 1e-8) is skipped;
- the next LP is the bound nearest to that point;
- the limit is 100 pivots per LP, as Ibex gives SoPlex;
- the LU of the basis comes from SoPlex (`SLUFactor`), with Forrest–Tomlin
  updates and a refactorization every 50 updates;
- x_B is recomputed every 32 pivots;
- the ratio test is Harris with a relative band of 1e-3;
- the pass stops as soon as the relaxation is proved empty.

This is the configuration called *SC base* in the paper
*Certificate-driven polytope-hull contraction for interval branch and bound*
(I. Araya). The research code it was extracted from is `examples/dfb`.

## Files

| file | content |
|---|---|
| `ibex_CtcSC.{h,cpp}` | the contractor (linearization, merging, the pass, certification) |
| `ibex_SCSimplex.{h,cpp}` | the revised dual/primal simplex |
| `ibex_SCBasis.{h,cpp}` | sparse LU of the basis (SoPlex `SLUFactor`, hidden behind a pimpl) |
| `sc_solve.cpp`, `sc_opt.cpp` | drivers: a system solver and an optimizer (not part of the contractor) |
| `Makefile` | builds the drivers against an installed Ibex (`pkg-config ibex`) |
| `Makefile.inplace` | the same against the in-place Ibex of the repository of the paper |
| `affine/`, `parse_args.h` | third-party code used only by the drivers |

The contractor uses only the Ibex core and the SoPlex headers. It has no
global or static mutable state and reads no environment variables.

## Interface

```cpp
#include "ibex_CtcSC.h"

CtcSC(Linearizer& lr, CtcSC::WarmStart ws = CtcSC::DUAL,
      double tol = 1e-9, int max_iter = -1 /* 100 */);
```

- `CtcSC::DUAL` relocates the nonbasic variables and then uses dual steps.
  It is the best choice on systems of equations.
- `CtcSC::PRIMAL` keeps the nonbasic variables in place and uses primal steps
  when needed. It is the best choice on optimization problems.
- `CtcSC::ADAPTIVE` alternates the two on the first calls, keeps the one with
  the lower CPU time per LP, and tries the other again every 25 calls. Primal
  steps use devex pricing. The choice depends on measured time, so the search
  tree is not reproducible from run to run. The paper uses it with the
  combined (affine + XTaylor) relaxation on systems.
- `tol` is the feasibility tolerance of the simplex. Use 1e-6 for relaxations
  with nearly parallel rows, such as `LinearizerCompo`.

It drops in where `CtcPolytopeHull` goes, for example inside a solver
strategy:

```cpp
LinearizerXTaylor lin(sys);                 // or any Linearizer
CtcSC sc(lin);                              // instead of CtcPolytopeHull hull(lin);
CtcHC4 hc4(sys.ctrs, 0.01);
CtcFixPoint fp(*new CtcCompo(sc, hc4), 0.2);
```

`add_property` forwards to the linearizer, like `CtcPolytopeHull` does.

## Building

`CtcSC` is a plugin: it is not part of Ibex, and it builds against an
installed Ibex, release 2.8 or develop. Ibex must use SoPlex as its LP
library, because `SCBasis` uses the LU factorization of SoPlex:

```
# Ibex (once)
cmake -DCMAKE_BUILD_TYPE=Release -DINTERVAL_LIB=gaol -DLP_LIB=soplex \
      -DCMAKE_INSTALL_PREFIX=$HOME/ibex ..
make && make install

# this folder
make PKG_CONFIG_PATH=$HOME/ibex/share/pkgconfig
./sc_solve <ibex>/benchs/solver/polynom/Bellido.bch --lin=art -t 60 --seed=1
./sc_opt <ibex>/benchs/optim/blowup/meanvar.bch --filtering=acidhc4 --lin=art --seed=1 -t 60
```

`make -f Makefile.inplace` builds instead against the in-place Ibex of the
repository of the paper.

To use `CtcSC` in your own program, compile the three files
`ibex_CtcSC.cpp`, `ibex_SCSimplex.cpp` and `ibex_SCBasis.cpp` with it, with the
flags of `pkg-config --cflags ibex`. The contractor uses only Ibex API that
exists in both 2.8 and develop:

- `Linearizer::linearize(box, LPSolver&, BoxProperties&)`;
- `LPSolver(n, LPSolver::Mode::NotCertified)`, `rows()`, `lhs_rhs()`,
  `clear_constraints()`;
- `ContractContext`, `BoxEvent`.

Third-party files bundled for the drivers only:

- `affine/` is the affine-arithmetic linearizer of Jordan Ninin (ENSTA
  Bretagne, LGPL), which Ibex develop does not ship. This copy is the one of
  the paper; it keeps a few research switches (`DFBH_LIN*`, `DFBH_AFNOPOOL`)
  that do nothing unless they are set.
- `parse_args.h` is the argument parser of Ibex (`src/bin`, MIT license).

The drivers use the same strategies as `examples/dfb/ibexsolve_dfb --lr=dfbhull`
and `examples/dfb/ibexopt --lr=dfbhull` of the paper:

- `sc_solve`: HC4, ACID and interval Newton, then `FixPoint(CtcSC ; HC4)` and
  LSmear. Ibex 2.8 and develop have no LSmear for a plain system, so with
  `Makefile` it bisects with SmearSumRelative; `Makefile.inplace` uses LSmear
  as in the paper. `--lin=art` uses `DUAL` with 1e-9; `--lin=compo` uses `ADAPTIVE`
  with 1e-6.
- `sc_opt`: HC4 and ACID, then `FixPoint(CtcSC ; HC4)` with ratio 0.2, then
  LSmear and beam search. It uses `PRIMAL` with 1e-9.

## Validation against the research code

Built against Ibex develop (commit 0a0f9796d), `sc_opt` gives the same tree as
the in-place build (ex2_1_9: same cells, passes, LPs and pivots); `sc_solve`
differs only by the bisector (see above).

The comparison is with the `examples/dfb` binaries in their SC base
configuration (rung R3 of `results/escalera_2026-11-06/escf.sh`), with
60 s runs and seeds 1 and 2.

| sample | identical trees | CPU time, plugin / research (quotient of sums) |
|---|---|---|
| 25 optimization problems, affine relaxation (`sc_opt`) | 50/50 | 0.96 |
| 24 systems, affine relaxation (`sc_solve`) | 44/48 | 0.99 |
| same 24 systems, plugin emulating the first-pass quirk below | 48/48 | 0.98 |
| 13 systems, `compo`, ADAPTIVE (not reproducible) | 13/26 | 0.96, with cells × 1.00 |

**One difference with the research code.** In the research code, the
system solver asks for its order of bounds through `DFBH_ORDEN=costo_o_cerca`,
which resolves to "nearest" in the revised engine. That choice is read from
the engine *before* the first relaxation is loaded, so the very first pass of
a search uses the *farthest* bound first. `CtcSC` always uses the nearest
bound, as the paper describes. On the sample above, this changes the tree in
2 of the 24 systems (BroydenTri-0100: 13 cells instead of 11; Geneig-icse:
419/397 instead of 285/393). A build that reproduces the quirk gives the same
trees in all 48 runs. The optimizer is not affected.
