# Case: non-finite multipliers in SCSimplex

`original.mbx` is a `func` instance of the ScriptDAG generator (n = 8, P = 12,
seed 3, domains `[-2, 2]`) with 6 solutions. `preprocesado_nsc.mbx` is an
equivalent rewriting produced by a CSE/NSC preprocessor: auxiliary variables
for the n-ary sums and the shared subexpressions, plus the redundant rows of a
Gauss–Jordan conditioner. `soluciones_completas.txt` has the 6 solutions with
all their variables (original and auxiliary), in the order of declaration.

## What happened

`sc_solve preprocesado_nsc.mbx` used to end INFEASIBLE (31 cells with
`--lin=art`, 1053 with `--lin=xt`): `CtcSC` emptied boxes that contain the
solutions. `caja_art.txt` and `caja_xt.txt` are those boxes, just before the
wrong contraction, in hexadecimal (`%a`), one variable per line. In both,
`_y2 = _w0^(-1)` with `_w0 ∈ [1.2e-308, 52.8]` (affine) or `[1e-153, 52.8]`
(XTaylor).

Three defects, all fixed:

- **Empty sum taken as a certificate.** `Interval(NaN)` is empty, so a NaN
  multiplier made `γ·z` empty, and `!sum.contains(0.0)` read it as a proof of
  emptiness. An empty sum no longer certifies anything.
- **Non-rigorous merging of nearly proportional rows** (`caja_xt.txt`). Rows
  with `b = c·a` within a relative tolerance of 1e-9 were merged as if exactly
  proportional; the residual `(b − c·a)·x` is now enclosed over the box.
- **The origin of the NaN** (`caja_art.txt`). The affine linearization of
  `_y2 = _w0^(-1)` has the slope `1/(lb·ub) ≈ 1.54e306`, next to the
  coefficient 1 of `_y2` and the −1 of the slack. When `_w0` enters the basis,
  the update `xB[i] -= alpha[i]*delta` overflows to +inf; that row leaves at the
  next step, and its BTRAN on LU factors of ~1e306 gives inf − inf = NaN. With
  a NaN pivot row the ratio test finds no entering variable, and the simplex
  declares infeasibility with the NaN row as Farkas ray. Such a row cannot be
  pivoted in double precision: `CtcSC` now drops every row with a coefficient
  above 1e100 (`nb_dropped_rows()`, `drop=` in `sc_solve`), which only
  enlarges the polytope.

## Reproduce

```
make PKG_CONFIG_PATH=<prefix>/share/pkgconfig
g++ -O3 -frounding-math -std=c++14 -I. -Iaffine casos/nan_lambda/sc_una.cpp \
    ibex_CtcSC.o ibex_SCSimplex.o ibex_SCBasis.o affine/*.o -o sc_una \
    $(pkg-config --cflags --libs ibex)
cd casos/nan_lambda
../../sc_una preprocesado_nsc.mbx caja_art.txt art sc
#   art sc: not empty (10 variables contracted), LPs=25, NaN certificates=0, dropped rows=1
../../sc_una preprocesado_nsc.mbx caja_xt.txt xt sc
#   xt sc: not empty (0 variables contracted), LPs=27, NaN certificates=0, dropped rows=1
../../sc_solve preprocesado_nsc.mbx -t 60
#   [sc] passes=5361 LPs=173330 pivots=1566729 nan=0 drop=681
#   2293 <time> 0 6 0 0 0
```

With `ph` (`CtcPolytopeHull`) instead of `sc`, neither box is emptied.
