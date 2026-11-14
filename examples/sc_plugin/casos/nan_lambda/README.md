# Caso: SCSimplex devuelve multiplicadores NaN

Origen: `ibex_cse` (preprocesamiento CSE/NSC a Minibex). `original.mbx` es una
instancia `func` del generador de ScriptDAG (n = 8, P = 12, semilla 3, dominios
`[-2, 2]`), con 6 soluciones. `preprocesado_nsc.mbx` es la reescritura
equivalente, con auxiliares para las sumas n-arias y las subexpresiones
compartidas, más las filas redundantes del conditioner de Gauss–Jordan.
`soluciones_completas.txt` tiene las 6 soluciones con todas sus variables (las
originales y las auxiliares), en el orden de declaración del archivo.

## Qué pasaba

Antes de `b11a701`, `sc_solve preprocesado_nsc.mbx` (`--lin=art`) terminaba
INFEASIBLE en 31 celdas, y `--lin=xt` en 1053. `CtcSC` vaciaba cajas que
contienen las soluciones:

- `caja_art.txt` (relajación afín): **`SCSimplex` devuelve λ con NaN**.
  `Interval(NaN)` es vacío, `γ·z` queda vacío, y `!sum.contains(0.0)` lo tomaba
  como prueba de vacío.
- `caja_xt.txt` (XTaylor): no hay NaN. El vacío venía de la fusión de filas
  proporcionales (tolerancia 1e-9) sin encerrar el residuo `(b − c·a)·x`.

Las dos cajas están en hexadecimal (`%a`), una variable por línea: son la caja
justo antes de la contracción que eliminaba las soluciones. En ambas,
`_y2 = _w0^(-1)` con `_w0 ∈ [1.2e-308, 52.8]` (afín) o `[1e-153, 52.8]` (XTaylor).
Por eso `_y2` llega a 8e307 o 1e153, y la relajación tiene coeficientes del orden
de `1/_w0²`.

## Estado

`b11a701` (CtcSC) y `b9d6787` (CtcDFBHull) impiden que un certificado con NaN
vacíe la caja y hacen rigurosa la fusión de filas. Con eso se recuperan las 6
soluciones (2293 celdas, 22 s, `--lin=art`).

**Origen de los NaN (2026-11-14).** Una sola fila de la relajación afín: la
linealización de `_y2 = _w0^(-1)`, con pendiente `1/(lb·ub) ≈ 1.54e306` junto al
coeficiente 1 de `_y2` y al −1 de la holgura. En el LP 2 entra `_w0`, su columna
en la base trae una entrada de 1.54e306, la actualización incremental
`xB[i] -= alpha_piv[i]*delta` se desborda (x_B = +inf), esa fila sale en el paso
siguiente y su BTRAN sobre factores de ~1e306 da inf − inf = NaN en ρ. Con la
fila pivote en NaN el ratio test no encuentra entrante y el simplex declara
infactible con λ = ρ: el «rayo» NaN. La fila no se puede pivotear en doble
precisión, así que `CtcSC` (y `CtcDFBHull`) descartan las filas con un
coeficiente mayor que 1e100 (`nb_dropped_rows()`, `drop=` en `sc_solve`).
Descartar una fila solo agranda el poliedro. Ahora `nan=0`.

## Reproducir

```
cd examples/sc_plugin && make PKG_CONFIG_PATH=<prefix>/share/pkgconfig
g++ -O3 -frounding-math -std=c++14 -I. -Iaffine casos/nan_lambda/sc_una.cpp \
    ibex_CtcSC.o ibex_SCSimplex.o ibex_SCBasis.o affine/*.o -o sc_una \
    $(pkg-config --cflags --libs ibex)
cd casos/nan_lambda
../../sc_una preprocesado_nsc.mbx caja_art.txt art sc
#   art sc: no vacia (10 variables contraidas), LPs=25, certificados con NaN=0
#   (antes del filtro: 1 variable contraida, LPs=19, NaN=2)
../../sc_una preprocesado_nsc.mbx caja_xt.txt xt sc
#   xt sc: no vacia (0 variables contraidas), LPs=27, certificados con NaN=0
../../sc_solve preprocesado_nsc.mbx -t 60
#   [sc] passes=5361 LPs=173330 pivots=1566729 nan=0 drop=681
#   2293 <tiempo> 0 6 0 0 0   (antes: nan=3; mismo arbol)
```

Antes de `b11a701`, las dos primeras líneas decían `VACIA`. `ph`
(`CtcPolytopeHull`) en lugar de `sc` no vacía ninguna de las dos cajas.

Con el Ibex de `develop` el `ibex.pc` pide `-llibsoplex`: hay que cambiarlo a
`-lsoplex`.
