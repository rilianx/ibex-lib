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

**Lo que falta es el origen de los NaN en `SCSimplex`.** Posibles causas: un
overflow en `x_B`, en el ratio test o en la actualización de λ, con cotas de
±8e307. `CtcSC::nb_nan_certificates()` cuenta los LPs cuyo certificado se
descartó por NaN, y `sc_solve` lo imprime en la línea `[sc] ... nan=`.

## Reproducir

```
cd examples/sc_plugin && make PKG_CONFIG_PATH=<prefix>/share/pkgconfig
g++ -O3 -frounding-math -std=c++14 -I. -Iaffine casos/nan_lambda/sc_una.cpp \
    ibex_CtcSC.o ibex_SCSimplex.o ibex_SCBasis.o affine/*.o -o sc_una \
    $(pkg-config --cflags --libs ibex)
cd casos/nan_lambda
../../sc_una preprocesado_nsc.mbx caja_art.txt art sc
#   art sc: no vacia (1 variables contraidas), LPs=19, certificados con NaN=2
../../sc_una preprocesado_nsc.mbx caja_xt.txt xt sc
#   xt sc: no vacia (1 variables contraidas), LPs=28, certificados con NaN=0
../../sc_solve preprocesado_nsc.mbx -t 60
#   [sc] passes=5360 LPs=172257 pivots=1583381 nan=3
#   2293 22.2 0 6 0 0 0
```

Antes de `b11a701`, las dos primeras líneas decían `VACIA`. `ph`
(`CtcPolytopeHull`) en lugar de `sc` no vacía ninguna de las dos cajas.

Con el Ibex de `develop` el `ibex.pc` pide `-llibsoplex`: hay que cambiarlo a
`-lsoplex`.
