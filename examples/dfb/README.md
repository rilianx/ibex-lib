# Ejemplo `example_dfb` — DFB vs. PolytopeHull sobre sistemas no lineales

Programa de comparación experimental entre el contractor **DFB** (*Dual Feasible
Bounding*, una variante de simplex sobre intervalos) y los contractores clásicos
de Ibex (HC4, `CtcPolytopeHull` con SoPlex). Carga un sistema de restricciones en
formato Minibex, construye **una copia independiente del sistema por método**,
aplica cada método sobre su propia copia y reporta la caja contraída junto con el
número de iteraciones consumidas.

Código fuente: [dfb_non_linear_system.cpp](dfb_non_linear_system.cpp) →
ejecutable `example_dfb` (así lo define el [Makefile](Makefile)).

El punto de este diseño es que **cada método parte exactamente de la misma caja
inicial**, sin contracciones previas de otros contractores que contaminen la
comparación. Por eso hay diez `System`, diez `ExtendedSystem` y diez
`LinearizerXTaylor` paralelos, uno por variante.

## Uso

```
./example_dfb <archivo_sistema>
```

- `<archivo_sistema>`: sistema en formato Minibex (`.bch` o `.txt`). En este
  directorio hay varios casos: [ex_oct25_01.bch](ex_oct25_01.bch),
  [ex_oct25_paper.bch](ex_oct25_paper.bch), [caprasse.txt](caprasse.txt),
  [bellido.txt](bellido.txt), [geneig.txt](geneig.txt), y más en
  [data_tests/](data_tests/).
- No hay más argumentos: a diferencia de [example_dfb.cpp](example_dfb.cpp), este
  programa **no** lee `argv[2]`; ejecuta el conjunto fijo de métodos activos y
  cualquier argumento extra se ignora silenciosamente.

Ejemplo:

```bash
./example_dfb ex_oct25_01.bch
```

Salida útil (al final, tras las trazas de depuración):

```
 xFinal DFB_ONLY =
([-9.0050016653701, 1.1586525026982] ; [-5.4383177468756, 0.92] ; ...)
ITERS DFB ONLY: 37

xFinal PH =
([-9.0050016653701, 1.1586525026982] ; [-5.4383177468752, 0.92] ; ...)
ITERS PH ONLY: 17
```

Es decir: DFB puro y PolytopeHull alcanzan (en este caso) prácticamente la misma
caja, y lo que se compara es el **costo** — 37 pivotes DFB contra 17 iteraciones
de simplex de SoPlex.

## Qué construye el programa

Para cada variante *v*:

1. `System sys_v(filename)` — copia propia del sistema.
2. `ExtendedSystem sys2_v(sys_v, 0.0)` — `eqeps = 0.0`, las igualdades se
   mantienen exactas.
3. `LinearizerXTaylor lr_v(sys2_v, RELAX, RANDOM, HANSEN)` — relajación X-Taylor
   con esquina aleatoria y matriz de Hansen. Es la fuente del sistema lineal
   `A·x = 0` que consumen tanto DFB como PolytopeHull.

Las diez variantes previstas son `hc4`, `dfb_hc4`, `dfb`, `dfb2`, `dfb_only`,
`ph`, `it_ph`, `it_ph2`, `acid_dfb` y `acid_ph`.

## Métodos: qué está activo y qué no

| Variante | Contractor | Estado | Se imprime |
|---|---|---|---|
| `dfb_only` | `CtcDFBPropag(…, 0.01, true, false, true)` — DFB puro (`only_dfb`), sin HC4 en la cola | **activo** | **sí** (caja + `count_dfb`) |
| `ph` | `CtcPolytopeHull(lr_ph)`, `n_soplex_iterations = 0` | **activo** | **sí** (caja + `n_soplex_iterations`) |
| `hc4` | `CtcDFBPropag(…, 0.001, true, true)` — HC4 puro (`only_hc4`) | activo | no (`cout` comentado) |
| `dfb_hc4` | `CtcDFBPropag(…, 0.001, true, true)` | activo | no |
| `dfb` | `CtcDFBPropag(…, 0.01)` — DFB + HC4 en propagación conjunta | activo | no |
| `dfb2` | `CtcDFBPropag(…, 0.01)` y luego `CtcDFBPropag(…, 0.001)` sobre la misma caja | activo | no |
| `it_ph` | PolytopeHull iterado con `optimizer()` (sin relinealizar) | **comentado** | no |
| `it_ph2` | PolytopeHull iterado con `contract()` (relineariza) | **comentado** | no |
| `acid_dfb` | `CtcDFBManager(dfb_acid, CtcAcid(dfb_acid), lr)` | **comentado** | no |
| `acid_ph` | `CtcAcid(hc4)` + `CtcPolytopeHull` (dos rondas) | **comentado** | no |

Consecuencia práctica: `hc4`, `dfb_hc4`, `dfb` y `dfb2` **sí se ejecutan** pero
sus resultados no se muestran, de modo que consumen tiempo sin aportar salida.
Si se quiere medir tiempos de `dfb_only` y `ph` limpiamente, hay que comentar
esas cuatro contracciones o descomentar sus `cout` correspondientes.

## Métricas reportadas

- `ITERS DFB ONLY` = `dfb_only.count_dfb`: suma de las iteraciones internas
  (`CtcDFB::iters`) de todas las llamadas a contractores DFB durante la
  propagación. Cada iteración es un **pivote**: `largestImpact()` elige la
  columna entrante, `calculateAlpha()` hace el *ratio test* y elige la fila
  saliente, `A[0] += alpha*A[j]` con `makeColumnIdentity()` aplica la
  eliminación, y `gaussSeidel()` extrae la cota de la variable `k`. En
  `CtcDFBPropag` los `CtcDFB` se construyen con `max_iters = 1`, es decir un
  pivote por invocación.
- `ITERS PH ONLY` = `ctc_ph.n_soplex_iterations`: iteraciones de simplex
  acumuladas que reporta SoPlex (`mysoplex->numIterations()`) al resolver los
  dos LP por variable.

Ambas cifras son comparables como *conteo de pivotes*, y son la base natural
para evaluar mejoras al pivoteo de DFB.

## Salida verbosa

`disableCout()` / `enableCout()` redirigen `std::cout` a `/dev/null` para
silenciar las trazas de los contractores. `enableCout()` se llama **antes** de
`dfb_only.contract(...)`, así que las dos variantes que interesan corren con
`cout` activo y su ruido sí aparece:

- la lista de números sueltos (`0`, `0`, …) viene del `cout << iters << endl` de
  [ibex_CtcDualFeasibleBounding.cpp:156](ibex_CtcDualFeasibleBounding.cpp#L156)
  al terminar cada llamada a `CtcDFB::contract`;
- los `{` `}` los imprime `CtcPolytopeHull`.

Al inicio de `main` hay un bloque comentado que redirige `std::cout` a un archivo
(`resultado_oct25_06.txt`); de ahí provienen los `resultado_oct25_*.txt` de este
directorio, que son salidas guardadas de corridas anteriores.

## Detalles a tener en cuenta

- `dfb_hc4` se anuncia en el comentario como «DFB + HC4», pero se construye con
  `only_hc4 = true`, o sea **sólo HC4**; es idéntico a `hc4`. La variante que
  realmente combina DFB y HC4 es `dfb`.
- `hc4_it_ph2` se construye sobre `sys2_hc4` / `lr_hc4` y no sobre `sys2_it_ph2`;
  sólo importa si se reactiva el bloque `it_ph2`.
- En el bloque comentado de `it_ph`/`it_ph2` se usa `ctc_ph`, que se declara
  bastante más abajo en la función: al descomentarlo no compilará sin mover esa
  declaración (o usar `ctc_it_ph`/`ctc_it_ph2`).
- `dfb2` aquí son dos pasadas con `ratio` 0.01 y 0.001; no activa
  `CtcDFBPropag::b_contraction` (a diferencia del modo `dfb2` de
  [example_dfb.cpp](example_dfb.cpp)).

## Parámetros de `CtcDFBPropag`

```cpp
CtcDFBPropag(ExtendedSystem& sys, Linearizer& lr, double ratio=0.1,
             bool stand_alone=true, bool only_hc4=false, bool only_dfb=false);
```

- `ratio`: umbral de reducción relativa (`ratiodelta`) por debajo del cual una
  contracción no vuelve a encolar la variable/restricción en la cola de
  prioridad de la propagación.
- `stand_alone`: si es `true`, el contractor lineariza e inicializa por sí mismo
  la matriz de referencia `refA`. Se pone en `false` cuando esa tarea la asume un
  `CtcDFBManager` (caso `acid_dfb`).
- `only_hc4`: no crea contractores DFB (queda HC4 puro).
- `only_dfb`: no crea contractores HC4 (queda DFB puro).
- `CtcDFBPropag::b_contraction` (estático): además de las variables `x`, contrae
  las variables `b` del sistema linealizado.

## Compilación

```bash
make example_dfb     # o simplemente: make
```

Requiere la biblioteca Ibex compilada en `../../src`, junto con Gaol y SoPlex;
las rutas de `-I`/`-L` están fijadas en el [Makefile](Makefile). Los fuentes del
ejecutable son:

```make
example_dfb_SRCS = dfb_non_linear_system.cpp ibex_CtcDualFeasibleBounding.cpp \
                   ibex_Optimizer05Config.cpp ibex_CtcDFBPropag.cpp \
                   ../../src/contractor/ibex_CtcPolytopeHull.cpp
```

Para compilar en su lugar la variante de un método por ejecución
([example_dfb.cpp](example_dfb.cpp), que sí lee `argv[2]`), reemplace
`dfb_non_linear_system.cpp` por `example_dfb.cpp` en esa variable.

## Archivos relacionados

| Archivo | Rol |
|---|---|
| [dfb_non_linear_system.cpp](dfb_non_linear_system.cpp) | **Este ejemplo**: todos los métodos activos en una corrida, sistemas no lineales. |
| [dfb_linear_system.cpp](dfb_linear_system.cpp) | Variante equivalente para sistemas lineales. |
| [example_dfb.cpp](example_dfb.cpp) | Variante que ejecuta un solo método, elegido por `argv[2]`, tras una contracción previa HC4. |
| [ibex_CtcDualFeasibleBounding.h](ibex_CtcDualFeasibleBounding.h) / [.cpp](ibex_CtcDualFeasibleBounding.cpp) | `CtcDFB`: el contractor DFB para una variable `k` en `A·x = 0` (pricing, ratio test, pivoteo, Gauss-Seidel). |
| [ibex_CtcDFBPropag.h](ibex_CtcDFBPropag.h) / [.cpp](ibex_CtcDFBPropag.cpp) | Propagación con cola de prioridad que combina contractores DFB y HC4. |
| [ibex_CtcDFBManager.h](ibex_CtcDFBManager.h) | Envoltorio que lineariza e inicializa los DFB para usarlos dentro de otro contractor (p. ej. ACID). |
| [ibex_Optimizer05Config.h](ibex_Optimizer05Config.h) / [.cpp](ibex_Optimizer05Config.cpp) | Configuración de optimizador que integra DFB, usada por `ibexopt`. |
| [ibexopt.cpp](ibexopt.cpp) | Optimizador global con DFB como opción de filtrado (`--filtering=dfb`). |
| [TestingDualFeasibleBounding.cpp](TestingDualFeasibleBounding.cpp) | Pruebas del contractor DFB. |
| `resultado_oct25_*.txt` | Salidas guardadas de corridas anteriores. |
