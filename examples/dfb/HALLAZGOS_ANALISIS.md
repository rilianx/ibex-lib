# Hallazgos del análisis de línea base de DFB

Resultados de los experimentos del §3.1 y §3.4 de [PLAN_MEJORA_DFB.md](PLAN_MEJORA_DFB.md),
sobre el código **tal como estaba antes de M0**, sin modificarlo.

> **Nota de estado.** Este documento describe la línea base. Las correcciones que
> siguieron (M0 y el recorte de las cotas `b`) están en
> [CORRECCION_M0.md](CORRECCION_M0.md), junto con la revisión del diagnóstico:
> la causa raíz del vaciado de cajas **no** era el ratio test, sino el recorte
> de las cotas de las variables `b` a ±1e50 en `CtcDFBPropag::linearize`.

- Instancias: las 153 de [data_tests/](data_tests/), una por proceso, límite de
  60 s y 8 GB de memoria virtual.
- Herramientas construidas para esto: [bench_dfb.cpp](bench_dfb.cpp) (una línea
  CSV por instancia) y [diag_dfb.cpp](diag_dfb.cpp) (diagnóstico por contractor
  y reproducción manual del primer pivote). Targets `make bench_dfb` y
  `make diag_dfb`.
- Datos crudos y scripts en [results/](results/).
- Referencia de comparación: `CtcPolytopeHull` con SoPlex, y `ibexsolve` como
  árbitro independiente de si una caja contiene soluciones.

---

## 1. Resultado principal: DFB vacía cajas que contienen soluciones

**De las 139 instancias que terminan, DFB declara la caja vacía en 89 (64 %).
PolytopeHull no vacía ninguna.**

Vaciar la caja es una afirmación fuerte: es una prueba de que no hay solución
ahí. La relajación X-Taylor contiene todas las soluciones de la caja, así que si
la caja contiene una solución, esa prueba es necesariamente falsa.

Se pasó `ibexsolve` (límite 30 s) por **las 89 instancias que DFB vacía**
([results/soundness_89_vaciadas.csv](results/soundness_89_vaciadas.csv)):

| resultado de `ibexsolve` | instancias |
|---|---|
| **caja con soluciones certificadas** → el vaciado de DFB es **erróneo** | **46** |
| caja sin soluciones → el vaciado sería legítimo | **0** |
| no terminó en 30 s → sin veredicto | 43 |

**46 casos confirmados de vaciado erróneo y ningún contraejemplo.** Los 43 sin
veredicto no son evidencia en contra: simplemente no se resolvieron a tiempo.

Algunas de las cajas vaciadas contienen muchísimas soluciones:

| instancia | cajas-solución certificadas | DFB la vacía |
|---|---|---|
| Brent-10 | **289** | sí |
| Brent-8 | 152 | sí |
| I5 / I5bis | 30 | sí |
| CountercurrentReactors2-12 | 28 | sí |
| Redeco9 | 16 | sí |
| Bellido, Redeco8 | 8 | sí |
| Brown-05 | 3 | sí |
| Brown-10 | 2 | sí |
| BroydenBanded-010, BroydenBanded-020, DiscreteBoundary-0020 | 1 | sí |

Es reproducible con los binarios del repositorio, sin herramientas nuevas:

```
$ ./example_dfb data_tests/BroydenBanded-010.bch
empty vector
ITERS DFB ONLY: 0
xFinal PH = ([-100, 66.99...] ; ... )     <- PH sí contrae, no vacía
```

**Conclusión: DFB no es sólido en su forma actual.** No es pérdida de precisión
ni de eficiencia: emite pruebas de infactibilidad falsas en la mayoría del banco.

## 2. Mecanismo: `i == -1` en el ratio test se interpreta como infactibilidad

Diagnóstico completo en [results/diag_BroydenBanded-010.txt](results/diag_BroydenBanded-010.txt).
En BroydenBanded-010, 18 de los 20 contractores `CtcDFB` vacían la caja **en su
primer pivote**. Reproduciendo ese pivote a mano:

```
largestImpact : fila entrante j=1, delta=[1.87e10, 1.87e10]
calculateAlpha: fila saliente i=-1        <- no encuentra fila de bloqueo

candidatos del ratio test (se exige alpha.ub() <= -1e-5):
  i=0   gamma_i=[-746.27, -746.27]  A[j][i]=[-1.119e8, -1.119e8]  alpha=[-6.6666e-06, -6.6666e-06]
  i=10  gamma_i=[ 0.00497, 0.00497] A[j][i]=[ 746.27,  746.27]    alpha=[-6.6666e-06, -6.6666e-06]
```

Hay candidatos válidos, con `alpha ≈ −6.67e−06`. La tolerancia dura de
`calculateAlpha` ([ibex_CtcDualFeasibleBounding.cpp:308](ibex_CtcDualFeasibleBounding.cpp#L308))
exige `alpha.ub() <= -1e-5`, así que los descarta todos, `i` queda en `-1`, y
`CtcDFB::contract` ([:166](ibex_CtcDualFeasibleBounding.cpp#L166)) responde con
`x_new.set_empty()`.

El problema de fondo es que **`i == -1` mezcla tres situaciones distintas** y las
trata todas como prueba de vacío:

1. no existe candidato del signo correcto — el único caso donde un argumento de
   infactibilidad tendría sentido;
2. existen candidatos pero los rechaza la **tolerancia absoluta** `1e-5`;
3. los rechaza la **guarda de cero** `Aj[i].lb() != 0 && Aj[i].ub() != 0` (que
   descarta también intervalos legítimos como `[0, 5]`), o el cociente
   `gamma_i / A[j][i]` es un intervalo que contiene el 0 y por lo tanto no
   cumple `alpha.ub() <= -tol`.

En los casos 2 y 3 la conclusión correcta es «no puedo concluir, me detengo»,
nunca «la caja es vacía».

Además, el `alpha` diminuto de este ejemplo es puro efecto de escalado:
`gamma_0 = −746` contra `A[j][0] = −1.1e8`. Una tolerancia **absoluta** sobre una
cantidad cuya escala depende del problema no puede funcionar en general.

## 3. Dos variantes experimentales para acotar la causa

16 instancias × 5 repeticiones por variante ([results/variantes_base_tol_safe.csv](results/variantes_base_tol_safe.csv)).
Los resultados son idénticos en las 5 repeticiones, o sea que el vaciado es
determinista pese a la esquina aleatoria del linealizador.

| variante | instancias vaciadas (de 16) |
|---|---|
| base (tolerancia `1e-5`) | 9 |
| tolerancia bajada a `1e-12` | 3 |
| no vaciar cuando `i == -1` | 0 |

Bajar la tolerancia arregla 6 de los 9 casos. Las 3 restantes (Brown-05,
Brown-10, I5 — con 3, 2 y 30 soluciones certificadas) siguen vaciando por la
causa 3 de arriba, así que la tolerancia no es la única culpable.

Con la variante que simplemente no vacía, DFB **sí contrae**, y por primera vez
hay números comparables:

| instancia | ganancia DFB | ganancia PH |
|---|---|---|
| Katsura-12 | 37.5 % | 45.8 % |
| Rose | 24.4 % | 24.4 % |
| Brown-05 | 16.0 % | 0 % |
| Redeco8 | 10.7 % | 0 % |
| DiscreteBoundary-0020 | 8.2 % | 16.0 % |

Esa variante es **diagnóstica, no una propuesta**: elimina por completo la
capacidad de probar infactibilidad, que es una función legítima del contractor.

## 4. Excepciones no atrapadas y explosión de memoria

De las 153 instancias, 14 no terminan:

- **3 lanzan excepción**. Dos de ellas
  (`DiscreteBoundary-0200`, `Neveu1`) con `std::runtime_error: No hay ninguna
  fila con valor no nulo en la columna k`, lanzada desde
  `makeColumnIdentity` ([:355](ibex_CtcDualFeasibleBounding.cpp#L355)). Nadie la
  atrapa: `CtcDFBPropag::init_dfb_contractors` no tiene manejo de excepciones, de
  modo que **`ibexopt --filtering=dfb` aborta** en esas instancias en lugar de
  degradar a otro contractor. La tercera (`Eiger-1000`) es `std::bad_alloc`.
- **11 agotan los 60 s o los 8 GB** (`BroydenBanded-200`, `BroydenBanded-1000`,
  `BroydenTri-1000`, `DiscreteBoundary-1000`, `ExtendedFreud-1000`,
  `Ex14-2-3`, y 5 `Discrete-Integralf2-*`), mientras PolytopeHull resuelve todas
  las que sí terminan en menos de 0.09 s.

La causa estructural está en `CtcDFBPropag`: crea `2n` contractores `CtcDFB`
([ibex_CtcDFBPropag.cpp:17](ibex_CtcDFBPropag.cpp#L17)) y **cada uno guarda su
propia copia completa de la matriz** de intervalos `m × (n+m)`
([:150](ibex_CtcDFBPropag.cpp#L150)), o sea O(n·m·(n+m)) por caja.

## 5. Costo: DFB usa más pivotes y cada pivote es más caro

Sobre las instancias donde DFB no vacía:

- **Pivotes**: 6975 en DFB contra 3001 iteraciones de simplex en SoPlex (2.3× más).
- **Tiempo**: `t_dfb / t_ph` tiene mediana 1.29 y llega a **49.6×**
  (ExtendedFreud-0100). Máximos: DFB 1.51 s, PH 0.083 s.

Esto invalida la métrica que usa hoy el ejemplo. El caso `ex_oct25_01` sugería
«37 pivotes DFB contra 17 iteraciones de SoPlex», una desventaja aparentemente
moderada; en el banco completo la desventaja es de 2.3× en pivotes **y** cada
pivote de DFB cuesta O(m·(n+m)) operaciones de intervalo sobre un tableau denso
contra un paso de simplex revisado disperso en punto flotante.

Casos donde el gasto no compra nada:

| instancia | pivotes DFB | ganancia DFB | iters PH | ganancia PH |
|---|---|---|---|---|
| Discrete-Integralf2-21 | 800 | 0.00 % | 276 | 16.5 % |
| Discrete-Integralf2-18 | 731 | 0.00 % | 236 | 5.5 % |
| Discrete-Integralf2-16 | 519 | 0.00 % | 162 | 9.3 % |
| ExtendedFreud-0100 | 100 | 0.00 % | 50 | 8.3 % |

Cientos de pivotes con ganancia exactamente nula: hay un problema de criterio de
terminación y de pricing, no sólo de costo por pivote.

## 6. Clasificación A/B/C

Con el código actual, sobre las 50 instancias que terminan y que DFB no vacía:

| grupo | criterio | instancias |
|---|---|---|
| A | DFB contrae ≥ que PolytopeHull | **2** (Rose, Synthesis) |
| B | DFB contrae < que PolytopeHull | **35** |
| C | ninguno contrae | 13 |

Las 89 instancias vaciadas quedan fuera de toda clasificación, porque su
resultado no es interpretable. Es decir: **el experimento del §3.1 del plan no
se puede completar con el código actual**, y esa es la razón por la que el
grupo A tiene sólo 2 instancias.

Para el objetivo declarado —que DFB **reemplace** a PolytopeHull— el grupo A es
la métrica que importa, y hoy es 2 de 50.

## 7. Densidad de la linealización: la vía de dispersión se justifica

Densidad de la matriz que produce `LinearizerXTaylor` (fracción de entradas no
nulas de `refA`):

| | mediana | p25 | p75 | máx |
|---|---|---|---|---|
| todas (139) | 0.190 | 0.080 | 0.344 | 0.750 |
| sólo n ≥ 50 (11) | **0.025** | | | 0.318 |

Las matrices grandes son dispersas (2.5 % de entradas no nulas en la mediana),
así que **M6 (dispersión) queda justificada** para instancias grandes — que son
justamente las que hoy agotan tiempo y memoria. En instancias chicas no aporta.

## 8. Nota metodológica: la esquina aleatoria (corregido)

`LinearizerXTaylor` usa `RANDOM` para elegir la esquina, así que DFB y
PolytopeHull reciben linealizaciones **distintas** de la misma caja, igual que
en [dfb_non_linear_system.cpp](dfb_non_linear_system.cpp).

Lo que se afirmaba aquí —que eso metía ruido entre corridas y obligaba a fijar
la semilla— **resultó falso**. Medido después sobre 5 repeticiones completas del
banco: 145 de 145 instancias dan exactamente el mismo resultado, porque cada
instancia corre en su propio proceso y el generador de Ibex se inicializa igual
en cada arranque. Ver [MEDICIONES_M2_M3.md](MEDICIONES_M2_M3.md) §1.

## 9. Qué implica para el plan

1. **Hay un M0 de corrección antes de todo lo demás.** Ninguna medición de
   pricing, ratio test o factorización significa nada mientras el contractor
   emita pruebas de vacío falsas en el 64 % del banco.
2. La separación de los tres casos de `i == -1` es, además, la primera pieza
   del ratio test tipo Harris de M3: el mismo código que decide «no puedo
   concluir» es el que debe elegir el candidato de mayor magnitud de pivote.
3. La clasificación A/B/C debe repetirse **después** de M0. La de la §6 sirve
   como registro del punto de partida, no como base para priorizar.
4. Los 2/50 del grupo A ponen el objetivo de reemplazar a PolytopeHull a mucha
   distancia: falta cerrar una brecha de contracción, no sólo una de costo.
