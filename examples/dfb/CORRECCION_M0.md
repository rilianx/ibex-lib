# Corrección M0: solidez del contractor DFB

Registro de la corrección de la mejora **M0** del [plan](PLAN_MEJORA_DFB.md) y de
la revisión del diagnóstico que la motivaba. Punto de partida:
[HALLAZGOS_ANALISIS.md](HALLAZGOS_ANALISIS.md), donde DFB declaraba vacías 89 de
139 instancias y en 46 de ellas `ibexsolve` certificaba soluciones dentro de la
caja.

**Resultado en una línea: la causa raíz no era el ratio test, sino el recorte de
las cotas de las variables `b` a ±1e50 en `CtcDFBPropag::linearize`.**

---

## 1. Lo que M0 arregló, y lo que no

M0 partía de un diagnóstico correcto en su descripción del mecanismo —
`calculateAlpha` devolvía `i == -1` por tres razones distintas y `contract`
trataba las tres como prueba de infactibilidad — pero incompleto en su
atribución de la causa.

Tras implementar M0, el barrido pasó de 89 a **10** instancias vaciadas. Las 10
restantes eran todas `Brown-*` y todas vaciaban por la vía de la cota. Y
Brown-10 tiene solución certificada en x = (1, …, 1). Es decir: **M0 no bastaba**.

La verificación que lo destapó fue el invariante del §3.2 del plan («la solución
conocida debe seguir en la caja contraída»), que hasta entonces no se había
ejecutado. Implementado en [diag_dfb.cpp](diag_dfb.cpp) con un testigo obtenido
de `ibexsolve -s`.

## 2. Cómo se localizó la causa

Cadena de descartes, toda reproducible con `./diag_dfb <instancia> <testigo>`:

1. **No son los contractores individuales.** Ninguno de los 20 de Brown-10
   pierde el testigo por sí solo; la propagación completa sí. El problema estaba
   en el estado acumulado.
2. **No es el vaciado.** Desactivando la prueba de vacío por la cota
   (`CtcDFB::prove_empty_by_bound = false`), la propagación seguía perdiendo el
   testigo, ahora por contracción. El vaciado era el síntoma visible, no la
   enfermedad.
3. **El paso exacto**: pivote 85, contractor `k=2` cota inferior, que contrae
   `x[2]` de `[-4e7, 4e7]` a `[1.5714…, 4e7]`, excluyendo el valor 1 del
   testigo. El ancho acumulado de `gamma` era 2.7e-14: **no** era degradación
   numérica.
4. **No es `gamma`.** `gamma · w = [-5.3e-14, 5.1e-14]`, contiene 0, así que la
   fila era una relación válida. Esto descartó la hipótesis de que el forzado de
   entradas a 0 y 1 exactos tras cada pivote fuera el paso no riguroso. De
   hecho, la matriz linealizada resulta **puntual**: 0 % de entradas con
   diámetro > 0, de modo que las operaciones de fila son exactas salvo redondeo.
   (El experimento que parecía confirmar la hipótesis del forzado estaba
   confundido: al quitar los forzados cambia la secuencia de pivotes.)
5. **Son las cotas de `b`.** Las filas 18 y 19 de `refA` valen ≈ −2e72 y −4e72
   en el testigo, mientras sus cotas `b` estaban recortadas a `[-1e50, 1e50]`.
   El linealizador no tiene culpa: sus 30 filas se satisfacen en el testigo, 0
   violaciones.

## 3. La causa

En `CtcDFBPropag::linearize`:

```cpp
x[nb_var+i] = lhs_rhs[nb_var+i];
if (x[nb_var+i].lb() < -1e50) x[nb_var+i] = Interval(-1e50, x[nb_var+i].ub());
if (x[nb_var+i].ub() >  1e50) x[nb_var+i] = Interval(x[nb_var+i].lb(), 1e50);
```

El recorte **aprieta** una cota no acotada hasta un valor finito arbitrario, que
es la dirección prohibida en aritmética de intervalos: una cota más estrecha que
la verdad excluye soluciones. Con las linealizaciones de `Brown-*` sobre cajas
de radio 1e9 los coeficientes son del orden de 1e72, así que el valor verdadero
de `b_i` cae muy por fuera de 1e50 y la solución se perdía.

## 4. El arreglo

`b_i` está definida por `b_i = fila_i · x`, así que su rango válido es la
evaluación por intervalos de la fila sobre la caja, intersectada con la cota que
reporta el solver LP. Las dos son válidas para toda solución, luego su
intersección también lo es — y es finita siempre que la caja lo sea, sin
constantes arbitrarias:

```cpp
Interval enclosure(0.0);
for (int j=0; j<nb_var; j++)
    enclosure += Interval(rows[nb_var+i][j]) * box[j];

Interval bi = lhs_rhs[nb_var+i] & enclosure;
x[nb_var+i] = bi.is_empty() ? enclosure : bi;
```

Si la intersección fuera vacía, eso probaría que la caja no contiene soluciones;
por prudencia se cae al enclosure (siempre válido) en vez de propagar vacuidad
desde la linealización. Aprovecharlo como detección de infactibilidad es una
mejora posible, no hecha.

El mismo recorte estaba duplicado en `CtcDFBManager::linearize`
([ibex_CtcDFBManager.h](ibex_CtcDFBManager.h)) y se corrigió igual.

Efecto colateral favorable: la intersección con el enclosure sólo puede
**apretar** las cotas de `b` respecto de lo que reporta el solver, así que además
de ser sólida puede contraer más.

## 5. Cambios en el código

### `ibex_CtcDualFeasibleBounding.h` / `.cpp` (M0)

1. **`calculateAlpha` clasifica su resultado** en `ratio_test_status`:
   `BLOCKING_ROW`, `NO_CANDIDATE` (ninguna fila puede bloquear: paso dual no
   acotado) o `INCONCLUSIVE` (hay filas cuyo signo no se puede determinar:
   `A[j][i]` contiene el 0, o el cociente contiene el 0). Se eliminó la
   tolerancia **absoluta** `alpha.ub() <= -1e-5`, que descartaba candidatos
   legítimos de una cantidad cuya escala depende del problema; la candidatura es
   ahora el signo riguroso `alpha.ub() < 0`, sin umbral de magnitud.
2. **Sólo `NO_CANDIDATE` vacía la caja.** Con `INCONCLUSIVE` se detiene el
   pivoteo (`state = FINAL`) y se conserva la caja.
3. **Detección de vacío por la vía rigurosa.** En `gaussSeidel`, si la cota
   deducida de `gamma·x = 0` no intersecta a `x[k]`, eso sí prueba que no hay
   solución. Esa rama antes no hacía nada. Hay que salir en cuanto la caja queda
   vacía: seguir pivoteando sobre una caja vacía leía cotas inexistentes (tres
   sitios de llamada). Se puede desactivar con
   `CtcDFB::prove_empty_by_bound = false`, útil para aislar efectos.
4. **La excepción de `makeColumnIdentity` ya no escapa.** `init()` y
   `regenerateA()` la atrapan y marcan el contractor inutilizable vía `init_ok`;
   `contract()` devuelve la caja intacta. Antes abortaba el programa (por
   ejemplo `ibexopt --filtering=dfb` en `DiscreteBoundary-0200` y `Neveu1`).
5. **La rama `INITIAL`** devolvía la caja con la dimensión extendida (x ∪ b) en
   lugar de la original.
6. **Contadores de diagnóstico** (`n_no_candidate`, `n_inconclusive`,
   `n_init_failed`, `n_empty_by_bound`, con `reset_counters()`) y el desempate
   tipo Harris (`harris_tie_break`), **apagado por defecto** para poder medir su
   efecto por separado. Corresponde a M3(a) y todavía no se evaluó.

### `ibex_CtcDFBPropag.cpp` y `ibex_CtcDFBManager.h`

El reemplazo del recorte a ±1e50 descrito en §4.

## 6. Verificación

### Suite de pruebas del proyecto

`TestingDualFeasibleBounding.cpp` trae 7 casos, de los que sólo
`standar_test()` está activo en `main()`. En [tests_dfb_all.cpp](tests_dfb_all.cpp)
(`make tests_dfb_all`) están los 7:

| caso | antes | después |
|---|---|---|
| STANDAR TEST 1 y 2 | ok | **misma caja** |
| NO SOLUTION TEST 1 | **cuelga** | **vacío detectado** |
| NO SOLUTION TEST 2 | no se llegaba a ejecutar | **vacío detectado** |
| UNBOUNDED TEST | no se llegaba | no vacía (correcto) |
| ILL-CONDITIONED, 15X20 | no se llegaban | corren |

El original se cuelga en `NO SOLUTION TEST 1` porque al vaciar no ponía
`state = FINAL` y el bucle `while (state != FINAL)` del test giraba para
siempre. La suite no pasaba antes tampoco.

### Campaña de testigos

Invariante: **un contractor no puede eliminar soluciones**, luego la caja
contraída debe seguir intersectando una caja-solución certificada.

Metodología: `ibexsolve -s` da una caja que **contiene** una solución; hay que
usar la caja, no su punto medio. Con el punto medio aparecían violaciones
espurias (Brown-05: 4 de 10 filas) porque el error de ~1e-16 del punto medio se
amplifica por coeficientes de ~1e6 en la linealización. Con la caja certificada:
0 de 10.

Resultados en [results/witness_campaign.csv](results/witness_campaign.csv):

| | |
|---|---|
| Testigos válidos obtenidos | 33 de 40 instancias |
| Propagación conserva el testigo | **33 de 33** |
| Filas con cotas `b` que excluyan el testigo | **0** |

Las 7 sin dato son limitaciones del arnés, no de DFB: 5 `ibexsolve` no terminó
en 60 s (Osborne1, Synthesis, Redeco10, Dietmaier, Fredtest), 1 sin soluciones
en la caja (Brown-15) y 1 error de parseo (brown5b).

## 7. Barrido completo del banco, antes y después

Mismo protocolo del §3.1 del plan, sobre las 153 instancias
([results/baseline_153_post_M0.csv](results/baseline_153_post_M0.csv), resumen
en [results/baseline_153_post_M0_resumen.txt](results/baseline_153_post_M0_resumen.txt)).

| | antes | después |
|---|---|---|
| Instancias que terminan | 139 | **147** |
| Timeout / excepción | 14 | 6 |
| **Cajas vaciadas** | **89** | **2** |
| Vaciados erróneos confirmados | 46 | **0** |
| Grupo A (DFB ≥ PolytopeHull) | 2 | **56** |
| Grupo B (DFB < PolytopeHull) | 35 | 41 |
| Grupo C (ninguno contrae) | 13 | 48 |
| Pivotes DFB / iteraciones SoPlex | 6975 / 3001 (2.3×) | 20474 / 12430 (1.6×) |
| `t_dfb / t_ph` mediana | 1.29 | 2.83 |

Los dos vaciados que quedan son **correctos**: `ibexsolve` reporta
`infeasible problem` para `Prolog` y `Prolog-icse`, o sea que esas cajas
efectivamente no contienen soluciones.

El grupo A pasó de 2 a 56 de 145: DFB ahora contrae al menos tanto como
PolytopeHull en el 39 % del banco. El caso más notorio es la familia
`Discrete-Integralf2-*`, que antes gastaba cientos de pivotes con ganancia
**0.00 %** y ahora obtiene **96.9 %** contra 5–50 % de PolytopeHull, con menos
pivotes que iteraciones de simplex (por ejemplo `Discrete-Integralf2-21`: 142
pivotes y 96.9 % frente a 276 iteraciones y 16.5 %).

El costo subió, y era esperable: antes DFB abandonaba temprano con un vaciado
falso. Ahora hace trabajo real, y por eso los pivotes totales y los tiempos
crecen mientras la razón pivotes/iteraciones **mejora** de 2.3× a 1.6×.

El grupo C creció de 13 a 48 porque muchas instancias que antes se vaciaban
ahora simplemente no contraen: pasan de un resultado erróneo a uno inútil pero
correcto. Ese grupo es el material de trabajo de M2 (pricing) y M3.
