# Dónde está realmente el tiempo en el optimizador, y qué techo tiene cada mejora

El plan tenía la **factorización de la base** como paso 1, con el argumento de
que la inicialización de los `2n` tableaus es el 40.8 % del tiempo *del
contractor*. Este documento mide el techo real de esa mejora —y de las otras— en
el tiempo *del optimizador*, que es donde está el criterio incumplido.

Instrumentación: `DFB_PHASES=1` en `ibexopt` imprime una línea `PHASES,...` con
el tiempo total, el del contractor y su desglose.

---

## 1. El contractor no es la mayor parte del tiempo, y varía muchísimo

| instancia | filtrado | t total | t contractor | % contractor | % init | % pricing |
|---|---|---|---|---|---|---|
| ex14_2_4 | `dfb` | 5.92 s | 1.12 s | 19.0 % | 36.6 % | 30.6 % |
| ex6_2_12 | `dfb` | 12.71 s | 2.36 s | 18.6 % | 31.8 % | 26.9 % |
| ex2_1_6 | `dfb` | 1.45 s | 0.61 s | 42.0 % | 19.8 % | 28.4 % |
| alkyl | `dfb` | 1.36 s | 1.11 s | **81.2 %** | 13.4 % | 49.6 % |
| ex14_2_4 | `acid_dfb` | 8.83 s | 1.59 s | 18.0 % | 22.2 % | 36.7 % |
| ex6_2_12 | `acid_dfb` | 20.00 s | 1.92 s | 9.6 % | **7.9 %** | 28.6 % |
| ex2_1_6 | `acid_dfb` | 20.00 s | 2.00 s | 10.0 % | **4.8 %** | 31.7 % |
| alkyl | `acid_dfb` | 9.48 s | 8.71 s | **91.9 %** | **0.4 %** | 48.0 % |

## 2. El techo de la factorización es chico, y depende de la configuración

Eliminar **toda** la inicialización ahorraría:

| | techo sobre el tiempo total |
|---|---|
| `dfb`, ex14_2_4 | 19.0 % × 36.6 % = **7.0 %** |
| `dfb`, ex6_2_12 | 18.6 % × 31.8 % = **5.9 %** |
| `dfb`, alkyl | 81.2 % × 13.4 % = 10.9 % |
| **`acid_dfb`** | 18.0 % × 22.2 % = 4.0 %, y en tres instancias **< 1 %** |

Y la brecha a cerrar es **×1.31** en `dfb` y **×3.21** en `acid_dfb`. O sea que
**la factorización no puede cerrarla**: su techo es de un dígito porcentual.

Peor: en `acid_dfb` —la configuración más cerca del objetivo en árbol— la
inicialización es del 0.4 % al 22 %, porque el `CtcDFBManager` inicializa una
vez por llamada de ACID mientras el contractor se aplica muchas veces dentro.
Ahí la factorización ahorraría **casi nada**.

**Conclusión: la justificación de la factorización es la memoria, no el tiempo.**
Es lo único que resuelve las 6 instancias que no terminan por almacenar `2n`
tableaus de `m × (nb_var+2m)`. Como mejora de tiempo, su techo medido es 4–11 %.

## 3. La brecha con producción es de árbol, no de tiempo por nodo

| | celdas × producción | tiempo × producción |
|---|---|---|
| `dfb` | 2.48 | 1.31 |
| `acid_dfb` | 1.48 | 3.21 |

En `dfb` el tiempo crece mucho menos que el árbol: **por nodo ya es más barato
que la configuración de producción** (×1.31 de tiempo con ×2.48 de nodos implica
alrededor de la mitad de costo por nodo). Si `dfb` igualara el árbol de
producción, sería aproximadamente el doble de rápido que producción.

Es decir: **el criterio que falta se cierra reduciendo el árbol, no acelerando el
contractor.**

## 4. Hipótesis probada y rechazada: priorizar la variable objetivo

En un optimizador lo que poda es la cota inferior de la variable objetivo, y
`CtcPolytopeHull` la ataca en forma dirigida (`optimizer()`,
`set_contracted_vars`), mientras DFB reparte el esfuerzo entre las `2n` cotas por
igual. Parecía la explicación del ×2.48.

Se implementó (`CtcDFBPropag::goal_priority`, `DFB_GOAL_PRIO=1`: ese contractor
entra primero en la cola) y **no funciona**:

| instancia | base | con prioridad |
|---|---|---|
| ex14_2_4 | 7012 celdas | **8880** |
| alkyl | 264 celdas | **874** |
| ex6_2_12 | 29 884 | 29 908 |
| ex2_1_6 | 1764 | 1738 |
| ex14_1_5 | 28 | 28 |

Empeora en dos y es neutro en tres. Adelantar un contractor retrasa a los demás
y la propagación llega a un punto fijo peor: el orden de la cola importa más que
la elección de a quién favorecer. Queda implementado y **apagado**.

También hay que notar que la medición sobre una caja no puede responder esta
pregunta: en la caja raíz de un problema de optimización la variable objetivo no
tiene cotas (diámetro 1e8) y **ni DFB ni PolytopeHull la contraen**, porque
todavía no hay `loup`. La poda por objetivo solo existe dentro del árbol.

## 5. Descomposición del árbol: DFB poda bien, lo que falla es el costo

Con las cuatro configuraciones sobre el mismo banco
([results/decomposicion_ibexopt.csv](results/decomposicion_ibexopt.csv)), y un
**duplicado exacto de producción como control de ruido**, que da ×1.00 celdas y
×0.99 tiempo — o sea que todo lo de abajo está muy por encima del ruido.

Tomando `ACID+HC4` sin relajación lineal como referencia 1.0:

| configuración | árbol relativo | instancias resueltas |
|---|---|---|
| ACID+HC4 | 1.000 | **60** |
| **ACID+DFB** | **0.490** | 98 |
| ACID+HC4+PolytopeHull (producción) | 0.380 | **116** |

**DFB corta el árbol a la mitad y cubre el 74 % de lo que logra PolytopeHull**
(en escala logarítmica). El ×2.48 que se reportaba antes medía mayormente otra
cosa: `dfb` corre **sin ACID**, y quitarle ACID a producción ya multiplica el
árbol por 2.63.

Comparaciones de una sola variable:

| comparación | celdas | tiempo |
|---|---|---|
| quitar PolytopeHull a producción | ×2.63 | ×1.28 |
| agregar DFB a ACID+HC4 | **×0.49** | ×2.24 |
| `dfb` solo vs ACID+HC4 | **×0.83** | ×0.97 |
| `acid_dfb` vs producción | ×1.53 | ×3.53 |

Dos lecturas:

- **PolytopeHull es indispensable para la robustez**: quitárselo a producción
  baja de 116 a 60 instancias resueltas. DFB recupera de 60 a 98–102.
- **`dfb` solo ya poda un 17 % mejor que ACID+HC4 sin costar más** (×0.83
  celdas, ×0.97 tiempo).

**El déficit es costo por unidad de poda, no poder de poda**: DFB cuesta ×2.24
donde PolytopeHull cuesta ×1.28, o sea que PolyHull es **1.75× más barato por
unidad de poda**.

## 6. Hipótesis probada y rechazada: aplicar DFB cada K nodos

Si la contracción de un nodo sirviera para su subárbol, aplicar DFB cada K nodos
conservaría la poda a 1/K del costo. Es la idea que en `ibex-ipopt` funciona con
Ipopt. Implementado (`CtcDFBPropag::dfb_period`, `DFB_PERIOD=K`; en las llamadas
intermedias no se linealiza ni se inicializan los `2n` contractores), **falla
rotundamente**:

| instancia | período 1 | período 5 | período 20 |
|---|---|---|---|
| ex14_2_4 | 7012 celdas | **51 832** | 95 232 |
| ex6_2_12 | 29 884 | **145 896** | 149 318 |
| alkyl | 264 | **49 852** | 61 964 |
| ex2_1_6 | 1764 | **117 990** | 118 462 |

Todas las corridas con período > 1 agotan el límite de 30 s.

**La razón es una distinción que conviene tener presente**: el truco de «aplicar
cada K nodos» sirve para lo que produce **información global durable** —Ipopt
encuentra un punto factible, que mejora el `loup` y sirve para todo el resto de
la búsqueda—, y **no sirve para un contractor**, cuyo beneficio es **local a la
caja**. Saltearse un nodo deja ese nodo sin contraer, sus hijos salen más
grandes, y el árbol explota. Queda apagado.

## 7. Qué queda como lever real

La brecha está identificada —costo por unidad de poda— y **tres hipótesis para
cerrarla ya se rechazaron con datos**: warm start entre nodos, priorizar la
variable objetivo en la cola, y aplicar DFB cada K nodos.

Lo que queda por probar, en orden de lo que yo intentaría:

- **Contraer solo un subconjunto de las `2n` cotas por nodo.** Es distinto de
  saltearse nodos (que falla, §6) y de reordenar la cola (que falla, §4):
  reduce el costo proporcionalmente **conservando** contracción en cada nodo.
  `CtcPolytopeHull` tiene `set_contracted_vars` justamente para eso, así que la
  configuración de producción probablemente ya lo aproveche. Cómo elegir el
  subconjunto (las variables que más se contrajeron en el padre, las que
  aparecen en más restricciones activas) es la pregunta abierta.
- **En `acid_dfb` el cuello es el pricing** (28–48 % del contractor, y el
  contractor llega al 92 % del total en `alkyl`), no la inicialización. Las
  mejoras baratas de pricing ya se hicieron; lo que queda ahí es la lista de
  candidatos, medida como ×0.969 de tiempo en esa configuración.
- **`alkyl` es un caso patológico** que conviene mirar aparte: `acid_dfb` tarda
  9.48 s contra 1.36 s de `dfb`, con el 91.9 % del tiempo en el contractor.

---

## 8. El mecanismo real de PolytopeHull: una base tibia para las `2n` cotas

La hipótesis del §7 era contraer solo un **subconjunto** de las cotas. Al mirar
el código de `CtcPolytopeHull` resultó que **no restringe nada**:
`contracted_vars = BitSet::all(nb_var)`, calcula las `2n`. Lo que hace distinto
es **cómo**: mantiene una sola base y cambia solo el objetivo
(`set_cost(i, 1.0)` / `minimize()` / `set_cost(i, 0.0)`), más los indicadores
`inf_bound`/`sup_bound` que saltean los LP que la solución actual ya certifica
(«call to simplex useless, cf Baharev»).

Medido: **costo por cota**, dividiendo el trabajo total por las `2n` cotas.

| | por cota (mediana, 146 instancias) |
|---|---|
| PolytopeHull, iteraciones de simplex | **1.08** |
| DFB, pivotes | **1.95** |
| razón | **1.75** |

Y de la descomposición del árbol (§5): DFB cuesta ×2.24 de tiempo donde
PolytopeHull cuesta ×1.28, o sea **1.75×**. Las dos mediciones, independientes,
dan el mismo número: **el déficit de costo por unidad de poda ES el costo por
cota**, y su causa es que DFB arma `2n` tableaus independientes y busca cada
base desde cero, mientras PolytopeHull calcula las `2n` cotas desde una base
tibia a ~1 iteración cada una.

### Prueba del mecanismo en DFB

[chain_dfb.cpp](chain_dfb.cpp) (`make chain_dfb`) compara dos formas de obtener
las **mismas** `n` cotas inferiores:

- **A, independiente**: un `CtcDFB` por variable, cada uno inicializado desde
  `refA` (copia + eliminación).
- **B, en cadena**: un solo tableau; para cada `k` se renormaliza la columna `k`
  en la fila 0 (`renormalize_to`, **solo una eliminación, sin copia**) y se
  pivotea desde el estado que dejó `k-1`.

Resultado sobre 148 instancias, 3807 cotas
([results/cadena_una_base.csv](results/cadena_una_base.csv)):

| | independiente | en cadena |
|---|---|---|
| **pivotes** | 16 276 | **11 357** (−30.2 %) |
| **copias de matriz** | 3864 | **1 por instancia** |
| cotas idénticas | — | 3049 (**80 %**) |
| cotas mejores / peores | — | 248 / **510 (13 %)** |
| razón de pivotes por instancia | — | mediana **0.90**, p25 0.33 |

**El mecanismo funciona**: −30 % de pivotes y se eliminan casi todas las copias
de matriz, que son el 40.8 % del tiempo del contractor (§5.3 del plan). Ataca
exactamente el ×1.75.

**Con un costo real de calidad**: el 13 % de las cotas sale peor. La cadena
arrastra el estado de la cota anterior, que no siempre es un buen punto de
partida —`Prolog` pasa de 146 a 497 pivotes con 17 de 21 cotas iguales—,
mientras en otras es dramáticamente mejor: `Brown-10` de 90 a 25 pivotes con las
10 cotas idénticas.

### Qué falta para llevarlo a producción

La cadena requiere aplicar DFB como **barrido por lotes** (las `2n` cotas en
secuencia) en vez de `2n` contractores independientes en la cola de prioridad.
Eso es una reestructuración de `CtcDFBPropag`, y cambia la semántica: se pierde
el intercalado fino con HC4 — que, según §3 de
[MEDICIONES_INTERACCION.md](MEDICIONES_INTERACCION.md), rinde ×1.03, o sea casi
nada. Así que el intercambio parece favorable, pero hay que medirlo en el
optimizador.

Los ingredientes ya están: `CtcDFB::renormalize_to(k)` hace el cambio de cota
sobre el tableau vigente conservando el invariante de `λ`, y la certificación
garantiza que **cualquier** `λ` da una cota válida, así que un barrido con
heurística imperfecta puede dar cotas débiles pero nunca falsas.

---

## 9. El barrido por lotes: implementado, y peor en el optimizador

El mecanismo del §8 se implementó en `CtcDFBPropag` (`DFB_SWEEP=1`): un solo
`CtcDFB` recorre las variables, renormaliza la columna `k` en la fila 0 y saca
**los dos lados** de cada cota desde la misma base —pivotea hacia la inferior y
luego niega la fila 0 completa (`flip_side`, con lo que `γ_k` pasa de +1 a −1 y
el invariante de `λ` se conserva) para pivotear hacia la superior. Se alterna
con el punto fijo de HC4 en rondas gruesas.

Resultado sobre 191 instancias con límite de 30 s
([results/barrido_lotes_ibexopt.csv](results/barrido_lotes_ibexopt.csv)):

| | celdas | tiempo | resueltas |
|---|---|---|---|
| barrido (8 pivotes por cota) vs cola | **×1.37** | **×1.37** | 106 vs 116 |
| barrido (16 pivotes) vs cola | ×1.39 | ×1.43 | 105 vs 116 |

0 óptimos incompatibles con la referencia: es sólido, pero **37 % peor**.

### Por qué el mecanismo no se transfiere

El −30 % de pivotes del §8 era real, pero medido **en aislamiento**: calcular las
`n` cotas una vez. En el optimizador el barrido reemplaza el **intercalado fino**
de la cola —donde cada contractor se reaplica sobre una caja recién apretada por
HC4— por rondas gruesas, y eso pierde más de lo que ahorra.

Y ahí está la tensión estructural, que es la conclusión de toda esta línea:

- **El intercalado fino exige estado por contractor**, o sea `2n` tableaus.
- **La cadena exige un solo tableau**, y cambiar de cota cuesta una eliminación
  completa, O(m·(n+2m)).

`CtcPolytopeHull` no tiene que elegir: su solver mantiene **una factorización**,
y re-resolver con otro objetivo cuesta ~1 iteración de simplex. Para él «cambiar
de cota» es barato, así que puede hacerlo tantas veces como haga falta, en
cualquier orden. Para un tableau denso no lo es.

**Conclusión: la factorización de la base no es una mejora de memoria ni de
velocidad bruta — es lo que vuelve barato «cambiar de objetivo», y sin eso la
cadena no se puede combinar con el intercalado fino.** Es el habilitador del
mecanismo que explica el 1.08 contra 1.95 iteraciones por cota.

### Nota metodológica: el árbol es caóticamente sensible

Al variar los parámetros del barrido en una sola instancia:

| ex14_2_4 | ROUNDS=1 | ROUNDS=3 | ROUNDS=6 | | PIVOTS=2 | PIVOTS=8 | PIVOTS=16 |
|---|---|---|---|---|---|---|---|
| celdas | 9934 | **4212** | 10 130 | | 11 482 | **4212** | 7648 |

2.4× de variación entre parámetros vecinos, sin monotonía. No es ruido de
medición —las corridas son deterministas— sino sensibilidad del
branch-and-bound: un cambio chico en la contracción mueve el punto de bisección
y el árbol entero diverge. **Para cambios que alteran el camino de búsqueda, las
comparaciones por instancia no significan nada**; solo el agregado sobre muchas
instancias mide algo.

### Un límite de la certificación que conviene tener presente

La primera versión de `flip_side` declaraba **infactibles problemas factibles**.
La causa: `contract_float` niega `work[k]` cuando `upper_contract` es true,
porque en el camino de la cola la convención de signo vive en la **columna `k`
de la matriz**; con `flip_side` vive en la **fila 0**, y negar además la caja
evalúa un `γ` correcto contra una caja transformada.

**La certificación protege los multiplicadores, no la evaluación.** Garantiza
que cualquier `λ` dé una cota válida, y no protege contra evaluar ese `γ` contra
la caja equivocada. Quedó explícito en el predicado `CtcDFB::box_sign_flip()`.

---

## 10. Intento de atajo a la factorización: el modo liviano, y por qué falla

El §9 concluye que el habilitador es una factorización, porque es lo que vuelve
barato «cambiar de objetivo». Antes de escribir LU con actualización
Forrest–Tomlin —que son miles de líneas y hecho a medias mide ruido— se probó un
atajo: **no transformar nada**.

**El modo liviano** (`CtcDFB::light_mode`, `DFB_LIGHT=1`): las direcciones de
pivoteo son las filas **originales** de `refA`, compartidas y de solo lectura
entre los `2n` contractores, y el estado de cada cota es solo dos vectores,

```
gam : la combinacion vigente, gam = lam^T . refA     (n+m doubles)
lam : los multiplicadores                             (m doubles)
```

Con eso cambiar de objetivo es O(n+m) —se elige otra fila de arranque— y la
memoria por contractor pasa de O(m·(n+2m)) a O(n+m): con n=m=100, de 234 KB a
2.3 KB, **100× menos**, lo que atacaría las 6 instancias que no terminan.

**No funciona.** Es unas 3× más rápido y contrae casi nada:

| instancia | tableau | liviano |
|---|---|---|
| Katsura-12 | 45.83 % (67 pivotes) | **0.00 %** (0 pivotes) |
| Brown-30 | 93.44 % (1710) | **0.00 %** (1692 pivotes) |
| Discrete-Integralf2-21 | 95.63 % (1205) | **0.00 %** (20) |
| BroydenBanded-140 | 15.87 % (838) | 11.39 % (302) |

Y los 7 casos de prueba pierden los dos vacíos que deben detectar.

### Lo que NO lo explica

La primera hipótesis fue que sin eliminación la `γ` certificada dejaría de ser
dispersa —y la cota es `−(Σ_{i≠k} γ_i x_i)/γ_k`, así que cuantos menos términos,
más ajustada. **Medido, es falso**: la fracción de entradas de `γ` que son ~0 es
similar en los dos modos (Brown-30: 65.6 % con tableau, 62.9 % liviano).

### Lo que queda sin explicar

Dos síntomas distintos, ninguno cerrado:

- **Brown-30 pivotea 1692 veces y no consigue nada.** Con dispersión de `γ`
  comparable, las direcciones sin transformar simplemente no hacen progreso: se
  pivotea en falso.
- **Katsura-12 no llega a certificar ninguna cota** (0 contractores pasan la
  guarda de `γ_k`), lo que sugiere además un problema de inicialización en
  `set_bound_light` que no se persiguió.

**Conclusión: el atajo no reemplaza a la factorización.** La transformación de
las filas no es un detalle de implementación: es lo que hace que un pivote
avance. Una factorización sirve porque calcula esas mismas filas transformadas
**a demanda** (BTRAN), en vez de renunciar a ellas.

El modo queda implementado y **apagado**, con el camino por defecto verificado
intacto (7 casos con sus 2 vacíos, testigos conservados, Katsura-12 45.83 % y
Brown-30 93.44 %). Sirve como punto de partida si alguien quiere retomar la
línea, sabiendo que el problema a resolver es el progreso por pivote, no la
memoria ni la dispersión.

---

## 11. El costo de certificar: la mitad del tiempo de PolytopeHull

**Corrección previa**: el modo certificado es de **Ibex**, no de SoPlex.
`LPSolver::Mode::Certified` dispara `neumaier_shcherbina_postprocessing()`
después de cada LP:

```cpp
Matrix A_trans = rows_transposed();                 // copia transpuesta, por LP
IntervalVector rest = A_trans * uncertified_dual_;  // gamma = lambda^T.A en intervalos
rest -= cost();
obj_ = uncertified_dual_*b - rest*ivec_bounds_;     // evaluado sobre la caja
```

Eso es **exactamente lo que hace la certificación de DFB**. Los dos resuelven en
punto flotante y certifican con Neumaier–Shcherbina; la diferencia está solo en
cómo se encuentran los multiplicadores (1.08 iteraciones de simplex desde una
base tibia contra 1.95 pivotes de DFB).

**Consecuencia**: «pivotear en flotantes y certificar la cota» **no es un
diferenciador de DFB frente a PolytopeHull**. Lo es frente al DFB anterior, que
pivoteaba en intervalos y por eso necesitaba `regenerateA`, pero no frente a
PolyHull.

### Cuánto cuesta

Primer intento inválido: poniendo el solver en `NotCertified`, PolytopeHull
**deja de contraer** (0.00 %), porque `optimizer()` solo usa la cota cuando el
estado es `OptimalProved`, que solo produce el camino certificado. Eso medía
«PolyHull sin certificar no hace nada».

Medición correcta ([results/sonda_ns_lpwrapper.cpp](results/sonda_ns_lpwrapper.cpp)):
una copia parcheada del envoltorio que **saltea el postprocesado pero reporta
`OptimalProved` igual** (`PH_SKIP_NS=1`), de modo que la contracción y las
iteraciones no cambian y la diferencia de tiempo es el costo puro de la
certificación. Es una sonda deliberadamente no sólida.

| instancia | con N-S | salteada | **certificación** |
|---|---|---|---|
| DiscreteBoundary-0100 | 0.0666 s | 0.0223 s | **67 %** |
| BroydenBanded-140 | 0.0821 s | 0.0344 s | **58 %** |
| BroydenBanded-120 | 0.0519 s | 0.0259 s | **50 %** |
| Katsura-12 | 0.0018 s | 0.0015 s | 17 % |
| Katsura-50 | 0.0396 s | 0.0343 s | 13 % |

Misma contracción y mismas iteraciones en todos: la diferencia es solo el
postprocesado.

**La certificación es el 50–67 % del tiempo de PolytopeHull en las instancias
grandes.** Y la misma tarea en DFB cuesta **4.8 %** de su contractor.

### Por qué la diferencia, y qué implica

DFB certifica recorriendo **solo los `λ` no nulos** —y `nnz(λ)/m` tiene mediana
**0.13** (§8)—, mientras la implementación de Ibex hace un producto denso
`A_trans * dual` sobre las `m` filas **y reconstruye una copia transpuesta de la
matriz en cada LP**.

Eso corta en dos direcciones y conviene no confundirlas:

1. **Es un diferenciador real de DFB hoy**: certifica ~10× más barato que
   PolyHull, y esa era la ventaja que faltaba encontrar.
2. **Pero es en buena medida una ineficiencia de implementación de Ibex, no una
   ventaja estructural.** El dual de un simplex también es disperso, así que la
   N-S de PolyHull podría explotarlo y evitar la transpuesta por llamada. Si
   alguien la optimiza, PolytopeHull se vuelve **1.5–2× más rápido** y la
   posición competitiva de DFB empeora en la misma medida.

**Recomendación**: antes de invertir en la factorización por etapas, optimizar la
N-S de Ibex (evitar `rows_transposed()` por llamada y explotar la dispersión del
dual) y **volver a medir la comparación**. Es un cambio chico en `ibex_LPSolver`,
mejora el solver de producción para todos sus usuarios, y define contra qué
vara hay que competir. Medir contra una vara que tiene un 50 % de grasa lleva a
conclusiones equivocadas sobre DFB.

---

## 12. La N-S de Ibex optimizada: PolytopeHull ×1.66, y la vara se movió

Hecho lo que recomendaba el §11. El cambio está en
`src/numeric/ibex_LPSolver.cpp::neumaier_shcherbina_postprocessing()` y es
exactamente el que se había diagnosticado: la misma fórmula

    obj_ = λᵀb − (Aᵀλ − c)ᵀ·[x]

acumulada **por filas, salteando los `λ` nulos**, en vez de construir `rows()`
(`m×n`) y su transpuesta (`n×m`) en cada resolución del LP. Una fila con
`λ_i = 0` no aporta ni a `λᵀb` ni a `Aᵀλ`, así que omitirla es exacto, no una
aproximación. El trabajo pasa de `m·n` a `nnz(λ)·n` y, sobre todo, desaparecen
las dos matrices densas que se asignaban y llenaban **una vez por cota, o sea
`2n` veces por caja**.

Se conservan dos interruptores: `IBEX_NS_LEGACY=1` vuelve al camino original
(para poder medir A/B sin relinkear) y `IBEX_NS_RIGOROUS=1` acumula el residuo
en intervalos (ver más abajo).

### Medición

Banco de una caja, 148 instancias comparables, tiempo total del contractor:

| | `t_ph` total | |
|---|---|---|
| N-S original | 1.5906 s | |
| **N-S optimizada** | **0.9580 s** | **×1.66 (−39.8 %)** |
| `t_dfb` (control) | 6.2713 → 6.2464 s | sin cambio, como debe ser |

**Contracción e iteraciones idénticas en las 148 instancias** (0 diferencias en
`per_ph` ni en `soplex_iters`): el cambio es puro ahorro, no un intercambio.

Por instancia, contra el piso que marcaba la sonda `PH_SKIP_NS` del §11:

| instancia | original | **optimizada** | piso (sin certificar) |
|---|---|---|---|
| DiscreteBoundary-0200 | 0.2833 s | **0.0601 s** | — |
| BroydenBanded-200 | 0.1853 s | **0.0619 s** | — |
| BroydenBanded-140 | 0.0821 s | **0.0353 s** | 0.0344 s |
| DiscreteBoundary-0100 | 0.0686 s | **0.0229 s** | 0.0223 s |
| BroydenBanded-120 | 0.0514 s | **0.0262 s** | 0.0259 s |
| Katsura-50 | 0.0396 s | **0.0351 s** | 0.0343 s |

La certificación **ya no se nota**: donde antes era el 50–67 % del tiempo de
PolytopeHull, ahora el tiempo con y sin certificar coincide dentro del ruido
(0.0353 contra 0.0344, 0.0229 contra 0.0223). Es decir, el 50–67 % era grasa de
implementación, no el costo de certificar.

### Consecuencia para DFB

La ventaja que el §11 había encontrado —«DFB certifica 10× más barato»—
**desapareció**, tal como se anticipaba ahí. Y la vara empeoró: la relación de
costo del contractor en el banco de una caja pasa de

    t_dfb / t_ph = 6.27 / 1.59 = ×3.9   →   6.25 / 0.96 = ×6.5

Conviene ser explícito sobre lo que esto significa y lo que no. **No** significa
que DFB sea 6.5× peor en el árbol: el banco de una caja corre DFB sin HC4 y
hasta su punto fijo, y en `ibexopt` la diferencia medida era ×2.24 contra ×1.28
(§5). Sí significa que **el margen que DFB tiene que recuperar es mayor que el
que se creía**, y que la comparación del §5 hay que releerla con la vara nueva.

### Un agujero de solidez que quedó a la vista

Al reescribir el postprocesado quedó claro que el residuo `Aᵀλ − c` **se acumula
en punto flotante, no en intervalos**: en el código original
`A_trans*uncertified_dual_` resuelve a `Vector operator*(const Matrix&, const
Vector&)`, o sea un producto de dobles, y solo después el resultado se convierte
a `IntervalVector` degenerado. El error de redondeo de ese residuo entra
multiplicado por el ancho de la caja, así que **la cota que Ibex reporta como
`OptimalProved` no es rigurosa en sentido estricto**.

El camino nuevo reproduce ese comportamiento a propósito, para que la medición
A/B sea limpia (solo cambia el orden de asociación). `IBEX_NS_RIGOROUS=1` activa
la variante que acumula el residuo en aritmética de intervalos. **Es un hallazgo
de paso, no algo que se haya decidido cambiar**: cerrarlo ensancharía un poco
todas las cotas de PolytopeHull y hay que medirlo antes de tocar el
comportamiento por omisión del solver de producción.

### La vara en el optimizador: casi no se movió

Re-medida la comparación completa con la N-S optimizada (191 instancias, 30 s,
mismo banco):

| config | resuelve | celdas (geom) | cpu (geom) |
|---|---|---|---|
| `acidhc4+ph` (producción) | **137** | 1.000 | 1.000 |
| `dfb` | 116 | ×2.755 | ×1.363 |
| `acid_dfb` | 105 | ×1.447 | ×3.405 |

Sobre las 99 instancias que las tres resuelven; 0 óptimos incompatibles.

Contra los valores previos a la optimización (×2.83 / ×1.40 para `dfb`;
×1.47 / ×3.75 para `acid_dfb`) **las razones no empeoran; si algo, mejoran un
pelo**, que es lo contrario de lo esperado si acelerar la referencia importara.
Y las tres cuentas de instancias resueltas quedan iguales (137 / 116 / 105).

O sea: **el ×1.66 del contractor no se traduce en el árbol**. Encaja con el §1
—el contractor es el 9–42 % del tiempo del optimizador— y con que en el árbol
las cajas son más chicas y PolyHull hace menos iteraciones, mientras el costo
de la N-S depende de `m·n` y no de las iteraciones. La diferencia queda dentro
del ruido de corrida.

Está corriendo un A/B directo (`acidhc4+ph` con y sin `IBEX_NS_LEGACY`, mismas
191 instancias) para poner un número a eso en vez de inferirlo de una
comparación entre corridas.
