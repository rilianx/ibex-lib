# DFB: estado, arquitectura y pasos a seguir

**Documento de traspaso.** Está escrito para que una sesión nueva pueda
retomar el trabajo sin releer todo el historial. Los detalles de cada medición
están en los cuatro documentos de resultados (§9); acá va lo necesario para
entender el código, no romper lo que funciona, y seguir.

---

## 1. El objetivo

**Reemplazar `CtcPolytopeHull` por DFB** como contractor lineal del filtrado de
Ibex. No es una línea exploratoria: DFB tiene que llegar a igualar la
contracción que hoy consigue PolytopeHull con SoPlex, no solamente a ser más
barato.

La forma prevista es un **DFB incremental** que interactúe con la propagación de
restricciones (HC4) o que viva dentro de ACID/CID, en lugar de recomputar desde
cero en cada llamada. El medio es tratar a DFB como lo que es —una variante de
simplex— e incorporarle las optimizaciones del simplex moderno.

### Criterios de aceptación y su estado

| criterio | estado |
|---|---|
| 1. **Solidez absoluta.** Un contractor que vacía cajas con soluciones no puede reemplazar a nada | **cumplido**: 2 vaciados en 147 instancias, ambos corroborados por `ibexsolve`; 32/32 testigos certificados conservados; 0 óptimos incompatibles en `ibexopt`; `valgrind` limpio |
| 2. **Igualar contracción** (grupo A del §5.2) | **grupo A 73 de 98 productivas** (empezó en 2 de 98) |
| 3quinquies. **El híbrido domina al simplex puro** | `DFB_SX_HIBRIDO=1`: regla de impacto primero, simplex solo donde aquella contrajo menos del 1 %. Misma poda que el simplex (×1.065 contra ×1.062), cpu de ×1.195 a **×1.129**, y 122 resueltas contra 120. Pero sigue sin superar a producción ([MEDICIONES_PODA.md](MEDICIONES_PODA.md) §12) |
| 3quater. **El warm start entre nodos no puede pagar, y es estructural** | `CtcDFBPropag::linearize` se llama en cada nodo y `LinearizerXTaylor` expande en una esquina de la caja actual: **los coeficientes cambian en cada nodo**, así que la base anterior es una conjetura sobre otro poliedro. Con base **por contractor** (`DFB_SX_PROPIO=1`) el warm start sí paga **dentro** de una caja —misma contracción, −39 % de pivotes— pero en el árbol no mueve nada ([MEDICIONES_PODA.md](MEDICIONES_PODA.md) §11) |
| 3ter. **El warm start: gana en la caja y pierde en el árbol** | `DFB_SX_WARM_MIN=1` (tibio solo para cotas inferiores, que es la mitad sana) da en el banco de una caja la **mejor clasificación** (A 96, B 2), la mejor contracción y **un tercio del tiempo**; en el árbol cuesta **28 % más** que el frío (cpu ×1.563 contra ×1.221) y resuelve 118 contra 121. El banco de una caja y el árbol premian configuraciones distintas ([MEDICIONES_PODA.md](MEDICIONES_PODA.md) §10) |
| 3bis. **Con la regla de pivote de simplex** (`DFB_SIMPLEX_RULE=1`) | **la poda se cierra**: celdas ×1.068 contra producción (eran ×1.321), o sea prácticamente su árbol. **El costo pasa a ser el cuello**: cpu ×1.227, porque cuesta ×21 pivotes. Resuelve 120 contra 137. 0 óptimos incompatibles ([MEDICIONES_PODA.md](MEDICIONES_PODA.md) §9) |
| 3. **Ser competitivo en el optimizador** | **pendiente, pero mucho más cerca de lo que se creía.** Con la sustitución directa (`--filtering=acidhc4 --lr=dfb`, DFB en el hueco exacto de PolyHull) el costo por nodo está **a la par**: ×1.069 de cpu geométrica y ×1.384 de celdas sobre 121 instancias comunes, y sin ACID DFB es **21 % más rápido** que PolyHull (×0.788) siendo más rápido en 77/121. Resuelve **123** contra 137. Lo que falta es **poda**, no velocidad ([MEDICIONES_FACTORIZACION.md](MEDICIONES_FACTORIZACION.md) §6) |
| 4. **Incrementalidad como arquitectura**, no como optimización | **medida y con veredicto** ([MEDICIONES_INCREMENTAL.md](MEDICIONES_INCREMENTAL.md)): reoptimizar cuesta **cero pivotes** en la instancia mediana con aprietos realistas. Pero el warm start **entre nodos** se implementó y se descartó: exige reusar la linealización, que es más floja y cuesta un 17 % más de celdas |

---

## 2. Cómo compilar y correr

Todo vive en `~/ibex2024/examples/dfb`. Requiere la biblioteca Ibex compilada en
`../../src` más Gaol y SoPlex; las rutas están fijadas en el `Makefile`.

```bash
make bench_dfb tests_dfb_all diag_dfb ibexopt
```

| target | qué es |
|---|---|
| `bench_dfb` | una línea CSV por instancia: contracción de DFB y de PolytopeHull sobre la misma caja, pivotes, desglose de tiempo por fase, contadores de diagnóstico |
| `diag_dfb` | diagnóstico por contractor: cuál vacía la caja y por qué, validez de `gamma`, validez de la linealización, y **verificación con testigos certificados** |
| `incr_dfb` | métrica incremental: pivotes al reoptimizar tras un cambio de cotas contra arrancar de cero |
| `chain_dfb` | las `n` cotas desde una base en cadena contra `n` contractores independientes |
| `tests_dfb_all` | los 7 casos de prueba del proyecto (copia de `TestingDualFeasibleBounding.cpp` con todos habilitados en `main`) |
| `ibexopt` | el optimizador, con `--filtering=hc4\|acidhc4\|dfb\|dfb_hc4\|acid_dfb` y `--lr=xn` (agrega PolytopeHull) o `--lr=no` |
| `example_dfb` | el ejemplo original, ver [README.md](README.md) |

Barridos (en `results/`, corren en paralelo con un proceso por instancia):

```bash
results/run_batch_par.sh salida.csv 60 12 4000000     # banco de una caja, 153 instancias
python3 results/run_opt.py salida.csv                  # ibexopt, 191 instancias
python3 results/analyze.py salida.csv                  # clasificación A/B/C y agregados
```

Con el banco ampliado la brecha de `dfb` **se ensancha**: ×2.49 → ×2.83 celdas y
×1.23 → ×1.43 tiempo respecto de producción. El banco liviano la subestimaba.

**El banco de `ibexopt` incluye los conjuntos pesados a propósito**: `easy`
(101), `medium` (31), `hard` (7), `blowup` (25) y `benchs-minlp` (27), en total
191. Se seleccionan con `OPT_SETS=easy,medium,...` y el límite de tiempo con
`OPT_TL=<segundos>`. Las mediciones anteriores a esto usaban solo `easy` y
`medium` con 20 s, que es un banco liviano: enmascara las mejoras de costo por
nodo (casi todo termina rápido de todos modos) y manda a *timeout* justamente
las instancias que distinguen las configuraciones.

**Las instancias que no terminan no se descartan del análisis**: cuántas
resuelve cada configuración es en sí mismo un resultado — en este trabajo fue la
métrica que más movió el diagnóstico (`dfb` pasó de 62 a 103 resueltas).

**Correr los barridos en paralelo usando todos los cores.** Está verificado que
los tiempos no varían mucho respecto de una corrida secuencial, así que no hace
falta una pasada aparte para medir tiempos. Al elegir el tope de memoria hay que
multiplicarlo por el número de jobs: las instancias grandes son las que explotan.

### Interruptores por variable de entorno

Los defaults son las decisiones tomadas; las variables sirven para volver atrás
y comparar. Se leen en los inicializadores estáticos de `CtcDFB`, así que valen
en los tres binarios.

| variable | efecto |
|---|---|
| `DFB_INTERVAL_PIVOT=1` | vuelve al pivoteo en intervalos (default: flotantes) |
| `DFB_PRICING_NORM=0` | vuelve al pricing tipo Dantzig (default: normalizado) |
| `DFB_NO_HARRIS=1` | desactiva el desempate de Harris en el ratio test |
| `DFB_NO_FUSED=1` | vuelve al pricing con `IntervalVector` intermedios |
| `DFB_NO_PAIR_INIT=1` | inicializa los `2n` contractores por separado |
| `DFB_PARTIAL=1`, `DFB_PARTIAL_SIZE`, `DFB_PARTIAL_REFRESH`, `DFB_PARTIAL_ADAPT` | lista de candidatos en el pricing (apagada) |
| `DFB_NO_BOUND_EMPTY=1` | desactiva la prueba de vacío por la cota |
| `DFB_NO_LIN_EMPTY=1` | no usa el `-1` del linealizador para vaciar la caja |
| `DFB_MAXITERS=k` | pivotes por visita de un contractor (default 1; `-1` = hasta su punto fijo) |
| `DFB_PIVOT_BUDGET=k`, `DFB_PIVOT_BUDGET_N=c` | tope **total** de pivotes DFB por caja, absoluto o `c·n` (default sin tope) |
| `BENCH_WITH_HC4=1` | en `bench_dfb`, crea también los HC4 para que la cola intercale |
| `IBEX_NS_LEGACY=1` | vuelve a la Neumaier–Shcherbina original de Ibex (con las dos matrices densas) |
| `IBEX_NS_RIGOROUS=1` | acumula el residuo `Aᵀλ − c` en intervalos en vez de dobles (ver §12 de MEDICIONES_TECHOS) |
| `DFB_DENSE=1` | vuelve al tableau denso (default: **disperso**, resultados idénticos y ×2.2) |
| `--lr=dfb` (no es variable de entorno) | pone DFB en el hueco exacto de `CtcPolytopeHull`; **`--filtering=acidhc4 --lr=dfb` es la configuración objetivo** |
| `DFB_SPARSE_CHECK=1` | mantiene las dos representaciones a la vez y compara tableau y decisiones en cada pivote (carísimo, solo diagnóstico) |
| `DFB_SIMPLEX_RULE=1` | reemplaza la **regla de pivote** por la de un simplex dual, que alcanza el óptimo del LP: grupo B de 33 a 3 (§7 de MEDICIONES_PODA). Alias: `DFB_SIMPLEX` |
| `DFB_SX_WARM=1` | reusa la base entre cotas: **1.5 pivotes por cota**, pero `γ_k` deriva y la cota certificada queda débil (§7.6 de MEDICIONES_PODA) |
| `SX_ESCALA=1` | en `test_dfbsimplex`, régimen de escalas calibrado al banco (filas hasta 1e6, cajas hasta 1e6) |
| `DFB_SX_SCALE=1` | activa el escalado de filas y columnas (apagado por omisión: **daña** en el banco real) |
| `DFB_SX_REQUEUE=t` | umbral de reencolado del simplex (default 0.01); medido sin efecto |
| `DFB_SX_CHECK=1` | mide el residuo de `d` sobre las columnas básicas (consistencia de `y` con la base) |
| `DFB_SX_MEMO_N=k` | saltea el simplex `k` nodos tras una invocación sin resultado. **Medido y peor** (§12.2 de MEDICIONES_PODA) |
| `DFB_SX_HIBRIDO=1`, `DFB_SX_HIB_TOL=t` | impacto primero y simplex solo donde contrajo menos de `t` (default 0.01). **La mejor combinación medida de poda y costo** |
| `DFB_SX_PROPIO=1` | un simplex **por contractor** en vez de uno compartido: −39 % de pivotes dentro de una caja, neutro en el árbol |
| `DFB_SX_NODOS=1` | conserva la base entre nodos cuando la linealización mantiene dimensiones (neutro: la matriz cambia igual) |
| `DFB_SX_WARM_MIN=1` | tibio solo para las cotas inferiores: la mitad sana del warm start (la superior tiene un defecto, §10 de MEDICIONES_PODA) |
| `DFB_SX_EXACT=1` | refresco exacto desde la factorización en cada iteración: arregla el warm start pero el grupo B vuelve a 30 (§7.9 de MEDICIONES_PODA) |
| `DFB_SX_TRACE=1` | imprime por cota los pivotes, la cota flotante, `γ_k` y si certificó |
| `DFB_SX_STATS=1` | imprime estados del simplex (óptimo / infactible / límite) por instancia |
| `DFB_DUAL=1` | pricing por violación de cota (simplex dual). **Implementado y descartado**: no termina, ver §5 de MEDICIONES_PODA |
| `DFB_STOP_TOL=t` | umbral de mejora estimada por debajo del cual DFB se detiene (default 1e-6) |
| `BENCH_CORNER=inf\|sup` | en `bench_dfb`, misma política de esquina determinista para DFB y PolyHull |
| `DFB_PRICE_STATS=1` | acumula la dispersión del pricing que mide `probe_price` (no medir tiempos con esto puesto) |
| `DFB_TRACE=1` | imprime cada pivote como `T k=.. up=.. j=.. i=.. alpha=..`, para diffear secuencias entre variantes |
| `DFB_REFACTOR_ALWAYS=1`, `DFB_REFACTOR_PERIOD=k` | en `DFBBasis`, refactorizar en cada cambio de base o cada `k` updates |

---

## 3. La arquitectura actual: pivoteo flotante con certificación

**Esto es lo primero que hay que entender antes de tocar el código.**

La linealización da `refA · x = 0` para toda solución de la caja (con `x`
incluyendo las variables `b`). Entonces para **cualquier** vector real `λ`, la
combinación `γ = λᵀ·refA` cumple `γ·x = 0`, y de ahí sale la cota

```
x_k  ≥  cota( −(Σ_{i≠k} γ_i·x_i) / γ_k )   evaluada sobre la caja
```

La validez **no depende de que `λ` sea bueno**, solo de que `γ` se calcule
rigurosamente a partir de `refA`. Un `λ` malo da una cota débil, **nunca una
cota falsa**. Es dualidad débil. Por eso el pivoteo —que es puramente la
búsqueda de `λ`— se hace en `double` sin ninguna garantía, y solo la evaluación
final es en intervalos. Es el mismo esquema de SoPlex en modo certificado.

### El invariante que hay que respetar

El tableau flotante se mantiene **aumentado con la identidad**:

```
Af = [ A | I ]        m × (nb_var+m + m)
```

y la parte aumentada de la fila `j` **es el `λ` de esa fila**. Invariante:
`λ_jᵀ·refA_k = parte A de la fila j`, con `refA_k` = `refA` con la columna `k`
negada cuando se contrae la cota superior (por eso la certificación usa el
`refA` compartido y solo niega la entrada `k`, sin almacenar una referencia por
contractor).

**El punto crítico**: el forzado de entradas a 0 exactos se aplica **solo a la
parte A, nunca a `λ`**. Así `λ` queda siempre veraz y la certificación absorbe
toda la suciedad numérica del pivoteo. Si alguien toca `λ` para "limpiarlo", se
rompe la solidez.

`certify_gamma()` calcula `γ = λᵀ·refA` en intervalos recorriendo solo los
`λ_j ≠ 0`, y descarta la cota si `γ_k` contiene el 0 (dividir por un intervalo
que contiene el 0 no acota nada; no es un problema de solidez, es descartar un
`λ` malo). Costo medido: `nnz(λ)/m` tiene mediana **0.13**, o sea que certificar
una cota cuesta ~13 % de **un solo** pivote de la vieja versión de intervalos.

### Qué desapareció con esta decisión

`get_Aerror()` devuelve 0, `regenerateA()` no se dispara (731 regeneraciones →
**0**), y con ellas la política de refactorización programada quedó sin objeto. También se resignó el vaciado por «ninguna fila bloquea», que con
pivoteo flotante no es prueba: costó nada, porque `n_no_candidate` fue 0 en todo
el banco y los 2 vaciados vienen de la cota.

### El camino de intervalos sigue existiendo

`DFB_INTERVAL_PIVOT=1` lo reactiva, y con él `regenerateA`, `get_Aerror`, el
pricing fusionado y la deducción del par de contractores lb/ub. Sirve para
comparar; no es el camino de producción.

---

## 4. Incrementalidad: la restricción transversal

El objetivo es que DFB se reoptimice cuando la propagación o ACID **modifican
las cotas de la caja**, sin rearmar nada.

**A favor**: un cambio de cotas **preserva la factibilidad dual**. Es la razón
por la que los solvers de LP usan simplex dual dentro de branch-and-bound, y DFB
la hereda. Mientras no se relinealice, `refA` no cambia y solo se mueve la caja.
Y ahora que no hay regeneraciones, **conservar la base entre llamadas es posible
de verdad**: antes cada regeneración era un reinicio forzado.

**Dos regímenes, y solo uno es el objetivo**:
- *Dentro de una linealización* (la propagación aprieta cotas): la base sigue
  válida, warm start, pocos pivotes. Es el caso frecuente y es el objetivo.
- *Al relinealizar* (nodo nuevo, esquina nueva): `refA` cambia entera y hay que
  reiniciar. Inevitable, y PolytopeHull paga lo mismo.

**Métrica propia, ya medida** con [incr_dfb.cpp](incr_dfb.cpp) (`make incr_dfb`,
ver [MEDICIONES_INCREMENTAL.md](MEDICIONES_INCREMENTAL.md)): pivotes tras un
cambio de cotas frente a arranque en frío. Resultado: mediana **0.032** por
instancia y 0 pivotes en varias familias enteras —la propiedad dual se cumple en
la práctica—, pero razón agregada **1.000**, porque tres instancias
(`Geneigbis` ×1.98, `Fourbar-icse` ×1.47) concentran los pivotes del banco y en
ellas el warm start cuesta **más**: la base heredada puede quedar lejos en la
geometría nueva. Y el 5 % de los contractores llega a una cota **peor** tibio
que frío.

**Consecuencia**: el warm start necesita una **salvaguarda** —si la
reoptimización supera el presupuesto de la corrida anterior, reiniciar en frío—,
que además es práctica estándar en los solvers de LP. Sin ella es neutro en el
banco completo.

**Nota**: `max_iters = 1` hace que cada llamada avance un solo pivote. Está sin
decidir si es diseño o valor heredado (§8).

---

## 5. Hechos medidos que condicionan las decisiones

### 5.1 La contracción extra NO se traslada al árbol

El hallazgo más importante para priorizar. Normalizar el pricing subió la
contracción media de 24.03 % a 30.03 % sobre una caja (+6 puntos, grupo A de 64 a 72) y en el
optimizador dio **×0.98 de celdas y ×0.97 de tiempo**. En cambio el cambio a
flotantes (−64 % de tiempo del contractor) movió `acid_dfb` de ×6.47 a ×3.44.

**Conclusión: contraer más una caja grande no es el factor que limita dentro del
optimizador.** La propagación y ACID recuperan buena parte de lo que el
contractor lineal deja sobre la mesa. Lo que rinde es **tiempo por nodo y
robustez**, no contracción pura. Cualquier propuesta futura debería justificarse
en esos términos.

### 5.2 Clasificación A/B/C, y por qué el grupo C no cuenta

Sobre las 145 instancias donde DFB no vacía:

| grupo | criterio | instancias |
|---|---|---|
| A | DFB contrae ≥ que PolytopeHull | **72** |
| B | DFB contrae < que PolytopeHull | 24 |
| C | ninguno de los dos contrae | 49 |

**El grupo C se excluye de los promedios**: ahí tampoco contrae PolytopeHull
(gasta hasta 257 iteraciones de simplex en el intento), así que no es una
debilidad de DFB — la relajación lineal simplemente no corta en esas cajas.

### 5.3 Perfil de tiempo en el camino flotante

| fase | % |
|---|---|
| **inicialización de los `2n` contractores** | **40.8 %** |
| pricing | ~37 % |
| certificación (`λᵀ·refA`) | 4.8 % |
| | *(el 4.8 % sigue siendo cierto para DFB; lo que cambió es que la N-S de Ibex ya no cuesta el 50–67 %, así que dejó de ser una ventaja comparativa — §8.0 A)* |
| eliminación | 3.4 % |
| cota, ratio test, resto | ~5 % |
| regeneración de `A` | 0 % |

(En el camino de intervalos era: pricing 54 %, eliminación 18 %, init 18 %,
regeneración 9 %.) Dos consecuencias: **la certificación es barata**, con lo que
la duda sobre su costo queda zanjada; y **la inicialización es lo que más pesa**,
no porque empeorara sino porque todo lo demás se abarató. Es lo que pone la
factorización de la base como prioridad.

### 5.4 Los resultados son deterministas

Se creía que la esquina `RANDOM` del linealizador metía ruido entre corridas.
**No lo hace**: 5 repeticiones completas del banco dan el mismo resultado en 145
de 145 instancias, porque cada instancia corre en su propio proceso y el
generador de Ibex se inicializa igual. **Una corrida basta** para comparar
variantes. Variar la semilla a propósito mediría otra cosa (sensibilidad a la
esquina), y sería un experimento legítimo pero distinto.

### 5.5 Densidad de la linealización

Mediana 0.190 en el banco, pero **0.025 para n ≥ 50** — las matrices grandes son
dispersas, y son las que agotan tiempo y memoria. Solo justifica explotar la
dispersión en instancias grandes.

### 5.6 Lo que aún no termina

6 instancias de 153 no terminan (`Eiger-1000` por `bad_alloc`,
`BroydenBanded-200/1000`, `BroydenTri-1000`, `DiscreteBoundary-1000`,
`ExtendedFreud-1000` por tiempo o memoria). La causa es estructural: `2n`
contractores × una matriz `m × (n+m+m)` cada uno. Solo lo resuelve dejar de
almacenar el tableau, es decir la factorización de la base (paso 4 del §8).

---

## 6. Protocolo de validación: correr esto antes de aceptar cualquier cambio

1. **Los 7 casos** de `./tests_dfb_all` deben pasar, con `NO SOLUTION TEST 1` y
   `2` detectados **vacíos**. Si el binario se cuelga en el tercer caso, hay una
   regresión: el original se colgaba ahí porque al vaciar no ponía
   `state = FINAL`.
2. **Campaña de testigos**: la caja contraída debe seguir **intersectando** una
   caja-solución certificada por `ibexsolve -s`.
   ```bash
   ./diag_dfb data_tests/Brown-10.bch "lb1:ub1,lb2:ub2,..."
   ```
   **Hay que usar la caja, no su punto medio.** `ibexsolve` certifica una caja
   que *contiene* una solución; con el punto medio aparecen violaciones espurias
   porque su error de ~1e-16 se amplifica por coeficientes de ~1e6. El script
   `results/witness_campaign.py` automatiza la campaña.
3. **Todo vaciado de caja debe ser corroborable.** Hoy son `Prolog` y
   `Prolog-icse`, ambos `infeasible problem` según `ibexsolve`.
4. **En `ibexopt`, ningún óptimo incompatible** con `--filtering=acidhc4 --lr=xn`.
5. **`valgrind` sin errores** sobre `ibexopt --filtering=dfb`.

---

## 7. Trampas conocidas (todas costaron tiempo al menos una vez)

1. **No recompilar mientras un barrido corre.** Relinquear `ibexopt` o
   `bench_dfb` mientras un barrido los ejecuta corta el barrido (permission
   denied) o mezcla dos binarios en los mismos datos. Pasó dos veces.
2. **Hay dos bucles de inicialización de los contractores DFB**, uno en
   `init_dfb_contractors` y otro dentro de `contract` (el que se usa en modo
   *stand_alone*). Están factorizados en `CtcDFBPropag::init_all_dfb`; parchear
   solo uno da mediciones sin efecto.
3. **El despacho al camino flotante va en `init()`**, no solo en la propagación:
   hay llamadores que usan `init()` directamente (los tests).
4. **El header del CSV de `run_batch_par.sh` tiene que coincidir** con las
   columnas que emite `bench_dfb`. Si se agregan columnas al `.cpp` hay que
   actualizar el `HDR` del script, o el análisis lee campos corridos.
5. **`Ctc::contract` no debe cambiar el tamaño de la caja del llamador.**
   `ibex::Optimizer` conserva referencias a sus componentes entre llamadas y
   `resize()` libera el arreglo interno. Se trabaja sobre el buffer `work`.
6. **ASan no reprodujo el segfault de memoria liberada; valgrind sí.** Para
   corrupción de memoria en este código, ir directo a valgrind.
7. **Los asserts de Ibex están activos** en estos TU (no hay `-DNDEBUG`), así
   que un `IntervalVector::operator[]` fuera de rango aborta con mensaje. Si en
   cambio hay segfault silencioso, el problema está en otro lado.
8. **`pkill -f patrón` mata la propia shell** si el patrón aparece en la línea
   de comandos (por ejemplo si en el mismo comando se menciona el archivo).
9. **`bench_dfb` daba a DFB y a PolytopeHull linealizadores distintos**, los dos
   con política de esquina `RANDOM`, así que **no resolvían el mismo poliedro** y
   la comparación de contracción mezclaba dos efectos. Con `BENCH_CORNER=inf`
   los dos usan la misma política determinista, y ahí DFB queda **peor**: grupo
   A 74 → 65, B 25 → 33. Todos los números de grupos A/B/C anteriores a esa
   medición hay que leerlos con ese ruido encima
   ([MEDICIONES_PODA.md](MEDICIONES_PODA.md) §4). La comparación en `ibexopt` no
   está afectada, porque cada corrida usa un solo contractor.
10. **Con redondeo dirigido, `-x` y `0 - x` no son la misma operación.** Ibex
   compila con `-frounding-math` y gaol deja el modo de redondeo hacia arriba,
   y ahí el redondeo **no es simétrico respecto del signo**:
   `round_up(-(f·b)) ≠ -round_up(f·b)`. Al escribir el tableau disperso puse
   `v = -f*b` donde el denso calcula `Af[i][c] -= f*Af[row][c]` con
   `Af[i][c] == 0`, o sea `0.0 - f*b`. La diferencia es **de 1 ulp**, pero se
   amplifica pivote a pivote y termina cambiando los desempates del ratio test:
   20 de 148 instancias daban otra secuencia de pivotes. Costó encontrarlo
   porque el test unitario de la representación pasaba con tolerancia 0 —el
   caso no aparecía— y porque la hipótesis natural (dispersión mal mantenida)
   era falsa. La lección general: **al reimplementar una operación numérica,
   hay que reproducir la expresión, no solo su valor matemático.**

---

## 8. Pasos a seguir, en orden

### 8.0 Las dos tareas pedidas: hechas y medidas

**A. Optimizar la Neumaier–Shcherbina de Ibex — hecha.**
`src/numeric/ibex_LPSolver.cpp::neumaier_shcherbina_postprocessing()` acumula
ahora `Aᵀλ − c` **por filas, salteando los `λ` nulos**, en vez de construir
`rows()` (`m×n`) y su transpuesta en cada resolución del LP. Resultado:
**PolytopeHull ×1.66 en el banco de una caja** (1.5906 s → 0.9580 s en 148
instancias) con **contracción e iteraciones idénticas**, y en las instancias
grandes el tiempo con certificación ya coincide con el piso que marcaba
`PH_SKIP_NS`. El 50–67 % que costaba certificar era grasa de implementación.
Detalle en [MEDICIONES_TECHOS.md](MEDICIONES_TECHOS.md) §12.

**En el optimizador de producción rinde poco, pero rinde gratis.** A/B directo
(`acidhc4 --lr=xn` con y sin `IBEX_NS_LEGACY`, 191 instancias, 30 s): **celdas
exactamente iguales** —187 388 en las dos, como debe ser si la contracción no
cambia— y cpu **338.81 s contra 363.95 s**, o sea **×0.93 en total y ×0.973 en
media geométrica** sobre las 137 instancias que resuelve. Las mismas 137
resueltas en las dos.

Es poco, pero es la clase de mejora más limpia que hay: **el árbol es idéntico**,
así que no hay nada que pueda salir mal por sensibilidad caótica. Y beneficia a
todo usuario de `CtcPolytopeHull`, no solo a este trabajo.

Consecuencias que hay que tener presentes al leer el resto del documento:

- La ventaja «DFB certifica 10× más barato» que había encontrado el §11 de
  MEDICIONES_TECHOS **desapareció**. Era el diferenciador que se había
  encontrado para DFB y hay que dejar de contarlo.
- La relación de costo del contractor en el banco de una caja pasa de ×3.9 a
  **×6.5**.
- En cambio la comparación en el árbol **no se movió**: re-medida con la N-S
  optimizada da 137 / 116 / 105 resueltas y ×2.755 celdas / ×1.363 cpu para
  `dfb`, contra ×2.83 / ×1.40 antes ([MEDICIONES_TECHOS.md](MEDICIONES_TECHOS.md)
  §12). O sea que las conclusiones del §5 sobre el árbol siguen en pie: lo que
  cambió es el costo del contractor en aislamiento, no su posición competitiva.
- Quedó a la vista un **agujero de solidez de Ibex**: el residuo `Aᵀλ − c` se
  acumulaba (y por omisión se sigue acumulando) en punto flotante, no en
  intervalos, así que la cota que el solver marca `OptimalProved` no es rigurosa
  en sentido estricto. `IBEX_NS_RIGOROUS=1` activa la variante en intervalos.
  **No se cambió el comportamiento por omisión**: cerrarlo ensancha todas las
  cotas de PolytopeHull y hay que medirlo antes de tocar el solver de
  producción. Es, probablemente, el hallazgo más importante de esta tanda.

**B. El DFB perezoso — implementado, medido, y el resultado va al revés de la
hipótesis.** Detalle en [MEDICIONES_PEREZOSO.md](MEDICIONES_PEREZOSO.md).

Lo primero que apareció: **ya estaba a medias implementado**. Los `2n`
contractores se construían con `max_iters = 1` y al reencolarse vuelven al final
de su grupo, así que el comportamiento vigente ya era un round-robin de un
pivote por cota con HC4 intercalado. Faltaban los mandos, que ahora existen
(`DFB_MAXITERS`, `DFB_PIVOT_BUDGET`/`_N`) y son neutros en sus valores por
omisión.

Los dos resultados, en el banco de una caja **con HC4** (`BENCH_WITH_HC4=1`,
que también es nuevo: antes el banco medía DFB en aislamiento):

1. **El presupuesto parcial es el peor punto de la curva.** La contracción es
   **convexa** en el presupuesto: con 4 n pivotes se ganan 3.42 puntos sobre
   HC4 solo, y al duplicarlos se ganan 13.25. Recortar no compra eficiencia, la
   destruye. La razón es el reparto: un tope de `c·n` con round-robin le da a
   cada cota `c/2` pivotes, y el primer pivote desde la base inicial es el que
   menos aporta.
   Se probaron las dos políticas de reparto y **concentrar es mejor que
   repartir**, pero ninguna le gana a no poner tope: con 4 n pivotes,
   concentrado da 29.95 % contra 28.27 % del reparto fino, y sin tope se llega a
   36.92–38.10 %. Puntos ganados por segundo adicional sobre HC4 solo: 2.59
   (concentrado 4 n), 3.93 (sin tope, 1 pivote/visita), **4.38** (sin tope,
   punto fijo). El máximo de eficiencia está en el extremo «sin tope».
2. **`max_iters = 1` no era el mejor default; el óptimo está en 5.** Esto es lo
   único que sobrevivió de la tarea, va en dirección contraria a la hipótesis, y
   **ya está aplicado**: validado con el protocolo del §6 (32/32 testigos,
   `valgrind` limpio) y con el default cambiado.
   En `ibexopt` (191 instancias, 30 s, 107 comunes, relativo a `max_iters = 1`):

   | `DFB_MAXITERS` | celdas (geom) | cpu (geom) | resuelve |
   |---|---|---|---|
   | 1 (default vigente) | 1.000 | 1.000 | 116 |
   | 2 | 0.961 | 0.944 | 117 |
   | **5** | **0.928** | **0.922** | 116 |
   | 10 | 0.944 | 0.979 | 116 |
   | 20 | 0.946 | 1.095 | 113 |
   | 50 | 0.942 | 1.254 | 108 |
   | −1 (punto fijo) | — | — | **74** (ver abajo) |

   Curva unimodal con máximo en 5 en las dos métricas, **0 óptimos
   incompatibles**, y las instancias resueltas no bajan. **−7 % celdas y −8 %
   tiempo** en media geométrica.
3. **Y el punto fijo por visita (`-1`) no es usable: 72 procesos hubo que
   matarlos.** Contra 1 en todas las demás variantes. Es la trampa de la no
   terminación del §7: sin tope por visita, si el pivoteo en flotantes cicla por
   degeneración `contract_float` no retorna, y el límite de tiempo de `ibexopt`
   se comprueba **entre nodos**, así que no puede actuar. El 38.10 % que el
   banco de una caja le atribuía se midió solo donde una contracción alcanzó a
   terminar. **El tope por visita no es solo una política de reparto: es lo que
   garantiza que el contractor retorne.**

### 8.1 Dónde estamos, en una tabla

| | estado |
|---|---|
| Configuración objetivo | **`--filtering=acidhc4 --lr=dfb`** (DFB en el hueco de `CtcPolytopeHull`, todo lo demás igual a producción) |
| Solidez | cumplida: 1 vaciado corroborado, 32/32 testigos, 0 óptimos incompatibles, `valgrind` limpio |
| Costo por nodo | **a la par**: ×1.069 de cpu contra PolyHull en el mismo hueco; **×0.788 sin ACID**, más rápido en 77/121 instancias |
| Poda | **es lo que falta**: ×1.384 celdas con ACID, ×1.268 sin ACID |
| Resueltas | **123** contra 137 de producción, de 190 |

**La brecha ya no es de velocidad, es de poda.** Eso invierte la prioridad que
tenía este documento durante todo el trabajo anterior, que apuntaba a abaratar
el contractor. Abaratarlo ya se hizo —N-S de Ibex, `max_iters`, tableau
disperso— y alcanzó para llegar a la paridad de tiempo. Lo que queda es que cada
llamada a DFB pode más.

**Por qué es de esperar, y por qué no es un techo duro.** `CtcPolytopeHull`
resuelve el LP **hasta la optimalidad** para cada una de las `2n` cotas, así que
obtiene la mejor cota que la relajación lineal puede dar. DFB con 5 pivotes
obtiene una cota válida pero subóptima. La brecha de poda es entonces función
del presupuesto de pivotes, y eso está medido: subir el tope de 1 a 5 dio −7 %
de celdas. **La curva del tope se midió para `--filtering=dfb`, no para
`--lr=dfb`**, y ahora que el costo está a la par hay margen para gastar más.

### 8.2 Pasos a seguir

Cada paso dice qué decide, para poder abandonarlo si la respuesta sale negativa.

1. **Validar la configuración objetivo** con el protocolo del §6. Los testigos y
   `valgrind` se corrieron sobre `--filtering=dfb`, no sobre `--lr=dfb`, que es
   un camino distinto (DFB con `only_dfb=true` dentro de un `CtcFixPoint`).
   Es media hora y es requisito de todo lo demás.
   **Decide**: si la configuración objetivo es sólida. Sin esto, ningún número
   de más abajo cuenta.

2. **DIAGNOSTICADO: la brecha de poda es la regla de pivote, y está medido con
   tres causas descartadas.** Detalle en [MEDICIONES_PODA.md](MEDICIONES_PODA.md).
   - **No es presupuesto**: con topes 5, 40 y 200 los grupos son idénticos
     (74/25/49) y **ninguna de las 25 instancias del grupo B hace un pivote
     más**; los 22× de pivotes extra se gastan en los grupos A y C.
   - **No es el ratio test**: `n_no_candidate = n_inconclusive = 0` en las 25.
   - **No es el umbral de parada**: bajar `DFB_STOP_TOL` de 1e-6 a 1e-14 mueve
     una instancia.
   - **La corrida que lo cierra**, con el mismo poliedro para los dos
     contractores, presupuesto 200 y umbral 1e-14: **27× más pivotes, los
     mismos grupos**, y DFB sigue por debajo del óptimo del LP en **33 de 148**
     instancias. `ExtendedWood-04` se detiene tras **3** pivotes teniendo 200
     disponibles; `Brent-10` no hace **ninguno** y PolyHull contrae 0.92 %.

   La explicación es que DFB maximiza una función **cóncava lineal a trozos** de
   `λ` cuyo máximo **es** la cota de PolyHull (dualidad LP), con un ascenso
   restringido a `m` direcciones por paso. El ascenso por coordenadas sobre una
   función no suave se detiene en puntos no óptimos.

3. **HECHO: el simplex dual, y la brecha de poda se cierra.** Está en
   [ibex_DFBSimplex.h](ibex_DFBSimplex.h), se activa con `DFB_SIMPLEX_RULE=1`, y el
   detalle está en [MEDICIONES_PODA.md](MEDICIONES_PODA.md) §7.

   No hizo falta fase 1: la linealización pone `A[i][nb_var+i] = −1` exacto, así
   que la base de las `b` da `B⁻¹ = −I` y, con `c = e_k`, costos reducidos
   `d = e_k ≥ 0` — **el punto de partida es dual-factible**.

   Validado contra SoPlex: **400/400 óptimos** con error relativo máximo
   1.5e-15, invariante `γ = yᵀĀ` en todos los casos (que es lo que hace
   certificable la cota) y dualidad fuerte donde hay cota, con 1.5 pivotes por
   cota (`test_dfbsimplex`).

   Resultado en el banco de una caja, **con el mismo poliedro** para los dos
   contractores:

   | | grupo A | grupo B | contracción | pivotes |
   |---|---|---|---|---|
   | regla actual | 64 | **33** | 29.55 % | 29 241 |
   | **simplex** | **95** | **3** | **31.61 %** | 625 277 |

   **De 33 instancias donde DFB contraía menos que PolytopeHull quedan 3.** Y
   en varias lo supera (`Virasoro-icse` 100.00 contra 77.78), porque DFB trabaja
   sobre la caja que las cotas anteriores ya contrajeron. **32/32 testigos
   conservados.**

4. **HECHO el simplex primal; y ahí aparece el motivo real de la
   factorización.** El bucle es ahora **compuesto** —pasos duales mientras haya
   infactibilidad primal, primales mientras un costo reducido viole su signo,
   con *bound flips*— y está validado contra SoPlex **también en el camino
   tibio**: 1 066 cotas en secuencia sobre una sola base, 1 066 coincidencias,
   **1.5 pivotes por resolución**.

   Con warm start (`DFB_SX_WARM=1`) los pivotes se derrumban —Virasoro-icse de
   1 407 a **39**, Brown-20 de 2 571 a **65**— **pero la contracción se cae**
   (Virasoro 100 % → 50 %). Y el motivo está medido: `γ_k` tiene que valer
   exactamente ±1 porque la columna `k` es básica, y en tibio **deriva hasta un
   2.5 %** (`−1.025008`) porque el tableau se mantiene por operaciones de fila
   acumuladas. Como **la certificación usa `λ` y no la cota flotante**, un `λ`
   derivado da una cota válida pero débil —los 32/32 testigos confirman que
   nunca es falsa— y ahí se pierde la ganancia.

5. **PARCIAL: `DFBBasis` cableado. Da el costo, no todavía la poda.** `λ` sale
   ahora de **un solo BTRAN** sobre la factorización, que se mantiene en
   paralelo con `change_basis` por pivote. El mecanismo es el correcto:
   `DiscreteBoundary-0040` con warm start y `DFB_REFACTOR_ALWAYS=1` da **15.63 %
   —exactamente la contracción en frío— con 57 pivotes contra 1 925**, 34 veces
   menos.
   Pero el agregado (147 instancias, mismo poliedro) muestra que falta algo:

   | | grupo A | B | contracción | pivotes |
   |---|---|---|---|---|
   | regla actual | 64 | 33 | 29.55 % | 29 241 |
   | **simplex frío** | **95** | **3** | **31.61 %** | 625 277 |
   | simplex tibio + refact. | 74 | 23 | 23.18 % | **32 130** |

   Los pivotes del tibio quedan **a la par de la regla actual** —el costo está
   resuelto— pero la contracción cae. Y no es solo deriva: refactorizar en cada
   pivote arregla `DiscreteBoundary-0040` del todo y no mueve `Virasoro-icse`
   ni `Brown-20`.

   **La hipótesis del `z_k` no básico se midió y es falsa**: en tibio `z_k` es
   básica *mucho más* seguido (113/120 contra 40/120), o sea que hay más cotas
   certificables, no menos. La pérdida era **debilidad** de las cotas.

   **La causa real, y el error era de arquitectura**: había dos estructuras en
   paralelo —el tableau para decidir pivotes y la factorización para `λ`— y
   refactorizar arreglaba solo la segunda. Pero `d` es el **test de
   optimalidad** y `x_B` el de factibilidad; con los dos derivados del tableau
   el bucle corta antes del óptimo y `λ` —exacto— certifica una cota débil sobre
   una base que no es óptima. En una cota de `Brown-15` devolvía **16**, la
   propia cota inferior, donde el máximo es 6.67e6.

6. **PARCIAL: el refresco exacto resuelve el warm start pero introduce una
   regresión.** Calculando `y`, `d`, la ubicación y `x_B` desde la
   factorización en cada iteración (un BTRAN, un FTRAN y dos productos
   dispersos), **frío y tibio dan resultados idénticos en las 147 instancias**
   con 4× menos pivotes. Pero el grupo B vuelve de 3 a 30, porque la versión
   implementada **rederiva la ubicación de las no básicas por el signo del costo
   reducido en cada iteración**, y eso es ambiguo cuando `d_i = 0` —el caso
   degenerado, que es omnipresente—: la ubicación de una degenerada no cambia el
   objetivo pero **sí** el test de factibilidad primal.

   **HECHA la corrección**: la ubicación de las no básicas se conserva como
   **estado** y se restaura solo al cambiar de objetivo, moviendo únicamente las
   columnas cuyo signo está violado. Sobre 140 instancias comparables:

   | | grupo A | grupo B | contracción | pivotes |
   |---|---|---|---|---|
   | regla actual | 62 | 28 | 29.12 % | 24 897 |
   | compuesto (default) | 88 | **3** | **30.65 %** | 279 431 |
   | exacto frío | 88 | **3** | **30.70 %** | 286 427 |
   | **exacto tibio** | **90** | **1** | 29.39 % | **112 555** |

   `exacto frío` iguala al compuesto, así que la corrección está confirmada. Y
   **`exacto tibio` tiene la mejor clasificación de todas —grupo B en 1— con
   2.5× menos pivotes**, aunque la contracción media baja por 10 instancias que
   pierden la contracción *extra* de las pasadas repetidas sobre la caja ya
   apretada.

   **Faltan dos cosas concretas para promoverlo**, y el default sigue siendo el
   compuesto:
   - **6 `unknown_error`** en `exacto frío` sobre las 153 instancias (el tibio,
     0): hay una excepción que se escapa y hay que encontrarla.
   - **las 10 instancias** donde el tibio pierde contracción, que es la de las
     pasadas repetidas, no la del óptimo del LP (el grupo B en 1 lo confirma).

   **Defecto anotado**: con warm start y la cadena Forrest–Tomlin por omisión
   (período 50), **63 de 153 instancias terminan en `unknown_error`** — la
   actualización de la factorización falla propagando una excepción. No afecta
   el camino por omisión, pero hay que arreglarlo antes de encender el tibio.

7. **PASO SIGUIENTE (era el 5): cerrar la brecha del tibio.** Es
   exactamente el argumento que este plan tenía para la factorización, ahora con
   la medición que lo respalda. El motor está escrito y verificado
   ([ibex_DFBBasis.h](ibex_DFBBasis.h)): `load`, BTRAN, FTRAN y `change()` tipo
   Forrest–Tomlin, equivalencia comprobada fila por fila. Recalcular las filas
   como `λᵀĀ` desde la factorización, en vez de arrastrar las eliminaciones,
   devuelve `γ_k = ±1` y hace que la cota certificada siga a la flotante.
   **Decide**: si DFB queda estrictamente mejor que producción. La poda ya está
   igualada (paso 3, en frío) y el costo por cota ya está resuelto (1.5 pivotes,
   en tibio); **falta tenerlos al mismo tiempo.**

8. **Después: medir la sustitución `--filtering=acidhc4 --lr=dfb` con el
   simplex**, que es la comparación que decide todo, y correrle el protocolo
   del §6 completo.

9. **Lo que quedó descartado del intento anterior: un simplex de variables
   acotadas «desde cero» no era el camino más caro.** El intento barato —usar el criterio de violación del
   simplex dual sobre la base que DFB ya mantiene— **está implementado
   (`DFB_DUAL=1`) y falla**, por un motivo que conviene tener claro: el tableau
   es `[Ā | I]` y la parte `I` no son variables del problema sino la
   contabilidad de `λ`, así que **la base de DFB es artificial en un 80–90 %**
   (entran 2–5 columnas reales de las `m`). No es una base del LP, sino fase 1,
   y el valor de una básica artificial no es la señal de ascenso. Medido: la
   regla ve violación en casi toda fila, no termina —1 326 pivotes contra 65 en
   Katsura-12, saturando el `step_cap`— y aterriza en un `λ` peor (13.33 %
   contra 54.17 % de contracción). Ver §5 de
   [MEDICIONES_PODA.md](MEDICIONES_PODA.md).

   **Lo que el intento sí dejó**, y son dos cosas útiles:
   - **Evidencia de que la brecha es recuperable**: con direcciones distintas,
     `ExtendedWood-04` y `-08` alcanzan **exactamente** la cota de PolyHull
     (12.50) y `Brent-10` pasa de 0.00 a 0.87 contra 0.92. El punto donde la
     regla actual se detiene no es óptimo y hay dirección de mejora.
   - **El seguimiento del mejor `λ` de la secuencia**, ya activado: 0
     instancias peores, 3 mejores, costo nulo. Es requisito de cualquier regla
     de pivote no monótona, porque la certificación se hace una sola vez al
     final.

   **Decide**: si vale la pena escribir un simplex especializado. A favor: DFB
   ya es 21 % más rápido por nodo, así que alcanzar la cota del LP lo dejaría
   estrictamente mejor, y la base factorizada, el tableau disperso y la
   certificación ya están y son independientes del pivoteo. En contra: es
   escribir un simplex de variables acotadas con fase 1, que es lo que hace
   SoPlex.

4. **Mirar las 14 instancias que producción resuelve y la sustitución no**, para
   confirmar que comparten la causa del grupo B y no otra.

5. **La factorización de la base**, con la justificación que quedó después de
   medir: **no abarata el pricing** —no puede, ver §2 de
   [MEDICIONES_FACTORIZACION.md](MEDICIONES_FACTORIZACION.md)— y el tableau
   disperso ya se llevó la ganancia de memoria y de tiempo que se le atribuía.
   Lo que sigue justificándola es **el cambio de objetivo barato**, y eso rinde
   exactamente donde el contractor domina el tiempo: `acid_dfb` (×3.49 de cpu) y
   las instancias grandes. El motor está escrito y verificado
   ([ibex_DFBBasis.h](ibex_DFBBasis.h)): `load`, BTRAN, FTRAN y `change()` tipo
   Forrest–Tomlin, con la equivalencia `Af = B⁻¹[Ā|I]` comprobada fila por fila
   (51 016 comparaciones, 0 fallos).
   **Decide**: si `acid_dfb` puede dejar de costar ×3.5. No toca la brecha de
   poda, así que va después de 2–4.

6. **Que `step_cap` cuente pivotes y no visitas de contractores.** Su
   presupuesto efectivo se multiplica por `max_iters`, así que no acota el
   trabajo de forma robusta, que es para lo que se agregó (§3.1 de
   [MEDICIONES_PEREZOSO.md](MEDICIONES_PEREZOSO.md)). Es un cambio chico pero
   altera los resultados medidos, así que hay que re-correr el banco después.
   **Decide**: nada de rendimiento; es deuda de corrección.

**Escalado: implementado, midió neutro, y la conclusión anterior fue
RETRACTADA.** Había escrito acá que el simplex no era robusto a la escala y que
el escalado era prerrequisito, sobre la base de un test que fallaba 217 de 400.
**El test estaba mal**: generaba coeficientes dispares *dentro de una misma
fila* (hasta 1e36), y medido con [results/spread.cpp](results/spread.cpp) la
dispersión intra-fila real es **8 en Brown-15 y 10.5 en Brown-20** —justo las
instancias donde el tibio falla, o sea que están bien condicionadas—. Ese
régimen inventado además **rompe a SoPlex**, que lanza `SPxStatusException` y
declara infactibles LPs que no lo son, así que la referencia tampoco servía.

Calibrado a lo medido —escalas **por fila** hasta 1e6 y cajas hasta 1e6— el
simplex da **0 discrepancias** con y sin escalado, en frío y en tibio. Y hay un
argumento algebraico que debí hacer antes de escribir el código: el escalado de
**columnas es un no-op** para los términos `a_ji·z_i`, porque
`(r_j a_ji c_i)·(width_i/c_i) = r_j a_ji width_i`, independiente de `c_i`.

El escalado queda implementado (equilibrado en norma infinito, escalas en
**potencias de 2**, columnas de `b` atadas a la fila para preservar `−I`) pero
**APAGADO por omisión** (`DFB_SX_SCALE=1`): midió neutro en el test sintético y
**daña** en el banco real, `Virasoro-icse` de **100 % a 66.67 %**. Misma lección
que con el régimen de escalas del test: lo que no se mide sobre el banco real no
se sabe.

**Sigue abierto**: la pérdida de contracción del tibio en 10 instancias (§8.4 de
[MEDICIONES_PODA.md](MEDICIONES_PODA.md)). Ya se descartaron **por medición** la
escala, la deriva de la factorización, el umbral de reencolado, la hipótesis del
`z_k` no básico y el criterio de optimalidad sobre datos derivados.
8. **`alkyl` aparte**: `acid_dfb` tarda 9.48 s contra 1.36 s de `dfb`, con el
   91.9 % del tiempo dentro del contractor.

**Lo que no priorizaría**, todo con medición detrás:
- **Topes totales de pivotes por caja** (`DFB_PIVOT_BUDGET`/`_N`), en cualquiera
  de sus dos políticas de reparto. La contracción es **convexa** en el
  presupuesto, así que recortar destruye eficiencia en vez de comprarla: en
  `ibexopt` los topes de 1 n, 2 n y 4 n resuelven 77, 80 y 88 instancias contra
  116, con ×1.95, ×1.42 y ×1.00 celdas ([MEDICIONES_PEREZOSO.md](MEDICIONES_PEREZOSO.md)
  §2.1). Queda implementado y apagado.
- **Quitar el tope por visita** (`DFB_MAXITERS=-1`): 72 procesos hubo que
  matarlos, es la trampa de la no terminación del §7.
- **El modo liviano** (abandonar la transformación para abaratar el cambio de
  objetivo): 3× más rápido y contracción ~0 (§10 de
  [MEDICIONES_TECHOS.md](MEDICIONES_TECHOS.md)).
- **El barrido por lotes solo**, sin factorización: ×1.37 celdas y ×1.37 tiempo,
  106 resueltas contra 116. El mecanismo da −30 % de pivotes en aislamiento pero
  no se transfiere, porque pierde el intercalado fino (§9 de
  [MEDICIONES_TECHOS.md](MEDICIONES_TECHOS.md)). Queda implementado tras
  `DFB_SWEEP=1`, listo para reevaluar **con** factorización.
- **Warm start entre nodos**: reusar la linealización cuesta ×1.167 celdas, y
  los pivotes que ahorraría ya son cero (§7 de
  [MEDICIONES_INCREMENTAL.md](MEDICIONES_INCREMENTAL.md)).
- **Priorizar la variable objetivo en la cola**: empeora en dos instancias y es
  neutro en tres (§4 de [MEDICIONES_TECHOS.md](MEDICIONES_TECHOS.md)).
- **La reactivación de DFB desde la propagación HC4**: es el comportamiento
  intencionado y quedó **arreglado y activado** (estaba muerto en el camino
  flotante), pero su efecto es ×1.03 celdas y ×1.03 tiempo — neutro, mejor en
  16 instancias y peor en 21 ([MEDICIONES_INTERACCION.md](MEDICIONES_INTERACCION.md)).
- **Aplicar DFB cada K nodos**: el árbol explota (ex14_2_4 de 7012 a 51 832
  celdas con K=5). El truco sirve para información global durable, como el punto
  factible de Ipopt, y no para un contractor, cuyo beneficio es local a la caja
  (§6 de [MEDICIONES_TECHOS.md](MEDICIONES_TECHOS.md)).
- **devex propiamente dicho** con pesos actualizables: la versión barata ya
  capturó la ganancia de contracción, y esa ganancia no llega al árbol.
- **Explotar la dispersión**: solo instancias grandes, y ahí el problema es
  memoria. Después de la factorización, si acaso.
- **Bound flipping**: el ratio test es el 0.2 % del tiempo.
- **Lista de candidatos**: con pivoteo flotante dejó de aportar en `dfb`
  (×1.40 vs ×1.38); en `acid_dfb` da ×0.969, que es lo único que aún podría
  justificarla.

### Decisiones abiertas

- **~~¿`dfb` o `acid_dfb`?~~ La pregunta estaba mal planteada**: ninguna de las
  dos es la sustitución que el objetivo pide, porque las dos cambian el
  componente **y** la composición. La configuración objetivo es
  **`--filtering=acidhc4 --lr=dfb`**, que deja producción intacta y pone DFB en
  el hueco de `CtcPolytopeHull`. Resuelve 123 contra 118 y 104, y su media
  geométrica de celdas es 1.38 contra 2.38 y 1.44
  ([MEDICIONES_FACTORIZACION.md](MEDICIONES_FACTORIZACION.md) §6). Las
  composiciones alternativas le estaban costando a DFB casi un factor 2 de
  árbol, y eso no era una propiedad del método.
- **¿Se cierra el agujero de solidez de la N-S de Ibex?** El residuo `Aᵀλ − c`
  se acumula en dobles (§8.0 A). `IBEX_NS_RIGOROUS=1` ya implementa la variante
  en intervalos; falta medir cuánto ensancha las cotas y cuánto cuesta antes de
  proponer un cambio en el solver de producción.

### Decisiones ya tomadas (no reabrir sin motivo)

- DFB **reemplaza** a PolytopeHull; no es una alternativa más barata.
- **Pivoteo en flotantes con certificación de la cota final** (§3).
- El vaciado por «ninguna fila bloquea» se **resigna** (era 0 casos).
- Banco de referencia: `data_tests/` completo (153) para el nivel de una caja, y
  para `ibexopt` los conjuntos `easy`, `medium`, `hard`, `blowup` y
  `benchs-minlp` (191), **incluyendo los pesados**.
- Los barridos se corren **en paralelo**, con todos los cores.
- **`max_iters` por visita es 5**, no 1: la curva en `ibexopt` es unimodal con
  máximo ahí (§8.0 B), y quitar el tope no es una opción porque rompe la
  terminación. Validado con el protocolo del §6 (32/32 testigos, `valgrind`
  limpio, 0 óptimos incompatibles) y ya cambiado el default.

---

## 9. Historial: qué se hizo y dónde está documentado

> **Sobre las etiquetas `M0`…`M7`**: venían de la numeración de mejoras del plan
> original, que este documento reemplaza. Sobreviven solo en los nombres y el
> texto de los documentos de resultados, así que se mantienen abajo como
> referencia cruzada para encontrar cada medición. En los pasos a seguir (§8) se
> nombran las cosas por lo que son.

### Defectos corregidos (7)

| defecto | efecto |
|---|---|
| Tolerancia absoluta `1e-5` en el ratio test | descartaba candidatos legítimos y `i == -1` se interpretaba como infactibilidad |
| **Recorte de las cotas `b` a ±1e50** en `CtcDFBPropag::linearize` | **la causa raíz**: apretaba una cota no acotada a un valor finito arbitrario y excluía soluciones |
| Excepción de `makeColumnIdentity` sin atrapar | abortaba `ibexopt` |
| No se chequeaba el `-1` de `lr.linearize` | `resize(-1,...)` → `bad_alloc` |
| La propagación HC4 no corría sin linealización útil | `--filtering=dfb` terminaba con «possibly unbounded objective» |
| `contract` cambiaba el tamaño de la caja del llamador | segfault por memoria liberada en `Optimizer` |
| `old_box[v]` actualizado dentro del bucle equivocado, más `ratiodelta` con cajas de 1e100 | la propagación no terminaba y el límite de tiempo no podía actuar |

### Mejoras aceptadas y activadas

| mejora | resultado |
|---|---|
| **Pivoteo flotante + certificación** | banco 15.62 s → **5.65 s**, regeneraciones 731 → **0**, grupo A 57 → 64 |
| **Pricing normalizado por la norma de la fila** (`M2(b)`) | contracción media 24.03 % → **30.03 %**, grupo A 64 → **72**, +870 pp netos |
| **Desempate de Harris en el ratio test** (`M3(a)`) | +63 pp netos, −3 % de pivotes. Reserva: no funciona por la razón que se había previsto (las regeneraciones subieron, en vez de bajar) |
| **Pricing fusionado** (`M2`) | −15 % de tiempo con resultado **idéntico bit a bit** (solo camino de intervalos) |
| **Una eliminación por variable** (`M5(a)`) | camino de intervalos: `t_init` −19.4 %, resultados idénticos. Camino flotante: `t_init` **−30.6 %**, `t_dfb` **−15.4 %**, equivalente en contracción (difiere ±1 pivote en 22 instancias, ninguna por más de 0.01 pp) |
| **Neumaier–Shcherbina de Ibex por filas, salteando los `λ` nulos** (en `libibex`, no en DFB) | `CtcPolytopeHull` **×1.66** en el banco de una caja con contracción e iteraciones **idénticas**; en `ibexopt` **×0.93 de cpu con el árbol exactamente igual** |
| **`max_iters` por visita de 1 a 5** | −7 % celdas y −8 % tiempo en media geométrica sobre 107 instancias de `ibexopt`, curva unimodal, 0 óptimos incompatibles |
| **Tableau disperso** (`ibex_DFBSparse`) | **resultados idénticos** (0 diferencias en 148 instancias del banco, celdas exactamente iguales en las 190 de `ibexopt`) y ×2.21 el contractor (×3.12 en n ≥ 80); `Eiger-1000` pasa de `bad_alloc` a terminar; en el árbol 3.4 % (`dfb`) y 11.6 % (`acid_dfb`) menos cpu, +2 y +1 instancias resueltas |

### Medidas y descartadas, con datos

| | por qué |
|---|---|
| Lista de candidatos en el pricing | sobre una caja hundía la contracción en `Brown-*`; en el optimizador solo daba 5–6 % de tiempo y con flotantes **dejó de aportar nada** (el pricing ya no es el cuello) |
| Refresco adaptativo de la lista | en `Brown-*` los pivotes de la lista sí mejoran localmente y aun así llevan a peor lugar: la señal nunca dispara |
| Compartir una base entre los `2n` contractores (la forma original de `M5(a)`) | se usan todos y divergen al primer pivote. Con una factorización sí es viable (paso 1 del §8) |
| **Warm start entre nodos** (reusar la linealización) | la relajación del padre es más floja: ×1.167 celdas. Y los pivotes que ahorraría son cero |
| **DFB perezoso con tope total de pivotes por caja** | la contracción es **convexa** en el presupuesto: recortar destruye eficiencia. En `ibexopt`, topes de 1 n / 2 n / 4 n resuelven 77 / 80 / 88 contra 116 |
| **Quitar el tope de pivotes por visita** (`max_iters = -1`) | contrae más en aislamiento pero **72 procesos hubo que matarlos**: sin tope, el pivoteo flotante puede ciclar y `contract_float` no retorna |

### Progreso del banco de una caja

| etapa | ok | vaciadas | A | B | C | contracción | tiempo |
|---|---|---|---|---|---|---|---|
| original | 139 | **89** (46 erróneas confirmadas) | 2 | 35 | 13 | 11.17 % | 6.35 s |
| tras la corrección de solidez | 147 | 2 | 56 | 41 | 48 | 24.13 % | 18.31 s |
| tras Harris, fusionado y par lb/ub | 147 | 2 | 57 | 40 | 48 | 24.56 % | 15.62 s |
| tras flotantes | 147 | 2 | 64 | 33 | 48 | 24.03 % | **5.65 s** |
| tras la N-S optimizada (no cambia DFB) | 148 | 2 (correctas) | 72 | 25 | 49 | 29.94 % | 6.25 s |
| tras `max_iters = 5` | 148 | 1 (correcta) | 73 | 25 | 49 | 30.29 % | 6.73 s |
| **final, con tableau disperso** | **149** | **1** (correcta) | **74** | 25 | 49 | 30.29 % | **3.02 s** |

La corrida final está en `results/baseline_153_sparse.csv`. El vaciado que se
pierde entre las dos últimas filas es `Prolog-icse`, y es pérdida de poda, no de
solidez: era un vaciado correcto y la detección de vacío por la cota depende de
la secuencia de pivotes (§3 de [MEDICIONES_PEREZOSO.md](MEDICIONES_PEREZOSO.md)).

### En `ibexopt` (132 instancias, límite de 20 s)

| config | resueltas | celdas × producción | tiempo × producción |
|---|---|---|---|
| `acidhc4+ph` (producción) | 116 | 1.00 | 1.00 |
| `dfb` | 103 | 2.48 | **1.31** |
| `acid_dfb` | 96 | **1.48** | 3.21 |
| `acidhc4` (sin relajación lineal) | 60 | 2.58 | 1.35 |

DFB supera claramente al escenario sin relajación lineal (103 contra 60
resueltas). Antes de los arreglos, `dfb` resolvía 62 con 15 abortos.

### Documentos

| documento | contenido |
|---|---|
| [HALLAZGOS_ANALISIS.md](HALLAZGOS_ANALISIS.md) | línea base del código original: 89 cajas vaciadas, 46 erróneas confirmadas con `ibexsolve` |
| [CORRECCION_M0.md](CORRECCION_M0.md) | la causa raíz de la falta de solidez y su corrección |
| [MEDICIONES_M2_M3.md](MEDICIONES_M2_M3.md) | perfil por fase, desempate de Harris, pricing fusionado, lista de candidatos, una eliminación por variable |
| [MEDICIONES_IBEXOPT.md](MEDICIONES_IBEXOPT.md) | los cinco defectos del optimizador y la comparación contra PolytopeHull |
| [MEDICIONES_FLOAT.md](MEDICIONES_FLOAT.md) | pivoteo flotante con certificación, y pricing normalizado |
| [MEDICIONES_INCREMENTAL.md](MEDICIONES_INCREMENTAL.md) | la métrica incremental (pivotes tibio vs frío) y por qué el warm start entre nodos se descartó |
| [MEDICIONES_TECHOS.md](MEDICIONES_TECHOS.md) | qué fracción del tiempo del optimizador es el contractor, el techo de cada mejora, y la descomposición del árbol |
| [MEDICIONES_INTERACCION.md](MEDICIONES_INTERACCION.md) | las dos vías de reactivación de DFB, el error que mató una en flotantes, y cuánto rinde la interacción con la propagación |
| [MEDICIONES_PODA.md](MEDICIONES_PODA.md) | por qué DFB poda menos: se traba en un punto no óptimo, con presupuesto, ratio test y umbral descartados por medición |
| [MEDICIONES_FACTORIZACION.md](MEDICIONES_FACTORIZACION.md) | la equivalencia `Af = B⁻¹[Ā\|I]` verificada, por qué el pricing impide calcular filas a demanda, el tableau disperso y la comparación en igualdad de condiciones |
| [MEDICIONES_PEREZOSO.md](MEDICIONES_PEREZOSO.md) | el reparto del presupuesto de pivotes: por qué el tope total fracasa, por qué quitar el tope por visita rompe la terminación, y por qué el óptimo está en 5 |
| [README.md](README.md) | qué hace el ejemplo `example_dfb` |

Datos crudos y scripts en [results/](results/).
