# La brecha de poda: DFB se traba, y no es por presupuesto

Con la sustitución directa medida ([MEDICIONES_FACTORIZACION.md](MEDICIONES_FACTORIZACION.md)
§6) la situación quedó así: **el costo por nodo está a la par o mejor** —DFB es
21 % más rápido que PolytopeHull en el mismo hueco sin ACID— y lo que separa a
DFB de producción es que **poda menos**: ×1.27–1.38 celdas, 123 instancias
resueltas contra 137.

Este documento contesta por qué.

## 1. El punto de partida teórico

DFB obtiene su cota de una combinación lineal de las filas de la relajación:
con multiplicadores `λ`, la fila `γ = λᵀĀ` da una relación válida y la cota sale
de evaluarla por intervalos sobre la caja. El valor de esa cota, como función de
`λ`, es

    f(λ) = − Σ_{i≠k} ub( γᵢ(λ) · [zᵢ] )     sujeto a γ_k(λ) = 1

que es **cóncava y lineal a trozos** (una suma de `−max` de funciones lineales).
Maximizarla es exactamente el dual del LP «minimizar `z_k` sobre la relajación
∩ caja», o sea que **su máximo es la cota que calcula `CtcPolytopeHull`**. No
hay un límite estructural: si DFB llegara al óptimo de su propia búsqueda,
igualaría a PolyHull.

El pivoteo de DFB es un ascenso sobre `f` restringido a **`m` direcciones por
paso**: agregar `α·fila_j` a la fila 0, con `α` del ratio test. Es decir, un
ascenso por coordenadas. Y el ascenso por coordenadas sobre una función **no
suave** se detiene en puntos no óptimos: es su fallo conocido.

## 2. Tres causas descartadas por medición

Banco de una caja, 148 instancias comparables. Clasificación habitual: grupo A
= DFB contrae ≥ que PolyHull, B = contrae menos, C = ninguno contrae.

### 2.1 No es el presupuesto de pivotes

| tope por visita | grupo A | B | C | contracción | pivotes |
|---|---|---|---|---|---|
| 5 | 74 | 25 | 49 | 30.26 % | 65 312 |
| 40 | 74 | 25 | 49 | 30.17 % | 314 712 |
| 200 | 74 | 25 | 49 | 30.17 % | 1 453 432 |

Los grupos son **idénticos** y la contracción no se mueve, con 22× más pivotes.
Y desglosado por grupo, los pivotes extra no van a donde falta:

| grupo | inst. | pivotes @5 | @200 | | instancias que pivotean más |
|---|---|---|---|---|---|
| A | 74 | 54 868 | 1 318 168 | ×24 | 18/74 |
| **B** | **25** | **5 446** | **5 403** | **×1.0** | **0/25** |
| C | 49 | 4 998 | 129 861 | ×26 | 2/49 |

**Ninguna de las 25 instancias del grupo B hace un solo pivote más** al subir el
tope de 5 a 200. Llegan a un punto fijo y se detienen. Los 22× se gastan en el
grupo A, donde DFB ya igualaba, y en el C, donde nadie contrae.

### 2.2 No es el ratio test

En las 25 instancias del grupo B: `n_no_candidate = 0` y `n_inconclusive = 0`.
DFB **nunca** se detiene por no encontrar fila de bloqueo; se detiene siempre
porque el pricing no encuentra ninguna fila con mejora estimada por encima del
umbral.

### 2.3 No es el umbral de parada

El umbral estaba fijo en `1e-6` (ahora es `DFB_STOP_TOL`):

| umbral | grupo A | B | C | contracción |
|---|---|---|---|---|
| 1e-6 | 74 | 25 | 49 | 30.17 % |
| 1e-10 | 74 | 25 | 49 | 30.20 % |
| 1e-14 | 75 | **24** | 49 | 30.20 % |

Bajarlo **ocho órdenes de magnitud** mueve una instancia. De las 25 del grupo B,
solo dos mejoran: `Geneig` cierra del todo (0.00 → 0.06, igual que PolyHull) y
`Virasoro-icse` pasa de 66.67 a 72.22 con PolyHull en 94.44.

## 3. La corrida que lo cierra

Todo controlado a la vez: **el mismo poliedro** para los dos contractores
(`BENCH_CORNER=inf`, ver §4), presupuesto 200 y umbral 1e-14.

| variante | grupo A | B | C | contracción | pivotes |
|---|---|---|---|---|---|
| tope 5, umbral 1e-6 | 65 | 33 | 50 | 29.52 % | 31 240 |
| tope 200, umbral 1e-14 | **65** | **33** | **50** | 29.80 % | **855 394** |

**27× más pivotes, los mismos grupos.** Con el mismo poliedro, presupuesto
prácticamente ilimitado y umbral en el ruido de la máquina, DFB sigue por debajo
del óptimo del LP en **33 de 148 instancias**.

Los casos extremos no dejan lugar a dudas:

| instancia | DFB | PolyHull | brecha | pivotes que hizo |
|---|---|---|---|---|
| Virasoro-icse | 72.22 | 77.78 | 5.56 | 256 |
| Brown-20 | 84.79 | 90.25 | 5.46 | 727 |
| Redeco9 | 4.86 | 9.72 | 4.86 | **36** |
| ExtendedWood-04 | 8.33 | 12.50 | 4.17 | **3** |
| ExtendedWood-08 | 8.33 | 12.50 | 4.17 | **6** |
| Discrete-Integralf2-12 | 6.37 | 9.50 | 3.14 | 127 |
| **Brent-10** | **0.00** | 0.92 | 0.92 | **0** |

`ExtendedWood-04` se detiene tras **3 pivotes** teniendo 200 disponibles.
`Brent-10` no hace **ninguno**: desde la base inicial, ninguna dirección de una
sola fila muestra mejora estimada por encima de `1e-14`, y sin embargo el LP
tiene una cota mejor. Es la demostración más limpia posible de que el criterio
de parada de DFB no es el criterio de optimalidad del LP.

## 4. Un defecto del arnés que conviene conocer

`bench_dfb` construía **dos linealizadores distintos** para DFB y para
PolytopeHull, los dos con política de esquina `RANDOM`. Con `RANDOM` cada objeto
sortea esquinas distintas, así que **los dos contractores no resolvían el mismo
poliedro** y la comparación de contracción mezclaba dos efectos.

Se agregó `BENCH_CORNER=inf` (o `sup`) para que los dos usen la misma política
determinista. El efecto sobre la clasificación no es despreciable:

| política | grupo A | B | C |
|---|---|---|---|
| `RANDOM` (dos poliedros) | 74 | 25 | 49 |
| `INF` (el mismo poliedro) | **65** | **33** | 50 |

O sea que **con el mismo poliedro DFB queda peor**, no mejor: parte del grupo A
que se venía reportando era suerte de linealización. Todos los números de
grupos A/B/C anteriores a esta medición hay que leerlos con ese ruido encima.
La comparación en `ibexopt` **no** está afectada, porque ahí cada corrida usa un
solo contractor.

## 5. El primer intento de arreglo, y por qué falló

La idea era directa: si el criterio de parada de DFB no es el del LP, ponerle el
del LP. En términos de simplex, `γ = λᵀĀ` **es** el vector de costos reducidos,
y la cota de DFB es el valor objetivo de la solución básica con las no básicas
en la cota que dicta el signo de su costo reducido. De ahí parecía seguir que el
criterio correcto era el del simplex dual: elegir la variable básica que más
viola sus cotas, y terminar cuando ninguna viola.

Se implementó (`DFB_DUAL=1`, con la asignación **única** de las no básicas en
vez de la esquina más favorable por fila) y **no funciona**. El motivo es un
error en la identificación, y conviene dejarlo escrito:

**La base de DFB es casi toda artificial.** El tableau es `[Ā | I]` y la parte
`I` no son variables del problema: es la contabilidad de `λ`. La base arranca
siendo esa identidad artificial, y como DFB hace 2–5 pivotes por cota, entran
apenas 2–5 columnas reales de las `m` de la base. O sea que **la base de DFB no
es una base del LP**, sino una situación de fase 1, y el valor de una básica
artificial no es la señal de ascenso para maximizar `f(λ)`: es el objetivo de la
fase 1.

Consecuencia medida: como casi toda fila tiene básica artificial y su valor
casi nunca es 0, la regla **ve una violación en casi todas las filas y no
termina nunca**. Sobre Katsura-12 con tope 1 por visita hace **1 326 pivotes**
contra 65 de la regla actual —o sea que satura el `step_cap`— y aterriza en un
`λ` malo: 13.33 % de contracción contra 54.17 %.

### 5.1 Lo que sí quedó del intento

**Primero, evidencia de que la brecha es recuperable.** Con direcciones
distintas, tres de las instancias que se trababan alcanzan la cota del LP:

| instancia | regla actual | regla dual | PolyHull |
|---|---|---|---|
| ExtendedWood-04 | 8.33 | **12.50** | 12.50 |
| ExtendedWood-08 | 8.33 | **12.50** | 12.50 |
| Brent-10 | 0.00 | **0.87** | 0.92 |

En esas, el punto donde la regla actual se detiene **no** es óptimo y hay una
dirección de mejora que otra regla encuentra. El estancamiento no es un techo.

**Segundo, una mejora que se queda: guardar el mejor `λ` de la secuencia.** La
certificación se hace una sola vez al final, con el `λ` que quedó, así que el
método tenía que ser monótono en la cota para no perder terreno. El pricing por
impacto lo es aproximadamente; cualquier regla que se mueva por optimalidad del
LP no lo es. Guardar el mejor `λ` visto lo vuelve monótono por construcción:

| | contracción media | instancias peores | mejores | tiempo |
|---|---|---|---|---|
| sin seguimiento | 30.72 % | — | — | 32.555 s |
| **con seguimiento** | **30.73 %** | **0** | 3 | 32.687 s (×1.00) |

Cero instancias empeoran, tres mejoran (`Brown-19`, `Brown-20-1`,
`Geneig-icse`) y el costo es nulo. Queda **activado**, y es requisito de
cualquier intento futuro de cambiar la regla de pivote.

## 6. Qué hacer con esto

El diagnóstico es específico y, con el intento del §5 descartado, el camino que
queda es más caro de lo que parecía: **un simplex de variables acotadas de
verdad sobre este LP**, con una base real de `Ā` y una fase 1 para llegar a
ella. No alcanza con cambiar el criterio de selección sobre la base artificial
que DFB mantiene hoy.

Lo que ya está y sirve para eso:

- **La base y su factorización**: verificado que el pivoteo de DFB es un cambio
  de base, `Af = B⁻¹[Ā|I]`, y el motor con BTRAN, FTRAN y `change()` está escrito
  y probado ([ibex_DFBBasis.h](ibex_DFBBasis.h)).
- **El tableau disperso**, que ya dio ×2.2 y arregló la memoria.
- **La certificación**, que es independiente del pivoteo.

- **El seguimiento del mejor `λ`** (§5.1), sin el cual cualquier regla no
  monótona pierde terreno.

Lo que falta es una **base real y una regla de pivote sobre ella**: hoy la base
es artificial en un 80–90 % y el criterio de optimalidad del LP no se puede
evaluar sobre ella.

Y conviene ser explícito sobre lo que eso implica y lo que no. Implica que DFB
se parezca más a un simplex especializado para esta familia de LPs. **No**
implica que pierda sentido: DFB ya es **21 % más rápido por nodo** que
PolytopeHull en el mismo hueco, así que obtener la misma cota lo volvería
estrictamente mejor. La pregunta pasa de «¿puede DFB competir?» a «¿se puede
arreglar la regla de pivote sin perder el costo por nodo?».


---

## 7. La regla de pivote de simplex: la brecha de poda se cierra

**Una aclaración de nombre antes de los números**, porque «DFB + simplex» se lee
como si fueran dos métodos pegados y no lo son. DFB siempre tuvo la **maquinaria**
de un simplex —un tableau que pivotea, una regla de pricing, un ratio test,
filas que se combinan por eliminación— pero no su **criterio**:

1. **Su base no es una base del LP.** El tableau es `[Ā | I]` y la parte `I` no
   son variables del problema sino la contabilidad de `λ`. Como hace 2–5 pivotes
   por cota, la base es **artificial en un 80–90 %**: vive en una situación de
   fase 1 permanente.
2. **Sacrifica una fila de restricción** para llevar `γ`, en vez de tener una
   fila objetivo aparte.
3. **Su criterio de parada es una heurística.** Elige la fila que maximiza un
   «impacto» calculado con evaluación de intervalos y se detiene cuando ninguna
   muestra mejora estimada. Eso es **ascenso por coordenadas**, y sobre una
   función cóncava lineal a trozos se detiene en puntos no óptimos.

Lo que sigue reemplaza **esa regla**, no el método: la certificación en
intervalos, el pivoteo flotante y la base compartida por las `2n` cotas se
conservan, y son justamente lo que hace que cueste ×0.754 por nodo contra
`CtcPolytopeHull`.



Escrito en [ibex_DFBSimplex.h](ibex_DFBSimplex.h), activado con `DFB_SIMPLEX_RULE=1`.

### 7.1 Por qué no hace falta fase 1

La linealización de DFB pone `A[i][nb_var+i] = -1` **exactamente**
([ibex_CtcDFBPropag.cpp](ibex_CtcDFBPropag.cpp), bloque de coeficientes), o sea
que la parte `b` de `Ā` es `−I`. Tomando como base las `m` columnas de `b`:

- `B = −I`, así que `B⁻¹ = −I` y el tableau inicial es `−[Ā | I]`: no hay nada
  que factorizar;
- el objetivo es `c = e_k` con `k < nb_var`, o sea `c_B = 0`, de donde `y = 0` y
  los costos reducidos son `d = e_k ≥ 0`. Con **todas** las no básicas en su cota
  inferior, el punto de partida es **dual-factible**.

O sea que se arranca dual-factible y se itera el simplex dual hasta factibilidad
primal, que es la condición de optimalidad. Sin fase 1.

El tableau tiene `m+1` filas: la 0 es la fila objetivo `[d | −y]`, aparte de las
`m` de restricción. Esa es la diferencia estructural con el DFB actual, que
sacrifica una fila de restricción para llevar `γ` y por eso su base es casi toda
artificial (§5).

### 7.2 Validado contra SoPlex

[test_dfbsimplex.cpp](test_dfbsimplex.cpp) genera relajaciones aleatorias con la
misma estructura, resuelve con las dos y compara, mezclando cotas inferiores y
superiores:

    optimo coincide con SoPlex   : 400 / 400   (peor error relativo 1.5e-15)
    discrepancias                : 0
    limite de iteraciones        : 0
    invariante gamma = y^T*Abar  : 400 ok, 0 mal
    dualidad fuerte donde hay cota: 103 ok, 0 mal
    pivotes: 598 en 400 resoluciones -> 1.5 por cota

El invariante `γ = yᵀĀ` es el que hace certificable la cota: con él, `λ` sirve
para recalcular `γ` **en intervalos** y evaluarlo sobre la caja, que es lo que
preserva la solidez. Un `λ` malo da una cota débil, nunca una falsa.

### 7.3 Resultado: el grupo B se derrumba

Banco de una caja, **con el mismo poliedro** para los dos contractores
(`BENCH_CORNER=inf`), 147 instancias comparables:

| | grupo A | grupo B | grupo C | contracción | `t_dfb` | pivotes |
|---|---|---|---|---|---|---|
| regla actual | 64 | **33** | 50 | 29.55 % | 2.19 s | 29 241 |
| **simplex** | **95** | **3** | 49 | **31.61 %** | 63.48 s | 625 277 |

**De 33 instancias donde DFB contraía menos que PolytopeHull quedan 3.** 35
mejoran, 1 empeora. Con esquina `RANDOM` el efecto es el mismo en dirección
aunque menor en magnitud (B de 25 a 17, A de 74 a 83).

Casos donde el simplex no solo iguala sino que **supera** a PolyHull, porque DFB
trabaja sobre la caja que las cotas anteriores ya contrajeron:

| instancia | regla actual | simplex | PolyHull |
|---|---|---|---|
| Virasoro-icse | 66.67 | **100.00** | 77.78 |
| Discrete-Integralf2-19 | 58.19 | **90.83** | 4.49 |
| Brown-30sp | 22.50 | **86.62** | 0.00 |

**Solidez: 32/32 testigos conservados**, y los 7 casos de `tests_dfb_all` pasan
con los dos `NO SOLUTION` detectados vacíos. Era la verificación imprescindible:
una contracción que de golpe es mucho mayor es justo cuando hay que comprobar
que no se perdió una solución.

### 7.4 Tres arreglos numéricos que costaron encontrar

Ninguno es del algoritmo, los tres son de tolerancias absolutas donde tenían que
ser relativas. Vale anotarlos porque son la clase de error que no se ve en un
test sintético:

1. **Tolerancia de pivote relativa a la fila.** Con coeficientes de ~1e72 —las
   linealizaciones de `Brown-*` sobre cajas de radio 1e9— un umbral absoluto de
   1e-9 acepta pivotes que son ruido, y el ratio test se queda sin candidatas y
   **declara infactible un LP factible**: pasaba en **40 de 40** resoluciones de
   `Brown-20`.
2. **Tolerancia de factibilidad relativa a la magnitud de la variable.** Con
   variables de magnitud 1e9, el error de redondeo relativo de 1e-16 da una
   violación absoluta de 1e-7, diez veces la tolerancia: el simplex pivotea
   sobre ruido hasta agotarse.
3. **Desempate de Harris en el ratio test dual.** Sin él cicla en instancias
   degeneradas: `ExtendedWood-04` bajó de **88 628 pivotes a 17**.

Y un cuarto, del lado del contractor: **hay que reencolarse mientras rinda.** El
óptimo del LP es el mejor posible *para esta caja*, pero la propagación sigue
apretando otras variables, así que volver a resolver sobre la caja nueva da una
cota mejor. Sin eso el simplex hacía una sola pasada donde la regla vieja hace
varias, y en las `Brown-*` eso le costaba más de lo que ganaba.

### 7.5 Lo que falta: el warm start, y por qué el obvio no sirve

El costo es **×21 pivotes**, y ahí está todo el trabajo que queda. La causa es
que cada cota se resuelve desde la base de las `b`, sin reusar nada, mientras
`CtcPolytopeHull` gasta **1.08 iteraciones por cota** reusando la base.

El intento obvio está implementado (`DFB_SX_WARM=1`) y es **peor**:

| instancia | frío | tibio |
|---|---|---|
| DiscreteBoundary-0040 | 1 925 | 8 864 |
| Katsura-12 | 156 | 286 |
| Virasoro-icse | 1 407 | 2 608 |

El motivo, que conviene tener claro: en un LP de variables acotadas la
factibilidad **dual** es solo una condición de signo, y uno es libre de elegir
en qué cota se apoya cada no básica, así que reubicarlas por el signo del costo
reducido nuevo da un punto dual-factible sin pivotear. Pero eso **destruye la
factibilidad primal**, y la base óptima para otro objetivo queda lejos de ser
primal-factible; la base de las `b`, en cambio, tiene una infactibilidad chica y
natural.

**El warm start que sirve conserva la base Y la ubicación de las no básicas**
—o sea que se queda primal-factible— **y corre un simplex primal** para el
objetivo nuevo. Eso está implementado: el bucle es ahora **compuesto**, da pasos
duales mientras haya infactibilidad primal y primales mientras algún costo
reducido viole su signo, y termina cuando no pasa ninguna de las dos, que es
optimalidad. Incluye *bound flips* (cuando lo que primero se agota es el rango
de la entrante, no hay pivote, solo cambia de cota).

Validado contra SoPlex **también en el camino tibio**, que es el que faltaba:
secuencias de las `2n` cotas sobre una sola base, **1 066 cotas comparadas,
1 066 coincidencias**, con 966 resoluciones tibias y **1.5 pivotes por
resolución**.

### 7.6 Y ahí aparece el motivo real por el que hace falta la factorización

Con el warm start los pivotes se derrumban, tal como se esperaba:

| instancia | pivotes en frío | tibio |
|---|---|---|
| Virasoro-icse | 1 407 | **39** |
| Brown-20 | 2 571 | **65** |
| DiscreteBoundary-0040 | 1 925 | **222** |
| Katsura-12 | 156 | **111** |

**Pero la contracción se cae**: Virasoro 100 % → 50 %, Brown-20 89.54 % →
11.17 %, Katsura-12 54.17 % → 5.00 %. Y no porque las cotas flotantes sean
peores —varias son *mejores*, `k=4 up=1` da 0.1098 tibio contra 0.3999 en frío—
sino por esto, que se ve en la traza (`DFB_SX_TRACE=1`):

    frio :  gamma_k = -1.000000   -1.000000   -1.000000   ...
    tibio:  gamma_k = -0.999999   -1.025008    0.000001   ...

`γ_k` **tiene que valer exactamente ±1**, porque la columna `k` es básica y su
costo reducido es 0. En frío lo vale; en tibio deriva hasta un **2.5 %**.

La causa es que el tableau se mantiene por operaciones de fila acumuladas, sin
refactorizar. Y el punto clave es **por qué se nota tanto**: la certificación no
usa la cota flotante, usa `λ`. Un `λ` derivado sigue siendo válido —la cota
certificada nunca es falsa, y los 32/32 testigos lo confirman— pero da una cota
**débil**, y ahí se pierde toda la ganancia.

**Esto es exactamente el argumento que el plan tenía para la factorización, y
ahora con la medición que lo respalda.** El motor está escrito y verificado
([ibex_DFBBasis.h](ibex_DFBBasis.h)): `load`, BTRAN, FTRAN y `change()` tipo
Forrest–Tomlin, con la equivalencia `Af = B⁻¹[Ā|I]` comprobada fila por fila.
Lo que falta es usarlo para **recalcular el tableau de la base tibia** en vez de
arrastrar las operaciones de fila, con lo que `γ_k` volvería a valer ±1 y la
cota certificada seguiría a la flotante.

### 7.7 La factorización cableada: da el costo, no todavía la poda

Se conectó `DFBBasis` dentro del simplex: la base se mantiene en paralelo con
`change_basis` por pivote y `λ` se obtiene con **un solo BTRAN** sobre la
factorización (`y = c_k ·` fila `r` de `B⁻¹`, con `r` la fila donde `k` es
básica; si `k` es no básica, `c_B = 0` y `y = 0`).

Un caso muestra que el mecanismo es el correcto: `DiscreteBoundary-0040` con
warm start y **refactorización en cada cambio de base** (`DFB_REFACTOR_ALWAYS=1`)
da **15.63 %, exactamente la contracción en frío, con 57 pivotes contra 1 925**
— 34× menos. Y con la cadena de updates Forrest–Tomlin por omisión (período 50)
da 15.09 %: la deriva es real y refactorizar la elimina.

Pero el agregado dice que falta algo más. Banco de una caja, mismo poliedro,
147 instancias:

| | grupo A | B | C | contracción | `t_dfb` | pivotes |
|---|---|---|---|---|---|---|
| regla actual | 64 | 33 | 50 | 29.55 % | 2.19 s | 29 241 |
| **simplex frío** | **95** | **3** | 49 | **31.61 %** | 63.48 s | 625 277 |
| simplex tibio + refact. | 74 | 23 | 50 | 23.18 % | 12.60 s | **32 130** |

**El warm start entrega el costo pero no la poda.** Los pivotes bajan 19× y
quedan **a la par de la regla actual** (32 130 contra 29 241), y el tiempo baja
5×; pero la contracción cae a 23.18 %, peor que frío y peor que la regla actual,
con 47 instancias por debajo del frío y 1 por encima.

Y la causa **no es solo la deriva**: refactorizar en cada pivote arregla
`DiscreteBoundary-0040` por completo pero no mueve `Virasoro-icse` (50 % contra
100 % en frío) ni `Brown-20` (11.17 % contra 89.54 %).

### 7.8 La hipótesis del `z_k` no básico: medida y refutada

La sospecha era que al conservar la ubicación de las no básicas el tibio
aterrizara en otra base óptima donde `z_k` quedara **no básica**, con lo que
`c_B = 0`, `y = 0` y no habría cota que certificar. Medido, es al revés:

| instancia | frío (`z_k` básica / no) | tibio |
|---|---|---|
| Virasoro-icse | 26 / 26 | **37 / 6** |
| Brown-20 | 38 / 39 | **42 / 2** |
| DiscreteBoundary-0040 | 40 / 80 | **113 / 7** |

En tibio `z_k` es básica **mucho más** seguido, o sea que hay *más* cotas
certificables, no menos. La pérdida no era falta de cotas: era **debilidad** de
las cotas.

### 7.9 La causa real, y una simplificación que resultó incorrecta

El error era mío y de arquitectura: había **dos** estructuras en paralelo —el
tableau para decidir pivotes y la factorización para `λ`— y refactorizar
arreglaba solo la segunda. Pero `d` no es únicamente un criterio de selección:
es el **test de optimalidad**, y `x_B` el de factibilidad. Con los dos derivados
del tableau el bucle corta antes del óptimo, y entonces `λ` —exacto— certifica
una cota válida pero **débil** sobre una base que no es óptima. En una cota de
`Brown-15` eso devolvía **16**, que es la propia cota inferior, donde el máximo
es 6.67e6.

Se implementó el refresco exacto en cada iteración (`y`, `d`, ubicación y `x_B`
desde la factorización: un BTRAN, un FTRAN y dos productos dispersos) y **eso
resolvió el warm start del todo**: frío y tibio dan resultados **idénticos en
las 147 instancias**, con 4× menos pivotes.

Pero el agregado dice que el algoritmo así es **peor**:

| | grupo A | B | C | contracción | pivotes |
|---|---|---|---|---|---|
| regla actual | 64 | 33 | 50 | 29.55 % | 29 241 |
| **simplex compuesto (tableau)** | **95** | **3** | 49 | **31.61 %** | 625 277 |
| simplex de refresco exacto, frío | 67 | 30 | 50 | 29.75 % | 145 136 |
| simplex de refresco exacto, tibio | 67 | 30 | 50 | 29.75 % | 135 188 |

El grupo B vuelve de 3 a 30. La causa es una simplificación que parecía elegante
y no lo es: el refresco **reubica todas las no básicas por el signo de su costo
reducido en cada iteración**, y eso es ambiguo cuando `d_i = 0`, que es el caso
**degenerado y omnipresente**. La ubicación de una columna degenerada no cambia
el objetivo pero **sí cambia el test de factibilidad primal**, así que el
criterio de optimalidad queda mal planteado: se declara óptimo un punto que con
otra ubicación de las degeneradas sería infactible, y al revés.

Un simplex dual de verdad **conserva la ubicación como estado** —la variable que
sale va a la cota que violaba— en vez de rederivarla de los signos. Eso es lo
que falta: calcular `x_B` y `d` exactos desde la factorización **manteniendo la
ubicación como estado**. Es una corrección acotada sobre lo que ya está.

Por eso el camino por omisión vuelve a ser el **compuesto sobre el tableau**,
que es el que mide mejor, y el refresco exacto queda tras `DFB_SX_EXACT=1`.

### 7.10 Un defecto anotado

Con warm start y la cadena Forrest–Tomlin por omisión (período 50, sin
`DFB_REFACTOR_ALWAYS=1`), **63 de 153 instancias terminan en
`unknown_error`**: la actualización de la factorización falla de un modo que
propaga una excepción. No afecta el camino por omisión —el warm start está
apagado— pero hay que arreglarlo antes de encender nada de eso.

### 7.11 La ubicación como estado: la corrección que faltaba

Un simplex dual **conserva la ubicación de las no básicas como estado** —la que
sale va a la cota que violaba— y la restaura solo al **cambiar de objetivo**,
moviendo únicamente las columnas cuyo signo está violado. El refresco exacto la
rederivaba de los signos en **cada iteración**, y de ahí venía la regresión.

Corregido eso (reubicación mínima, solo al entrar, dejando las degeneradas donde
están —que además es lo que preserva la factibilidad primal heredada por el warm
start—), sobre 140 instancias comparables:

| | grupo A | grupo B | grupo C | contracción | pivotes |
|---|---|---|---|---|---|
| regla actual | 62 | 28 | 50 | 29.12 % | 24 897 |
| compuesto (default) | 88 | **3** | 49 | **30.65 %** | 279 431 |
| exacto frío | 88 | **3** | 49 | **30.70 %** | 286 427 |
| **exacto tibio** | **90** | **1** | 49 | 29.39 % | **112 555** |

Dos cosas:

1. **`exacto frío` iguala al compuesto** (88/3, 30.70 contra 30.65). O sea que
   la corrección está confirmada en el agregado, no solo en las cuatro
   instancias de prueba: el grupo B en 30 venía de reubicar demasiado.
2. **`exacto tibio` tiene la mejor clasificación de todas —grupo B en 1— con
   2.5× menos pivotes.** Pero la contracción media baja a 29.39 %, arrastrada
   por **10 instancias** que pierden (`Brown-15` 53.11 contra 93.56,
   `Virasoro-icse` 88.89 contra 100.00), contra 2 que mejoran.

Que el grupo B baje a 1 y la media baje a la vez no es contradictorio: el grupo
mide *cuántas veces DFB alcanza a PolyHull* y la media mide *cuánto contrae*.
El tibio alcanza a PolyHull casi siempre pero pierde la contracción **extra**
—la que venía de las pasadas repetidas sobre la caja ya apretada— en esas 10.

**Bloqueante para promoverlo**: `exacto frío` produce **6 `unknown_error`** en
las 153 instancias (el tibio, 0). Hay una excepción que se escapa y hay que
encontrarla antes de mover el default.

Resumen del estado: **el default sigue siendo el compuesto** —mejor contracción
media y sin fallos—; el refresco exacto con warm start está a un paso (grupo B
en 1, 2.5× menos pivotes) y le faltan dos cosas concretas: las 6 excepciones y
las 10 instancias donde pierde la contracción de las pasadas repetidas.


---

## 8. Una conclusión que hubo que retractar, y lo que quedó en su lugar

**Lo que afirmé y era falso.** Al agregar «escalas extremas» al test, el simplex
fallaba contra SoPlex en 217 de 400 casos, y escribí que no era robusto al
régimen real y que **el escalado era prerrequisito**. Las dos cosas son
incorrectas, y el error estaba en el test.

### 8.1 El régimen que construí no existe, y además rompe a SoPlex

El test generaba coeficientes con escalas dispares **dentro de una misma fila**
(hasta 1e36 junto al `−1` de la columna `b`). Medido sobre las instancias reales
con [spread.cpp](spread.cpp), la dispersión de los términos `a·z` **dentro de
una fila** es:

| instancia | mediana | máximo |
|---|---|---|
| Brown-15 | **8** | 8.01 |
| Brown-20 | **10.5** | 10.6 |
| DiscreteBoundary-0040 | 11.2 | 11.3 |
| Virasoro-icse | 1e8 | 1.5e8 |
| Katsura-12 | 7e5 | 7e10 |

O sea que `Brown-*` —justo donde el camino tibio falla— está **bien
condicionada**: dispersión intra-fila de 8. Lo que varía entre instancias es la
escala global, no la mezcla dentro de una fila.

Y el régimen inventado no solo es irreal: **rompe también a SoPlex**, que lanza
`SPxStatusException` y declara infactibles LPs que no lo son. Comparar contra
una referencia rota no mide nada.

### 8.2 Con el régimen calibrado, el simplex es robusto

Recalibrado a lo medido —escalas **por fila** hasta 1e6 y cajas hasta 1e6, con
dispersión intra-fila chica— el resultado es:

| | coinciden | discrepancias |
|---|---|---|
| frío, sin escalado | 400 / 400 | **0** |
| frío, con escalado | 400 / 400 | **0** |
| tibio (804 cotas), sin escalado | 804 / 804 | **0** |
| tibio, con escalado | 804 / 804 | **0** |

**El simplex es robusto en el régimen del banco, y el escalado no cambia nada.**

### 8.3 Por qué el escalado diagonal no podía ser la respuesta

Hay además un argumento algebraico que debí haber hecho antes de escribir el
código: el escalado de **columnas es un no-op** para los términos `a_ji·z_i`.
Con `z = C z_s` el coeficiente pasa a `r_j a_ji c_i` y el ancho a
`width_i / c_i`, así que el producto vale `r_j a_ji width_i`, **independiente de
`c_i`**. Solo la escala de fila mueve los términos, y en bloque, sin cambiar la
dispersión relativa dentro de la fila. Si el problema fuera la dispersión
intra-fila, ningún escalado diagonal lo arreglaría; y como no lo es, no hay nada
que arreglar.

El escalado queda implementado (equilibrado en norma infinito, escalas
redondeadas a **potencias de 2** para que multiplicar por ellas sea exacto, y
las columnas de `b` atadas a la fila con `C[nx+j] = 1/R[j]` para que la parte
`b` siga siendo exactamente `−I` y la base inicial siga siendo trivial), pero
**APAGADO por omisión** (`DFB_SX_SCALE=1` lo enciende).

Y acá la misma lección otra vez, que por eso conviene subrayarla: midió
**neutro** en el test sintético y **daña** en el banco real —`Virasoro-icse` cae
de **100 % a 66.67 %** de contracción—. Lo que no se mide sobre el banco real no
se sabe, por bien que se vea el test.

### 8.4 Lo que sigue abierto

La pérdida de contracción del camino tibio en **10 instancias** —y el caso
concreto de `Brown-15`, donde una cota superior devuelve **16**, la cota
*inferior* de la variable, con base consistente (residuo `d[básica]` ~1e-14),
factibilidad dual y factibilidad primal satisfechas— **sigue sin explicación**.
No es escala, no es deriva de la factorización, no es el umbral de reencolado
(medido: sin efecto), no es `z_k` no básica (medido: al revés) y no es el
criterio de optimalidad sobre datos derivados (corregido).

Lo que sí quedó: **el camino frío es sólido y cierra la brecha de poda** (grupo
B de 33 a 3), y **el tibio tiene la mejor clasificación de todas** (grupo B en 1)
con 2.5× menos pivotes, a cambio de esas 10.

### 8.5 La lección de método

Un test sintético vale lo que vale su régimen. Este pasó 400/400 y 1 066/1 066
sin cubrir el caso que importaba, y después «falló» 217/400 en un régimen que no
existe y que rompe a la propia referencia. Las dos veces el número parecía
concluyente. Lo que lo resolvió fue **medir el régimen real primero**
([spread.cpp](spread.cpp)) y calibrar el test contra esa medición.


---

## 9. En el árbol: la poda se cierra, el costo pasa a ser el cuello

La medición que contesta la pregunta del proyecto: la sustitución directa
(`--filtering=acidhc4 --lr=dfb`, DFB en el hueco exacto de `CtcPolytopeHull`)
con la regla de pivote de simplex. 191 instancias, 30 s.

Sobre las **113 que las cuatro resuelven**, relativo a producción:

| config | resuelve | celdas (geom) | cpu (geom) |
|---|---|---|---|
| `acidhc4+ph` (producción) | **137** | 1.000 | 1.000 |
| `acidhc4+dfb` (regla de impacto) | 123 | 1.321 | 1.098 |
| **`acidhc4+dfb` + regla simplex** | 120 | **1.068** | **1.227** |
| `hc4+dfb` + regla simplex | 119 | 1.261 | 1.374 |

0 óptimos incompatibles en las cuatro.

**La poda está esencialmente resuelta**: las celdas bajan de ×1.321 a **×1.068**,
o sea prácticamente el árbol de producción. Es la confirmación en el árbol de lo
que el banco de una caja mostraba (grupo B de 33 a 3).

**Y el costo pasa a ser el cuello**: el cpu sube de ×1.098 a **×1.227**, porque
el simplex cuesta ×21 pivotes al resolver cada cota desde cero. Resuelve 120
contra 137.

### 9.1 Dos correcciones sobre lo que se reportó antes

Esta medición hubo que hacerla **tres veces**, y las dos primeras estaban
contaminadas. Vale dejarlo escrito porque explica números que circularon:

1. **La primera corrida daba ×0.873 de cpu** —«12.7 % más rápido que
   producción»— y era **falso**: corrió con el escalado activado. Y el escalado
   no ayudaba, **enmascaraba**: al evitar las bases numéricamente difíciles
   reducía los crashes de SoPlex de 47 a 23.
2. **La segunda, sin escalado, daba peor** (108 resueltas) porque las
   excepciones de SoPlex se escapaban: `SPxException` **no deriva de
   `std::exception`**, así que atravesaba cualquier `catch (std::exception&)` y
   mataba el proceso. Atajadas dentro de `DFBBasis` y degradando a «la
   factorización no sirve», el banco de una caja pasa de 6 fallos a **0**.

### 9.2 Lo que esto ordena

El warm start deja de ser una optimización opcional y pasa a ser **la ruta
crítica**: baja los pivotes 2.5×, que es del orden de lo que falta para cerrar
el ×1.227. Y su bloqueante son las 10 instancias del §8.4.


---

## 10. El warm start: gana en la caja, pierde en el árbol

El defecto del camino tibio quedó **localizado**: sobre la **misma caja**, las
cotas inferiores (`up=0`) coinciden con el frío y las superiores (`up=1`)
devuelven la respuesta del minimizar —`Brown-15` da **16** en vez de 6.67e6—
con base **dual-factible y primal-factible verificadas las dos**. O sea que el
criterio de parada no está mal planteado: lo que falla es cómo se establece el
objetivo nuevo sobre una base heredada, que es lo único asimétrico entre las dos
direcciones.

**No reproduce en el test unitario**, y se probaron las tres dimensiones que
faltaban: caja fija, caja que se aprieta entre cotas (`SX_ENCOGE=1`) y cotas de
`b` anchas como en el uso real (`SX_BANCHAS=1`). Las tres pasan 528/528. El
reproductor necesita algo más del contexto —probablemente los reencolados o los
`2n` contractores compartiendo la base— y queda anotado.

### 10.1 La mitad sana del warm start ya paga… en el banco de una caja

`DFB_SX_WARM_MIN=1` usa tibio solo para las cotas inferiores. Sobre 147
instancias, mismo poliedro:

| | grupo A | grupo B | contracción | `t_dfb` | pivotes |
|---|---|---|---|---|---|
| regla de impacto | 64 | 33 | 29.55 % | 2.19 s | 29 241 |
| simplex frío | 95 | 3 | 31.58 % | 95.04 s | 622 449 |
| **simplex tibio-min** | **96** | **2** | **31.64 %** | **28.30 s** | **295 418** |

Mejor clasificación, mejor contracción, **la mitad de los pivotes y un tercio
del tiempo**. Una instancia apenas peor, cuatro mejores.

### 10.2 Pero en el árbol es al revés

191 instancias, 30 s, sobre las 115 que las cuatro resuelven:

| config | resuelve | celdas (geom) | cpu (geom) |
|---|---|---|---|
| `acidhc4+ph` (producción) | **137** | 1.000 | 1.000 |
| `acidhc4+dfb` (regla de impacto) | 123 | 1.319 | 1.061 |
| `acidhc4+dfb` + regla simplex, **frío** | 121 | **1.065** | **1.221** |
| `acidhc4+dfb` + regla simplex, **tibio-min** | 118 | 1.097 | 1.563 |

**El tibio-min, que en la caja cuesta un tercio, en el árbol cuesta un 28 %
más.** Y resuelve 118 contra 121.

La explicación probable, y es medible: la configuración del tibio incluye
`DFB_REFACTOR_ALWAYS=1` y el refresco exacto en cada iteración, y **los dos
cuestan por caja**. En el banco de una caja se pagan una vez y se amortizan; en
el árbol se pagan en cada uno de los miles de nodos. Además, entre nodos la caja
cambia por completo, así que la base heredada del nodo anterior vale menos que
dentro de un nodo.

**La lección, que es la misma del §5.1 del plan aplicada al costo**: el banco de
una caja y el árbol **premian configuraciones distintas**. Una mejora de ×3 en
el banco puede ser una pérdida de 28 % en el árbol, y no hay forma de saberlo
sin medir en el árbol.

### 10.3 Dónde queda el objetivo del proyecto

La poda está cerrada —celdas ×1.065, prácticamente el árbol de producción— y lo
único que separa a DFB de ser preferible es el **costo del simplex**: ×1.221 de
cpu, que viene de resolver cada cota desde cero. El warm start es la vía, pero
la versión que funciona en la caja no es la que sirve en el árbol. Lo que falta
es **barrer las combinaciones en el árbol** —refactorización periódica en vez de
siempre, refresco exacto solo al cambiar de objetivo, tibio dentro del nodo y
frío entre nodos— y elegir ahí, no en el banco de una caja.


---

## 11. Por qué el warm start no paga en el árbol: la linealización cambia en cada nodo

Las pruebas del §10 dejaban una duda razonable: el warm start se probó con **un
solo simplex compartido** por los `2n` contractores, o sea que cada resolución
heredaba la base de **otro objetivo** —el caso donde menos vale, y encima con
`d` dado vuelta entero—. Lo que correspondía era que **cada contractor guardara
su propia base** y arrancara tibio desde *su* llamada anterior: mismo objetivo,
caja apenas cambiada.

Implementado (`DFB_SX_PROPIO=1`), y en el banco de una caja funciona:

| instancia | base compartida | **base propia** | contracción |
|---|---|---|---|
| Brown-15 | 1 371 pivotes | **837** | 93.56 % en las dos |
| Virasoro-icse | 1 409 | **1 004** | 100.00 % en las dos |
| Katsura-12 | 156 | **123** | 54.17 % en las dos |
| DiscreteBoundary-0040 | 1 925 | 1 925 | 15.63 % (una sola visita por cota) |

**Misma contracción, hasta 39 % menos pivotes.** Es la incrementalidad que el
proyecto tenía como premisa, y funciona *dentro* de una caja, cuando el
contractor se reencola.

### 11.1 Pero en el árbol no cambia nada

| config | resuelve | celdas (geom) | cpu (geom) |
|---|---|---|---|
| `acidhc4+ph` (producción) | **137** | 1.000 | 1.000 |
| regla de impacto | 123 | 1.321 | **1.059** |
| simplex, base compartida | 121 | **1.064** | 1.190 |
| simplex, **base propia** | 121 | 1.091 | 1.233 |
| simplex, base propia + refresco exacto | 118 | 1.049 | 1.362 |

Y conservar la base **entre nodos** (`DFB_SX_NODOS=1`, que evita que `load`
resetee la base cuando la linealización mantiene las dimensiones) tampoco
cambia nada medible.

### 11.2 La razón es estructural

`CtcDFBPropag::linearize` se llama **en cada nodo** con la caja actual, y
`LinearizerXTaylor` expande en una esquina de esa caja: **los coeficientes de
`Ā` cambian en cada nodo**. Así que la base del nodo anterior es una conjetura
sobre un **poliedro distinto**, no sobre el mismo LP con las cotas cambiadas.

Un *branch and bound* clásico warm-startea bien porque la matriz de la
relajación es fija y entre nodos solo cambian las cotas de las variables. Acá
cambia la matriz entera, y por eso:

- el warm start **entre nodos** no puede pagar, y no es un defecto de
  implementación;
- el warm start **dentro de una caja** sí paga (−39 % de pivotes), pero en el
  árbol hay pocos reencolados por nodo, así que el efecto se diluye.

Eso también explica, retrospectivamente, por qué `CtcPolytopeHull` gasta 1.08
iteraciones por cota: **no es warm start entre nodos**, es la cadena de las `2n`
cotas sobre la misma base **dentro** del nodo, que es lo único reutilizable
cuando la matriz cambia en cada nodo.

### 11.3 Dónde deja esto al proyecto

La mejor configuración medida sigue siendo el **simplex con base compartida en
frío**: celdas ×1.064 —prácticamente el árbol de producción— a cpu ×1.190, con
121 instancias contra 137.

Y queda un dato que hay que mirar de frente: **la regla de impacto es mejor en
cpu total** (×1.059) pese a podar mucho peor (×1.321 celdas), porque es
muchísimo más barata. El simplex cierra la brecha de poda pero cuesta más de lo
que esa poda vale en el árbol.

**Lo que sigue con más chance es un híbrido**: la regla de impacto es barata y
falla solo donde se traba, y el estancamiento es **detectable** —es exactamente
cuando el pricing no encuentra fila con la caja todavía grande—. Correr impacto
por omisión e invocar el simplex solo ahí daría la poda del simplex pagándolo en
las pocas cotas que lo necesitan.


---

## 12. El híbrido: la poda del simplex a dos tercios de su costo

La idea sale de una asimetría medida: la regla de impacto se traba en un punto
no óptimo, **pero cuando no se traba llega a la misma cota que el simplex** —el
grupo A son 64 de 147 instancias—, así que pagar el simplex ahí es puro costo.
Y el estancamiento deja una firma barata: **contrajo poco o nada teniendo la
variable todavía ancha**.

`DFB_SX_HIBRIDO=1` corre primero la regla de impacto, mira cuánto contrajo esa
cota, y solo si fue menos que `DFB_SX_HIB_TOL` (1 % por omisión) paga el simplex
sobre la caja ya contraída.

En el banco de una caja:

| instancia | impacto | simplex | **híbrido** |
|---|---|---|---|
| Virasoro-icse | 72.22 % / 58 piv | 100.00 % / 1 409 | **100.00 % / 937** |
| Brown-15 | 93.56 % / 248 | 93.56 % / 1 371 | **93.56 % / 440** |
| Katsura-12 | 54.17 % / 67 | 54.17 % / 156 | **54.17 % / 131** |
| Discrete-Integralf2-12 | 6.37 % / 115 | 90.92 % / 4 036 | 63.32 % / 2 560 |

Y en el árbol, 191 instancias, 30 s, sobre las 115 que todas resuelven:

| config | resuelve | celdas (geom) | cpu (geom) |
|---|---|---|---|
| `acidhc4+ph` (producción) | **137** | 1.000 | 1.000 |
| regla de impacto | 123 | 1.296 | **1.020** |
| regla de simplex | 120 | **1.062** | 1.195 |
| **híbrido (1 %)** | **122** | 1.065 | **1.129** |
| híbrido (0.1 %) | 120 | 1.071 | 1.120 |

**El híbrido domina al simplex puro** en las tres columnas: misma poda, cpu de
×1.195 a ×1.129, y dos instancias más. 0 óptimos incompatibles.

### 12.1 El estado honesto del objetivo

Ninguna configuración de DFB supera a producción todavía:

- por **poda** el simplex y el híbrido ya están ahí (×1.06 contra ×1.00);
- por **cpu** la regla de impacto es la mejor (×1.020) y poda mal (×1.296);
- por **instancias resueltas** ninguna pasa de 123 contra 137.

O sea que DFB tiene ahora dos configuraciones buenas en dimensiones distintas y
el híbrido las combina parcialmente, pero el costo del simplex donde hace falta
sigue siendo demasiado. Lo que falta no es un mecanismo más: es **abaratar el
simplex por cota**, y las dos vías obvias —warm start entre nodos y factorización
compartida— están cerradas por el §11: la linealización cambia en cada nodo.

Lo que queda sin explorar, en orden de promesa:

1. **Un criterio de invocación mejor que «contrajo poco»**. Hoy el híbrido paga
   el simplex también en el grupo C, donde nadie contrae y no hay nada que
   ganar. Distinguir «se trabó» de «no hay nada» ahorraría esa fracción.
2. **Menos cotas**: invocar el simplex solo en la variable objetivo, que es la
   que poda en un optimizador, en vez de en las `2n`.
3. **Devex o steepest edge** en la regla de impacto, para que se trabe menos y
   el híbrido tenga que invocar el simplex menos veces.


### 12.2 Dos ideas para abaratar el híbrido, las dos cerradas

**(a) Un criterio de invocación que distinga «se trabó» de «no hay nada».** El
híbrido paga el simplex también en el grupo C —50 de 147 instancias donde nadie
contrae—, y a priori los dos casos se ven igual: la regla de impacto no contrajo.
A posteriori sí se distinguen, así que se probó **memorizar**: si el simplex ya
corrió para esa cota y no dio nada, saltearlo.

- **Dentro de una caja no sirve**: cada contractor corre el simplex una sola vez
  por caja, así que la memoria no llega a dispararse nunca (`memo=0` en las
  cuatro instancias de prueba).
- **Entre nodos es claramente peor**: saltear los 8 nodos siguientes lleva
  `alkyl` de **82 celdas y 1.44 s a 1 712 celdas y timeout**. En retrospectiva
  es razonable: «no dio nada en este nodo» no predice el siguiente, porque la
  bisección parte una variable al medio y lo que no tenía nada que sacar puede
  tenerlo de golpe.

Queda apagada (`DFB_SX_MEMO_N=k` la enciende con `k` nodos de salteo).

**(b) Devex o steepest edge en la regla de impacto.** No tiene margen, y por una
razón que conviene dejar escrita: el pricing **ya calcula la norma de fila
exacta** (`nrm = Σ a²`) **en la misma pasada** que el impacto, o sea gratis. Y
devex existe precisamente para *aproximar* esa norma cuando recalcularla es
caro. Además el plan ya registra que devex con pesos actualizables se midió y se
descartó: la versión barata capturó la ganancia de contracción y esa ganancia no
llega al árbol.

O sea que el margen no está en **cuándo** invocar el simplex ni en **cómo
pricea** la regla barata, sino en **cuánto cuesta el simplex por cota** — y las
vías para eso (warm start entre nodos, factorización compartida) están cerradas
por el §11.
