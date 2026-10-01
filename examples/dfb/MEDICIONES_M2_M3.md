# Mediciones: desglose por fase y evaluación de M3(a)

Continuación de [CORRECCION_M0.md](CORRECCION_M0.md), sobre el código ya sólido.
Datos crudos en [results/](results/); protocolo del §3 de
[PLAN_MEJORA_DFB.md](PLAN_MEJORA_DFB.md).

---

## 1. Nota metodológica: los resultados son deterministas

El plan advertía que la esquina `RANDOM` de `LinearizerXTaylor` introducía
varianza entre corridas y que habría que fijar la semilla para comparar reglas
de pivoteo. **Esa advertencia era incorrecta.** Medido sobre 5 repeticiones
completas del banco: **145 de 145 instancias dan exactamente el mismo
resultado**, rango máximo-mínimo de 0.000 puntos porcentuales.

La razón es que cada instancia corre en su propio proceso y el generador de
Ibex se inicializa igual en cada arranque, de modo que la secuencia de esquinas
"aleatorias" es la misma. La varianza aparecería sólo al cambiar la semilla
explícitamente.

Consecuencia práctica: **una corrida basta** para comparar variantes, y las
diferencias observadas son el efecto de la variante, no ruido. (Lo que sí falta
para robustez estadística es variar la semilla a propósito, que es un
experimento distinto: mediría la sensibilidad a la esquina, no el ruido.)

## 2. Dónde se va el tiempo de DFB

Suma sobre las 147 instancias que terminan, 18.35 s de `t_dfb` total
([results/fases_base.csv](results/fases_base.csv)):

| fase | tiempo | % |
|---|---|---|
| **pricing** (`largestImpact` + `calculateImpacts`) | 9.92 s | **54.0 %** |
| eliminación (`makeColumnIdentity` y `A[0] += α·A[j]`) | 3.34 s | 18.2 % |
| resto (init de los `2n` contractores, linealización) | 3.34 s | 18.2 % |
| regeneración de `A` (`regenerateA`) | 1.69 s | 9.2 % |
| cota (`gaussSeidel`) | 0.05 s | 0.3 % |
| **ratio test** (`calculateAlpha`) | 0.02 s | **0.1 %** |

En las instancias más caras el pricing sube al 62 % (`BroydenBanded-100/120/140`).

### Qué implica para el orden del plan

1. **M2 (pricing) es la prioridad, por encima de M1.** El 54 % del tiempo está
   en recalcular el impacto de las `m` filas en cada pivote. La lista de
   candidatos (*partial pricing*) ataca ese 54 % con un cambio local, mientras
   M1 (simplex revisado) es rediseño y apunta al 18 % de eliminación más parte
   del 18 % de «resto».
2. **El ratio test no se justifica por costo**: 0.1 %. M3 sólo se justifica por
   calidad de pivote — que es justo lo que mide el §3.
3. **La regeneración es un 9 % promedio pero muy concentrada**: 345
   regeneraciones en el banco, de las cuales **186 son de `Katsura-50` sola**,
   donde pasa a dominar el perfil. Ese caso aislado es el argumento de M4(b)
   (política de refactorización programada), no un problema general.

## 3. M3(a): desempate tipo Harris en el ratio test

Entre los candidatos cuyo cociente está dentro de una banda relativa del mínimo
(`harris_rel_band`, 1e-3), se elige el de **mayor |A[j][i]|**. El mínimo del
ratio test nunca se excede, así que la elección sigue siendo válida.

Resultado sobre las 153 instancias
([results/m3a_base.csv](results/m3a_base.csv) vs
[results/m3a_harris.csv](results/m3a_harris.csv)):

| | base | con M3(a) |
|---|---|---|
| Contracción media | 24.29 % | **24.73 %** |
| Instancias que mejoran / empeoran / iguales | — | **8 / 6 / 130** |
| Suma de deltas de contracción | — | **+63.5 pp** (ganan +68.1, pierden −4.7) |
| Pivotes por corrida | 20898 | **20237** (−3.2 %) |
| Tiempo total | 18.19 s | 18.42 s (+1.2 %) |
| Regeneraciones de `A` | 345 | 353 (+2.3 %) |
| Ratio tests indeterminados | 19 | **12** |
| Cajas vaciadas | 2 | 2 |
| Grupo A / B | 56 / 41 | **57 / 40** |

Los cambios son concentrados y grandes en pocas instancias:

| instancia | base | con M3(a) | delta |
|---|---|---|---|
| Brown-11 | 24.79 % | 49.59 % | **+24.79** |
| Brown-10-5 | 47.67 % | 67.50 % | **+19.83** |
| Brown-10 | 47.67 % | 61.83 % | **+14.17** |
| Katsura-12 | 37.50 % | 41.67 % | +4.17 |
| Katsura-15 | 46.25 % | 49.37 % | +3.12 |
| Brown-20 | 70.58 % | 72.15 % | +1.58 |
| Brown-15 | 70.10 % | 68.18 % | −1.92 |
| Katsura-40 | 47.64 % | 46.53 % | −1.11 |
| Katsura-25 | 41.67 % | 40.64 % | −1.03 |

**Veredicto: se acepta.** El neto es +63 puntos porcentuales de contracción y
−3 % de pivotes a cambio de +1 % de tiempo, sin perder solidez (mismas 2 cajas
vaciadas, ambas correctas). Queda **activado por defecto**
(`CtcDFB::harris_tie_break = true`); ponerlo en `false`, o `DFB_NO_HARRIS=1` en
el banco, reproduce la comparación.

Dos reservas honestas:

- **No hace lo que M3 predecía.** El desempate se propuso para mejorar la
  estabilidad numérica y por lo tanto reducir las regeneraciones de `A`; las
  regeneraciones **subieron** un 2.3 %. El beneficio real vino por otro lado:
  mejores cotas en las familias `Brown-*` y `Katsura-*` y menos ratio tests
  indeterminados (19 → 12). La hipótesis del plan era incorrecta aunque la
  mejora sea real.
- Hay 6 instancias donde empeora. Son pérdidas chicas (máximo −1.92 pp) frente a
  ganancias de hasta +24.79, pero indican que el desempate no domina: una regla
  de pricing/ratio test adaptativa podría quedarse con lo mejor de las dos.

## 4. El grupo C no cuenta para reemplazar a PolytopeHull

De las 48 instancias donde DFB no contrae nada, **en las 48 tampoco contrae
PolytopeHull** (y gasta hasta 257 iteraciones de simplex en el intento). No es
una debilidad de DFB frente a PH: la relajación lineal simplemente no corta en
esas cajas. La comparación relevante las excluye:

| | instancias |
|---|---|
| A — DFB ≥ PolytopeHull | 57 |
| B — DFB < PolytopeHull | 40 |
| **DFB gana o empata** | **59 %** de las productivas |
| C — ninguno contrae (excluido) | 47 |

**La brecha del grupo B es chica**: mediana de **3.94 puntos porcentuales**, con
25 de 42 instancias por debajo de 5 puntos. La cola es corta y con estructura
propia: `Brown-10-5` (33.3 pp, que M3(a) reduce a 13.5), `Virasoro-icse`
(27.8), `brown5b` (20.0), `Virasoro-prodmin` (19.3) y los `ExtendedWood`
(12.5).

Del grupo C sale además un objetivo separado: DFB gasta **1924 pivotes** ahí sin
conseguir nada (`Geneig-icse`: 330 pivotes, 0 %). 47 de las 48 pivotean en
falso. Eso es criterio de terminación (M7), no pricing.

---

## 5. M2: pricing

El pricing es el 54 % del tiempo de DFB (§2). Se probaron dos cosas distintas:
una optimización de implementación que **no cambia el resultado**, y una técnica
de simplex que **sí cambia las filas elegidas**.

### 5.1 Pricing fusionado — **aceptado, activado por defecto**

`calculateImpacts` construía dos `IntervalVector` por fila y hacía el producto
punto en una pasada aparte: `2m` asignaciones y dos recorridos por pivote.
Además reevaluaba `m` veces la elección de esquina, que **depende sólo de
`gamma`, no de la fila `j`** — sólo el tercer caso (cuando `gamma[i]` está en la
banda muerta) mira `A[j][i]`.

La versión fusionada (`CtcDFB::pricing_fused`, `row_impact` + `prepare_columns`)
clasifica las columnas una vez por pivote y acumula el producto en la misma
pasada, en el mismo orden de índices que `mulVV` (`y += v1[i]*v2[i]` desde
`y=0`), de modo que el resultado es **idéntico**:

| | original | fusionado |
|---|---|---|
| Tiempo total | 18.33 s | **15.55 s** (−15 %) |
| Fase de pricing | 9.64 s | **6.83 s** (−29 %) |
| Instancias con resultado distinto | — | **0 de 146** |

Cero riesgo: mismos pivotes, mismas cotas, misma solidez.

### 5.2 Lista de candidatos (partial pricing) — **medido, no activado**

Se evalúa sólo un subconjunto de filas (`partial_pricing_size`, por defecto
`⌈√m⌉`), refrescándolo cada `partial_pricing_refresh` pivotes o cuando ninguna
candidata mejora (con pasada completa antes de concluir que no hay fila que
mejore, para no confundir «la lista no tiene nada» con «no hay nada»).

| config | tiempo | pricing | pivotes | regen | contracción media |
|---|---|---|---|---|---|
| fusionado (referencia) | 15.55 s | 6.83 s | 20237 | 353 | **24.73 %** |
| candidatos, refresco auto (√m) | **11.10 s** | 4.68 s | 15171 | **52** | 20.57 % |
| candidatos, refresco cada 1 | 13.38 s | 5.54 s | 19873 | 226 | **24.80 %** |
| candidatos, refresco cada 2 | 12.39 s | 5.08 s | 17052 | 170 | 22.31 % |

Con refresco automático el ahorro es grande (−29 % de tiempo, −85 % de
regeneraciones) pero **la contracción se hunde en la familia `Brown-*`**:
Brown-30 pasa de 73.99 % a **0.00 %**, Brown-20 de 78.50 % a 5.00 %. Pérdida
acumulada 627 pp contra 27 pp de ganancia. Inaceptable para un contractor.

**El parámetro que gobierna la calidad no es el tamaño de la lista sino la
frecuencia de las pasadas completas.** Con listas más grandes Brown-30 apenas
recupera (27 % con 20 filas de 60, 29 % con 40), mientras que alternando lista y
pasada completa vuelve a 73.99 %.

La configuración `refresco = 1` (alterna una pasada completa y una por lista) es
la única que no degrada:

- tiempo −14 %, pricing −19 %, **regeneraciones −36 %** (353 → 226)
- contracción media apenas mejor (24.73 % → 24.80 %): mejora en 8 instancias
  (+21.7 pp) y empeora en 4 (−12.0 pp)
- mismas 2 cajas vaciadas, ambas correctas
- pero el grupo A baja de 57 a 56 (`Katsura-18` pasa a B) y `Brown-18` cae
  10 pp, aunque sigue muy por encima de PolytopeHull (66.2 % contra 36.7 %)

Gana donde más importaba el costo: `Katsura-50` 45.20 % → 47.16 % con 1.49 s en
vez de 2.34 s y 155 regeneraciones en vez de 249; `Katsura-25` +7.05 pp;
`Katsura-30` +3.23 pp.

**Decisión: queda implementado y apagado** (`CtcDFB::partial_pricing = false`).
El motivo es el objetivo: mientras la restricción que limita a DFB sea igualar
la contracción de PolytopeHull, cambiar contracción por tiempo va en contra, y
`refresco = 1` cuesta una instancia del grupo A. **La medición que debe decidirlo
es la del nivel 3 del §3.3 del plan**: dentro de `ibexopt --filtering=dfb`, donde
el contractor se llama en cada nodo y el tiempo sí es la restricción; ahí un
−14 % de tiempo por −0.3 pp de contracción probablemente sea buen negocio.
Activarlo: `partial_pricing = true; partial_pricing_refresh = 1`, o
`DFB_PARTIAL=1 DFB_PARTIAL_REFRESH=1` en el banco.

### 5.3 Refresco adaptativo — **probado y descartado**

Un refresco cada K pivotes es ciego. La idea natural era invalidar la lista en
cuanto un pivote deja de mejorar la cota («si el partial pricing se estanca,
ampliar»). Implementado en `partial_pricing_adaptive`, **no funciona**:

| instancia | fusionado | refresco 1 | adaptativo |
|---|---|---|---|
| Brown-30 | 73.99 % | 73.99 % | **0.00 %** |
| Brown-18 | 76.23 % | 66.22 % | **25.88 %** |
| Katsura-50 | 45.20 % | 47.16 % | 48.40 % |

La razón es informativa: en la familia `Brown-*` los pivotes elegidos por la
lista **sí mejoran** la cota en cada paso, así que la señal nunca se dispara, y
sin embargo llevan la búsqueda a un punto peor que el que alcanza el pricing
completo. Es decir, el problema no es el estancamiento local sino que la lista
toma decisiones localmente buenas y globalmente peores. Queda apagado.

## 6. Estado de las mejoras

| mejora | estado |
|---|---|
| M0 solidez | **hecho**, ver [CORRECCION_M0.md](CORRECCION_M0.md) |
| M3(a) desempate Harris | **aceptado y activado** (§3) |
| M2 pricing fusionado | **aceptado y activado** (§5.1) |
| M2 lista de candidatos | medido, **apagado**; lo decide la medición dentro de `ibexopt` (§5.2) |
| M2 refresco adaptativo | probado y **descartado** (§5.3) |
| M2 devex (normalizar por norma de fila) | sin probar |
| M4(a) escalado | sin probar |
| M4(b) refactorización programada | sin probar; el partial pricing ya baja las regeneraciones un 36–85 %, lo que resta urgencia |
| M1, M5, M6, M7 | sin implementar |

Acumulado hasta acá, sobre el mismo banco de 146 instancias comparables:
tiempo **18.33 s → 15.55 s** y contracción media **24.29 % → 24.73 %**, sin
perder solidez.

---

## 7. M5(a): una eliminación por variable en vez de dos — **aceptado**

El plan proponía «compartir una sola base entre los `2n` contractores». La
medición previa **descartó ese diseño**: los `2n` contractores se usan todos y
muchas veces (Katsura-50 hace 1953 aplicaciones sobre 102 contractores;
BroydenBanded-140, 1119 sobre 280), así que no hay nada que diferir, y comparten
poco porque divergen en cuanto pivotean.

Lo que sí se puede: **deducir el contractor de cota superior del de cota
inferior**. Los dos parten de la misma matriz de referencia y difieren solo en
que el superior niega antes la columna `k`. Desarrollando `makeColumnIdentity`
con esa columna negada:

- la **fila 0** queda igual que en el caso inferior con las entradas `i ≠ k`
  negadas (la entrada `k` se fuerza a 1 en ambos);
- las **filas `jj ≥ 1` quedan idénticas**, porque el factor de eliminación
  también cambia de signo y `(−factor)·(−fila 0) = factor·fila 0`.

En aritmética de intervalos la igualdad es **exacta**: negar es simétrico bajo
redondeo dirigido, luego `−(a/b)` y `(−a)/b` dan el mismo intervalo.

Resultado sobre las 147 instancias que terminan
([results/m5a_off.csv](results/m5a_off.csv) vs [results/m5a_on.csv](results/m5a_on.csv)):

| | sin pares | con pares |
|---|---|---|
| Instancias con resultado distinto | — | **0** |
| `t_init` | 2.904 s | **2.342 s** (−19.4 %) |
| — de la cual copia / eliminación | 1.847 / 1.024 | 0.902 / 0.510 |
| `t_dfb` total | 16.20 s | 15.62 s (−3.6 %) |

Activado por defecto (`CtcDFB::pair_init`); `DFB_NO_PAIR_INIT=1` vuelve al camino
anterior y sirve para comprobar la identidad de los resultados.

**Lo que este experimento dejó como dato para M1**: la inicialización pesa el
**19.4 % del tiempo en la mediana y hasta el 45 %**, y lo que queda tras el
arreglo es la **copia** de la matriz, no la eliminación. Copiar `2n` matrices de
`m × (n+m)` por caja no se puede optimizar más sin dejar de almacenarlas, que es
exactamente lo que hace M1 al pasar a una factorización de tamaño O(m + nnz).
Es un argumento nuevo y cuantificado a favor de M1.

Un detalle de implementación que costó una medición en falso: había **dos**
bucles de inicialización, uno en `init_dfb_contractors` y otro dentro de
`contract` (el que se usa en modo *stand_alone*). Parchear uno solo daba cero
mejora. Quedaron factorizados en `CtcDFBPropag::init_all_dfb`.
