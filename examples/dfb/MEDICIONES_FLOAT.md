# Pivoteo en punto flotante con certificación de la cota

Resuelve la decisión que el plan tenía abierta: **se pivotea en `double` y se
certifica solo la cota final en intervalos.** Datos en
[results/float_box.csv](results/float_box.csv) y
[results/float_ibexopt.csv](results/float_ibexopt.csv).

---

## 1. Por qué es sólido

La linealización da `refA · x = 0` para toda solución de la caja (con `x`
incluyendo las variables `b`). Entonces para **cualquier** vector real `λ`, la
combinación `γ = λᵀ·refA` cumple `γ·x = 0`, y de ahí sale la cota

```
x_k  ≥  cota( −(Σ_{i≠k} γ_i·x_i) / γ_k )   evaluada sobre la caja
```

La validez **no depende de que `λ` sea bueno**, solo de que `γ` se calcule
rigurosamente a partir de `refA`. Un `λ` malo da una cota débil, nunca una cota
falsa. Es dualidad débil: todo punto dual factible acota. Por eso el pivoteo
—que es puramente la búsqueda de `λ`— puede hacerse en flotantes sin ninguna
garantía. Es el mismo esquema de SoPlex en modo certificado y de
Neumaier–Shcherbina, en el que ya se apoya `CtcPolytopeHull`.

## 2. Cómo está implementado

El tableau flotante se mantiene **aumentado con la identidad**:

```
Af = [ A | I ]        m × (nb_var+m + m)
```

y la parte aumentada de la fila `j` **es el `λ` de esa fila**. Invariante:
`λ_jᵀ·refA_k = parte A de la fila j`, con `refA_k` = `refA` con la columna `k`
negada cuando se contrae la cota superior (así la certificación usa el `refA`
compartido y solo niega la entrada `k`, sin almacenar una referencia por
contractor).

**El punto fino**: el forzado de entradas a 0 exactos se aplica **solo a la
parte A**, nunca a `λ`. Así `λ` queda siempre veraz y la certificación absorbe
toda la suciedad numérica del pivoteo. Eso es lo que hace que pivotear sucio sea
seguro.

La certificación es una pasada al final de cada llamada: `γ = λᵀ·refA` en
intervalos (solo sobre los `λ_j ≠ 0`) y luego `gaussSeidel` sobre la caja. Con
una guarda: si `γ_k` contiene el 0 la cota no acota nada y se descarta — no es
un problema de solidez, es descartar un `λ` malo.

**Costo**: `nnz(λ)` es el número de pivotes acumulados, y medido contra `m` da
mediana **0.13** (p90 0.40). O sea que certificar una cota cuesta ~13 % de **un
solo** pivote de la versión de intervalos.

Lo que desaparece: `get_Aerror()` devuelve 0, `regenerateA()` no se dispara, y
con ellos **M4(b) queda sin objeto**.

Lo que se resigna: el vaciado por «ninguna fila bloquea» deja de ser prueba,
porque se apoyaría en datos flotantes. Cuesta poco: en todo el banco
`n_no_candidate` fue **0**, los 2 vaciados vienen de la cota y esa vía sigue
siendo válida.

## 3. Resultados: banco de una caja

| | intervalos | flotantes |
|---|---|---|
| **Tiempo total** | 15.62 s | **5.65 s** (−63.8 %) |
| **Regeneraciones de `A`** | 731 | **0** |
| Grupo A (DFB ≥ PolytopeHull) | 57 | **64** |
| Grupo B | 40 | **33** |
| Grupo C | 48 | 48 |
| Cajas vaciadas | 2 (correctas) | 2 (correctas) |
| Contracción media | 24.56 % | 24.03 % |
| Pivotes | 20878 | 25646 |
| Instancias que no terminan | 6 | 6 |

Los pivotes **aumentan** un 23 % y el tiempo cae dos tercios: cada pivote cuesta
~3.5× menos, y el método se puede permitir pivotear más. DFB gana o empata ahora
en el **66 %** de las instancias productivas (era 59 %).

La contracción media queda igual porque hay movimiento en las dos direcciones
(contrae más en 27 instancias, menos en 13), y los cambios grandes se concentran
en la familia `Brown-*`:

| instancia | intervalos | flotantes |
|---|---|---|
| Brown-07sp | 46.47 % | **92.73 %** |
| Brown-10sp | 46.02 % | **91.02 %** |
| brown5b | 80.00 % | **100.00 %** |
| Brown-30 / Brown-30-1 | 73.99 % | 43.99 % |
| Brown-20-1 | 78.50 % | 47.50 % |
| Brown-19 | 70.74 % | 40.73 % |

Son secuencias de pivotes distintas por el redondeo, no un sesgo del método. Si
una regla de pricing mejor estabiliza esa familia hay ~30 puntos porcentuales
por recuperar; es el argumento más concreto que queda a favor de M2(b) devex.

## 4. Resultados: dentro de `ibexopt`

| config | resueltas (int) | resueltas (float) |
|---|---|---|
| `acidhc4+ph` (producción) | 116 | 116 |
| `dfb` | 103 | 103 |
| `acid_dfb` | 86 | **96** |

Comparación directa del mismo filtrado, intervalos → flotantes:

| | celdas | tiempo | tiempo total |
|---|---|---|---|
| `dfb` (101 inst.) | ×0.914 | **×0.645** | 187.0 s → 101.2 s |
| `acid_dfb` (84 inst.) | ×0.961 | **×0.439** | 205.8 s → 58.9 s |

Contra la configuración de producción, sobre las 91 instancias comunes:

| config | celdas × ref | tiempo × ref | (antes, en intervalos) |
|---|---|---|---|
| `dfb` | 2.48 | **1.38** | ×1.80 |
| `acid_dfb` | 1.51 | **3.44** | ×6.47 |

**0 óptimos incompatibles** con la referencia en ambas configuraciones.

### La lista de candidatos deja de aportar

Con pivoteo flotante, `dfb/partial` da ×1.40 contra ×1.38 y
`acid_dfb/partial` ×3.45 contra ×3.44: **el efecto desaparece**. Tiene sentido —
la lista de candidatos servía para abaratar un pricing que costaba operaciones
de intervalo; en flotantes el pricing ya no es el cuello. Queda apagada, ahora
por una razón más clara que antes.

## 5. Estado del objetivo

| criterio | antes de flotantes | ahora |
|---|---|---|
| Solidez | cumplido | cumplido |
| Igualar contracción | 59 % de las productivas | **66 %** |
| Competitivo en el optimizador | `acid_dfb` ×6.47 | **×3.44** |

Ninguna de las dos configuraciones domina todavía: `dfb` está cerca en tiempo
(×1.38) y lejos en árbol (×2.48 celdas); `acid_dfb` está cerca en árbol (×1.51)
y lejos en tiempo (×3.44). La brecha de tiempo se redujo casi a la mitad con un
solo cambio.

Lo que sigue, por orden: **M2(b) devex** (la familia `Brown-*` dice que hay 30
puntos de contracción en juego), la **métrica incremental** (ahora que no hay
regeneraciones, conservar la base entre llamadas es posible de verdad) y
**M5(b) warm start entre nodos**, que con `λ` explícito ya tiene el objeto que
necesitaba.

---

## 6. M2(b): normalizar el pricing por la norma de la fila — **aceptado**

La otra mitad de M2, que estaba sin probar. El pricing era tipo Dantzig: elegía
la fila de mayor mejora *estimada*, **sin mirar cuánto se puede avanzar** en esa
dirección. Pero el paso `alpha` que permite el ratio test es aproximadamente
inversamente proporcional a la magnitud de la fila, así que `delta / ‖A_j‖`
aproxima la mejora **real** (`delta·alpha`) en vez de la mejora por unidad de
paso. Es el fundamento del steepest edge, en su forma más barata.

**En el camino flotante sale casi gratis**: el pricing ya recorre todas las
columnas de cada fila, así que la norma se acumula en la misma pasada
(`nrm += a*a` y una raíz por fila).

Un detalle que hay que cuidar: con normalización `delta` queda escalado y **no
sirve como criterio de parada** (el umbral de `1e-6` dejaría de significar «no
hay mejora»). La normalización cambia solo el **orden de elección**; la parada
se evalúa sobre la mejora sin normalizar de la fila elegida.

### Banco de una caja

| | Dantzig | normalizado |
|---|---|---|
| **Contracción media** | 24.03 % | **30.03 %** |
| **Grupo A** (DFB ≥ PolytopeHull) | 64 | **72** |
| **Grupo B** | 33 | **24** |
| Cajas vaciadas | 2 (correctas) | 2 (correctas) |
| Tiempo | 5.65 s | 6.02 s (+6.5 %) |
| Pivotes | 25646 | 33176 (+29 %) |

Mejora en 31 instancias y empeora en 14, con **+870 puntos porcentuales netos**
(ganan +901, pierden −31): una relación de 30 a 1.

| instancia | Dantzig | normalizado | pivotes |
|---|---|---|---|
| Brown-07 | 10.20 % | **73.47 %** | 59 → 83 |
| Brown-05 / brown5a | 2.00 % | **64.00 %** | 31 → 38 |
| Brown-11 | 24.79 % | **82.64 %** | 134 → 209 |
| Discrete-Integralf2-19 | 0.88 % | **58.19 %** | 182 → 465 |
| Brown-19 | 40.73 % | **91.18 %** | 301 → 536 |
| Brown-30 | 43.99 % | **93.44 %** | 714 → 1710 |
| Katsura-50 | 45.93 % | 50.29 % | **2155 → 716** |

Las pérdidas son chicas (máximo −6.19) y se concentran en la familia
`Discrete-Integralf2-*`, que parte de 96.9 % y por lo tanto tiene poco que
perder, aunque gasta muchos más pivotes para llegar. La única pérdida fuera de
esa familia es `EQCombustion` (−5.00).

DFB gana o empata ahora en el **75 %** de las instancias productivas (era 66 %).

### Dentro de `ibexopt`

Y acá está el resultado que conviene no maquillar: **la ganancia no se traslada**.

| | celdas | tiempo | resueltas |
|---|---|---|---|
| `dfb` normalizado vs Dantzig (98 inst.) | ×1.009 | ×0.988 | 103 → 100 |
| `acid_dfb` normalizado vs Dantzig (94 inst.) | ×0.982 | ×0.969 | 96 → 96 |

Contra la configuración de producción:

| config | celdas × ref | tiempo × ref | (con Dantzig) |
|---|---|---|---|
| `dfb` | 2.48 | **1.31** | ×1.38 |
| `acid_dfb` | 1.48 | **3.21** | ×3.44 |

0 óptimos incompatibles. Las 3 instancias que `dfb` deja de resolver **no son
fallas**: se verificó a mano que corren bien y solo son más lentas
(`ex14_2_4`: 5432 celdas y 4.74 s con Dantzig, 7970 y 6.63 s normalizado), y en
el barrido paralelo chocaron con el límite de pared.

**Conclusión que importa más que el veredicto**: seis puntos porcentuales de
contracción media sobre una caja se traducen en un 2 % de celdas y un 3 % de
tiempo en el árbol. Es decir que **contraer más una caja grande no es el factor
que limita dentro del optimizador** — la propagación y ACID recuperan buena
parte de lo que el contractor lineal deja sobre la mesa. Eso reordena las
expectativas de todo lo que queda: las mejoras de contracción pura rinden poco
en el árbol, y lo que rinde es tiempo por nodo y robustez.

**Aceptado y activado** (`CtcDFB::pricing_norm = 1`); `DFB_PRICING_NORM=0`
vuelve a la regla tipo Dantzig.

## 7. Verificación de la configuración final

Con pivoteo flotante y pricing normalizado, los dos activados por defecto:

| invariante | resultado |
|---|---|
| Los 7 casos de `tests_dfb_all` | pasan, con `NO SOLUTION 1` y `2` vacíos |
| **Campaña de testigos** ([results/witness_campaign_final.csv](results/witness_campaign_final.csv)) | **32 de 32 conservan el testigo**, 0 filas con cotas `b` que lo excluyan |
| Cajas vaciadas en el banco | 2, ambas corroboradas por `ibexsolve` |
| Óptimos en `ibexopt` | 0 incompatibles con la referencia, en `dfb` y en `acid_dfb` |
| `valgrind` | sin errores |

Nótese que la solidez se sostiene **aunque todo el pivoteo sea flotante**, que es
exactamente lo que la certificación garantiza: un `λ` malo da una cota débil,
nunca una falsa.

---

## 8. El perfil por fase en flotantes es otro

Volver a medirlo era necesario: el reparto conocido se había medido en el camino
de intervalos y ya no aplica. Sobre las 147 instancias que terminan
([results/fases_float.csv](results/fases_float.csv)), `t_dfb` = 6.04 s:

| fase | flotante | (intervalos) |
|---|---|---|
| **inicialización de los `2n` contractores** | **49.4 %** | 18 % |
| pricing | 37.3 % | 54 % |
| **certificación (`λᵀ·refA`)** | **4.8 %** | — |
| eliminación | 3.4 % | 18 % |
| cota (`gaussSeidel`, en intervalos) | 1.4 % | 0.3 % |
| ratio test | 0.2 % | 0.1 % |
| regeneración de `A` | **0.0 %** | 9 % |
| resto (linealización, propagación) | 3.5 % | — |

Dos lecturas:

1. **La certificación cuesta 4.8 %.** La preocupación por el costo de certificar
   queda zanjada: es barata, tal como anticipaba el `nnz(λ)/m` de 0.13.
2. **La inicialización pasó a ser la mitad del tiempo**, no porque empeorara
   sino porque todo lo demás se abarató (la eliminación cayó de 18 % a 3.4 %).
   Eso convierte la factorización de la base en la prioridad clara: es lo único
   que elimina el armado de `2n` tableaus de `m × (nb_var+2m)`.

## 9. El truco del par lb/ub, ahora en el tableau flotante

Estaba implementado solo para el camino de intervalos. La identidad **se extiende
al tableau aumentado**: con la columna `k` negada, el pivote de la fila 0 pasa de
`R[0][k]` a `−R[0][k]`, así que la fila 0 completa —parte `A` **y parte `λ`**—
queda negada, salvo la entrada `k` de la parte `A` que se fuerza a 1 en ambos; y
las filas `jj ≥ 1` quedan idénticas, porque el factor de eliminación también
cambia de signo:

```
e_jj − (−R[jj][k])·(−(e_0/R[0][k]))  =  e_jj − R[jj][k]·e_0/R[0][k]
```

Resultado sobre las 147 instancias
([results/pairf_off.csv](results/pairf_off.csv) vs [results/pairf_on.csv](results/pairf_on.csv)):

| | sin par | con par |
|---|---|---|
| `t_init` | 2.983 s | **2.069 s** (−30.6 %) |
| `t_dfb` | 5.99 s | **5.07 s** (−15.4 %) |
| init como % de `t_dfb` | 49.8 % | 40.8 % |

**Sobre la identidad de los resultados**: 22 instancias difieren, pero **ninguna
en contracción por más de 0.01 pp** (suma de deltas: −0.0001 pp). Lo que cambia
es el conteo de pivotes en ±1 o ±2, casi todo en instancias del grupo C que
contraen 0 % de todos modos. La causa es que negar la fila y reconstruirla
producen atajos distintos en la eliminación (`if (f == 0.0) continue`) para
valores en el borde de las tolerancias. Así que la afirmación correcta es
**equivalente en contracción**, no bit a bit idéntico como en el camino de
intervalos.

Activado por defecto (comparte el interruptor `DFB_NO_PAIR_INIT=1`).
