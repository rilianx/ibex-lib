# DFB en Ibex — estado del trabajo

Documento de traspaso: lo que quedó en pie, con las cifras que lo sostienen.
El registro cronológico está en `MEDICIONES_DUELO.md` y, para la sesión del
2026-09-29/30, en `MEDICIONES_SESION_2026-09-30.md`. Versiones largas de este
documento en `respaldo/`.

---

## 1. Qué es

`CtcDFBHull` reemplaza a `CtcPolytopeHull` en el hueco exacto del contractor
lineal: misma interfaz, misma composición, `--filtering=acidhc4 --lr=dfbhull`.

La relajación del nodo se escribe `Ā z = 0` con `z = (x, b)`, `na = nb_var + m`,
y la parte `b` de `Ā` es exactamente `−I`. De ahí:

- **Arranque dual-factible gratis**: con las `m` columnas de `b` como base,
  `B = −I` y `c_B = 0`. Sin fase 1.
- **Certificado**: de los multiplicadores `λ` se arma `γ = λᵀĀ` **en
  intervalos**; como `Āz = 0` para toda solución, `0 ∈ γ·z` y
  `z_k ∈ (−Σ_{i≠k} γ_i z_i) / γ_k`.

Tres consecuencias que definen el método:

1. **Se pivotea en punto flotante.** Un `λ` malo da una cota débil, nunca falsa.
2. **`γ` no depende de la caja**: vale para cualquier sub-caja de aquella en que
   se linealizó.
3. **Prueba de vacío por Farkas**: si `0 ∉ γ·z`, la caja no contiene solución de
   la relajación. Es el test de Neumaier–Shcherbina verificado en intervalos.

### Tres ejes que se llaman «dual»

| eje | qué decide | hoy |
|---|---|---|
| **el acotamiento** | de dónde sale la cota | dual siempre, por dualidad débil |
| **la selección** (`DFBH_ORDEN`) | qué cota se resuelve después | por costo dual, §3 |
| **el pivoteo** (`DFB_SX_PIVOTEO`) | cómo se re-optimiza al cambiar de objetivo | **`dual` (puro, omisión)**, `primal` o `mixta` |

### Las estrategias de pivoteo

Por nodo se linealiza, `load` reinicia la base, y se recorre una pasada de hasta
`2n` cotas con **un solo** simplex compartido. La cota 0 arranca en frío con
`B = −I`, dual-factible y primal-infactible. Las estrategias difieren en las
cotas `1 … 2n−1`, donde cambia el objetivo y la caja no se mueve:

- **`por_motor`** (omisión 2026-10-08…10; desde 2026-10-10 la omisión volvió a `dual`, §4): `dual` en los nodos que van al
  motor revisado (que no tiene pasos primales) y `primal` en los de tableau,
  donde el primal es más rápido (ablación, §9). El certificado se calcula igual
  en los dos; se pierde sólo la cota anytime monótona en los nodos de tableau.
- **`dual`** (omisión 2026-09-30…10-08 y de nuevo desde 2026-10-10): reubica **todas** las no básicas por
  el signo de su costo reducido, `k` incluida. La base es dual-factible en cada
  iteración y la cota certificada sube monótonamente durante el LP. Si `k`
  termina no básica, el óptimo es su cota y no había nada que contraer. Medido:
  0.1–1 % de los LPs dan algún paso primal (ver §9).
- **`primal`**: no reubica nada; conserva la factibilidad primal heredada y
  re-optimiza con pasos primales, el arranque tibio de SoPlex. El 78–92 % de los
  LPs dan pasos primales.
- **`mixta`**: la que hasta el 2026-09-30 se llamaba `dual`. Excluye `k` de la
  reubicación, y en los LPs que empiezan con `k` no básica (38–83 %) termina
  pivoteando en primal. Se conserva para reproducir mediciones.

**La poda es la misma** en las tres: el vacío se decide en la cota 0, común a
todas. La diferencia es costo, y con la tolerancia a escala de la caja, el paso
largo y devex (§2) el dual queda a la par del primal: cpu ×1.014 / ×1.002 /
×0.959 (§9).

Cuando lo que cambia es la **caja** y no el objetivo, el que hereda es el dual.
Es el escenario de la propagación propia archivada (§8).

---

## 2. Las optimizaciones que quedaron

Activas por omisión, sin variables de entorno (en DFB: `--lr=dfbhull`; la
relajación por omisión de `ibexopt` sigue siendo `xn`, producción). Además, en
`ibexopt` y la biblioteca: `eps_x` = 1e-7 real y `LargestFirst` con precisiones
nulas (§8).

| | efecto medido |
|---|---|
| **Prueba de vacío por Farkas** | 87 → 128 instancias resueltas |
| **Filtro de estabilidad en el test de razón primal** | celdas ×1.139 → ×0.935 contra producción (> 1 s) |
| **Orden de cotas por infactibilidades duales** | pivotes/LP ×0.78; cpu ×0.90 / ×0.92 / ×0.96 |
| **Pivoteo dual puro** (§9) | reubica todas las no básicas, `k` incluida: base dual-factible en cada iteración, cota certificada monótona; < 1 % de LPs con algún paso primal. Con paso largo, τ y devex queda a ~×1.03 del primal (sin `bearing`) |
| **Ratio test de paso largo**, con el arreglo del último quiebre (§9) | cpu ×0.97 contra el dual sin él |
| **Tolerancia primal a escala de la caja**, τ = 1e-3 (§9) | neutra para el primal; sin ella el dual puro declaraba factibles poliedros vacíos en cajas angostas |
| **Devex en la fila que sale** (§9) | cpu ×0.98–0.99 con el dual puro |
| **Estrategia por nodo y umbral 80** (`por_motor`: dual en los nodos revisados, primal en los de tableau) (§4, §9) | solver cpu ×0.97 por la estrategia y ×0.97–0.98 por el umbral, celdas iguales; `ibexopt` neutro |
| **Motor híbrido**: simplex revisado sin tableau en relajaciones de más de 30 filas (§9) | `ibexsolve` cpu ×0.75–0.93 con las mismas celdas; neutro en `ibexopt` |
| **LU sin el reloj de SoPlex** (`SLUFactor` leía `times()` en cada FTRAN/BTRAN) (§9) | `ibexsolve` cpu ×0.92–0.93 con el afín y ×0.83 con `compo`, mismas celdas; `launch` ×0.85 |
| **Sondas leídas una vez y `mid()` directo** (más de 20 `getenv` por LP en `una_pasada`; el `mid()` de gaol guarda y restaura el redondeo aun en intervalos de un punto, en `Af = mid(A)` cada pasada; también los `getenv` de `CtcPolytopeHull`, `LSmear` y el optimizador) | mismas celdas, DFB ~×0.97 en `ex6_1_3` y `ex7_2_3` (corridas sueltas, sin banco) |
| **Un solo simplex secuencial** | 6.8–12.3 % más rápido que `2n` ejemplares sembrados |
| **`x_B` sólo si está sucio en `primal_step`** | cpu ×0.99, resultados **idénticos bit a bit** |
| **LU diferida** (§9) | cpu por celda ×0.94; celdas dentro del ruido de redondeo |
| **`λ`, `d`, `x_B` del tableau con guardia** (§9) | cpu ×0.89 contra la diferida; la LU sólo cuando la guardia detecta deriva |
| **Dos costos fijos del linearizador afín** | cpu ×0.942 (> 1 s), 58 más rápidas y ninguna más lenta, **idénticos bit a bit** |

**El filtro de estabilidad — leer antes de tocar `primal_step`.** El test de
razón primal aceptaba cualquier coeficiente no nulo, `pivot` lo rechazaba por
`|p| < PIV_TOL`, y el rechazo llegaba al bucle indistinguible de «soy óptimo».
El 48 % de los `OPTIMAL` eran con cero pivotes, y el 83 % de ésos
dual-infactibles. El paso dual ya tenía el filtro de Harris; el primal no.

**Los dos costos del linearizador afín** los encontró un perfilado: `pow(2,-50)`
dentro del bucle de variables, `nb_ctr × nb_var` llamadas por linealización, y
una asignación al heap por nodo del árbol de expresión, 400–2200 por
linealización, hoy recicladas en un pool por tamaño. Cuidado: es código
vendorizado en el camino caliente, y las asignaciones de `ibex_AffineVar.cpp`
tienen que ir al mismo asignador que las de `ibex_Affine2_fAF2.cpp`, porque las
libera el mismo destructor.

**`compute_basics` se llamaba dos veces por iteración**: `dual_step` la
recalcula al entrar, y si devuelve falso no tocó nada, pero `primal_step` la
rehacía. Con la guarda `if (xb_sucio)` las llamadas bajan 25–34 %. No es la
actualización incremental, que sigue apagada por deriva numérica.

**Arquitectura de estados.** Hubo `2n` ejemplares sembrados desde la base óptima
de la cota 0: ×0.75 pivotes, pero sembrar cuesta lo que ahorra, 6.8–12.3 % más
lento. Hoy hay un solo `DFBSimplex` secuencial.

---

## 3. Selección de cota

Regla por omisión: **salteo de Achterberg** (si `x*_j` ya está en el borde, esa
cota está resuelta) más **orden por infactibilidades duales**
(`DFBH_ORDEN=costo`): para `x_j` básica en la fila `r`, se cuentan las no
básicas cuyo costo reducido `−ck·ā_rc` tiene el signo violado. Es el número de
columnas que el primal tendría que hacer entrar, `O(nnz(fila r))` sin pivotear.
Se elige el conteo menor; desempate por la más lejana.

Contra Achterberg + «lejana», 137 instancias, 4 semillas:

| corte | celdas | cpu | +rápido / −lento | piv/LP | cpu/LP |
|---|---|---|---|---|---|
| > 0.1 s | ×0.974 | **×0.903** | **128 / 21** | **×0.777** | ×0.922 |
| > 1 s | ×0.971 | **×0.924** | 41 / 9 | | |
| > 5 s | ×0.979 | ×0.960 | 17 / 7 | | |

Pivotes por LP, mediana por instancia: 3.35 → 2.39. La distancia al borde es un
sustituto del costo; el conteo de signos violados **es** el costo.

| variante | cpu (> 0.1 s / > 1 s / > 5 s) |
|---|---|
| **orden por costo** (omisión) | **×0.903 / ×0.924 / ×0.960** |
| hull + orden por costo | ×0.918 / ×0.943 / ×0.965 |
| hull solo (`DFBH_HULL=1`) | ×0.943 / ×0.956 / ×0.964 |
| orden más cercana / natural | ×1.015 / ×1.038 |
| sin salteo (`DFBH_NOSKIP=1`) | ×1.133 |

**El hull de puntos primales** (`DFBH_HULL=1`): la caja mínima `H` de los `x*`
vistos da una cota superior exacta de la ganancia de cada cota, `H.lb_j −
X.lb_j`. Se saltea bajo `τ` relativo (`DFBH_HULLTAU`, omisión `1e-6`; el
resultado no depende de `τ`). No suma con el orden por costo: el hull compra
celdas pagando LPs por nodo (×1.03–1.05).

---

## 4. Resultados contra producción

Producción = `--filtering=acidhc4` + `CtcPolytopeHull` + linearizador afín.

**Comparación controlada: el mismo contractor con dos motores (2026-10-09).**
Para el paper, la comparación justa no es contra `CtcPolytopeHull` de Ibex sino
el MISMO contractor (misma pasada, mismo salto, misma política ante un LP sin
certificado: seguir con la cota siguiente) con dos motores de LP:
DFB, o SoPlex como caja negra (`DFBH_MOTOR=caja`), cuyos multiplicadores (dual
o rayo de Farkas) pasan por el MISMO certificado gamma en intervalos
(Neumaier–Shcherbina en la forma homogénea). La tolerancia de SoPlex es un
parámetro de la caja negra (`DFBH_CAJA_TOL`): con 1e-9 (la de Ibex) SoPlex
declara infactibles LPs factibles en los `Brown` (casi-igualdades de ancho
1e-14 en cajas de ±1e8; verificado resolviendo el LP aparte), y en `ibexopt`
1e-6 es demasiado laxa. Mejor tolerancia: solver 1e-6 (≈ 1e-9), `ibexopt` 1e-9.
`results/controlado_caja_negra_2026-10-09.csv` (264 sistemas × 2 semillas,
139 problemas × 4; en `ibexopt` sin `alkyl`, `alkylbis`, `bearing`):

| comparación | > 0.1 s | > 1 s | > 5 s | terminan |
|---|---|---|---|---|
| solver, DFB (orden propio) contra caja negra 1e-6 | celdas ×0.89, cpu ×1.03 | ×0.86, ×1.08 | ×0.81, ×0.98 | 354 contra 356 |
| `ibexopt`, DFB (orden propio) contra caja negra 1e-9 | celdas ×1.00, **cpu ×0.77** | ×1.05, **×0.83** | ×1.04, **×0.78** | 129 = 129 (4/4) |
| solver, mismo orden (lejos): caja negra 1e-6 contra DFB | celdas ×1.13, cpu ×0.85 | ×1.17, ×0.81 | ×1.23, ×0.88 | 356 contra 354 |
| `ibexopt`, mismo orden: caja negra 1e-9 contra DFB | celdas ×0.97, cpu ×1.19 | ×0.94, ×1.15 | ×0.96, ×1.14 | |
| DFB, orden propio (por costo, usa el tableau) contra lejos | solver ×1.00, cpu ×0.87 | ×1.00, ×0.88 | ×1.00, ×0.91 | igual |

Lectura: en el optimizador (relajaciones chicas) DFB es 17–23 % más rápido con
el mismo árbol. En el solver, cpu a la par, DFB poda 11–19 % más y es más
robusto en los `Brown` (la caja negra deja de terminar `Brown-30`; ×15–20 en
`Brown-18…20`) y en `CountercurrentReactors2-9` (×25); SoPlex gana en las
relajaciones grandes (`Bratu-0175…0200`, que DFB no termina;
`Discrete-Integralf2-19…21`, ×3–5), donde su LU es mejor que el motor revisado.
Por qué DFB poda más en el solver con el mismo certificado (sonda
`DFBH_ABCAJA`: los dos motores sobre la misma caja y la misma relajación):
(1) **fallas numéricas de SoPlex** en relajaciones con casi-igualdades
angostas: en Geneig-icse declara infactibles 18 838 LPs factibles (el
certificado prueba vacío en 148; las demás cotas se pierden), más ~6500
desconocidos y ~1000 topes de iteraciones; DFB deja sin óptimo 211 de 75 446.
Contracción media 0.303 contra 0.248. (2) **precisión del dual** con
tolerancia 1e-6: diferencias chicas y sistemáticas (CountercurrentReactors2-9:
DFB gana 253 nodos contra 0; con 1e-9, 47 contra 1). Se probó darle a SoPlex
su mejor configuración: presolve, escalado y representación no cambian; el
primal baja los infactibles espurios (18 838 → 2549) pero en el banco
empeora (solver: termina 345–346 contra 356, pierde `Geneig-icse` y
`Brown-15…20`; `ibexopt` igual). Mejor caja negra: la de fábrica con
tolerancia 1e-6 (solver) o 1e-9 (`ibexopt`). La diferencia
de poda contra `CtcPolytopeHull` de Ibex venía sobre todo del marco: abandonar
la pasada y certificar a posteriori con un rayo que puede faltar. En
DiscreteBoundary-0200 DFB da 3 cajas solución y la caja negra 4: se solapan
en cadena y cada una interseca alguna del otro brazo, no es pérdida de
soluciones.

**Motor revisado más barato (2026-10-09).** `perf` en Bratu-0150: el
recálculo de `x_B` en cada pivote (producto por `Abar` y FTRAN densa) era el
16 % del total, y la columna entrante se resolvía dos veces (devex y
Forrest–Tomlin). Ahora `preparar_cambio`/`aplicar_cambio` reutilizan el eta de
`solveRight4update`, y `x_B` se actualiza de forma incremental en los nodos
revisados con recálculo cada 32 pivotes (`DFB_SX_XBINC_REV=0`,
`DFB_SX_XBREFRESCO`). Solver (`results/revisado_ftran_xbinc_2026-10-09.csv`):
cpu ×0.97 / ×0.96 / ×0.95 contra sin incremental (Bratu-0210 ×0.40, 0120
×0.65, Discrete-Integralf2-21 ×0.83), celdas ×1.00; el incremental redondea
distinto y cambia algunos árboles (BroydenBanded-200: 55 → 89 celdas, al tope;
no es deriva acumulada, refrescar cada 8 pivotes da lo mismo). `ibexopt`
neutro. Contra la caja negra con el mismo binario: solver celdas ×0.81–0.89,
cpu ×1.01 / ×1.04 / ×0.93 (antes ×1.03 / ×1.08 / ×0.98); `ibexopt` ×0.77 /
×0.82 / ×0.77. Lo que queda en el solver son las relajaciones más grandes
(Bratu-0175…0200, Discrete-Integralf2-18…19 ×3): DFB hace ~2 veces más
pivotes por LP que SoPlex, que reoptimiza en primal; el revisado no tiene pasos
primales (opción B).

**Pasos primales en el motor revisado (2026-10-10, `DFB_SX_REVPRIMAL=1`).**
`primal_step_revisado`: Dantzig sobre `d`, FTRAN preparada de la entrante
(reutilizada por el pivote inmediato), test de razón sobre `x_B`, BTRAN de la
fila y `pivot_revisado`; con `por_motor` los nodos revisados se resuelven en
primal (el dual queda para la primera cota, en frío). Dos defectos corregidos
en el camino: la LU tiene que asegurarse antes de pedir resoluciones a SoPlex
(con `fact_sucia` se corrompía el heap), y una columna preparada que quedó sin
usar (la entrante sólo cambió de cota) no puede reutilizarse en un pivote
posterior (segfault en `SLUFactor::change`, Yamamura). Pivotes por LP:
BroydenTri-0090 33 → 6, Discrete-Integralf2-21 78 → 39, Bratu-0150 51 → 35.
Banco (`results/revisado_primal_2026-10-10.csv`), contra la omisión dual:

| umbral | solver cpu > 0.1 / > 1 / > 5 s | terminan | `ibexopt` cpu |
|---|---|---|---|
| 80 | ×0.98 / ×0.98 / ×0.99 | 354 = 354 | ×1.00 |
| **40** | **×0.95 / ×0.94 / ×0.97** | 354 = 354 | ×0.99–1.00 |
| 20 | ×0.95 / ×0.93 / ×0.96 | 353 contra 354 | ×0.98–0.99 |

Celdas iguales en todos. **Por omisión desde 2026-10-10 con umbral 40.** Gana BroydenTri-0090…0100 (×0.45–0.64),
Discrete-Integralf2-17…19 (×0.57–0.62); pierde Bratu-0120 (×1.25). Contra la
caja negra, umbral 40: solver celdas ×0.81–0.89, cpu ×0.96 / ×0.97 / ×0.90;
`ibexopt` cpu ×0.77 / ×0.82 / ×0.78. Siguen sin terminar Bratu-0175…0200
(costo por pivote en relajaciones de ese tamaño).

**Umbral 30 desde 2026-10-13.** Las dos pérdidas de SC en `all` a 3600 s
(ex5_4_4 ×3.72, hhfair ×1.80; 28 variables, relajaciones de 39 y ~40–41 filas)
eran el motor: la sonda `DFBH_ABCAJA` da la misma contracción que SoPlex, pero
cada LP en tableau cuesta ×4.4 / ×2.6; el revisado hace los mismos pivotes más
baratos. Barrido en el banco de la ablación (132 sistemas + 139 opt, 4 semillas,
60 s): 35 y 30 empatan con 40 (±0.5 %, terminan 362 = 362, 519 = 519); ganan
Discrete-Integralf2-18/20 (×0.72–0.78) y ex5_3_2 (×0.93). En las dos instancias
(120 s, 4 semillas, orden lejos):

| | umbral 40 | 35 | 30 | caja |
|---|---|---|---|---|
| ex5_4_4, celdas en 120 s | 4760 | 9660 | 9560 | 14400 |
| hhfair, tiempo total | no termina (10100 celdas) | 107 s | 109 s | 106 s |

hhfair empata con la caja (20 % menos celdas); en ex5_4_4 queda ×1.5 por el
costo por LP del revisado (≈20 pivotes por LP contra 12 de SoPlex). Protocolo:
elegido en el banco del paper, se valida sin cambios en los 162 de `all`.
**Validación (2026-10-13, `results/all60_umbral30_2026-10-13.csv`)**: los 162
de `all` menos los 58 excluidos (104 opt, 4 semillas, 60 s), SC umbral 30
contra caja: cpu ×0.95 / ×0.95 / ×1.07 (> 0.1 / 1 / 5 s), resueltas 309
contra 305; contra umbral 40: ×0.98, 309 = 309, gana hhfair1 ×0.47, ninguna
pierde más de ×1.01. Caja repetida contra all60: ×1.00 (piso de ruido).
Siguen perdiendo con ambos umbrales hs099 ×2.4, linear ×1.8 (más celdas que
la caja) y ex8_4_4 ×1.5–1.8 (+8 % celdas, costo por LP).

**NA a 3600 s con 12 GB (2026-10-13, corridas de 600 s con `/usr/bin/time`)**:
loadbal crece ~12 MB/s (7.1 GB a 600 s, ~42 GB a 3600) y orthrege con SC
~5 MB/s (3.2 GB; la caja 1.8 GB porque hace la mitad de celdas): excluidas por
memoria como avion2. hs089 y hs091 usan < 300 MB pero `ibexopt` se pasa del
`-t` un 4–8 % (623 y 649 s con `-t600`) y a 3600 s el `timeout 3750` los mata;
no encuentran punto factible (`inf`), quedan como tope. matrix2 aborta en los
dos brazos (semillas 1–3, a los ~120 s) por el `assert(is_actif())` de
`AffineMain<AF_fAF2>::val` (`affine/ibex_Affine2_fAF2.cpp:225`): falla de la
linealización afín de Ibex, no de SC; con la semilla 4 termina en 0.1 s.

**Pérdidas que quedan en `all` (sonda `DFBH_ABCAJA`, 60 s, semilla 1)**: la
contracción en la misma caja es igual o mejor con DFB (hs099 0.341 contra
0.336, gana DFB 3934 nodos contra 70; linear 0.333 contra 0.326; ex8_4_4 0.330
contra 0.327); las celdas de más en la corrida completa son de trayectoria. Es
costo por LP: DFB hace ~2× pivotes por LP (hs099 20 contra 12, linear 31
contra 15, ex8_4_4 10 contra 5.5). hs099 (29 filas) y ex8_4_4 (25) van al
tableau con umbral 30, pero el revisado es más barato también ahí: µs/LP 277 →
150 y 120 → 78 (SoPlex 98 y 54). El umbral sólo por filas no captura el
cruce.

**Por qué DFB hace ~2× pivotes por LP (2026-10-13, `DFBH_PORPASO` por motor,
`DFBH_CRONO`, 60 s, semillas 1–2).** No es la cota 0 (en frío: hs099 39
contra 35 de SoPlex) sino las cotas tibias 1…2n−1 (hs099 20–25 contra 6–16;
ex8_4_4 8–9 contra 4–5; linear 27 contra 13). Al cambiar el objetivo, el dual
puro reubica todas las no básicas por el signo del costo reducido y tira la
factibilidad primal del vértice anterior; SoPlex conserva ese vértice y
re-optimiza en primal. Con `DFB_SX_PIVOTEO=primal` DFB baja a los pivotes de
SoPlex. Ojo: sin `DFB_SX_REVPRIMAL=1` el primal manda todos los nodos al
tableau (linear, 41 filas, no ganaba nada por eso). Corridas completas (cpu s,
media de 2 semillas):

| | dual (omisión) | primal + revisado > 30 filas | primal, todo revisado | caja |
|---|---|---|---|---|
| hs099 | 43.7 (19.8 piv/LP, 260 µs/LP) | 30.2 | **18.0** (10.9, 66) | 20.8 (11.9, 96) |
| linear | 11.5 (30.5, 245) | **5.0** (10.2, 76) | 5.0 | 6.7 (9.6, 85) |
| ex8_4_4 | 31.2 (10.0, 109) | 24.1 | 22.3 (5.8, 32; celdas ×1.3) | 18.5 (5.4, 48) |
| hhfair1 | **6.4** (946 celdas) | 9.2 (1680 celdas) | 9.2 | 10.4 |

Con primal y motor revisado el µs/LP queda por debajo de SoPlex en las
cuatro. La ablación del banco del paper (primal en ambos motores, umbral 40)
daba solver ×1.04 y `ibexopt` ×0.98: el dual gana en sistemas y el primal en
optimización.

**Primal en `ibexopt` (2026-10-13; `results/primal_banco_opt_2026-10-13.csv`,
`results/primal_all_2026-10-13.csv`; 4 semillas, 60 s, orden lejos).** Brazos:
dual (omisión), primal (`DFB_SX_PIVOTEO=primal DFB_SX_REVPRIMAL=1`) con
umbral 30 / 15 / siempre revisado, caja. **Elegido en el banco del paper**
(139 opt), contra el dual: p30 ×0.94 / ×0.96 / ×0.97, **p15 ×0.94 / ×0.96 /
×0.97** (peor ×1.03; launch ×0.75, himmel16 ×0.77, ex5_3_2 ×0.80), p0 ×0.95 /
×0.97 / ×0.98 (peor ×1.11); contra la caja p15 ×0.79 (dual ×0.83).
**Validado en `all`** (104 opt): p15 contra dual ×0.91 / ×0.90 / ×0.88, 309 =
309; contra la caja ×0.85 / ×0.85 / ×0.91 (dual ×0.93 / ×0.93 / ×1.04).
Ganan hs099 ×0.40, linear ×0.43, srcpm ×0.53, ex8_4_4-1/-2 ×0.54; pierden
hhfair1 ×1.19 y ex7_2_9 ×1.25 (el dual hace menos celdas). Falla corregida en
el camino (commit bd15955b): `pivot_revisado` salía por pivote chico sin
limpiar `col_primal_lista`; tras refactorizar, un pivote dual sobre esa
columna aplicaba la eta de la LU vieja (segfault en ramsey, 2 de 4
semillas; hallado con ASan). Los p15 repetidos tras el arreglo dan ×1.00
contra los anteriores. **Adoptado como omisión de `ibexopt` (2026-10-13)**; la cota
certificada monótona durante el LP queda sólo en `ibexsolve_dfb`.

**Por qué el dual gana en sistemas y el primal en optimización (2026-10-14,
`DFBH_BASECOMP=1`, semilla 1, 30 s).** La sonda compara el óptimo de cada cota
tibia con el de la anterior: columnas básicas nuevas y no básicas que siguen
siéndolo pero cambian de cota. La base cambia parecido en los dos programas
(sistemas 3–15 % de la base, optimización 5–19 %); **lo que difiere son los
lados**: en sistemas cambia de cota el 8–37 % de las no básicas (Bratu-0150 17
por LP, Discrete-Integralf2-19 12, Geneig-icse 3.5), en optimización el 1.4–7 %
(0.3–1.3 por LP). El dual reubica lados sin pivotear; el primal paga cada
cambio de lado con una iteración: pivotes primal ≈ dual + cambios de lado
(Bratu-0120 35.5 ≈ 9.9 + 19.4; Discrete-Integralf2-19 25 ≈ 7.5 + 12.4;
Geneig-icse 6.9 ≈ 2.2 + 3.5). En optimización casi nada cambia de lado y la
reubicación del dual por el signo del costo reducido mueve columnas que el
óptimo no mueve (linear 28.6 piv/LP en dual contra 7.6 en primal). La hipótesis
anterior («en sistemas la base es la misma») era falsa: es la base casi igual
con muchos lados opuestos. No está aislado por qué los lados cambian más en
sistemas.

**Banco final del paper (2026-10-14, `results/final_2026-10-14.csv`, resumen
en `results/final_2026-10-14_resumen.txt`; sonda misma caja en
`results/costo_2026-10-14/`).** Binario con todas las omisiones (solver: dual,
umbral 30; `ibexopt`: primal, umbral 15). 242 sistemas (sin los 22 excluidos)
× 2 semillas y 139 opt × 4, 60 s, 11 brazos. SC contra caja: opt ×0.75 / ×0.77
/ ×0.72, 131 / 130 resueltos; sistemas ×0.78 / ×0.76 / ×0.81, celdas
×0.88–0.92, 364 / 356 corridas. Ablación: arranque tibio del otro programa
sistemas ×1.04–1.05 (361), opt ×1.01–1.05; sin fusión ×1.27–1.43 (354); sólo
tableau ×1.08–1.15 (356); sin paso largo ×1.08–1.10; la perturbación no actúa
en `ibexopt` (búsqueda idéntica). Costo por pasada: sistemas ×0.50, opt ×0.39;
pivotes opt ×0.98 (con el dual era ×0.80), µs/pivote opt ×0.49. Paper al día
(Tablas 2 y 4, figura, resumen); la Tabla 3 no cambia.

**Sonda de lados sobre el banco del paper (2026-10-14, semilla 1, 20 s,
`DFBH_BASECOMP=200`, ≥ 1000 LPs tibios).** 162 sistemas: base cambiada 10 %,
cambios de lado 1.9 por LP (15 % de las no básicas), pivotes por LP tibio dual
3.1 / primal 4.9, primal − dual = 1.90 ≈ lados 1.87. 60 opt: base 17 %, lados
0.5 por LP (8 %), pivotes dual 1.8 / primal 1.9, primal − dual = 0.42 ≈ lados
0.49. **Corrige la nota anterior**: en el banco del paper el primal pivotea
~30 % MÁS que el dual en opt hasta 40 filas (menos arriba); lo de linear/hs099
(dual 2–4× más pivotes) es de relajaciones grandes de `all`. La ventaja del
primal en opt viene del costo por LP: la reubicación cambia las no básicas y
obliga a recalcular `x_B` desde la factorización (`compute_basics`, umbral 30
en ambos: ex5_3_2 61963 → 23365, launch 28760 → 11753; cpu −10–20 %;
ex8_4_4bis 0.67 → 1.26 piv/LP con cpu igual). Sin cronómetro por componente
todavía. Tabla `tab:sides` en el paper.

**Mecanismo completo de la regla por programa (2026-10-14).** Costo por LP
(`DFB_SX_PERF=1`, umbral 30 en ambos, semilla 1): el primal paga una iteración
por cambio de lado; el dual, un costo fijo por LP (reubicar cambia las no
básicas y obliga a recalcular `x_B`). Sistemas: Bratu-0120 222.8 → 422.2 µs/LP
con el primal (224 µs en pasos primales), Discrete-Integralf2-19 154 → 282.
Opt: launch 158 → 143, ex5_3_2 168 → 128 (basics 16–17 → 6–7 µs/LP, menos
pivotes), ex2_1_9 33.6 → 31.6 con 25 % más pivotes. **Por qué cambian más los
lados en sistemas** (`[basecomp2]`, banco del paper, omisión SC): las cotas de
caja cambian igual en los dos programas (14 % / 17 % de las no básicas); las
holguras de fila, 14 % en sistemas y 0.4 % en opt. La franja de cada fila
(ancho de la holgura / rango de la fila en la caja) es 0.13 en sistemas y 0.48
en opt; las filas de opt que cambian de lado son las angostas (0.04). Todas
las no básicas tienen las dos cotas finitas en ambos. Paper: §4 y Tabla
`tab:sides`; §7 (discusión) escrita; cotas infinitas y máquina declaradas.

**Rigor rehecho (2026-10-14, `DFBH_AB=1 DFBH_ABHULL=1`, semilla 1).** 3090
cotas donde SC contrae más que `CtcPolytopeHull` (antes 973: el resumen se
imprimía cada 100 nodos y contaba de menos): SC en el hull 2980 (96 %), lo
pasa en 27 (Discrete-Integralf2-13 22, Bellido 5), todas con cajas < 2e-5 de
ancho y exceso ≤ 3.4e-9 absoluto (2.8e-10 relativo), bajo la tolerancia 1e-9
del oráculo SoPlex, cuyo x* está fuera de la caja en 21 de las 27. Contador
`mas alla de la tolerancia del oraculo` = 0 en las cuatro. Diagnóstico
`[abhull-pasa]` en la sonda. Perfiles de desempeño (`fig:profiles`,
`paper/perfil_*.dat`): SC el más rápido en 96 % (sistemas) y 95 % (opt).

**Arranque tibio adaptativo (2026-10-14, `DFBH_ADAPT=u`).** Promedio móvil
de cambios de lado por LP tibio; al cargar cada nodo, dual (umbral de motor 30)
si supera `u`, primal (15) si no. Etapa de diagnóstico (7 sistemas + 9 opt, 2
semillas, 60 s), contra la omisión por programa: u = 1 sol ×1.02, opt ×1.16;
**u = 2 sol ×1.05, opt ×1.01**; u = 4 sol ×1.02, opt ×1.02. No gana en
ninguna: la elección por programa ya clasifica casi todo bien (en la ablación,
2 de 127 sistemas ganan con el primal; las pérdidas del primal en opt, hhfair1
y ex7_2_9, son de celdas y tienen pocos cambios de lado, 0.2), y unos pocos
nodos en la estrategia contraria cambian el camino (Geneig-icse con u = 2: 9
nodos primales, celdas 295 → 469). launch tiene 3.5 cambios de lado y aun así
gana con el primal: el predictor no es perfecto. No se pasó al banco completo;
la opción queda apagada.

**Orden de cotas para SoPlex (2026-10-14, `results/ordcaja_2026-10-14/`).**
Etapa 1 (109 instancias con caja de 1–60 s, semilla 1, 60 s) contra «lejos»:
`cerca` sistemas ×0.95 en tiempo con celdas ×1.02 (26 ganan > 10 %), opt ×1.00;
natural sistemas ×1.04, opt ×1.00. Etapa 2 (242 sistemas × 2 semillas, `cerca`):
×0.957 / ×0.949 / ×0.941 (> 0.1 / 1 / 5 s), celdas ×1.00–1.02, pero termina
352 contra 356: pierde Brown-19 y Brown-20-1 (celdas ×6, el orden cambia el
camino). Gana en los Bratu grandes (×0.62–0.81 con las mismas celdas: menos
costo por LP con el arranque en caliente de SoPlex). Con SoPlex en `cerca`, SC
(su orden) queda ×0.84 / ×0.85 / ×0.83 (antes ×0.78 / ×0.76 / ×0.81). En opt
SoPlex sigue en «lejos» (empate). **Paper pasado a SoPlex `cerca` en sistemas
(fila best of each)**: ×0.84 / ×0.85 / ×0.83, celdas 0.87–0.95, 364 / 352;
resumen 15–23 % y 5–13 %; Brown ×0.15 (la caja termina 4 de 12), Geneig-icse
×0.39, Bratu grandes ×0.92 (×0.72 con lejos en los dos), pierde Katsurabis-14
×1.20 (mismas celdas) y dos dentro del ruido. Perfiles: 271 corridas, SC el más
rápido en 92 %, la caja no termina 4 %. Interruptor nuevo `DFBH_ORDEN_CAJA`
(orden propio de la caja; sin él, la caja usa lejos cuando DFB usa el orden por
costo: la sonda vieja daba lo mismo que lejos, verificado en Bratu-0120,
ex5_3_2 y launch). En Bratu-0120 la caja con `cerca` hace 8.75 piv/LP contra
30.1 con lejos (165 contra 364 µs/LP): la figura de costo de sistemas se rehace
con `DFBH_ORDEN_CAJA=cerca` (en cola, `$M/costo2.sh`, binarios `*_oc`).

**Banco a 3600 s (2026-10-14, `results/l3600_2026-10-14/`).** SC con orden
lejos (binario final, 6 GB) contra la caja «lejos» de `l3600`/`all3600`
(2026-10-11/13, 5 GB), sin las excluidas; se emparejan las corridas que
terminan las dos. Sistemas del paper (48, 96 corridas): tiempo ×0.74, celdas
×0.76, 94 / 94 (CountercurrentReactors-30 no termina con ninguno). Opt del
paper (7, 28): ×0.74, celdas ×0.94, 28 / 28. Opt de `all` (16, 36): ×0.77,
celdas ×0.97, 33 / 25: SC termina ex6_2_13, hs101, hs112 y pollut (2400–3140 s)
y la caja no; matrix2 no encuentra punto factible con ninguno en las semillas
1–3 (NA). **SC con su propio orden** (best of each): sistemas contra SoPlex
`cerca` ×0.91, celdas ×1.01, 92 / 89 (contra lejos era ×0.72 / 0.76: Brown-30,
SC 1–2 s, SoPlex lejos 80 s y `cerca` no termina, pesaba ×0.85 en la media);
opt paper ×0.74, 28 / 28; opt all ×0.77, 33 / 25. SoPlex `cerca` contra lejos
a 3600 s: ×0.93, celdas ×0.95 (CountercurrentReactors2-24 ×0.03).

**Vacíos espurios en caliente (2026-10-14).** Con su orden propio, SC no
termina CountercurrentReactors2-26 (s2), -28 (s2) y -30 (s1) en 3600 s (600–800 mil
celdas; con lejos 13–38 mil). Traza de soluciones (`DFBH_SOLTRACE=1`, en
`ibexsolve_dfb`): mismo ritmo hasta ~147 s y 2200 soluciones, después nada.
Estados: infeasible 10.9 % de los LPs (lejos 0.4 %), pero vacíos probados 0.6 %
de los nodos; `DFBH_MARGEN`: 73311 rayos que no certifican, nnz 3.0 (los que
certifican 15.6), el cero en medio de la suma en el 99 %: el simplex en caliente
declara vacíos espurios, el certificado los rechaza (la cota es válida) y el
nodo no se poda. La caja del fallo 20000 (`DFBH_DUMPFALLO=N`) resuelta aparte
muere en 3 celdas con los dos órdenes: depende del estado heredado, no de la
región (anchos de 1e-16 a 1e6). **`DFBH_REINTENTA=1`**: un INFEASIBLE que no
certifica se resuelve otra vez en frío. Arregla los tres: 260 s / 18451
celdas, 373 s / 25739, 780 s / 39645 (lejos: 191, 336, 728 s). Banco de 60 s
(1040 corridas, contra el brazo sc de `results/final_2026-10-14.csv`):
sistemas t ×0.999/0.998/1.004, celdas ×1.000, terminan 361/364 (las 3 que pierde
estaban en 59.1–59.7 s); opt ×1.006/1.009/1.004, terminan 518/518; búsqueda
idéntica en el 96 %. **Desde 2026-10-15 está encendido por omisión** con el nombre
`DFBH_FRIO_SI_NO_CERTIFICA` (`=0` lo apaga); los binarios finales anteriores
quedan en `$M/ibexopt.fin` e `$M/ibexsolve_dfb.fin`. Las tablas de 60 s del paper
quedan con los binarios anteriores (el cambio es neutro ahí); el brazo propio
de sistemas a 3600 s se rehace con el binario nuevo (`$M/l3600r`).
Para comparar en igualdad, la caja negra recibe el mismo reintento (binarios
`*_oc`, 2026-10-15): si SoPlex declara infactible y el rayo no certifica, o no
entrega rayo, se borra su base (`clearBasis`) y se resuelve otra vez. Primeras
pruebas (cerca, semilla 1, 60 s): Brown-18 idéntico (sus pérdidas son estado
`Unknown`, 61730 de 71 mil LPs, no infactibles: el reintento no las toca);
Geneig-icse 549 → 625 celdas. Etapa de 60 s en sistemas en `$M/fcaja`. Etapa de
60 s (484 corridas, cerca): t ×1.002/1.003/1.003, celdas ×1.000, terminan 352/352,
celdas idénticas en 347. Neutro; se rehace la caja «cerca» a 3600 s (`$M/l3600rc`).

**SC propio a 3600 s con el reintento** (`$M/l3600r`, 96 corridas): arregla los
tres estancamientos (CCR2-26 s2 257 s, CCR2-28 s2 372 s, CCR2-30 s1 777 s), pero
CCR2-28 s1 pasa de 321 s / 21559 celdas a 1640 s / 337457 (reproducible), y
CCR-30 s1 deja de terminar (antes 1138 s; s2 no termina con nadie). Mecanismo de
CCR2-28 s1: **no es DFB, es el ajuste adaptativo de ACID**. Sonda
`DFBH_DUMPNODO=N1,N2` (vuelca la caja de entrada de esas llamadas a CtcDFBHull) y
`DFBH_SOLO=1` en `ibexsolve_dfb` (aplica cada contractor una vez a la caja): las
cajas del estancamiento (llamadas 25k, 100k, 200k) las vacía ACID solo, no DFB
(ni una pasada ni el punto fijo). `DFBH_ACIDTRACE=N` (subclase `AcidSonda`, misma
conducta: 21559 celdas): sin reintento `nbcidvar` se mueve entre 18 y 103; con
reintento cae a 0 en la llamada ~18000 y se queda en 0–3 hasta ~36000. ACID
reajusta cada 1000 llamadas con 50 de ajuste sobre 2·nbcidvar variables (mínimo
2): si una ventana de 50 cae en cajas donde recortar no rinde (ganancia < 0.002),
deja de recortar y tarda en recuperarse. El reintento solo cambia el camino; la
misma inestabilidad le puede pasar a cualquier configuración.

**Tabla de 3600 s en igualdad (paper, `tab:long`, 2026-10-15;
`results/frio_2026-10-15/`).** Sistemas: SC propio con reintento contra SoPlex
`cerca` con reintento: tiempo ×0.96, celdas ×1.04, terminan 94 / 89 (SC termina
Brown-30 y Brown-30-1 en 1–2 s y CCR2-24 s2); sin CCR2-28 s1 (×4.2, ACID) ×0.94 y
celdas ×1.01. Opt (sin reintento en ninguno): paper ×0.74, celdas ×0.97, 28 / 28;
`all` ×0.77, celdas ×0.97, 33 / 25. La caja con reintento a 3600 s cambia poco
(CCR2-28 s2 466 → 349 s). Los binarios finales se recompilaron desde el fuente
(reintento en los dos motores; respaldo `$M/*.fin2`).

**3600 s con reintento en los dos motores (2026-10-16,
`results/frio_2026-10-15/l3600z/`, `tab3600z.py`, 320 corridas).** Sistemas,
orden igual (SC lejos / SoPlex lejos): t ×0.76, celdas ×0.76, 94 / 94; sin
Brown (SoPlex 72–78 s, 15 mil celdas) ×0.89; sin Brown ni CCR2-24 s1 (SoPlex
3253 s, 633 mil celdas) ×0.92, celdas ×1.00. Opt paper (SC propio / SoPlex
lejos): ×0.60, celdas ×0.79, 28 / 28; lo baja bearing: el reintento cambia la
búsqueda de SC (s2/s3 112 → 4.8 s, 40132 → 1910 celdas, comprobado con
`DFBH_FRIO_SI_NO_CERTIFICA=0`; menos de 1000 reintentos), s1/s4 empeoran
72–81 → 93 s; sin bearing ×0.78, celdas ×1.04. Opt all: ×0.77, celdas ×0.96,
33 / 25 (igual que sin reintento). SoPlex con reintento en opt ×0.99 contra
su corrida vieja. matrix2 NA en 3 semillas con los dos (memoria, igual que
antes). Paper `tab:long` al día con las cuatro filas.

**CORRECCIÓN: la sombra comparaba otro LP en sistemas (2026-10-19).** La caja solo
recargaba las filas fusionadas con `DFBH_MOTOR=caja` o `DFBH_ABCAJA`; con
`DFBH_SOMBRA` SoPlex resolvía el LP SIN fusionar (mismo poliedro, el doble de filas).
Corregido (la sombra también recarga). Con el LP correcto, en sistemas SoPlex NO
pivotea mucho menos: Katsurabis-14 10.0 vs 9.1, Bratu-0060 7.9 vs 10.9 (SoPlex más),
Bellido 4.5 vs 3.5, Geneig-icse 3.5 vs 5.1. El «vértice de caja» de Katsurabis-14
(10.4 vs 1.2, 22 vs 6.5 estructurales) era artefacto; el banco (`DFBH_ABCAJA`) estaba
bien. En Katsurabis-14 SoPlex gana por costo por pivote a ~10 pivotes/LP (110 vs 172
µs por solve). En opt (también hay fusión con igualdades) el patrón se mantiene:
SoPlex 13–40 % de LPs sin pivotes, ex6_2_10 1.5 vs 2.5, ex7_2_3 1.3 vs 1.9; **SC
arrancando desde la base de SoPlex** (`DFBH_SOMBRA_BASE=1`, `instalar_base` por
pivotes de tableau, sin LU: la LU de SoPlex rechazaba hasta la base propia de SC)
pivotea como SoPlex en ex2_1_9 (1.68 vs 1.87), ex6_2_10 (1.92 vs 1.89), chem (3.33
vs 3.34); no en ex7_2_3 (2.03 vs 1.27). La diferencia está en la cadena de bases.
Paper corregido (Where it loses, Limits, Why SC wins).

**Banco de pasadas y cascada de atribución (2026-10-19, `results/banco_2026-10-19/`).**
71 instancias (36 sistemas, una por familia; 35 opt por cuantiles de tiempo), 1000
cajas capturadas por instancia (`banco_cap.sh`), reproducidas con `ibexreplay` y
`DFBH_ABCAJA=1`. Pasada SC/caja: sistemas ×0.55 (mediana 215 vs 412 µs, SC más
barato en 33/36), opt ×0.47 (121 vs 261, 35/35); contracción idéntica. Excepciones:
DiscreteBoundary-0200 ×2.03 (en la sombra, mismos pivotes: la diferencia son los
LPs que resuelve cada uno), Katsurabis-14 ×1.63, Bratu-0120 ×1.13. Cascada
(`cascada.sh`, `ibexreplay_sx` con SoPlex parchado `SOPLEX_SIN_VERIF`): agregar a
SC factorizar por LP 45/41 % de la brecha log (sist/opt), recalcular vectores
7/6 %, verificar 19/16 %, copiar 2/1 % (total 73/65 %); quitar a SoPlex verificar
18/10 %, refactorizar al óptimo 1/3 %, presolve 22/23 %, escalado −22/+4 %. Tablas
`tab:bank` y `tab:cascade` en el paper; el resumen usa ×0.47 / ×0.55.
**Caminata en la cara óptima (`DFB_SX_CAMINA=K`):** no sirve (Katsurabis-14
estructurales básicas 22.6 → 22.1, pivotes 10.0 → 11.6 con K=10; DiscreteBoundary
y Bratu-0120: cero pasos posibles). Apagada; resultado negativo en Limits.

**Por qué SoPlex pivotea menos: el vértice en que termina (2026-10-19).**
Con la sombra en `ibexreplay` (`DFBH_SOMBRA=1`, mide también la base al empezar y
al terminar cada LP): SoPlex termina en bases con las estructurales en sus cotas
(Katsurabis-14 ~6.5 estructurales básicas de 30 filas, Bratu-0060 ~25) y la
cota siguiente le sale en 0–1 pivotes (x_k básica al empezar 6–31 %); SC termina
con muchas estructurales básicas (22.6 y 50.4; x_k básica 32–79 %) y parte a ~10
infactibilidades (dual: básicas fuera de cota; primal: costos reducidos de signo
malo). En los LPs que la caja hace con ≤1 pivote, SC hace 39 (Katsurabis-14).
SoPlex en frío también los resuelve con 0.7 pivotes: son LPs casi triviales.
Descartado como causa o como remedio: warm start primal (24.6 piv/LP), devex
primal (`DFB_SX_DEVEXPT`, ×0.96 en opt), DSE, perturbación, motor, cotas de
holgura sin intersectar (`DFBH_HOLGURA=lados`: idéntico en sistemas, contrae
menos en opt), limpieza de básicas en cota (solo 10–20 % lo están), desempate
hacia las esquinas (`DFB_SX_ESQUINA=eps`: la base final no cambia). Todos
apagados. Abierto: cómo hacer que SC termine en vértices "de caja".

**Sombra LP a LP (2026-10-19, `DFBH_SOMBRA=1` en `ibexreplay`).** Cada LP que
resuelve SC se resuelve también con la caja (mismo k, sentido y caja, en la misma
cadena). En la misma secuencia SoPlex pivotea menos sobre todo porque termina en
vértices que ya son óptimos para la cota siguiente: LPs sin pivotes de la caja
ex6_2_10 51 %, ex7_2_3 50 %, Bratu-0060 56 % (SC ~1.4 %). Piv/LP SC/caja:
ex6_2_10 2.49/0.96, ex7_2_3 1.63/0.95, chem 3.92/3.46, Bellido 4.46/2.37,
**Katsurabis-14 10.43/1.22** (µs/LP 181 vs 34: explica su pérdida), Geneig-icse
3.08/4.43 (SC menos). En los LPs sin pivotes de la caja, los pivotes de SC no son
degenerados (2 de 4809 en ex2_1_9). En Katsurabis-14 no lo explican: warm start
primal de SC (24.6 piv/LP), devex primal (22.8), sin perturbación (21.2), DSE
(15.0); ni la representación o el algoritmo de SoPlex (1.2–1.4 en todos).
Hipótesis abierta: las cotas de holgura de SC (lado ∩ evaluación por intervalos)
contra los lados del solver en SoPlex; **la prueba hecha fue inválida** (las
cotas no calzaban con las filas: SoPlex declaró infactibles los 5823 LPs).

**Banco de pasadas y error en la sonda (2026-10-18).** Captura de cajas
(`DFBH_CAPTURA=archivo`, `DFBH_CAPTURA_N`, `DFBH_CAPTURA_CADA`, binarios `*_cap`) y
reproductor `ibexreplay problema.bch opt|sol cajas.txt` (mismo contractor que
`ibexopt`/`ibexsolve_dfb`, mismas omisiones del optimizador): determinista, y con
`DFBH_ABCAJA=1` mide SC y la caja sobre las mismas pasadas sin ruido de búsqueda.
**Error corregido en `DFBH_ABCAJA`:** «gana DFB» contaba los nodos que vacían los
dos motores. Corregido: Geneig-icse 128 vs 9 (antes 224), CCR2-9 66 vs 2 (105),
Brown-18 1 vs 1, Bellido 4 vs 4 (antes 138 vs 4). La «precisión del dual» se
sostiene en Geneig-icse y CCR2-9, no en Bellido. Tabla samebox corregida.
**Devex primal en el motor tableau (`DFB_SX_DEVEXPT=1`):** pivotes por LP ×0.96,
costo por LP ×1.00 en 20 opt; el ×0.87 de tiempo es ruido de búsqueda. Apagado.
**Pivotes:** comparar piv/LP exagera la brecha (la caja resuelve más LPs triviales);
por pasada, ex2_1_9 370 k vs 366 k.

**Emulación de SoPlex dentro de SC (2026-10-18, `DFB_EMU_FACT/INIT/VERIF/COPIA`,
binarios `*_emu`; datos en `results/directa_2026-10-17/emu/`).** Se le agrega a
SC, en cada LP y sobre vectores aparte (la búsqueda no cambia), el trabajo por
llamada de SoPlex. µs/LP de SC base → con todo (caja): ex2_1_9 7.7 → 13.3
(22.6), dipigri 6.8 → 12.2 (22.6), Geneig 11.2 → 18.2 (26.9), Bellido
23.0 → 33.0 (35.6), chem 24.1 → 31.8 (37.2). Factorizar por LP +25–40 %,
verificar +10–25 %, init y copia poco. Explica ~40 % de la brecha en los LPs
chicos; el resto es presolve/escalado por relajación y maquinaria no emulada.
En Katsurabis-14 (203 vs 105), Bratu-0060 (110 vs 80) y launch (93 vs 79) el LP
de SoPlex ya es más barato sin emular: pivotea menos y más barato. Párrafo en
la sección de costo. **De paso:** SC con `DFBH_ORDEN=cerca` pivotea mucho menos
en algunos opt (ex8_4_4 5.8 → 2.8 piv/LP, 31 → 21 s; himmel16 6.0 → 3.6,
20.8 → 16.1 s; ex5_3_2 9.9 → 5.2) y gana en 6 de 8: probar en todo el banco.

**SoPlex 7.1.6 (2026-10-17).** Compilado con el namespace renombrado
(`-Dsoplex=soplex7`, `$M/sx7ns`) para convivir con 4.0.2. (1) Caja con SoPlex 7
(`DFBH_CAJA_SOPLEX7=1`, `ibex_Caja7.{h,cpp}`, binarios `*_s7`): en 8 problemas
es ×1.2–1.9 más lenta que la caja directa 4.0.2 (ex2_1_9 14.8 vs 7.6 s mismas
celdas; sin presolve ni escalado 12.0 s); Brown-18 85 celdas vs 67. La caja
del paper sigue con 4.0.2. (2) SC con la LU de SoPlex 7 (`-DDFB_LU7`,
`ibex_DFBBasis7.cpp`, objetos en `$M/obj7`, binarios `*_lu7`): mismas celdas
en todos los sistemas que terminan; Bratu-0120/0150/0200 ×1.08–1.14 más lento,
el resto igual; en opt con tope, sambal procesa 34 % menos celdas. Frase en §5.

**Brown era un lastre de la caja: una fila libre con coeficientes de 1e108
(2026-10-17).** Volcando los LPs fallidos (`DFBH_CAJA_DUMP=dir`, formato
`.txtlp` exacto; el MPS no se puede escribir, `XMPSWR02`) y resolviéndolos
fuera de Ibex (`$M/lptest/resolver2.cpp`), SoPlex 4.0.2 **y 7.1.6** fallan en
los 40. En Brown-18 la fila 17 (la ecuación producto, forma afín sobre una caja
de ancho 5.5e6) tiene lados ±1e100 (libre para SoPlex) y coeficientes de
6e107 a 1.1e109: no restringe nada pero quiebra a SoPlex. Con
`DFBH_CAJA_SINLIBRES=1` (la fila libre se le pasa vacía) la caja resuelve
Brown-18/20-1/30 en 67/65/71 celdas, igual que SC, en 0.31/0.38/2.4 s. Geneig-icse
no cambia (767 vs 295 celdas): ahí el mecanismo es otro. El redondeo hacia
arriba de GAOL (SoPlex corre con FE_UPWARD dentro de Ibex) no es la causa.
**Orden de la caja directa en opt:** cerca es ×0.93–0.96 contra lejos con las
mismas celdas (antes «empatan», medido con la caja vieja en 19 problemas).
Barrido corregido en curso (`d60s.sh`, binarios `_dir2`, caja sin filas libres,
los dos órdenes en opt y sistemas). Las corridas de 3600 s (`d3600`) y la sonda
de costo (`costo5`) se hicieron con la caja sin corregir: rehacer.

**Por qué SoPlex pierde Brown (2026-10-17, `DFBH_CAJA_DBG=1`, `DFBH_CAJA_VERB=1`,
`DFBH_CAJA_FRIO=1` en la caja directa).** En Brown-18 el 90 % de los LPs de
SoPlex terminan con estado RUNNING: lanza `SPxInternalCodeException("XLEAVE02
This should never happen")` en `getLeaveVals` (estado de base inconsistente
en el paso de salida del dual); `_evaluateSolutionReal` lee `m_status` y queda
RUNNING (Ibex lo ve como Unknown: 51 k de 60 k). Pasa igual en frío
(85 %), con representación por columnas, tolerancia 1e-9 y sin
presolve/escalado (99 %). Geneig-icse: 1.3 % no óptimos (tope de iteraciones
e infactibles). La tolerancia escalada (`sintau`) NO explica Brown ni Geneig:
celdas idénticas sin ella. Paper corregido (§samebox, §5 y contribución 1:
la factibilidad dual por reubicación es estándar en LPs acotados, Koberstein).
Pendiente: no se aisló la causa del XLEAVE02 (¿filas casi dependientes de
Brown?); con SoPlex más nuevo podría no pasar.

**Caja negra directa (2026-10-17, `DFBH_CAJA_DIRECTA=1`, binarios
`ibexopt_dir` / `ibexsolve_dfb_dir`).** El `LPSolver` de Ibex guarda las cotas
de las n variables también como n filas (duales para Neumaier–Shcherbina): en
ex2_1_9 SoPlex ve 14 filas contra 3 de SC. La caja directa usa su propio
SoPlex con los mismos parámetros (copiados DESPUÉS de `clearLPReal`, que
devuelve el sentido a maximizar), las mismas filas y lados, cotas solo de
columna e `ignoreUnscaledViolations` como Ibex. Validada con
`DFBH_CAJA_CMP=1`: mismos estados y óptimos en 160 000 LPs de ex2_1_9. Costo
por LP −16 % en ex2_1_9 y hs117, igual en launch. Callgrind lado a lado sobre
los mismos LPs (`pasada_ab_caja` / `pasada_ab_sc`, ex2_1_9): caja 84 k
instrucciones por LP, SC 26.5 k; el pivoteo cuesta lo mismo (~10 k); el resto
de la caja es maquinaria por llamada (`init` 2.1 veces por LP 18 k, factorizar
9 k, presolve y recarga por relajación ~10 k, guardar/verificar ~10 k,
interfaz ~14 k). Barrido de 60 s con la caja directa (`d60.sh`, `results/directa_2026-10-17/`):
caja directa ×0.96–0.98 opt, ×0.95–0.97 sistemas; SC contra ella opt ×0.79/0.80/0.82,
sistemas best ×0.91, same order ×0.88–0.89; termina 364/352 (hayes1 rehecho:
tenía el mismo assert afín). **Sonda de costo con la caja directa (`costo5.sh`):
SC NO pivotea menos** (sistemas 3.0 vs 2.0 piv/LP, ×1.33; opt 2.7 vs 2.5,
×1.12): el «×0.73 de pivotes» era artefacto de las filas de cotas, que le
agregaban pivotes a SoPlex. Toda la ganancia es costo por LP (×0.69 sistemas,
×0.41 opt); pasada ×0.62 / ×0.44. Callgrind con la caja directa (ex2_1_9):
caja 74 k instr/LP, SC 26 k; pivoteo 9 k vs 10 k; SoPlex `init` 16 %,
guardar/desescalar/verificar 20 %, presolve 11 %, factorizar 7 %. Brown-18 con
la caja directa: 250 de 1862 LPs sin óptimo (antes 176).
**Ojo al compilar:** cambiar `ibex_CtcDFBHull.h` exige recompilar también
`ibex_LSmearDFB.o`, `ibex_Optimizer05Config.o` e `ibexsolve_dfb.o` (si no,
core dump por tamaño de clase distinto).

**Por qué un LP de la caja cuesta más (2026-10-16, `DFBH_CAJA_RELOJ=1
IBEX_PH_STATS=1`, `DFBH_CAJA_STATS=k` imprime la estadística de SoPlex cada k
LPs; `DFBH_CAJA_COTAS=cambio` carga solo las cotas que cambiaron).** En opt SC
pivotea casi igual que la caja (×0.94); la ganancia es por LP (×0.39). En
ex2_1_9 un LP de la caja son 31 µs: `solve` de SoPlex 24, cargar cotas y
objetivo 2.5, copiar la solución 1.1, leer multiplicadores 1.8; un LP entero de
SC son 8 µs. SoPlex factoriza ~1.2 veces por LP (1.16–1.36 en ex2_1_9, chem, hs117 y
launch), **pero no por el cambio de objetivo** (`unInit()` no descarta la
factorización; los LPs sin pivotes no factorizan): refactoriza al llegar al
óptimo si hubo alguna actualización (`spxsolve.cpp:403/717`, «eliminate
accumulated errors»). **Y no es el costo:** con un SoPlex parchado
(`SOPLEX_REFACT_FIN=k`, `$M/sxp/libsoplex.a`, binarios `ibexopt_sx`,
`ibexsolve_dfb_sx`) bajan a 0.25 LU/LP y el µs/LP no cambia (27.7 → 27.0).
Callgrind dentro de `resolver_caja` (ex2_1_9, instrucciones): pivotes
(`leave`) 12 %, factorizar 10 %, `init()` al inicio de cada solve 20 %,
guardar/desescalar/verificar la solución 20 %, presolve en LPs de una
relajación nueva 11 %, `set_bounds` de Ibex 9 %. Ninguna parte domina: es el
costo por llamada de un solver de propósito general. La primera versión del
paper (commit e2b7fe28/56f5114d) atribuía todo a la factorización; corregido. Cargar solo las cotas cambiadas ahorra 3 µs y no
evita la factorización: no es un lastre de nuestro arnés. En launch (10
pivotes por LP) la factorización se amortiza y el LP cuesta igual (95 vs
101 µs): es donde SC más pierde. Panel derecho de la figura de costo ahora es
tiempo por LP (sistemas ×0.59, opt ×0.39). Matiz en «Limits»: el costo por llamada es el de
SoPlex 4.0.2 vía Ibex; otros solvers tendrán otro, no medido.

**Banco único de opt de 240 problemas (2026-10-16, `results/valid_2026-10-16/`,
60 s en `results/final_2026-10-16.csv`).** A los 139 del paper se sumaron los
102 de `ibex2026/benchs/optim/all` que no estaban ni excluidos; robot se
excluye (`cuelga`: no respeta `-t`, también con el Ibex sin modificar), quedan
101. Corridos a 3600 s en sc/caja/lejos × 4 semillas (`h3600.sh`) y los 8
brazos de ablación a 60 s (`fix.sh`). Ignacio decidió un solo banco (no
validación aparte); en el paper se dice que SC se ajustó en 139 y que en los
otros 101 gana lo mismo (×0.79 en los dos a 3600 s).
- La caja no terminaba a 3600 s ex6_2_13, hs101, hs112 y pollut (SC sí, en
  2400–3100 s): con 7200 s (`c7200.sh`) termina en 3850–4390 s con las mismas
  celdas, SC ×0.56–0.80. Entran como pares; no se cuentan como victorias.
- 60 s, 240 problemas: best ×0.77 / 0.78 / 0.80, celdas 0.97–0.99, 217 / 216;
  same order ×0.79 / 0.80 / 0.81. Pierden >5 %: launch ×1.43 y ex8_4_4 ×1.12
  (3722–4130 celdas vs 2638–3374, como launch). Ablación opt: tibio ×1.07,
  tableau ×1.04, revisado ×1.03. Perfiles: SC más rápido en 92 % (489 corridas).
- 3600 s, los 24 con alguna corrida en el tope de 60 s: ×0.73, celdas 0.91,
  96 / 80 (sin bearing ×0.78, celdas 0.99). Reemplaza las filas «7 del banco» e
  «Ibex hard» de `tab:long`.
- **Bug del linearizador afín (commit `32446e22`):** una forma afín inactiva
  (vacía o no acotada) conserva `size()==nb_var` y `LinearizerAffine2` la leía;
  sin `-DNDEBUG` abortaba en el `assert` de `AffineMain::val`. Se agregó
  `af2.is_actif()` a los cuatro chequeos. Afectaba wall (los NA de siempre),
  matrix2, concon, eigmaxc, ex8_4_5(-1), haldmads, makela3, polak1–5; todos
  terminan con `ibexopt_fix`, y donde no abortaba da las mismas celdas.
  Pendiente: los barridos viejos usaron `./ibexopt` sin el arreglo; reemplazar
  `./ibexopt` por `ibexopt_fix` antes del próximo barrido.
- Todo el paper usa este banco a 60 s salvo `tab:long` (3600 s). Sonda de
  costo extendida (`costo4.sh`, 101 nuevos con `ibexopt_fix`): opt pasada ×0.43
  (116 vs 270 µs, n=175), pivotes ×0.94, µs/pivote ×0.59; sistemas igual.
  Sonda de lados a 60 s (`lados60.sh`, antes 20 s; `lados_ana.py`): 165
  sistemas casi idénticos; 140 opt: base 16 %, lados 0.49 por LP (6 %), pivotes
  dual 1.9 / primal 2.0, holguras cambian 0.3 %, franja media 0.48.
- En §5: en los 101 agregados SC es 14–22 % más rápido (×0.78/0.79/0.86),
  contra 23–29 % en los 139 de ajuste; el banco de ajuste favorecía a SC.

**Sonda de costo con binarios finales (2026-10-16, `results/frio_2026-10-15/costo3/`).**
Sistemas: pasada ×0.53 (mediana 422 vs 826 µs), pivotes ×0.73, µs/pivote ×0.75,
igual que antes. Opt: ahora la caja usa su orden (`DFBH_ORDEN_CAJA=lejos`; la
sonda vieja de opt no lo fijaba): pasada ×0.39, pivotes ×0.87 (2.3 vs 2.8 por LP;
antes ×0.98), µs/pivote ×0.57 (antes ×0.49). Figura y texto del paper al día.
Scatter en `results/frio_2026-10-15/scatter2.png`.

**Banco de 60 s rehecho con el reintento en los dos motores (2026-10-15,
`results/final_2026-10-15.csv`, 11 924 corridas, 12 brazos con `cajacerca`;
resumen en `results/final_2026-10-15_resumen.txt`).** El reintento es neutro:
contra `final_2026-10-14` celdas idénticas en 96–99 % de las corridas, tiempo
×0.996–1.013 por brazo; mismos NA (wall y los -1000/0400+ por memoria). Paper:
opt best ×0.75 / 0.77 / 0.71, 131 / 131 (la caja ahora termina ex8_4_4bis s2);
sistemas best ×0.84 / 0.85 / 0.82, 362 / 352 (SC pierde CCR2-20 s2 y
Yamamua1-0030sp s2, ambas a <1 s del tope y con menos celdas). Pierde >5 %:
Katsurabis-14 ×1.19, CCR2-18 ×1.10, Bratu-0180 ×1.06 (51 vs 47 celdas); launch
opt ×1.43. Ablación: `sintau` en opt pasa de ×1.02 a ×1.00 (sin efecto en
ningún banco; Ignacio decidió dejarla, es protección en cajas angostas).
Perfiles (modo `bb`, `perfil2.py`): SC más rápido en 93 % / 94 %; SC no
termina 1 corrida opt (ship-1 s3, tope).

**Brecha en las relajaciones grandes (2026-10-09).** x_B incremental también en
los pasos primales (commit 41d22ab8): Bratu-0150 36.9 → 31.2 s, Bratu-0120
25 → 21 s, Discrete-Integralf2-19 2.38 → 1.88 s, celdas iguales. Devex primal
en el motor revisado (`DFB_SX_DEVEXP=1`, **apagado por omisión**, falta banco):
Bratu-0150 pivotes 897k → 811k (34.8 → 33.1 por LP), 30.8 → 29.5 s;
Discrete-Integralf2-19/21 ×0.88; BroydenTri-0100 ×1.02; Bratu-0175 sigue en
el límite. Test unitario sin discrepancias.

*Mecanismo de lo que queda* (Bratu-0150, caja negra 18.2 s, 721k pivotes;
DFB+devex 29.4 s, 811k): la diferencia es sobre todo de costo por pivote, no
de pivotes. `vSolveUright` (FTRAN de la entrante, todo desde
`preparar_cambio` del paso primal) cuesta 2.46 s en DFB contra 0.45 s en
SoPlex. La relajación compo tiene m = 300 filas para n = 150 (dos
linealizaciones por ecuación), la base de DFB es de dimensión m y la FTRAN
sale densa (164 de 300 no ceros de media). SoPlex elige **solo** la
representación por filas (base de dimensión n): forzándole columnas
(`DFBH_CAJA_REP=col`, base de dimensión n + m con las filas de cota) sube a
24.8 s; con `DFBH_CAJA_REP=row`, 18.2 s, idéntico a la omisión. O sea, ~2/3
de la brecha es la representación. El período de refactorización no influye
(10/25/50/100: 30.4/29.0/29.3/31.3 s). Cerrarla pediría que el motor
revisado trabaje con la base de dimensión ≤ n (núcleo de columnas
estructurales básicas, o representación por filas) cuando m > n.

**Fusión de filas gemelas (2026-10-09, `DFBH_FUSIONA`, por omisión).** La
causa de m = 2n no era la linealización compo: el linealizador afín
(`linearize`, camino fresco) parte cada ecuación en dos filas unilateras con
la MISMA normal. En la forma homogénea basta una fila con rango
[lhs, rhs]; `linearize_fresca` fusiona filas consecutivas idénticas (si los
lados se cruzan las deja, para que el LP pruebe el vacío) y la caja negra
recibe la misma relajación fusionada (se recarga `mylinearsolver`). Prueba
puntual (`--lin=art`, `-t 100`), tiempo en s (celdas):

| instancia | DFB sin fusión | DFB | DFB devex | caja sin fusión | caja |
|---|---|---|---|---|---|
| Bratu-0120 | 21.0 (33) | 9.9 (33) | 9.6 | 12.4 (35) | 10.7 (35) |
| Bratu-0175 | 88.1 (55) | 44.0 (55) | 42.4 | 47.0 (59) | 41.8 (59) |
| Bratu-0200 | 71.6 (31) | 34.4 (31) | 32.3 | 43.3 (37) | 37.7 (37) |
| Discrete-Integralf2-21 | 3.14 (11) | 1.30 (11) | 1.20 | 2.03 (15) | 1.74 (13) |
| DiscreteBoundary-0200 | 17.3 (5) | 16.4 (7) | 11.0 (5) | 19.2 (7) | 18.0 (7) |
| Brown-15 | 0.28 (73) | 0.24 (65) | 0.24 | 0.94 (175) | 0.27 (65) |
| CountercurrentReactors2-9 | 1.39 (955) | 1.22 (957) | 1.24 | 38.0 (20 083) | 1.43 (985) |
| Geneig-icse | 18.6 (559) | 8.2 (431) | 8.2 | 35.9 (1779) | 17.3 (455) |

Dos consecuencias. (1) La brecha de Bratu/Discrete se cierra: con base de
dimensión m = n el motor revisado queda a la par o mejor que SoPlex por
filas. (2) **Las filas gemelas eran también la causa de buena parte de los
vacíos espurios de SoPlex** (casi-igualdades: dos semiespacios con la misma
normal y lados a ~1e-14): CountercurrentReactors2-9 20 083 → 985 celdas,
Brown-15 175 → 65. La ventaja de DFB por certificación en esas familias era
en parte un artefacto de la relajación; con la fusión la comparación es más
justa. Banco con los defaults finales
(`results/defaults_finales_fusion_2026-10-09.csv`, 2 semillas solver, 4
`ibexopt` sin alkyl/alkylbis/bearing; caja con tolerancia 1e-6 / 1e-9):

| comparación | solver cpu > 0.1 / > 1 / > 5 s | celdas | terminan | `ibexopt` cpu |
|---|---|---|---|---|
| **DFB orden lejos contra caja (controlada)** | **×0.88 / ×0.87 / ×0.93** | ×0.93 / ×0.90 / ×0.95 | 360 contra 356 | **×0.78 / ×0.79 / ×0.70** |
| DFB orden propio (costo) contra caja | ×0.84 / ×0.83 / ×0.89 | ×0.93 / ×0.90 / ×0.95 | 360 contra 356 | ×0.75 / ×0.77 / ×0.71 |
| DFB fusión contra sin fusión | ×0.80 / ×0.81 / ×0.84 | ×0.99 | 360 contra 354 | ×0.97 |
| caja fusión contra sin fusión | ×0.89 / ×0.90 / ×0.87 | ×0.95 / ×0.95 / ×0.91 | 356 = 356 | ×1.02 / ×1.04 / ×1.08 |
| devex primal contra Dantzig | ×0.99 / ×0.98 / ×0.98 | ×1.00 | 360 = 360 | ×1.00 |
| (antes) DFB sin fusión contra caja sin fusión | ×0.92 / ×0.93 / ×0.87 | ×0.89 / ×0.85 / ×0.81 | 354 contra 356 | ×0.77 / ×0.82 / ×0.78 |

Sin soluciones distintas contra la caja. Por familia (orden lejos contra
caja): Bratu queda neutro (ya no hace falta apartarla), gana Brown ×0.31,
Geneig-icse ×0.61, Bellido ×0.68, Geneig ×0.71; pierden sólo
Discrete-Integralf2 ×1.12 y Katsurabis-14 ×1.43. Con el orden propio
Discrete-Integralf2 pasa a ×0.86. Devex: gana sólo Bratu-0210 ×0.62 y
DiscreteBoundary-0200 ×0.67, nada pierde; queda opcional (no cambia la
historia y obligaría a remedir el brazo lejos).

Sonda misma caja (`DFBH_ABCAJA`) con filas fusionadas (2026-10-10; la sonda
ahora también recarga la caja negra con la relajación fusionada). LPs sin
óptimo DFB / caja: Geneig-icse 187 / 743 (antes ~26 000 en la caja), Brown-18
34 / 286, CountercurrentReactors2-9 116 / 151, Bellido 151 / 0. Nodos donde
cada uno contrae más: Geneig-icse 470 / 0, CCR2-9 117 / 5, Bellido 138 / 5,
Brown-18 32 / 10. Con la caja en 1e-9 persiste (Geneig-icse 544 / 1, CCR2-9
49 / 5): la ventaja de precisión del dual no es la tolerancia, sin aislar.

**Ablación con los defaults finales y filas fusionadas (2026-10-10,
`results/ablacion_fusion_2026-10-10.csv`).** Contra el completo (orden por
costo, umbral 40, primal en ambos motores, dual sólo en frío):

| brazo | solver cpu > 0.1 / > 1 / > 5 s | terminan | `ibexopt` cpu | `ibexopt` resueltas |
|---|---|---|---|---|
| dual puro (`DFB_SX_PIVOTEO=dual DFB_SX_REVPRIMAL=0`) | **×0.96 / ×0.96 / ×0.98** | 361 / 360 | ×1.02 / ×1.03 / ×1.05 (celdas ×1.01–1.06) | **131** / 129 (chembis, ship-1) |
| revisado en dual, tableau en primal (`DFB_SX_REVPRIMAL=0`) | ×0.98 / ×0.97 / ×0.98 | 360 | ×1.00 | 129 |
| umbral 20 | ×0.97 / ×0.98 / ×0.99 | 360 | ×0.99 | 129 |
| sólo revisado | ×0.97 / ×0.97 / ×0.97 | 361 | ×1.04 / ×1.02 / ×1.00 | 129 |
| sólo tableau | ×1.08 / ×1.11 / ×1.13 | 352 | ×1.01 | 129 |
| sin fusión | ×1.25 / ×1.23 / ×1.19 | 354 | ×1.02–1.04 | 129 |
| sin paso largo | ×1.00 | 360 | ×1.02–1.03 | 129 |
| sin devex | ×1.01 | 360 | ×1.02–1.06 | 129 |
| sin τ | ×1.00 | 360 | ×1.00–1.02 | 129 |

Celdas del solver iguales en todos (±1 %). **Con la fusión cambió el
balance**: el primal revisado, que ganaba con m = 2n, ahora pierde contra el
dual en los nodos revisados (Bratu-0150…0210 ×0.47–0.81 en dual); el dual
puro gana en el solver y resuelve dos problemas más en `ibexopt`, pagando
×1.02–1.05 de cpu allí. Pierde BroydenTri-0100 ×1.36 y Katsurabis-14 ×1.25.
Pendiente decidir la omisión (dual puro, o dual en revisado + umbral 20) y
remedir el brazo controlado.

**Decidido: dual puro como omisión** (commit 8918195f). Comparación
controlada con la omisión nueva (`results/controlado_dual_puro_2026-10-10.csv`):

| comparación | solver cpu > 0.1 / > 1 / > 5 s | celdas | terminan | `ibexopt` cpu | resueltas |
|---|---|---|---|---|---|
| **SC orden lejos contra caja** | **×0.83 / ×0.82 / ×0.88** | ×0.93 / ×0.90 / ×0.95 | 360 / 356 | **×0.80 / ×0.80 / ×0.72** | 129 / 129 |
| SC orden propio contra caja | ×0.81 / ×0.81 / ×0.88 | ×0.93 / ×0.90 / ×0.96 | 361 / 356 | ×0.75 / ×0.79 / ×0.73 | 131 / 129 |
| umbral 20 contra 40 (orden propio) | ×0.98 / ×0.99 / ×0.98 | ×1.00 | 362 / 361 | ×0.99 | 131 = 131 |

Por familia (lejos contra caja): Brown ×0.29, Geneig-icse ×0.41, Bellido
×0.70, Bratu ×0.80; pierden BroydenTri ×1.11 y Katsurabis-14 ×1.55. El
umbral 20 gana ~2 %: queda en 40 (no compensa remedir). Ablación contra la
omisión nueva (`results/ablacion_dual_puro_2026-10-10.csv`), contra SC
completo (dual puro, umbral 40, fusión):

| brazo | solver cpu > 0.1 / > 1 / > 5 s | terminan (361) | `ibexopt` cpu > 0.1 s | resueltas (131) |
|---|---|---|---|---|
| arranque tibio primal | ×1.04 / ×1.04 / ×1.02 | 360 | ×0.98 | 129 |
| sin paso largo | ×1.08 / ×1.08 / ×1.05 | 360 | ×1.00 | 130 |
| sin devex | ×1.02 / ×1.03 / ×1.02 | 362 | ×1.00 | 131 |
| sin τ | ×1.01 / ×1.02 / ×0.99 | 360 | ×0.99 | 130 |
| sólo tableau | ×1.12 / ×1.18 / ×1.20 | 353 | ×1.01 | 131 |
| sólo revisado | ×0.98 / ×0.99 / ×0.96 | 363 | ×1.05 | 129 |
| sin fusión | ×1.46 / ×1.47 / ×1.35 | 354 | ×1.01 | 129 |

Celdas del solver iguales (±0.5 %). Con el dual, el paso largo pasa a
importar (era neutro con la omisión primal). chembis y ship-1 están al borde
del límite: se pierden con varios brazos.

**Costo por LP sobre los mismos LPs (2026-10-10, sonda `DFBH_ABCAJA` con
tiempo y pivotes por motor; una semilla, instancias con ≥ 100 LPs).** SC
contra caja: solver µs/LP ×0.56 (220 sistemas, mediana 27 contra 53 µs),
piv/LP ×0.74, µs/piv ×0.78, más barato en el 84 %; `ibexopt` µs/LP ×0.38
(94 problemas, mediana 11 contra 30 µs), piv/LP ×0.80, µs/piv ×0.50, más
barato en el 95 %. Peores: Virasoro-prod ×9 (87 contra 18 piv/LP con el
dual puro), Katsura-30/40 ×2.5 (67 contra 17). Figura en el paper
(`paper/costo_lp_*.dat`).

**Virasoro: degeneración dual del dual puro (2026-10-10).** Virasoro-prod
hace 84–90 piv/LP con el dual puro contra 18 de SoPlex y 22 del primal; el
95 % de los pivotes duales son degenerados. Es el empate ya descrito en
`perturbar_costos`: con objetivo ±e_k y x_k básica en la fila r, tras
reubicar `d = −ā_r`, nulo en toda columna que no toca la fila r, y el ratio
test dual las hace entrar con paso cero. La perturbación de costos
(`DFB_SX_PERTURB=eps`), medida antes sólo en ex5_3_2, lo arregla (20 s):
Virasoro-prod 107 → 323 celdas, 84 → 21 piv/LP, degenerados 95 → 29 %;
Katsura-40 125 → 163 celdas, 47 → 24 piv/LP; Bellido 4.9 → 4.4 piv/LP. ε
entre 1e-7 y 1e-5 da lo mismo. Banco en curso (`re`: dfb, pert, pertl, caja).

Banco `re` (completo): perturbación contra dual sin perturbar, solver
×0.985 / ×0.98 / ×0.96, 360 / 361; `ibexopt` ×1.02 / ×1.02 / ×0.99, 130 /
131 (ship-1). Contra la caja, orden lejos: ×0.81 / ×0.80 / ×0.85, 361 / 356,
`ibexopt` ×0.81 / ×0.83 / ×0.72; orden propio: ninguna instancia pierde más
de ×1.1. En las que llegan al límite, celdas en 60 s ×2 (Virasoro, Katsura).

Barridos sobre las 26 instancias donde la perturbación cambia celdas/s más de
15 % (4 semillas, celdas/s contra el dual sin perturbar): **ε** 1e-9…1e-5
×1.19–1.21, 1e-4 ×1.15 (plano; queda 1e-7). **Cuándo perturbar**
(`DFB_SX_PERTURB_TRAS=K`, sólo tras K pivotes duales degenerados seguidos,
una vez por LP): K = 2 / 3 / 5 / 10 ×1.18 / 1.19 / 1.21 / 1.13 contra ×1.21
perturbando siempre. K = 2 gana más en Katsura (×2.4 contra ×2.15 en
Katsura-30) pero pierde en las CountercurrentReactors (ruido por semilla) y no
cambia las pérdidas sistemáticas: BroydenBanded-200 ×0.73 y Brown-30
×0.5–0.8 con cualquier K y cualquier ε (efecto de trayectoria: en Brown-30
basta perturbar 136 LPs). Conclusión: perturbar siempre, ε = 1e-7.
**Adoptada como omisión** (commit ab0c6e3e; `DFB_SX_PERTURB=0` la apaga).
Paper actualizado con el banco `re` (Tabla 2, abstract, §4 «Dual degeneracy»,
nota en la proposición anytime: el objetivo perturbado difiere de la cota
certificada en ≤ 2ε Σ mag z_j) y figura de costo por pasada (tres paneles,
sonda con perturbación). Pendiente: ablación contra la omisión nueva (`abl4`)
y sonda misma caja (Tabla 3).

Ablación contra la omisión con perturbación
(`results/ablacion_perturbacion_2026-10-11.csv`; controlado en
`results/controlado_perturbacion_2026-10-11.csv`): sin perturbación ×1.03 /
×1.04 / ×1.08 (360 / 362), `ibexopt` ×0.97 y 131 / 129; primal ×1.07–1.10,
`ibexopt` ×0.95; sin paso largo ×1.10–1.11; sin devex ×1.01–1.04; sin τ
×1.00; sólo tableau ×1.13–1.25 (356); sólo revisado ×0.98–0.99 (364),
`ibexopt` ×1.06; sin fusión ×1.40–1.45 (354). Celdas ±2 % (hasta 5 % en > 5 s).
Sonda misma caja con la omisión final: LPs sin óptimo SC / caja Geneig-icse
70 / 499, Brown-18 2 / 176; nodos donde contrae más Geneig-icse 277 / 0,
CCR2-9 126 / 4, Bellido 141 / 3. Paper al día (Tablas 2–4, figura).

**Comparación contra producción fuerte (2026-10-08).** DFB con la
estrategia por nodo (`DFB_SX_PIVOTEO=por_motor`: dual en los nodos del motor
revisado, primal en los de tableau) contra la producción más fuerte de cada
programa: en el solver, SoPlex con representación por columnas
(`IBEX_LP_REP=col`); en `ibexopt`, la de fábrica (allí las columnas son
×1.03 de cpu). Todo sin reloj de SoPlex. En `ibexopt` se excluyen `alkyl`,
`alkylbis` y `bearing` (caóticas por semilla; resueltas 4/4 iguales con y
sin ellas). `results/produccion_fuerte_2026-10-08.csv`:

| programa | > 0.1 s | > 1 s | > 5 s | terminan |
|---|---|---|---|---|
| solver, contra producción con columnas | celdas ×0.77, **cpu ×0.92** | ×0.73, **×0.95** | ×0.70, **×0.89** | **354 contra 323** |
| solver, contra producción de fábrica (como upstream) | ×0.53, ×0.76 | ×0.48, ×0.74 | ×0.46, ×0.71 | 354 contra 320 |
| `ibexopt`, contra producción de fábrica | celdas ×0.99, **cpu ×0.82** | ×1.03, **×0.86** | ×1.02, **×0.79** | 131 contra 131 |
| `ibexopt`, contra producción con columnas | ×1.02, ×0.79 | ×1.07, ×0.83 | ×1.03, ×0.77 | 131 contra 131 |

La estrategia por nodo contra la omisión dual: solver cpu ×0.97–0.98 con
celdas iguales (Bellido ×0.77, `Katsura-16…18` ×0.78–0.83, ninguna pierde);
`ibexopt` neutro. Umbral del motor híbrido con la estrategia por nodo
(`results/umbral_por_motor_2semillas.csv`, solver), contra 40 filas: 60 →
cpu ×0.98; 80 → ×0.97–0.98; 120 → igual que 80 (celdas iguales en todos).
Con el primal el tableau es más barato, así que conviene mandarle más nodos,
pero la ganancia se aplana desde 60–80.

**Línea de base vigente (2026-10-06, tarde): SoPlex de fábrica y sin reloj.**
SoPlex 4.0.2 leía `times()`, una llamada al sistema, en cada FTRAN/BTRAN de su
LU y en cada `optimize`, aun con `TIMER_OFF`, porque el parámetro no llegaba a
la LU ni a las estadísticas. Se parcheó el SoPlex incluido
(`lp_lib_wrapper/soplex/3rd/soplex-4.0.2.timer.patch`, aplicado por CMake) y el
`LPSolver` pone `TIMER_OFF` (`IBEX_LP_RELOJ=1` lo restaura). 138 instancias, 4
semillas, 15 procesos (`results/soplex_sinreloj_4semillas.csv`):

| comparación | > 0.1 s | > 1 s | > 5 s | resueltas 4/4 |
|---|---|---|---|---|
| **DFB contra producción** | celdas ×1.040, **cpu ×0.887**, 161/22 | ×1.016, **×0.866**, 56/13 | ×1.017, **×0.808**, 29/5 | 131 contra 131 |
| producción, sin reloj contra con reloj | celdas ×1.000, cpu ×0.687, 232/0 | ×1.000, ×0.677 | ×1.000, ×0.694 | |
| DFB, sin reloj contra con reloj | celdas ×1.000, cpu ×0.942, 132/0 | ×1.000, ×0.940 | ×1.000, ×0.951 | |

**La mayor parte de la ventaja de DFB en cpu era el reloj de SoPlex en
producción**: producción gastaba ~30 % de su tiempo leyendo el reloj y DFB,
que sólo usa SoPlex para el `loup`, ~6 %. Contra una producción sin ese
impuesto DFB queda en ×0.81–0.89, con las mismas celdas. Todas las razones
DFB/producción anteriores en `ibexopt` y en el solver están infladas por esto;
las comparaciones dentro de DFB (DFB contra DFB) no. 0 óptimos incompatibles;
DFB resuelve además `ship-1`, producción `bearing`. Lo que más pierde DFB:
`bearing` (celdas ×2.2), `ex7_2_3` (×1.46), los `chem` (×1.3); lo que más
gana: `alkyl` (cpu ×0.57), `ex14_2_4` (×0.78).

**Línea de base intermedia (2026-10-06): SoPlex de fábrica, con reloj.** Hasta aquí el
`LPSolver` de Ibex forzaba el simplex primal por columnas (puesto para comparar
con las primeras versiones de DFB); ahora usa la configuración de fábrica
(`IBEX_LP_PRIMAL=1` restaura el primal) y sin el reloj de SoPlex en lo que el
parámetro alcanza. 139 instancias, 4 semillas, 15 procesos, cuatro brazos en el
mismo barrido (`results/soplex_fabrica_4semillas.csv`), `analiza_lu.py`:

| comparación | > 0.1 s | > 1 s | > 5 s | resueltas 4/4 |
|---|---|---|---|---|
| DFB contra producción, los dos de fábrica | celdas ×1.038, **cpu ×0.646**, 214/11 | ×1.016, **×0.612**, 84/4 | ×1.027, **×0.608**, 45/1 | 131 contra 130 |
| DFB contra producción, los dos con el primal (reproduce la de abajo) | ×0.983, ×0.623, 215/9 | ×0.978, ×0.623 | ×0.982, ×0.609 | 131 contra 131 |
| producción, fábrica contra primal | ×1.056, ×1.040 | ×1.015, ×1.016 | ×0.930, ×0.958 | 130 contra 131 |
| DFB, fábrica contra primal | ×1.159, ×1.122 | ×1.098, ×1.057 | ×1.052, ×1.035 | 131 contra 131 |

0 óptimos incompatibles en todos los pares. DFB resuelve además `ex8_4_4bis`
y `ship-1`; producción, `bearing`. Los dos brazos con el primal reproducen el
barrido del 2026-10-02 celda por celda en 525 de 552 pares (el resto, al
tope). El cambio a fábrica pesa más en DFB porque sólo toca su buscador de
`loup`, y ahí está el daño: `alkyl`/`alkylbis` explotan en algunas semillas
(62 → 4 542–32 416 celdas, en los dos brazos) y `launch` pasa de ~64 a
~200–300 celdas en los dos. Es el LP de `LoupFinderXTaylor`, no la
contracción: con el primal el vértice que devuelve SoPlex servía mejor de
candidato.

**Línea de base anterior (2026-10-02, con el primal forzado)**, con `eps_x` = 1e-7 de verdad (§8) y la
omisión de DFB de entonces, incluido el arreglo de devex (§9). 138 instancias, 4
semillas, 8 cores, los brazos en el mismo barrido
(`results/devex_arreglo_4semillas.csv`):

| corte | celdas | cpu | +rápido / −lento |
|---|---|---|---|
| > 0.1 s (n=245) | ×0.985 | **×0.601** | **221 / 10** |
| > 1 s (n=91) | ×0.977 | **×0.591** | 83 / 7 |
| > 5 s (n=45) | ×0.982 | **×0.584** | **44 / 1** |

cpu por celda ×0.583 (> 1 s), ×0.511 en las 21 que llegan al tope en los dos.
0 óptimos incompatibles. Resueltas en 4/4: 131 y 131; DFB resuelve además
`ex8_4_4bis`, producción `bearing`, y `ship-1` se decide por semilla. La línea
anterior (2026-09-30, sin el arreglo) daba cpu ×0.617 / ×0.591 / ×0.575. Bajar `eps_x` de 1e-6 a 1e-7, por
sí solo, ayuda a los dos: producción celdas ×0.99 / ×0.96 / ×0.93, DFB ×0.995 /
×0.978 / ×0.960.

> **PROBLEMA ABIERTO — `bearing` con `eps_x` = 1e-7.** La cota inferior de DFB
> es buena; falta un `loup` que cierre. `bearing` tiene 8 igualdades, relajadas a
> `|h| ≤ 1e-8`, y el buscador de `loup` sólo acierta en cajas de 1e-4 a 1e-3.
> Tras el `loup` de DFB toda la búsqueda está en la zona del óptimo y el buscador
> falla 18 633 veces en cajas < 1e-6, que la contracción empuja contra el borde
> de la holgura (tres igualdades en ±1e-8). Con `eps_x` = 1e-6 el `loup` bueno
> salía de una caja de 5.95e-7. También le cuesta al primal con este `eps_x` (una
> semilla al tope): es caos del `loup`, no sólo de DFB. Análisis completo en
> `MEDICIONES_SESION_2026-09-30.md` §A. Verificación pendiente: llamar al
> buscador antes y después de contraer cada caja.

**El linearizador afín es el default** (`ibex-affine` vendorizado en `affine/`;
`DFBH_LR=xn` restaura xtaylor). Cambiar sólo el linearizador vale más que todo
lo hecho sobre el contractor: con `CtcPolytopeHull` intacto da cpu ×0.294 en
las de más de 1 s. La relajación afín es más apretada **y más chica**: 34 → 23 filas,
pivotes por LP 10.9 → 6.2, costo por LP 168 → 56 µs.

Evolución de DFB contra producción con `eps_x` = 1e-6 (cpu > 0.1 / > 1 / > 5 s):
LU actualizada ×0.684 / ×0.657 / ×0.621; tableau con guardia ×0.573 / ×0.558 /
×0.536. Tablas en `MEDICIONES_SESION_2026-09-30.md` §B.

| | celdas | cpu |
|---|---|---|
| PolyHull + compo / PolyHull + afín | ×0.858 | ×1.110 |
| DFB + compo / PolyHull + afín | ×0.831 | ×0.922 |
| producción sin fix / con fix | ×1.389 | ×1.008 |
| DFB sin fix / con fix | ×1.380 | ×1.077 |
| DFB / prod, ambos sin fix | ×0.979 | ×0.780 |

El punto fijo conviene a los dos y a DFB más: la ventaja no es un artefacto de
una composición sintonizada para `PolytopeHull`.

---

## 5. El régimen de DFB

La ventaja se sostiene con relajaciones de 23 y 34 filas y **desaparece con
54–81**. Sobre los mismos LPs, con `m` creciendo ×2.53: µs por iteración ×1.25
en SoPlex contra ×2.36 en DFB; iteraciones por LP ×1.30 contra ×2.52. El costo
por pivote crece con `m` porque `pivot` actualiza las `m` filas del tableau
explícito, y la ventaja del arranque tibio cae de 2.60× a 1.34×. **La
certificación no es el cuello de botella**: 2–3 % del tiempo y plana en `m`.

La ventaja no viene de pivotear menos (7.4 contra 3.1 en `ex7_3_4`) sino de que
el ciclo completo del LP, pivoteo más certificación, es 3–4 veces más barato que
SoPlex más Neumaier–Shcherbina: 16 contra 59 µs/LP en `ex7_2_3`, 29 contra 116
en `ex14_2_7`. Y DFB emite 17–26 % menos LPs por celda.

---

## 6. Lo que NO funciona — no reintentar sin leer esto

Todas medidas con el banco completo.

| idea | resultado |
|---|---|
| HC4 intercalado en la pasada, 4 disparadores | cpu ×1.020–1.034, celdas planas |
| Cotas colaterales (`0 ∉ γ_j·z` para j≠k) | no contraen ni ahorran pivotes |
| Rayos alternos de `B⁻¹` (`DFBH_OTRORAYO=1`) | recuperan el 100 % de los certificados perdidos y **cero** efecto en celdas |
| Cosecha de filas de `B⁻¹` (`DFBH_FILAS=1`) | contrae de verdad y baja LPs/nodo ×0.835, pero el `γ` por básica sube el costo por LP ×1.254: cpu ×1.06–1.09 |
| Cachear esos `γ` (`DFBH_FILASCADA=k`) | un certificado viejo sigue **válido** pero no **fuerte** |
| Saltear el LP si la cota de cero pivotes basta (`DFBH_LIBRE=1`) | iguala al óptimo justo cuando el LP costaba 0 o 1 pivote |
| Generación perezosa de filas | 34–56 % peor |
| Reuso de linealización (`DFBH_RELIN=1`) | catastrófico: 472 → 22 512 celdas |
| Congelar `A` por nodo y actualizar sólo las cotas de `b` (`DFBH_CONG=1`) | celdas ×1.20 en > 1 s, cpu ×1.14; nodo a nodo, contracción ×0.85 y pruebas de vacío ×0.66 |
| Re-linealizar sólo lo que se movió (`DFBH_CONGVAR=τ`) | parecía cpu ×0.70 y era **incorrecto** (38–53 óptimos incompatibles); corregido, cpu ×1.15 |
| Afinar el umbral del punto fijo (`DFBH_RATIO`) | compra poda y la paga cara |
| Ciclo interno propio con re-linealización | ×0.97 en las pesadas, pierde 4–5 % en el banco |
| Apretar la parte `b` antes de certificar | neutro |
| Perturbar los costos nulos (`DFB_SX_PERTURB`) | 23.0 contra 23.7 pivotes/LP en `ex5_3_2` |
| Ratio test de paso largo (`DFB_SX_BFRT=1`) | `alkylbis` de 62 a 19 982 celdas. **Era un defecto, arreglado el 2026-09-30 (§9)**: con el dual, cpu ×0.972 contra el dual solo, y el dual con paso largo queda a ×1.04 del primal (×0.98 en > 5 s) |
| Regla de la primera cota del nodo (`DFBH_PRIMERA`) | neutro en el banco, ×0.99–1.00 |
| Devex en la fila que sale | mixto con la estrategia que en los hechos era primal; **con el dual puro, cpu ×0.98–0.99 y por omisión** (§9) |
| Presupuestos / truncamiento anytime, 11 variantes sobre el camino primal | ninguna paga |
| Corte anytime por brecha sobre el pivoteo dual (`DFBH_ANYTAU=τ`) | ahorra pivotes (×0.84) pero celdas ×1.39–1.60 y cpu ×1.36–1.49 |
| Herencia de base entre nodos (`DFB_SX_BASENODOS=1`) | con el dual destapó el defecto de devex (§9); arreglado, baja pivotes ×0.78–0.96 pero reconstruir el tableau cuesta lo que ahorra: tiempo ×1.00–1.10 |
| Test de `γ` en las rodajas de ACID (`DFBH_ENACID=1`), re-medido tras el filtro de estabilidad | con `γ` propios acierta ≈ 0 %; con los de un ancestro mata rodajas que el lineal del nodo mata igual (§9) |
| Cortes en `x*` (`DFBH_CORTES=1\|2\|3`, tope 0.25, §9), solver | en el nodo celdas ×0.82–0.84, cpu ×1.18–1.24 (vueltas extra del punto fijo y más pivotes); heredados ×0.80, ×1.26–1.30; sólo en las hijas ×0.84–0.88, ×1.25–1.31 y ninguna instancia gana |
| Cota lagrangiana con monotonía (`DFBH_LAGRC=1\|2`, §9) | mejora la cota del LP en el 70–97 % de las aplicaciones y mata nodos que el LP no mata; en el objetivo solo, celdas iguales y cpu ×1.08–1.24; en las `2n` cotas, celdas ×0.81–0.96 y cpu ×1.6–3.5 |

**El patrón.** Adelantar trabajo dentro del nodo no mueve el árbol. Lo único
que movió la aguja fue cambiar **qué** se calcula: la relajación afín, el filtro
de estabilidad, la regla de selección. Corolario: lo que se ahorra son los LPs
baratos, porque la calidad de la cota intermedia y el costo del LP están
correlacionados.

**Qué aporta el punto fijo.** Tres brazos con pivoteo dual, celdas en > 1 s:
re-linealizando 1.00, con `A` congelada y sólo cotas actualizadas 1.20, sin
punto fijo 1.53. El 62 % del beneficio se consigue sin re-linealizar. El alcance
del congelamiento es el nodo, detectado por el umbral del `CtcFixPoint`
(`rel_distance > 0.2`); con un contador fijo el mismo brazo medía ×1.4.

**Sobre el anytime.** Con la garantía puesta (pivoteo dual, corte por brecha
contra el techo del hull) falla igual: dispara en el 24–69 % de los LPs y todo el
daño está en las celdas. En esta búsqueda la contracción vale mucho más que los
pivotes.

**Sobre el vacío.** Es propiedad del nodo: si una cota llega a `OPTIMAL`, ningún
`λ` posterior puede probar vacío. DFB mata entre el 1 % y el 17 % de los nodos
que ve; quitarle el test multiplica el árbol por 2.5 a 10.

---

## 7. Advertencias de método

- **Una semilla no es una medición.** `LinearizerXTaylor` usa esquinas
  aleatorias (`--seed`, omisión 42); `ship-1` oscila ×5.2 entre semillas y
  reescribir el mismo poliedro mueve `ex14_2_7` ±40 %. Cuatro semillas, y
  mirar el loup final antes de explicar un efecto grande.
- **Los totales del banco engañan**: una instancia de 28 000 celdas pesa como 300
  de 90. Media geométrica de razones, o ventaja acotada `(A−B)/max`.
- **Controles chicos no filtran**: el HC4 intercalado daba ×0.57 en el control y
  ×1.03 en el banco. Incluir las pesadas antes de concluir; sobre 6 fáciles el
  congelamiento empataba y con 5 pesadas perdía ×1.143.
- **La caja raíz no representa los nodos profundos**, ni las chicas a las grandes.
- **Reusar información entre cajas exige verificar contención.** Costó tres
  veces; una caja puede encoger y a la vez desplazarse.
- **El perfilador propio distorsiona**: `clock_gettime` por fase infla 22–37 % y
  crea baldes falsos. Preferir `perf` por muestreo.
- **Un control que no ejercita el término en duda no verifica nada.** Para el
  acotamiento afín de `c` hubo que separar identidad e implicación y muestrear
  donde el término es no nulo (`DFBH_CHKC=k`: 0 violaciones en 3.4 M de
  muestras). Al muestrear, no evaluar `A.x − a'.x` como resta de intervalos
  correlacionados: acumular término a término en precisión extendida.
- **Para diagnosticar**, lo contrario que para evaluar: la instancia más chica
  con mayor contraste y sin sensibilidad a la semilla, y después un nodo. Correr
  las dos configuraciones **sobre la misma caja** (`DFBH_ABDUAL`, `DFBH_ABCONG`),
  porque dos corridas completas divergen en el primer nodo distinto.

---

## 8. Estado del código

```
ibex_CtcDFBHull.{h,cpp}     el contractor
ibex_DFBSimplex.{h,cpp}     el simplex: tableau explícito disperso
ibex_DFBBasis.{h,cpp}       la factorización LU (SLUFactor de SoPlex)
ibex_Optimizer05Config.cpp  el cableado de --lr y DFBH_LR
affine/                     el plugin ibex-affine vendorizado
```

Archivadas pero **no borrar**: `ibex_CtcDualFeasibleBounding.cpp`,
`ibex_CtcDFBPropag.cpp`, modo `--lr=dfb`. Arneses obsoletos en `obsoleto/`.
Barridos con 4 semillas en `results/*_4semillas.csv` con su `analiza_*.py`.

### Interruptores vivos

| variable | |
|---|---|
| `DFBH_LR=xn\|art\|compo` | linearizador. Omisión `art` |
| `DFBH_ORDEN=costo\|lejos\|cerca\|nat` | orden de la próxima cota. Omisión `costo` |
| `DFB_SX_PIVOTEO=primal\|dual\|mixta`, `DFB_SX_ENTRAK=1`, `DFB_SX_MONO=1` | re-optimización al cambiar de objetivo: `dual` es el dual puro (reubica también `k`), `mixta` la anterior; entrada explícita de `k`; sonda de monotonía (§9). **Omisión: `dual` en `ibexsolve_dfb`; en `ibexopt` primal + revisado > 15 filas (2026-10-13)** |
| `DFB_SX_LU=actualizada\|diferida\|tableau` | la LU se actualiza en cada pivote, se refactoriza al leerla, o no se usa salvo que la guardia de deriva falle (§9). Omisión `tableau`; `DFB_SX_GUARDIA=tol` |
| `DFBH_HULL=1`, `DFBH_HULLTAU=τ` | selección por hull (§3) |
| `DFBH_NOFIX=1` | sin `CtcFixPoint` en el contractor lineal |
| `DFB_SX_REUBK=1`, `DFB_SX_NOXBSKIP=1`, `DFBH_LINSLOW=1`, `DFBH_AFNOPOOL=1` | restauran comportamientos previos, para reproducir mediciones |
| `DFBH_PRIMERA`, `DFB_SX_BFRT`, `DFB_SX_PERTURB`, `DFBH_FILAS`, `DFBH_FILASCADA`, `DFBH_LIBRE`, `DFBH_ANYTAU`, `DFBH_CONG`, `DFBH_CONGVAR` | variantes medidas en el §6 |
| `DFBH_PERFIL`, `DFBH_LINPERF`, `DFBH_AFALLOC`, `DFBH_CRONO`, `DFB_SX_PERF`, `IBEX_PH_STATS` | perfilado propio |
| `IBEX_TRACE=1`, `IBEX_LOUPSTATS=1`, `IBEX_SMEAR_STATS=1`, `IBEX_ARGS=1` | traza de `loup`/`uplo` con celda, tiempo y ancho; intentos y éxitos del buscador de `loup` por ancho; respaldos del bisector; parámetros recibidos (§4, §8) |
| `DFB_SX_VERIF=1`, `DFB_SX_FEASTOL=t`, `DFB_SX_TOLANCHO=τ` | verifica cada `OPTIMAL` desde la factorización; tolerancia primal, omisión `1e-9`; tope de la tolerancia en `τ·ancho` de la básica, **omisión 1e-3** (§9) |
| `DFBH_SINCOTA`, `DFBH_PREDICE`, `DFBH_GOALPROBE`, `DFBH_CHKC=k`, `DFB_SX_DERIVA=1`, `DFBH_LAGR=1` | sondas (las dos últimas: deriva del tableau contra la LU; cota lagrangiana contra la del LP, §9) |
| `ibexsolve_dfb archivo.bch [--lr=xn\|dfbhull] [--lin=art\|xt\|compo] [-b lsmear\|smearsumrel\|lsmeardfb] [-e eps] [-t s] [--seed=s] [--bfs] [--sols]` | `ibexsolve` con DFB en el hueco de `CtcPolytopeHull` linearizador y bisector a elección, **afín y `lsmear` por omisión**; `--lr=xn --lin=xt -b smearsumrel` reproduce al oficial (§9) |
| `DFB_SX_RELOJ=1`, `IBEX_LP_RELOJ=1` | restauran el reloj de SoPlex en la LU de DFB y en el `LPSolver` de Ibex (§9) |
| `IBEX_LP_PRIMAL=1` | el `LPSolver` vuelve al simplex primal por columnas (antes del 2026-10-06; §4) |
| `IBEX_LP_REP=col` | SoPlex con representación por columnas (fábrica elige filas en muchos LPs del solver; §9) |
| `DFBH_CORTES=1\|2\|3`, `DFBH_CORTES_TOPE=f` | cortes en `x*`: en el nodo, heredados, sólo en las hijas; tope del pozo en veces las filas afines (§9, cerrado) |
| `DFBH_SALTEO=ancho`, `IBEX_PH_SALTEO=ancho` | salto de Achterberg relativo al ancho de la caja en DFB / producción (§9) |
| `IBEX_PH_NOACH=1`, `DFBH_AB=1`, `DFBH_ABHULL=1` | sin salto en producción; misma caja con los dos contractores; cotas de DFB contra el hull exacto (§9) |
| `DFBH_FUSIONA=0` | no fusiona las filas gemelas (misma normal, un lado cada una) de la linealización fresca (por omisión las fusiona desde 2026-10-09) (§4) |
| `DFB_SX_PERTURB=eps`, `DFB_SX_PERTURB_TRAS=K` | perturbación de costos al reubicar, **omisión ε = 1e-7** desde 2026-10-11 (0 la apaga); con K > 0, sólo tras K pivotes duales degenerados seguidos (§4) |
| `DFB_SX_DEVEXP=1` | precio devex en los pasos primales del motor revisado (por omisión Dantzig; sin banco aún) (§4) |
| `DFB_SX_REVPRIMAL=1` | pasos primales en el motor revisado (omisión sólo el 2026-10-10 antes de la fusión; desde 2026-10-13 encendidos por omisión en `ibexopt`) (§4) |
| `DFB_SX_PIVOTEO=por_motor\|dual\|primal\|mixta` | **`ibexsolve_dfb`: omisión `dual` (puro) desde 2026-10-10. `ibexopt`: omisión primal con `REVPRIMAL=1` y umbral 15 desde 2026-10-13, fijada en `main` si `DFB_SX_PIVOTEO` no está definida (`DFB_SX_PIVOTEO=dual` vuelve al dual puro, umbral 30)**; `por_motor` (dual en revisados, primal en tableau) fue omisión 2026-10-08…10 (§1, §4, §9) |
| `DFBH_MOTOR=caja`, `DFBH_CAJA_TOL=t` | mismo contractor con SoPlex como caja negra y el mismo certificado gamma; tolerancia de SoPlex (omisión 1e-6) (§4) |
| `DFBH_CUENTAS=1` | `ibexsolve_dfb`: celdas, llamadas a DFB, LPs y pivotes al final |
| `IBEX_PH_STATS=1` | `CtcPolytopeHull`: estados de SoPlex por LP (óptimo, vacío probado / sin prueba, desconocido) |
| `DFBH_TESTIGO=archivo`, `DFBH_INFSTATS=1` | verifica un punto solución conocido contra la relajación y las cotas; cuenta filas donde el ±1e100 de SoPlex habría cortado (§9) |
| `DFBH_LATIDO=1`, `DFBH_LPLARGO=1` | latido del contractor cada 200 llamadas; LPs de más de 5000 vueltas y progreso de la pasada |
| `DFB_SX_MOTOR=tableau\|revisado\|auto`, `DFB_SX_MOTOR_M=m` | motor del simplex: tableau explícito, revisado sin tableau, o revisado por encima de `m` filas (**umbral 30** desde 2026-10-13; 40 entre el 10 y el 13, 80 entre el 08 y el 10) (§9). **Omisión `auto`** |
| `DFBH_LUSTATS=1` | cuenta factorizaciones, updates, BTRAN y FTRAN de la LU (§9) |
| `-b lsmeardfb`, `DFBH_SENS=redcost\|abs\|jac`, `DFBH_BISEC_CMP=1`, `IBEX_LSMEAR_DUMP=1` | bisector `LSmear` con sensibilidades de DFB guardadas en la celda; comparación con `LSmear` sobre el mismo nodo; volcado de puntajes (§9). Fuera de la omisión |
| `DFB_SX_DEVEXVIEJO=1`, `DFB_SX_DSE=1`, `DFB_SX_DSEABS=1`, `DFB_SX_HARRISABS=1` | restaura el defecto de devex (§9); dual steepest edge exacto y su variante con violación absoluta; Harris de dos pasadas. Apagados |
| `DFB_SX_BFRT=0`, `DFB_SX_BFRTINF=1` | ratio test de paso largo, **encendido por omisión** (`=0` lo apaga); la segunda restaura el defecto de la infactibilidad espuria (§9) |
| `DFBH_ADAPT=u`, `DFBH_ADAPT_TRAZA=1`, `DFBH_BASECOMP=1` | arranque tibio adaptativo por cambios de lado (apagado, no gana); sonda de base y lados entre cotas consecutivas (§4) |
| `DFBH_LAGRC=1\|2`, `DFBH_LAGRR=r` | cota lagrangiana con monotonía como contracción, en el objetivo o en las `2n`; rondas de monotonía (omisión 2). Apagada (§9) |
| `DFBH_ABDUAL=τ`, `DFBH_ABNODO=n`, `DFBH_ABPRIM=R`, `DFBH_ABCONG=1`, `DFBH_ABPURO=1` | comparan dos configuraciones sobre el mismo nodo (la última: primal contra dual puro) |

Queda una veintena de interruptores de líneas descartadas (`DFBH_HC4*`,
`DFBH_COLAT`, `DFBH_LAZY`, `DFBH_ADAPT`, `DFBH_TECHO`); varias sondas envuelven
código que hay que conservar, así que no sirve un borrado por rangos.

**`LargestFirst` con precisiones nulas y `eps_x` real (2026-09-30).** En otra
rama, el respaldo de `LSmear` (`OptimLargestFirst`) quedaba viciado con
`eps_x = 0`: `diam/prec` daba +∞ para todas las variables y se bisecaba siempre
la primera. Arreglado en `ibex_LargestFirst.cpp` e `ibex_OptimLargestFirst.cpp`:
primero las variables con precisión 0, por diámetro; si ninguna es bisecable,
las demás por `diam/prec` (`IBEX_LF_VIEJO=1` restaura el criterio anterior). En
esta rama no cambia nada (celdas idénticas en las 18 instancias que usan el
respaldo, ≈10 % de las bisecciones, siempre por jacobiana infinita; `ship-1` el
99 %; contador `IBEX_SMEAR_STATS=1`), porque `eps_x` no llegaba como 0 sino como
**1e-6**: `ibexopt.cpp` lo pasaba con `std::to_string`, y con el redondeo hacia
arriba que deja gaol `1e-7` se formateaba `"0.000001"`. **Arreglado**: los
`double` se pasan con `%.17g` y `eps_x` = 1e-7 queda menor que `eps` = 1e-6
(`IBEX_ARGS=1` imprime lo que llega); cambió la línea de base (§4). Defecto menor
sin tocar: en `OptimLargestFirst` la condición para bisecar el objetivo compara
`l = diam/prec` con el diámetro crudo del objetivo, así que el respaldo casi
nunca lo elige.

Arreglos en el núcleo: `CtcPolytopeHull` no inicializaba `n_soplex_iterations`
ni `n_soplex_calls`; instrumentación de costo por LP en `LPSolver::minimize`
tras `IBEX_PH_STATS` (los agregados sí son de fiar, el reparto por rubros no).

---

## 9. Lo que queda

### Dónde está el pivoteo

El dual puro con paso largo, tolerancia a escala de la caja y devex es la
omisión y está a la par del primal (banco con `eps_x` = 1e-7,
`results/dse_epsx1e-7_4semillas.csv`):

| dual (omisión) contra primal | cpu > 0.1 / > 1 / > 5 s | celdas | resueltas 4/4 |
|---|---|---|---|
| | **×1.014 / ×1.002 / ×0.959** | ×0.999 / ×0.997 / ×0.969 | 131 / 131 |

Cuatro defectos que hubo que arreglar para llegar acá, con su mecanismo en
`MEDICIONES_SESION_2026-09-30.md`:

1. **El «dual» anterior no era dual**: excluía `k` de la reubicación y el
   38–83 % de los LPs terminaba en primal. Hoy es `mixta`.
2. **Paso largo**: declaraba `INFEASIBLE` cuando los cambios de cota no cubrían
   la violación, con 100 % de LPs sin cota en `ex7_2_6`. Ahora pivotea en el
   último quiebre (`DFB_SX_BFRTINF=1` restaura el defecto).
3. **Tolerancia primal**: `1e-9·(1+|cota|)` era ~13 % del ancho de las cajas
   angostas, y el dual puro declaraba factibles poliedros que el primal probaba
   vacíos. Ahora se acota por `τ·ancho`, τ = 1e-3.
4. **`eps_x`** llegaba como 1e-6 (§8).

Lo que no sumó: entrada explícita de `k` (`DFB_SX_ENTRAK`, más pivotes; su
mejora en `bearing` era redondeo), Harris absoluto (`DFB_SX_HARRISABS`, elimina
casi todos los pasos primales pero cuesta pivotes), DSE exacto (`DFB_SX_DSE`,
pesos `‖e_rᵀB⁻¹‖²` desde la parte λ del tableau; mixto y ×1.03 por celda).

**Los pasos primales que quedan (0.1–1 % de los LPs)** rompen
momentáneamente la factibilidad dual —y con ella la monotonía en ese LP— pero no
la validez del certificado. La causa principal es la banda de Harris relativa al
cociente; queda otra fuente menor sin identificar en `bearing`.

**Defecto de devex, arreglado (2026-10-01).** El paso dual elegía la fila que
sale con `if (v > worst)` partiendo de `worst = 0`, con `v = viol²/peso`. Los
pesos de devex sólo crecen; cuando uno se disparaba, `v` daba 0, la fila violada
no se elegía nunca, y si pasaba con todas el paso concluía «factible primal» y
el LP se declaraba óptimo siendo infactible. Apareció al heredar la base entre
nodos (`DFB_SX_BASENODOS=1`), donde los pesos ya no se reinician en el arranque
en frío: **71.7 % de óptimos falsos** en `ex5_3_2` (`DFB_SX_VERIF=1`), 168 → 568
celdas; con el primal la herencia funcionaba (1.7 %), y sin devex también
(0.1 %). Arreglo: `worst = −1`, para que la decisión de factibilidad no dependa
de los pesos, y reinicio del marco de referencia al heredar la base o cuando un
peso pasa de 1e8 (`DFB_SX_DEVEXVIEJO=1` restaura el defecto). En la omisión
(banco, 138 instancias, 4 semillas, `results/devex_arreglo_4semillas.csv`): cpu
×0.960 / ×0.906 / ×0.982, celdas idénticas en 380 de 552 pares, cpu por celda
×1.007, 0 óptimos incompatibles. Buena parte del ×0.906 en > 1 s es
`ex7_2_2` (= `ex8_1_8bis`), cuya semilla 4 pasa de 2348 a 54 celdas: el caso
que antes se atribuyó a caos del `loup` pudo ser este defecto, pero con
`DFB_SX_VERIF` encendido las dos versiones dan 50 celdas y no se puede aislar.

**La herencia de base entre nodos, ya correcta, no paga** (8 instancias, 2
semillas): pivotes por LP ×0.78–0.96 (`ex7_2_3` 1.6 → 1.23, `ex6_1_3` 4.35 →
3.84) pero tiempo ×1.00–1.10 en las estables y peor en `launch` y `himmel16`.
Heredar exige refactorizar con la linealización nueva y reconstruir el tableau
completo (una BTRAN por fila y un producto por `Ā`), del orden de `m` pivotes,
que es lo que costaba la cota 0 en frío. Con un simplex sin tableau explícito
(opción 3) heredar costaría una factorización, y ahí sí podría pagar.

**La brecha que queda** está en las relajaciones grandes: `launch` y `ex5_3_2`
(×1.18–1.20) y `ex14_2_7` (×1.16). El dual reubica y reparte la infactibilidad
por muchas de las ~50 filas; el orden de cotas no es el problema (el predictor
por columnas acierta el costo del dual).

### Los duales de DFB en la bisección: no sirven como están (2026-10-02)

La idea: `LSmear` resuelve en cada nodo un LP con SoPlex sólo para obtener los
multiplicadores del objetivo, y DFB ya los tiene. Se construyó `LSmearDFB`
(`-b lsmeardfb`): la regla de `LSmear` con sensibilidades de DFB, guardadas
**en la celda** (`BxpDFBSens`), porque el optimizador biseca la celda al
sacarla de la cola y recién después contrae a las hijas, así que el último nodo
que pasó por el contractor es otro. Cuando la pasada no resuelve la cota del
objetivo, o termina con `y` no básica (`λ = 0`), se resuelve `min y` con la
cota inferior de `y` suelta, como hace `LSmear` (`LPs extra`, ~0.4 por nodo).
Tres sensibilidades (`DFBH_SENS=redcost|abs|jac`), 8 instancias, 2 semillas,
celdas contra `lsmearmg` (`results/lsmeardfb_sueltas.csv`):

| sensibilidad | coincide con `LSmear` | mejora | empeora |
|---|---|---|---|
| `redcost`: `\|γ_j/γ_goal\|·diam` | 13 % | `ex14_2_7` ×0.5 | `ex2_1_9` ×22, `ex7_2_3` ×12, `himmel16` sin `loup`, `ex6_1_3` y `ex8_4_4bis` al tope |
| `abs`: `Σ_i \|λ_i\|\|Ā_ij\|·diam` | 35 % | `ex2_1_9` ×0.55, `ex7_2_3` ×0.6, `ex14_2_7` ×0.5 | `himmel16` ×25, `launch` ×10, `ex6_1_3` ×2.5, `ex8_4_4bis` ×2.2 |
| `jac`: `\|Σ_c λ_c J_cj(mid)\|·diam`, la jacobiana de `LSmear` | 26 % | `ex14_2_7` | `ex2_1_9` ×22, `ex7_2_3` ×12, `himmel16` ×14 |

Mecanismo. El puntaje de `LSmear` resulta ser `\|Σ_i y_i J_ij\|·diam_j`, con
`y` los duales de las filas de su LP (cotas de las variables como filas aparte):
una *smear* ponderada por los duales, grande donde la variable pesa en las
restricciones activas. El costo reducido de DFB es **cero en toda variable
básica** por definición, aunque bisecarla achique mucho la relajación, así que
`redcost` sólo ordena las no básicas y el árbol explota. Y `jac`, que usa la
misma matriz que `LSmear`, también falla: la diferencia está en los `λ`. Los
del vértice óptimo de la relajación afín son pocos y dispersos (2–4 filas
activas), y concentran el puntaje en pocas variables; los de `LSmear`, de un LP
más flojo con las cotas como filas, reparten el peso entre más restricciones.
Para bisecar sirve un peso ancho, no la sensibilidad exacta de un vértice.

Banco de la variante `abs` contra `lsmearmg` (138 instancias, 4 semillas;
`results/lsmeardfb_4semillas.csv`): celdas ×1.53 / ×1.88 / ×0.62 y cpu ×1.28 /
×1.58 / ×0.73 en los cortes > 0.1 / > 1 / > 5 s; **resueltas en 4/4: 131 →
116**; de las instancias de > 0.1 s, mejora más del 10 % en 13 y empeora en
23, con colas enormes (`ex7_3_4` ×67, `ex6_2_12` ×58, `house` ×18) frente a
mejoras acotadas (`ex14_2_4` ×0.35, `ex2_1_9` ×0.55, `ex7_2_3` ×0.65). El
×0.62 de > 5 s es sólo `bearing` y `ship-1`. Cpu por celda ×0.77: lo que ahorra
es el LP de SoPlex de `LSmear`, que no compensa el árbol.
Cerrada como sustitución directa; el bisector y la propiedad quedan detrás de
`-b lsmeardfb` para quien quiera retomarla con otra regla. La idea del paso
dual de anticipación (puntuar las básicas por `viol·θ_r` del test de razón)
queda sin probar y con la expectativa bajada por esto mismo.

### Las opciones

0. **`bearing` con `eps_x` = 1e-7** (§4): verificar si el buscador de `loup`
   acierta antes de contraer y falla después.
1. **El `x*` de DFB como candidato de `loup`.** Hoy `LoupFinderXTaylor`
   resuelve otro LP con SoPlex por nodo, y el `loup` es lo que falla en
   `bearing`. Techo ~7 % de cpu. (Los multiplicadores para `LSmear` ya se
   probaron: arriba.)
2. **`ibexsolve`: hecho y medido (2026-10-05).** `ibexsolve_dfb` arma la
   composición de `DefaultSolver` (HC4, ACID, Newton si es cuadrado, punto fijo
   del lineal con HC4; mismo bisector y buffer) y sólo cambia el contractor
   lineal (`--lr=xn|dfbhull`) y, aparte, el linearizador (`--lin=xt|art`). Con
   `--lr=xn --lin=xt` reproduce al `ibexsolve` oficial celda por celda.

   **Defecto de corrección encontrado y arreglado.** SoPlex representa el lado
   abierto de una fila como ±1e100 (`soplex::infinity`), y `linearize_fresca`
   lo leía como cota finita de la holgura `b`. Con cajas de 1e8 y productos de
   muchas variables (`Brown-15`), XTaylor da filas de coeficientes ~1e113 y ese
   falso «−1e100» cortaba soluciones reales: DFB terminaba con **éxito** y 2
   soluciones de 3 (infactible con la estrategia primal). Lo destapó una
   discrepancia entre linearizadores; la sonda `DFBH_TESTIGO=archivo` (un punto
   solución conocido, verificado contra cada fila y cada cota) lo localizó. Ahora
   ±1e100 se lee como ±∞. En el banco de `ibexopt` no se activó nunca
   (`DFBH_INFSTATS=1`: 0 filas en 138 instancias), así que las mediciones de
   `ibexopt` siguen valiendo; en el solver sí, y las cifras que siguen son las
   corregidas.

   Banco `benchs/solver` (231 sistemas), 2 semillas, 60 s internos y 150 s
   externos, cuatro brazos (`results/ibexsolve_4brazos_2semillas.csv`). Celdas
   / cpu en los cortes > 0.1 / > 1 / > 5 s, sobre pares que terminan los dos:

   | comparación | > 0.1 s | > 1 s | > 5 s | terminan 2/2 |
   |---|---|---|---|---|
   | DFB/XTaylor contra prod/XTaylor | ×0.52 / ×0.96 | ×0.42 / ×0.83 | ×0.34 / ×0.86 | 142 contra 138 |
   | **DFB/afín contra prod/afín** (efecto del contractor) | ×0.76 / **×0.89** | ×0.71 / **×0.89** | ×0.74 / **×0.87** | 160 contra 148 |
   | prod/afín contra prod/XTaylor (efecto del linearizador) | ×0.86 / ×0.65 | ×0.72 / ×0.57 | ×0.67 / ×0.52 | 148 contra 138 |
   | **DFB/afín contra prod/XTaylor** (los dos cambios) | ×0.64 / **×0.56** | ×0.49 / **×0.50** | ×0.45 / **×0.46** | **160 contra 138** |

   0 discrepancias en cajas solución entre todos los brazos. Con XTaylor, DFB
   gana 13 sistemas y pierde 9, todos grandes (`BroydenBanded-100…140`,
   `Katsura-18`, `Eiger-0100`…): XTaylor da hasta 4 filas por igualdad (dos
   esquinas × dos lados) y el tableau explícito se vuelve caro. El afín da la
   mitad de filas y elimina esas pérdidas: DFB/afín contra producción gana 24 y
   pierde 2 (`Geneig-icse`, `hayes1`). El afín le sirve también a producción
   (×0.52–0.65), pero DFB le suma un ×0.87–0.89 encima y 12 sistemas más.

3. **El rediseño para relajaciones grandes** (BTRAN/FTRAN sin tableau
   explícito): la brecha que queda, la condición para que `compo` (celdas
   ×0.83) sea rentable, y lo que haría barata la herencia de base entre nodos.
   **Modelo de costo (2026-10-02, `perf` con pila; µs por operación inflados
   por igual por `perf`).** Pivote del tableau contra la suma de lo que haría
   un pivote revisado (Forrest–Tomlin + BTRAN + FTRAN, medidos con
   `DFB_SX_LU=actualizada` y `DFBH_LUSTATS=1`):

   | instancia, linealizador | filas `m` | pivote tableau | pivote revisado (estimado) | `pivot` + `compute_basics` en el tiempo total |
   |---|---|---|---|---|
   | `ex2_1_9`, afín | 4 | 0.46 | ~5.1 (×11) | 7 % |
   | `ex6_1_3`, afín | 20 | 3.99 | ~5.9 (×1.5) | 15 % |
   | `ex5_3_2`, afín | 32 | 8.73 | ~6.9 (×0.8) | 38 % |
   | `launch`, afín | 50 | 15.6 | ~7.2 (×0.46) | 26 % |
   | `ex6_1_3`, `compo` | 53 | 14.4 | ~9.7 (×0.67) | 59 % |
   | `ex5_3_2`, `compo` | 81 | 31.7 | ~15.2 (×0.48) | 73 % |

   Las curvas se cruzan hacia `m` ≈ 25–30, y con el linearizador afín **sólo 3
   de 134 instancias** pasan de 25 filas (`ex5_3_2`, `ex5_4_3`, `launch`); 93
   tienen 10 o menos, donde el revisado sería varias veces más caro por pivote.
   Para la omisión no paga. Donde sí podría pagar es con `compo` (53–81 filas),
   donde el pivoteo es el 47–61 % del tiempo y el revisado lo bajaría a la mitad
   o dos tercios: es lo que haría competitivo a `compo`, que poda más (celdas
   ×0.83). O sea que la opción 3 es, en los hechos, «hacer viable `compo`», con
   un motor que elija revisado o tableau según `m`.

   **Prototipo (2026-10-05), `DFB_SX_MOTOR=revisado`**, apagado. Simplex dual
   revisado sin tableau: la fila del pivote sale de una BTRAN y un producto por
   `Ā`, el test de razón es el mismo (paso largo, Harris), `d` e `y` se
   actualizan con la fila, la base cambia en la LU con Forrest–Tomlin, `x_B` se
   recalcula con un producto y una FTRAN, y devex usa una FTRAN de la columna
   entrante. Sin pasos primales. Test unitario 528/528 y 200/200. Medido con el
   orden `lejos` en los dos motores (el orden por costo lee el tableau),
   semilla 1:

   | sistema (filas) | tableau | revisado |
   |---|---|---|
   | `Eiger-0100` | 47.3 s | **0.56 s** |
   | `Bratu-0060` | 45.4 s | **4.9 s** |
   | `BroydenTri-0100` | tope (155 s) | **3.5 s**, resuelto |
   | `BroydenBanded-100` (253) | no termina en 200 s | 94 s |
   | `BroydenBanded-120` (282) | 114 s, 5 celdas | 187 s, 7 celdas |
   | `Katsura-16` (43) / `Katsura-14` (56) | 43.0 / 16.4 s | 35.5 / 14.8 s, mismas celdas |
   | `ex5_3_2` (32, `ibexopt`) | 2.91 s | 2.55 s |
   | `ex6_1_3` (~20, `ibexopt`) | 26.1 s | 34.3 s (×1.31) |
   | `ex2_1_9` (~4, `ibexopt`) | 5.0 s | 7.9 s (×1.57) |

   Confirma el modelo de costo: el revisado gana desde ~30 filas, por uno o dos
   órdenes de magnitud con 100 variables, y pierde con relajaciones chicas por
   el costo fijo de las llamadas a la LU. Sin devex hacía ×1.5–1.9 más pivotes
   por LP.

   **Híbrido por nodo (2026-10-05), `DFB_SX_MOTOR=auto`**: revisado si la
   relajación tiene más de `DFB_SX_MOTOR_M` filas (umbral 30 desde 2026-10-13), tableau si no.
   En los nodos revisados el refresco y el certificado usan la LU y el orden de
   cotas pasa a `lejos`. Con umbral 30 `ex5_3_2` (32 filas) empeoraba 1.75 →
   2.62 s por perder el orden por costo; con 40 el revisado entra sólo donde su
   ventaja lo compensa. Bancos:

   | comparación (celdas / cpu) | > 0.1 s | > 1 s | > 5 s | terminan |
   |---|---|---|---|---|
   | solver, auto contra tableau, XTaylor | ×0.99 / **×0.87** | ×1.00 / **×0.81** | ×1.00 / **×0.75** | 147 contra 142 |
   | solver, auto contra tableau, afín | ×1.00 / ×0.93 | ×1.00 / ×0.88 | ×1.00 / ×0.80 | 161 contra 160 |
   | **solver, DFB auto/afín contra prod/XTaylor** | ×0.64 / **×0.54** | ×0.49 / **×0.47** | ×0.44 / **×0.43** | **161 contra 138** |
   | `ibexopt`, auto contra omisión (4 semillas) | ×1.00 / ×1.01 | ×1.00 / ×1.01 | ×1.00 / ×1.01 | 131 y 131 |

   0 discrepancias en cajas solución y 0 óptimos incompatibles. En el solver el
   motor cambia el costo y no el árbol (celdas ×1.00); con XTaylor gana 6
   sistemas y pierde 1 (`DiscreteBoundary-0200`). En `ibexopt` es neutro: casi
   ninguna relajación pasa de 40 filas (525 de 552 pares con celdas idénticas).
   `results/ibexsolve_motor_auto_2semillas.csv`,
   `results/ibexopt_motor_auto_4semillas.csv`. Quedan `BroydenBanded-100…140`,
   con LPs de más de 200 pivotes: ahí el límite ya no es el motor.

   **Por omisión desde 2026-10-05:** el motor `auto` con umbral 40 en DFB
   (`DFB_SX_MOTOR=tableau` lo apaga) y el linearizador afín en
   `ibexsolve_dfb` (`--lin=xt` vuelve a XTaylor).

   **El reloj de SoPlex (2026-10-06).** `perf` sobre `compo` mostró que el
   17–32 % del tiempo era `times()`: el `SLUFactor` arranca un `UserTimer` en
   cada FTRAN/BTRAN, y esa LU es la del motor revisado. `DFBBasis` usa ahora
   una subclase con el reloj apagado (`DFB_SX_RELOJ=1` lo restaura). Semilla 1
   del banco de abajo, con contra sin reloj, celdas idénticas: cpu ×0.92–0.93
   con el afín y ×0.83 con `compo` (`results/ibexsolve_reloj_semilla1.csv`; la
   corrida nueva iba con 15 procesos y la vieja con 8, así que es una cota
   conservadora). En `ibexopt` sólo pesa donde hay relajaciones grandes
   (`launch` ×0.85; `ex6_1_3` neutro). Las cifras del motor híbrido de arriba
   tenían este impuesto. En el `LPSolver` de Ibex se puso `TIMER_OFF`
   (`IBEX_LP_RELOJ=1`), pero en SoPlex 4.0.2 ese parámetro no llega a la LU ni
   a las estadísticas: producción en `ex6_1_3` ×0.964 y siguen ~240 000
   lecturas de `times()` por segundo. Apagarlas exige parchear el SoPlex
   incluido. Ojo: el wrapper fuerza el simplex primal por columnas
   (`ibex_LPLibWrapper.cpp`, commit `b59a201f`), así que el brazo de
   producción no usa la configuración de fábrica de SoPlex.

   **Ablación (2026-10-08).** Sobre la omisión vigente, un brazo por pieza
   quitada, solver (264 sistemas, 2 semillas) e `ibexopt` (139 instancias, 4
   semillas) (`results/ablacion_2026-10-08.csv`). En `ibexopt` se excluyen
   `alkyl`, `alkylbis` y `bearing`, que explotan o no según la semilla y el
   brazo (`alkylbis` de 148 a 32 416 celdas en el mismo brazo) y dominan el
   corte > 5 s.

   | brazo contra DFB completo | solver: celdas, cpu > 0.1 / > 1 / > 5 s | `ibexopt`: cpu > 0.1 / > 1 / > 5 s |
   |---|---|---|
   | primal con el motor híbrido (defecto, ver abajo) | ×1.59–2.15, ×1.07 / ×1.14 / ×1.22; 325 contra 354 terminan | ×1.01 / ×1.04 / ×1.00 |
   | primal sobre tableau | ×1.00, ×0.96 / ×0.97 / ×1.01; 350 contra 354 | ×1.00 / ×1.00 / ×1.00 |
   | sin paso largo | ×1.01, ×0.99 / ×0.99 / ×1.01 | ×1.02 / ×1.00 / ×0.97 |
   | sin tolerancia a escala de la caja | ×1.00, ×1.00 / ×1.00 / ×1.00 | ×1.00 / ×0.99 / ×1.00 |
   | sin devex | ×1.00, ×1.03 / ×1.03 / ×1.05 | ×1.01 / ×1.01 / ×1.00 |
   | sólo tableau (sin motor híbrido) | ×1.00, **×1.13 / ×1.19 / ×1.25**; 350 contra 354 | ×1.01 / ×1.02 / ×1.00 |
   | primal contra dual, los dos sobre tableau | ×1.01, **×0.85 / ×0.82 / ×0.80** | ×0.99 / ×0.98 / ×1.00 |

   Lectura: ninguna pieza cambia la poda (celdas ×1.00); todas son costo. La
   que pesa es el **motor híbrido** (solver ×0.80–0.88), y el motor revisado
   exige el dual. Con el mismo motor (tableau) el **primal es más rápido que
   el dual** en el solver (×0.80–0.85): el valor del dual es habilitar el
   revisado en las relajaciones grandes, no la poda ni el costo por pivote en
   las chicas. Paso largo, τ y devex son de segundo orden ahora (±5 %); en el
   solver devex ayuda un 3–5 %. En `ibexopt` todo es neutro (±3 %).

   **Defecto encontrado y arreglado.** `primal_step` lee el tableau, así que en
   los nodos revisados una estrategia que no reubica (`primal`, `mixta`) no
   daba ningún paso primal: terminaba con una base dual-infactible declarada
   OPTIMAL y un certificado válido pero flojo (sonda `DFBH_ABHULL`: en 8000
   cotas donde producción ganaba, producción en el hull exacto siempre, DFB
   nunca). Ahora esas estrategias se quedan en el tableau
   (`rev = pivoteo_dual && ...`). La nota de §1 «la poda es la misma en las
   tres» se midió en `ibexopt`, donde casi no hay nodos revisados.

   **Por qué DFB contrae más que `CtcPolytopeHull` con la misma relajación
   (2026-10-07).** Sonda `DFBH_AB` (misma caja de entrada, mismo poliedro
   afín, determinista) más `IBEX_PH_STATS` y `DFBH_ABHULL`. **Rigor:** en 973
   cotas donde DFB contrae más, se resolvió el LP exacto del poliedro sobre la
   caja de entrada: DFB nunca lo pasa y llega a él en el 97 %; producción en el
   3 %. **Mecanismos**, por sistema: (1) **el salto de Achterberg de
   producción** da por hecha una cota si el `x*` anterior está a menos de
   1e-8·|cota|, relativo a la magnitud y no al ancho: en cajas angostas lejos
   del cero deja sin resolver cotas contraíbles (22 % de las cotas en
   Discrete-Integralf2-13 y Bellido; sin el salto producción supera a DFB en
   esos dos); (2) **la pasada se corta** cuando SoPlex da vacío sin prueba o
   desconocido (9–11 % de las cotas sin resolver en Bratu-0060 y Geneig), más
   los nodos que no se podan por esos vacíos; (3) **la holgura de
   Neumaier–Shcherbina** es menor: 2–6 % de la contracción posible por LP,
   20–28 % en Discrete-Integralf2-13. DFB usa la misma regla de salto, y el
   salto relativo al ancho (`DFBH_SALTEO=ancho`, `IBEX_PH_SALTEO=ancho`) sube
   la contracción de los dos sobre la misma caja sin cambiar la razón.

   Banco del salto relativo al ancho, cada contractor con y sin
   (`results/salto_achterberg_ancho.csv`), cpu > 0.1 / > 1 / > 5 s: solver,
   los dos ×1.02 / ×1.02 / ×1.03 con celdas iguales; `ibexopt`, DFB ×0.97 /
   ×0.84 / ×0.66 y producción ×0.98 / ×0.92 / ×0.87, casi todo por
   `alkyl`/`alkylbis` (las semillas que explotaban con SoPlex de fábrica
   vuelven a decenas de celdas: no era sólo el buscador de `loup`; hipótesis
   sin verificar). DFB contra producción con la misma regla: solver igual
   (celdas ×0.48–0.53, cpu ×0.72–0.78) con las dos; `ibexopt` ×0.86 / ×0.84 /
   ×0.78 con Achterberg y ×0.87 / ×0.89 / ×0.86 con el ancho. Conclusión: el
   salto explica parte de la contracción extra sobre la misma caja, pero no
   la ventaja en el árbol del solver, que viene de (2): el certificado en
   intervalos de DFB siempre es utilizable, no depende del rayo ni del dual de
   SoPlex, y la pasada no se abandona. Omisión: Achterberg, como upstream.

   **Cortes en `x*` (2026-10-06, cerrado como agregado).** Tras cada LP
   óptimo, para cada restricción que `x*` viola, una fila XTaylor en la esquina
   de la caja más cercana a `x*` (jacobiano en intervalos sobre la caja de
   linealización, lado derecho en intervalos: riguroso); viajan con la celda
   (`BxpDFBCortes`) y entran en la relajación afín de la siguiente llamada del
   punto fijo (`DFBH_CORTES=1`) o también en las hijas (`=2`), hasta
   `DFBH_CORTES_TOPE` veces las filas afines. Banco del solver, 264 sistemas, 2
   semillas, tope 0.25 (`results/ibexsolve_cortes_xstar_2semillas.csv`):

   | brazo contra sin cortes | > 0.1 s | > 1 s | > 5 s | terminan |
   |---|---|---|---|---|
   | en el nodo | celdas ×0.84, cpu ×1.24 | ×0.82, ×1.20 | ×0.84, ×1.18 | 350 contra 354 |
   | heredados | ×0.82, ×1.30 | ×0.80, ×1.28 | ×0.80, ×1.26 | 344 contra 354 |

   | sólo en las hijas (`=3`, contra un brazo sin cortes del mismo binario; control nulo ×1.000 / ×1.005) | ×0.88, ×1.31 | ×0.86, ×1.30 | ×0.84, ×1.25 | 349 contra 354 |

   **Mecanismo** (`DFBH_CUENTAS=1`: celdas, llamadas, LPs y pivotes). El costo
   por celda se descompone en llamadas/celda × LPs/llamada × pivotes/LP ×
   µs/pivote, y domina uno distinto según el sistema: (1) **vueltas extra del
   punto fijo** en el nodo: los cortes contraen más, `CtcFixPoint` supera su
   umbral y vuelve a llamar a DFB, que re-linealiza y resuelve otros `2n` LPs
   (Bellido 3.10 → 4.11 llamadas/celda, Redeco8 1.00 → 1.11); (2) **más
   pivotes por LP** por las filas extra (Katsura-16 15.0 → 19.6); (3) en
   sistemas grandes, generarlos (Trigexp1-100, µs/pivote ×3.8; no se arregló
   con gradientes por restricción). Aplicarlos sólo en las hijas elimina (1)
   pero los cortes heredados separan los `x*` del padre, que son extremos de su
   caja y caen fuera de la hija en la mitad de los casos: valen pero no
   aprietan (el precedente del caché de `γ`), y (2) se paga igual. El modo 3
   necesitó tope también al heredar: sin él el pozo acumula los cortes de todos
   los ancestros. Variante sin probar: **reemplazar** la fila afín de la
   restricción por el corte cuando es más fuerte en `x*`, sin agrandar el LP.

   **DFB contra producción en el solver, sin reloj (2026-10-06).** Las
   cifras de arriba contra producción pagaban el reloj de SoPlex y no valen.
   264 sistemas, 2 semillas, `lsmear` en los tres primeros brazos
   (`results/ibexsolve_prod_sinreloj_2semillas.csv`). Celdas / cpu:

   | comparación | > 0.1 s | > 1 s | > 5 s | terminan |
   |---|---|---|---|---|
   | DFB/afín contra producción/afín (SoPlex de fábrica, como upstream) | ×0.53 / ×0.79 | ×0.48 / ×0.77 | ×0.46 / ×0.73 | 354 contra 320 |
   | DFB/afín contra producción/afín con `IBEX_LP_REP=col` | ×0.76 / **×0.96** | ×0.73 / **×0.99** | ×0.71 / **×0.92** | 354 contra 323 |
   | DFB/afín contra el oficial (XTaylor, `smearsumrel`, fábrica) | ×0.58 / ×0.69 | ×0.50 / ×0.64 | ×0.49 / ×0.63 | 354 contra 318 |
   | producción, columnas contra fábrica | ×0.70 / ×0.90 | ×0.69 / ×0.88 | ×0.68 / ×0.88 | 323 contra 320 |

   **SoPlex de fábrica debilita a `CtcPolytopeHull` en el solver.** Upstream
   (ibex-team/ibex-lib, develop y master) no fija algoritmo ni representación,
   y en muchos LPs del solver SoPlex elige la representación por filas; ahí
   aparecen cientos de LPs declarados vacíos sin prueba de Neumaier–Shcherbina
   (la prueba es idéntica a la de upstream), que cortan la pasada sin vaciar la
   caja (`IBEX_PH_STATS=1` cuenta los estados). Discrete-Integralf2-15: 4641
   celdas y 2519 vacíos sin prueba de fábrica, 19 y 2 con columnas. No está
   aislado si el rayo falla la verificación o si el LP no era vacío. Las
   columnas no arreglan todo: producción deja de terminar `Brown-15…20`. En
   `ibexopt` casi no hay vacíos sin prueba con fábrica (0–1 por instancia) y
   las columnas sólo mueven el árbol en las dos direcciones. Omisión: la de
   upstream. **Mérito honesto del contractor en el solver: contra la
   producción más fuerte, celdas ×0.71–0.76 y cpu ×0.92–0.99, 31 sistemas
   más terminados.**

   **Bisector y `compo` en el solver (2026-10-06).** 264 sistemas, 2
   semillas, todo con `--lr=dfbhull` y sin el reloj
   (`results/ibexsolve_bisector_compo_2semillas.csv`). Celdas / cpu en
   > 0.1 / > 1 / > 5 s:

   | comparación | > 0.1 s | > 1 s | > 5 s | terminan |
   |---|---|---|---|---|
   | `lsmear` contra `smearsumrel` | ×0.95 / ×1.00 | ×0.96 / ×1.02 | ×0.95 / ×1.00 | 354 contra 346 |
   | `lsmeardfb` contra `smearsumrel` | ×1.13 / ×1.10 | ×1.16 / ×1.14 | ×1.15 / ×1.16 | 338 contra 346 |
   | `compo` contra afín (los dos `lsmear`) | ×0.60 / ×1.62 | ×0.56 / ×1.58 | ×0.53 / ×1.53 | 341 contra 354 |

   `lsmear` (omisión desde 2026-10-06) baja las celdas un 5 % y lo paga con su
   LP: neutro en tiempo, 8 terminadas más (`DiscreteBoundary-0200`,
   `Geneig-icse`, `Yamamura`); pierde en `exnewton` ×3.2 y
   `CountercurrentReactors2`. `lsmeardfb` en el solver suma las
   sensibilidades `abs` de todos los LPs de la pasada (no hay objetivo):
   gana toda la familia `BroydenBanded` (×0.52, termina la 200) y
   `Discrete-Integralf2`, pero deja de terminar `Redeco9`, `Eco10`,
   `Fredtest`, `Fourbarbis-1`. `compo` parte el árbol a la mitad, mucho más
   que el ×0.83 de `ibexopt`, y cuesta ×2.8 por celda aun sin el reloj; gana
   los `Bratu` grandes. Es el argumento para los cortes dirigidos en `x*`:
   las filas extra compran celdas, el problema es agregarlas a ciegas.

   **Nunca se implementó** (hasta este prototipo): lo que se midió como «camino sin tableau»,
   `DFB_SX_EXACT`, mantenía el tableau y además recalculaba `d` y `x_B` desde
   la LU en cada iteración (§21 de `MEDICIONES_DUELO.md`: no reprodujo, 31.6
   contra 29.9 piv/LP, con la estrategia que en los hechos era primal).
4. ~~**Corte del LP del objetivo contra el `loup`**~~ **Cerrado
   (2026-10-02): techo nulo por construcción.** El optimizador pone el `loup`
   como cota superior de `y` en la caja, así que cuando la relajación implica
   `y > loup − eps` el LP no termina `OPTIMAL` con una cota que lo supera: el
   poliedro con la caja es vacío, el dual lo detecta como `INFEASIBLE` y el nodo
   muere por Farkas en la primera cota, sea cual sea `k` (`DFBH_GOALPROBE=1`:
   0 nodos muertos por la cota del objetivo en 8 instancias; la sonda tenía un
   defecto, contaba después de la rama que vacía la caja, ya corregido). La
   variante «comprobar Farkas en cada pivote» tampoco paga
   (`DFBH_VACIO_PROBE=1`): el primer certificado que alcanza llega al 91–97 %
   del LP, y lo ahorrable es 0.2–0.6 % de los pivotes.
5. **Dual estricto**: eliminar los pasos primales que quedan sin que el Harris
   absoluto cueste pivotes.

No reabrir salvo evidencia nueva: cortes anytime por presupuesto o brecha,
caché de `γ`, HC4 intercalado, congelar la relajación, el test de `γ` en las
rodajas de ACID (re-medido el 2026-09-30: con los `γ` propios acierta ≈ 0 %;
con los de un ancestro mata rodajas que el lineal del nodo mata igual) y la cota
lagrangiana con monotonía (contrae de verdad, a precio de linealización). Todas
con la misma explicación (§6).

**Caja negra sin presolve y otros solvers (2026-10-19).**
- SoPlex sin presolve (`DFBH_CAJA_SIMP=0`) es la caja negra desde ahora: a
  60 s, ×0.97 de tiempo en ambos bancos (×0.94 en opt > 1 s), mismas celdas.
  Tabla principal contra ella: opt celdas 0.99–1.02, cpu 0.85/0.85/0.87,
  216/216; sistemas 0.95–0.98, 0.93/0.93/0.92, 364/362; sistemas con «cerca»
  en ambos 0.99–1.00, 0.94/0.95/0.95, 362/362 (`results/final_2026-10-19.csv`,
  brazo `cajasp`).
- Banco de pasadas, SC/caja por pasada (sistemas / opt), misma contracción:
  SoPlex sin presolve 0.62 / 0.56, Gurobi 12.0.3 0.26 / 0.20, HiGHS 1.15.1
  0.18 / 0.12. SoPlex es el más barato por llamada de los tres.
- Cascada contra SoPlex sin presolve: refactorizar 58/53 %, vectores 9/8 %,
  verificar 23/21 %, copiar 2/2 %, resto 8/16 %.
- Pendiente: tabla de 3600 s y sonda en la búsqueda (figura de costo) siguen
  con la caja con presolve; caja con SPxSolver directo (mantener la LU entre
  cotas) sin probar.
- Resultados y scripts en `results/solvers_2026-10-19/`.

**SPxSolver directo como caja negra (2026-10-19).** El núcleo de SoPlex 4.0.2
sin la clase SoPlex (`DFBH_CAJA_SPX=1`): LP cargado, entre cotas solo cambian
objetivo y cotas, relojes apagados (también los de la LU: sin eso, 41 mil
llamadas a `times()` lo hacían más lento que la clase).
- Refactoriza igual ~1 vez por LP: la hace el propio simplex de SoPlex antes
  de declarar el óptimo (spxsolve.cpp), no la envoltura.
- Banco de pasadas: SPx ×0.755 de la clase en las 71 instancias; SC/SPx 0.80
  (sistemas, 27/36 más barato) y 0.72 (opt, 25/35). Contra la clase eran 0.62
  y 0.56: la envoltura explicaba cerca de la mitad de la ventaja.
- 60 s, banco completo (brazo `cajaspx`): opt celdas 0.99–1.03, cpu
  0.91/0.91/0.93, 872/868; sistemas 0.94–0.98, 0.97/0.95/0.95, 364/362;
  sistemas con «cerca» en ambos 0.98–0.99, 0.98/0.97/0.98, 362/362. SPx
  contra la clase: cpu ×0.93 opt, ×0.96 sistemas, mismas celdas.
- SC pierde donde pivotea mucho más (Bratu-0120 10.3 contra 3.3 piv/LP,
  BroydenTri, launch, Katsurabis-14); gana grande por robustez
  (DiscreteBoundary-0200 ×0.19, Geneig-icse ×0.45).

**Reorganización del paper (2026-10-19).** Encuadre: un simplex diseñado
alrededor del certificado (cualquier base certifica → sin verificación, sin
refactorizar antes del óptimo, sin fase 1, un tableau para las 2n cotas).
Resumen y contribuciones nuevos; párrafo «What the certificate removes» al
inicio de §4; §4.2 sin los µs viejos; §6 reordenado: rigor → principal (con
filas de SPxSolver) → largas → costo de una pasada (banco con cuatro cajas +
figura escalera `fig:ladder`) → «Where the gain comes from» (ablación, luego
emulación, luego callgrind sin presolve) → dónde pierde cotas la caja (con la
explicación de los «no óptimos» de SC: infactibilidades cuyo rayo no
certifica, 1.0–1.1 % de sus LPs, nunca tope de iteraciones) → relación con el
PHull de Ibex; discusión y conclusión nuevas. Se quitó la sonda en la búsqueda
(`fig:cost`, medida con presolve). Escalera de la emulación contra SPxSolver:
SPx directo quita 42–54 % de la brecha en log; una refactorización por LP sola
explica 53–58 %. Callgrind sin presolve (ex2_1_9, 200 cajas): clase SoPlex 52
mil instr/LP (32 simplex: 10 antes del primer pivote, 8 pivotes; 10 guardar y
evaluar, 7 verificar; resto interfaz), SPxSolver 37 mil, SC 14 mil (23 con
orden y certificación).
Pendiente en cola (`cadena_larga.sh`): brazo de la caja a 3600 s sin presolve
(`d3600p`, 192) + 7200 s (12) y PHull de producción a 60 s (`prod`, 1444).

**Largas sin presolve y PHull de producción (2026-10-20,
`results/largas_2026-10-20/`).** `tab:long` contra SoPlex sin presolve
(`d3600p`, 7200 en `d7200p`): sistemas (SC propio de `frio_2026-10-15/sc3600`)
celdas ×1.00, cpu ×0.98, 94/93 (SC termina CCR2-28 s1 en 1640 s; la caja no);
opt (SC de `l3600c`) celdas ×1.01, cpu ×0.89, 96/84 (96 con 7200: ex6_2_13,
hs101, hs112 en 3620–4053 s); sin bearing ×0.95 / ×0.83. Fila «same bound
order» a 3600 s pendiente: SC cerca en sistemas corriendo (`l3600sc`, 96).
PHull de producción a 60 s (`prod`, brazo en `final_2026-10-19.csv`; solver
con `IBEX_LP_REP=col`, ibexopt de fábrica): SC/prod sistemas celdas 0.70–0.76,
cpu 0.66–0.68, 364/325; opt celdas 0.92–0.97, cpu 0.80–0.82, 872/865. Caja
sin presolve/prod: sistemas celdas 0.71–0.77, cpu 0.71–0.72; opt cpu
0.92–0.96. Párrafo del paper con números.
- Fila «same bound order» a 3600 s (SC cerca `l3600sc` contra `d3600p`):
  celdas ×1.00, cpu ×0.98, 93/93 (SC termina CCR2-28 s1 en 339 s; la caja
  termina CCR2-30 s2 en 799 s y SC no). `tab:long` completa, sin pendientes.

**Relajación compo a 60 s (2026-10-20, `results/compo_2026-10-20/`, brazos
`sccompo` y `cajacompo` en `final_2026-10-19.csv`).** SC/caja con compo: opt
celdas 0.95–1.00, cpu 0.97–0.98, 875/872 (la ventaja de 13–15 % casi
desaparece); sistemas celdas 0.95–0.98, **cpu ×1.38–1.46**, 355/368 (Bratu
×4–4.7, DiscreteBoundary-0200 ×15.8). Compo contra afín: caja en sistemas
celdas ×0.52–0.57 con cpu ×0.93–1.03 y 368 terminadas (6 más) → en sistemas
**la mejor configuración es caja+compo**: SC+afín es ×1.12–1.15 más lento con
×1.7 celdas. En opt SC+afín sigue siendo la más rápida (×0.74–0.76 de
caja+compo, ×1.26 celdas).
Mecanismo (sonda en la búsqueda y sombra, Katsurabis-14, Bratu-0120): no es
el costo por pivote (SC 21 µs contra 15; en Bratu 54 contra 60) sino los
pivotes: con compo SC 64 piv/LP; SoPlex sobre los mismos LPs 36; SoPlex en su
propia pasada 17. Dos efectos: la cadena de bases de SC deja LPs más difíciles
(73 básicas fuera de cota al empezar, contra 5 con afín) y sobre el mismo LP
SC necesita ×1.8 (degeneración 24 % → 53 %). El warm start primal no lo
arregla (71 piv/LP). Ninguna regla aislada lo explica (sin long-step 49 en
Katsurabis pero 159 en Bratu; sin perturbación 180; sin devex 86; sin fusión
103). `instalar_base` falla en el motor revisado (no se pudo probar SC desde
la base de SoPlex con m=146).
Confundido: tope de iteraciones SoPlex 100 contra SC 20m+200; con compo en
Bratu-0120 SoPlex abandona el 18 % de sus LPs (`optimos caja 2057/2501` en la
sombra). Nueva opción `DFBH_CAJA_TOPE=igual` (caja con 20m+200); sonda en
curso.
Fusión de filas paralelas (`DFBH_FILAS=paralelas|quita`, en el contractor,
para los dos motores): implementada, resultados inconsistentes (DiscreteBoundary
SC 174 → 6.5 piv/LP pero 705 de 736 óptimos sin certificado y la caja ×27 más
lenta por pasada; Brown-30 SC 2 → 17 piv/LP). Hay algo mal o mal entendido:
NO usar hasta depurar.
- Sonda con `DFBH_CAJA_TOPE=igual` (12 instancias, afín y compo): el costo por
  pasada SC/caja cambia < 5 %; los LPs sin óptimo de la caja caen a 0 en casi
  todas (Geneig-icse afín 377 → 8; compo 3429 → 2000; Bratu-0120 compo 93 → 0)
  y los nodos que SC gana bajan (Geneig-icse afín 129 → 26 de 1710; Bratu afín
  28 → 6). **La robustez de `tab:samebox` era en buena parte el tope de 100.**
  Barrido de 60 s en curso (`d60t`, 2888): `cajatope` (caja con tope igual) y
  `scx` (SC), ambos con los binarios `_x` del mismo fuente.
- `tab:samebox` rehecha con el tope igualado (`samebox2`, 60 s, seed 1): solo
  cambia Geneig-icse (caja sin óptimo 377 → 8; gana SC 129 → 26 contra 9);
  CCR2-9 66/2, Brown-18 1/1, Bellido 2/6 iguales con los dos topes. Barrido
  `d60t`: la caja con tope igual explora exactamente las mismas celdas que con
  100 en los dos bancos (×1.000, mismo tiempo): `tab:main` no cambia.
- Paper: §6 nueva subsección «A larger relaxation» (`tab:compo`, mecanismo,
  tope), robustez con tono bajo en samebox/discusión/conclusión, frase sobre
  compo en resumen y contribuciones.

**Barrido `d60t` completo (2026-10-20, `results/tope_2026-10-20/`): caja con
tope igual (`cajatope`) y SC (`scx`), binarios `_x`.** La caja con el tope
de SC explora las mismas celdas que con 100 en los dos bancos (×1.000) con el
mismo tiempo; SC `_x` = SC viejo (celdas ×1.000, cpu ×1.01). `tab:main` no
cambia. Diagnósticos etapa 2 en compo (`diag2`, Katsurabis-14 y Bratu-0120):
SoPlex sin escalado igual; `DFB_SX_SCALE=1` destruye la contracción (ya
sabido); perturbación 1e-5/1e-3 igual; **`DFBH_ORDEN=cerca` en SC: Bratu-0120
61.5 → 10.9 piv/LP (pasada 8.2 → 2.3 s, SoPlex 1.07), Katsurabis-14 63.6 →
47.6**. El orden por conteo engaña con compo (no ve las básicas que la
reubicación deja fuera de cota: 73 por LP). Brazo `sccompocerca` (sistemas,
compo, cota más cercana) en curso.
- Brazo `sccompocerca` (sistemas, compo, SC con la cota más cercana): contra la
  caja con compo cpu ×1.17–1.18 (era ×1.38–1.46), celdas 0.93–0.98, 354/368;
  contra SC conteo con compo ×0.77–0.80, mismas celdas (Bratu ×0.3,
  DiscreteBoundary-0200 ×0.21; pierde solo Yamamua1-0018sp ×1.74). Queda ×2–3
  en DiscreteBoundary-0200, Katsurabis-14, Discrete-Integralf2-21.
- **Causa de fondo del orden (2026-10-20):** el orden por conteo (omisión en
  sistemas) cae a «la más lejana» cuando el motor es el revisado (el tableau
  no está al día), o sea en TODAS las relajaciones grandes (afín > 30 filas:
  26 sistemas; compo: 104 de 126 medidos). La más lejana es el peor orden para
  el warm start dual en m grande. Nueva opción `DFBH_ORDEN=costo_o_cerca`
  (conteo con tableau, más cercana con el revisado). Bratu-0120 afín: 10.3 →
  7.5 piv/LP, 7.5 → 7.0 s, mismas celdas. Barrido `d60coc` (sistemas, afín y
  compo, 968 corridas) en curso.
- **Barrido `d60coc` completo (sistemas, `DFBH_ORDEN=costo_o_cerca`).** Afín:
  contra SC actual celdas ×1.00–1.02, cpu ×0.99–1.01, 363/364 (Bratu ×0.93–0.95;
  DiscreteBoundary-0200 ×4.6: con la más lejana se resolvía en la raíz, efecto
  de árbol); contra la caja ×0.93–0.94 (era ×0.92–0.93). Compo: contra SC
  actual ×0.77–0.80 con las mismas celdas; contra la caja con compo ×1.18–1.21
  (era ×1.38–1.46), 355/368.
- **Descomposición de la pérdida con compo** (`perdida/`, seis peores, orden
  nuevo): pasada SC/SoPlex = LPs × algoritmo × cadena × µs/pivote. La cadena
  ya es ~1 (salvo DiscreteBoundary ×4.9) y el costo por pivote ~1; queda el
  **algoritmo**: SC ×1.6–2.5 los pivotes de SoPlex sobre el mismo LP (Trigo
  ×5). Descartado: representación de SoPlex (por columnas: mismos pivotes,
  pivotes más baratos), tolerancia (SC a 1e-6: igual). Hipótesis: SoPlex
  reoptimiza en primal con devex/steepest (SPxSolver pasa a ENTER tras el
  primer LP); SC repara en dual tras reubicar (decenas de básicas fuera de
  cota) y su primal usa Dantzig.
- **SoPlex reoptimiza en primal** (contadores de SPxSolver, compo):
  Katsurabis-14 91 % de sus pivotes primales, Trigo1-0020sp 87 %, aun
  forzando el dual al empezar cada LP (la base queda dual infactible tras
  cambiar el objetivo). SC con warm start primal + devex (`DFB_SX_PIVOTEO=primal
  DFB_SX_REVPRIMAL=1 DFB_SX_DEVEXP=1 DFB_SX_DEVEXPT=1`): CCR2-12 pasada ×2.22 →
  ×1.36, CCR2-14 ×1.96 → ×1.11; Katsurabis-14, Trigo1-0020sp,
  Discrete-Integralf2-21 no mejoran (SC 2.3–3.7× los pivotes de SoPlex con
  cualquier warm start). Con SoPlex a 1e-9 la ventaja de contracción de SC con
  compo desaparece (Trigo 0.050 → 0.118 contra 0.111 de SC): era tolerancia.
  Pendiente para esos tres: degeneración primal (SoPlex desplaza cotas; SC
  solo perturba costos) o calidad del pricing (steepest edge).
- **`costo_o_cerca` es la omisión de `ibexsolve_dfb` desde 2026-10-20.** En
  cola (`cadena_sis.sh`): SC a 3600 s en sistemas (`l3600coc`, 96) y ablación
  de sistemas con el orden nuevo (`abls`, 11 brazos, 5324 corridas).
- **Causa de las infactibilidades falsas con compo (2026-10-20,
  `DFB_SX_INFDBG=k`, `results/compo_2026-10-20/infana.py`).** Trigo1-0020sp
  con compo: 13 832 declaraciones de vacío en 15 s, 96 % en holguras; el 98 %
  en filas cuyo intervalo de holgura mide < 1e-10 relativo (mediana 4e-13 en
  |x| ≈ 20: casi un punto); violación mediana 14 ulps de |x| (p90 22): es
  redondeo; tolerancia en el piso 1e-15(1+|x|) ≈ 5 ulps en el 74 %. La fila
  violada no tiene ninguna no básica de buen signo → INFEASIBLE → el
  certificado falla (el LP no es vacío) → reintento en frío → lo mismo.
  Discrete-Integralf2-21 con compo igual (ancho 2.6e-14, violación mediana 75
  ulps). Con afín, y en Katsurabis-14 con cualquiera, pocas declaraciones y con
  violaciones grandes (genuinas). La tolerancia escalada por caja baja hasta el
  piso de redondeo en filas casi puntuales, que compo produce y el afín no.
  Arreglo candidato (sin probar): piso de la tolerancia primal a la escala del
  error de calcular x_B, no 1e-15 relativo.
- **Piso de tolerancia al error de la fila (`DFB_SX_PISO=error`, K=100,
  apagado por omisión; `results/compo_2026-10-20/piso/`).** Compo, 60 s,
  seed 1:

  | | infactibles | piv/LP SC | pasada SC/SoPlex | tiempo, celdas |
  |---|---|---|---|---|
  | Trigo1-0020sp, base → piso | 59.6 % → 0.3 % | 22.0 → 8.2 | 6.19 → 1.23 | 10.3 → 6.3 s, 653 → 655 |
  | Discrete-Integralf2-21 | 18.4 % → 0.7 % | 21.4 → 12.8 | 2.34 → 1.34 | 2.90 → 2.04 s, 5 → 5 |
  | Katsurabis-14 | igual | 46.2 | 3.74 | igual (no era esto) |

  bearing (opt, el caso que motivó la tolerancia por ancho): 60 s en ambas
  semillas con y sin piso, mismo óptimo, 21.8–22.8 mil celdas: no conclusivo
  (instancia caótica, excluida de las medias). Falta el barrido de
  confirmación (afín y compo, ambos programas).
- **Peores instancias con el piso (`peores/`, 34 instancias × semillas, 60 s).**
  Cierra la brecha donde la causa era el redondeo: sistemas compo Trigo
  ×2.02 → ×1.06, Discrete-Integralf2-21 ×2.19 → ×1.17; opt compo haldmads
  ×1.80 → ×1.36, concon ×1.34 → ×1.18, hs099 ×1.41 → ×1.33. Afín sin
  cambios. Cambios de árbol: disc2 (opt compo) celdas ×1.9; ex7_2_3 seed 1
  (caótica) 81 → 198 mil celdas. Quedan dos mecanismos más:
  (a) **bases mal condicionadas** — BroydenBanded-100/120/140 con compo
  (×2–4.6, m ≈ 309): aun con el piso 13 % de declaraciones de vacío, ahora con
  violaciones reales; en el 81 % la fila tiene una entrante de buen signo pero
  |a| bajo el umbral relativo amax·1e-11, con max |a| de la fila ~5e8. SC
  pivotea 81/LP contra 7 de SoPlex (que escala). (b) **punto de partida en
  caras degeneradas** — Katsurabis-14, CCR2-14/16, Yamamua1-0016,
  DiscreteBoundary-0200: en los LPs que SoPlex resuelve sin pivotear, SC hace
  20+ pivotes en buena parte.
- **Filas antiparalelas de compo (2026-10-20).** XTaylor escribe cada ecuación
  como dos filas de un lado con normales opuestas (a·x ≤ u, −a·x ≤ −l), a
  distancia m/3 o similar: la fusión de gemelas (consecutivas e idénticas) no
  las ve. BroydenBanded-100: 100 pares de 300 filas; Katsurabis-14 34,
  CCR2-14 10, Trigo1-0020sp 33. La relajación en sí está bien condicionada
  (cond(A) ≈ 6); el mal condicionamiento es de las bases que incluyen ambas
  filas (s_i + s_j = 0).
- **ERROR de nombre:** la opción de fusión de paralelas se llamaba
  `DFBH_FILAS`, que ya existía (cosecha de cotas de otras filas, cara):
  encenderla encendía la cosecha. Eso explica los resultados «inconsistentes»
  anteriores (DiscreteBoundary, Brown). Renombrada a
  `DFBH_FILAS_PARALELAS=paralelas|quita`.
- Con la fusión bien (BroydenBanded-100 compo, piso, seed 1, 30 s): mismas
  celdas (7) en 28.8 s contra 50.1 s; pasada SC 39.8 → 17.1 s, SC/SoPlex 6.75
  → 3.10; piv/LP SC 92.7 → 39.4 (SoPlex 6.9); infactibles 12.1 % → 3.8 %.
  SoPlex casi no cambia con la fusión. Queda ×3 por pivotes (degeneración 66 %).
- **Peores 34 con piso + fusión de paralelas (`peores2/`).** SC/SoPlex, media
  geométrica: sistemas compo ×2.22 → ×1.96 (piso) → ×1.41 (ambos); opt compo
  ×1.39 → ×1.03; afín sin cambios (sistemas ×1.08, opt ×1.02). Por instancia:
  Discrete-Integralf2-21 2.19 → 0.66, DiscreteBoundary-0200 1.85 → 0.71,
  Trigo 2.02 → 0.95, BroydenBanded-100 4.62 → 2.07, concon (opt compo) 1.34 →
  0.68, haldmads 1.80 → 1.02, hydro 1.58 → 1.03. En opt compo la fusión
  cambia el árbol (celdas ×0.6–0.7 en concon, haldmads, hs114). Siguen sin
  cambio CCR2-14/16 (×1.7), Katsurabis-14 compo (×2.2), BroydenBanded-120/140
  (×1.8–1.9), Yamamua1-0016 (×1.5): el punto de partida.
- **Caja negra también con la fusión de paralelas (`peores_caja/`).** La fusión
  ayuda a SoPlex: caja fusionada / caja ×0.90 en sistemas compo, ×0.97 en opt
  compo. En igualdad (los dos con fusión): SC/SoPlex sistemas compo ×1.57
  (era ×2.22 sin arreglos y sin fusión en ninguno), opt compo ×1.06.
  DiscreteBoundary-0200 1.01, Discrete-Integralf2-21 0.82, Trigo 0.98;
  BroydenBanded-120/140 ×2.2–2.4, Katsurabis-14 ×2.3, CCR ×1.5–1.6.
- **`instalar_base` en el motor revisado (2026-10-20)** y sombra con la base
  final de SoPlex instalada en SC antes de cada LP (compo, piso + fusión,
  `sombra_base_rev/`). Pivotes por LP, SC / SoPlex sobre los mismos LPs:

  | | base propia | base de SoPlex, dual | base de SoPlex, primal+devex |
  |---|---|---|---|
  | Katsurabis-14 | 43.4 / 19.0 | 28.9 / 18.8 | 25.1 / 17.5 |
  | CCR2-14 | 21.3 / 13.0 | 22.0 / 11.8 | **12.9 / 12.9** |
  | BroydenBanded-120 | 39.9 / 8.4 | 55.6 / 8.8 | 35.4 / 7.8 |

  CCR: todo es la estrategia (la reubicación del dual destruye el vértice;
  desde la misma base con primal+devex, SC = SoPlex). Katsurabis: mitad base
  de partida (43 → 29), resto algoritmo (×1.4). BroydenBanded: algoritmo (×4.5
  aun desde la misma base y en primal; 8 % de LPs sin óptimo).
- **Tolerancia y BroydenBanded (2026-10-20, `tolerancia/`, compo, piso +
  fusión, seed 1).** Desde la base de SoPlex, el 76 % de las violaciones que
  SC ve al empezar son < 1e-8 relativo (`DFB_SX_INFHIST=1`): la base óptima
  para SoPlex a 1e-6 es infactible para SC. Pasada SC/SoPlex en
  BroydenBanded-120 (dual): base ×3.74 (16 % infactibles) → sin tolerancia por
  ancho ×2.25 → sin ancho y 1e-6 ×1.31 (0 % infactibles, piv/LP 46 → 18),
  mismas celdas. En primal es pésimo (×14). Katsurabis-14 no cambia con la
  tolerancia (×3.2–3.4 dual, ×2.6 primal): queda pendiente.
- **Misma tolerancia que SoPlex** (`DFB_SX_FEASTOL=1e-6 DFB_SX_DUALTOL=1e-6
  DFB_SX_TOLANCHO=0` en sistemas; `DFB_SX_DUALTOL` es nueva) **y qué distingue
  el warm start** (`ws/`, compo, piso + fusión, seed 1, 40 s). Pasada
  SC/SoPlex, dual / primal+devex: BroydenBanded-100 1.52 / 9.88, -120 1.47 /
  7.12; CCR2-14 1.62 / 1.03, CCR2-16 1.89 / 1.03; Katsurabis-14 3.23 / 2.69;
  Yamamua1-0016 1.83 / 1.70 (27 % de LPs infactibles en los dos: por ver).
  NO separan: cambios de lado por LP (BroydenBanded 2–3 %, Katsurabis 12 %, al
  revés de la regla del afín) ni básicas en una cota al empezar
  (`DFB_SX_DEGHIST=1`: BroydenBanded primal 4.7 %, CCR 5.4 %). SÍ separa la
  degeneración de los pasos primales (BroydenBanded 85–88 %, CCR 9 %), pero es
  un resultado, solo se ve probando.
- **ERROR del certificado de vacío de SC (2026-10-20).** Al declarar
  INFEASIBLE, SC entregaba como λ la fila de B⁻¹ de `k` (el multiplicador del
  objetivo), no la fila de la básica que hizo declarar la infactibilidad: solo
  certificaba cuando coincidían (por eso «siempre hay otra fila de B⁻¹ que
  certifica», `DFBH_OTRORAYO`). Opción `DFB_SX_RAYO=fila` (fila `r_inf`).
  Mismas cajas (`DFBH_ABCAJA`, 20 s, seed 1), SC sin óptimo viejo → fila:
  Yamamua1-0016 compo 12 251 → 0 (vacíos 390/745 → 149/149, gana 0/356 →
  7/10); afín 823 → 0; Bellido afín 152 → 0; CCR2-9 afín 110 → 2;
  Katsurabis-14 afín 22 → 5, compo 56 → 34. Explica los «LPs not optimal» de
  SC en `tab:samebox` (Bellido 152, CCR2-9 110): eran este error.
- **Katsurabis-14 con compo y todas las correcciones (`katsurabis/`).** Pasada
  SC/SoPlex: dual ×3.24, primal+devex ×2.70. Se descompone en tres factores
  que se multiplican (primal): (1) punto de partida ×1.5 — desde la base final
  de SoPlex, SC baja de 33 a 22 piv/LP, y en los LPs que SoPlex resuelve sin
  pivotear SC tampoco pivotea (100 %, contra 76 % desde su propia base);
  (2) reglas de pivoteo desde la misma base ×1.25 (22 contra 17.7);
  (3) costo por pivote ×1.4 (19.4 contra 14.0 µs): `DFB_SX_PERF=1` da el
  test de razón dual 26–40 % (BTRAN + fila ρᵀĀ), recálculo de x_B desde cero
  14–18 % (0.5–0.6 por pivote; SoPlex lo actualiza), pivote 28–37 %.
  1.5 × 1.25 × 1.4 ≈ 2.6. Remedios: (1) cómo termina SC dentro de la cara
  óptima; (3) x_B incremental (existía, apagado: cambiaba el árbol por
  redondeo) y fila del pivote hiperdispersa.

**Correcciones por omisión (2026-10-20, binarios `_y`).** `ibexsolve_dfb`,
`ibexopt` e `ibexreplay` fijan en `main`, si no vienen de fuera:
`DFB_SX_RAYO=fila`, `DFB_SX_PISO=error`, `DFBH_FILAS_PARALELAS=paralelas`,
`DFB_SX_FEASTOL`/`DFB_SX_DUALTOL` = 1e-6 (sistemas) o 1e-9 (opt),
`DFB_SX_TOLANCHO=0`, `DFBH_CAJA_TOPE=igual`; y `DFBH_ORDEN=costo_o_cerca`
(sistemas) / `cerca` (opt). Las tolerancias de SC se leen ahora al primer uso
(antes al cargar el programa, y el `setenv` de `main` no tenía efecto).
Verificado: mismas celdas que con las opciones explícitas (Katsurabis-14,
Bratu-0120, Trigo1-0020sp con compo). El warm start no cambia.
- **Punto 3: tope de tolerancia por ancho en el optimizador** (`p3/`, 11
  problemas × 4 semillas, `ibexopt_y` con y sin `DFB_SX_TOLANCHO=1e-3`). Sin el
  tope: mismos óptimos (ex14_2_7 ≈ 0 dentro de 1e-6), hhfair1 cpu ×0.45 y
  celdas ×0.41, concon ×0.91, launch ×1.23 (celdas ×1.37), el resto ×0.97–1.01;
  bearing termina 2 de 4 contra 1 de 4 (caótica). Se mantiene sin tope.
- **Punto 4a en curso: bancos a 60 s con las correcciones** (`d60y`, 5776
  corridas: `sc_art`, `caja_art`, `sc_compo`, `caja_compo`; binarios `_y`;
  la caja con las mismas correcciones del contractor, tolerancia de su
  programa, sin presolve, tope igual).
- **Tolerancia de SC en sistemas: 1e-9 sin tope (2026-10-21).** Con 1e-6 (la
  de la caja) SC perdía contracción con el afín: Geneig-icse 295 → 577 celdas,
  CCR2-9 955 → 1313, Bratu-0210 ×3 tiempo (`tolv/`). Lo que arreglaba
  BroydenBanded con compo era quitar el tope por ancho, no 1e-6 (1e-9 sin tope:
  7 celdas, 30.5 s; 1e-6: 28.6 s). 1e-9 sin tope es igual o mejor en las siete
  instancias probadas. Omisión de `ibexsolve_dfb` e `ibexreplay` cambiada; la
  caja sigue en su mejor valor (1e-6 en sistemas). Brazos de SC en sistemas del
  barrido `d60y` medidos con 1e-6 guardados en `d60y_sol_t6`; se repiten.
- **Punto 4a terminado (`d60y`, `results/final_2026-10-21.csv`, brazos
  `y_*`).** SC/caja en igualdad, celdas / cpu (> 0.1, > 1, > 5 s), terminadas:
  opt afín 0.99–1.00 / 0.84 0.84 0.86, 875/865; opt compo 0.99 / 0.92 0.94
  0.95, 875/871; sistemas afín 0.98–0.99 / 0.94 0.94 0.95, 362/362; sistemas
  compo 0.96–0.98 / 1.10 1.14 1.20, 360/367. Compo / afín: SC opt ×1.12–1.19
  cpu con celdas ×0.67–0.75; SC sistemas ×1.07–1.14 con celdas ×0.45–0.54;
  caja sistemas ×0.86–0.97. Mejor de cada programa: opt SC afín (×0.75–0.81 de
  la caja con compo); sistemas, SC afín / caja compo ×0.90–0.95 de tiempo con
  ×1.7 celdas, 362/367. SC sistemas compo con 1e-9 es ×1.03–1.04 más lento que
  con 1e-6 (afín al revés: ×0.98).
- **Remedios 1 y 2 (2026-10-21, `r12/`, `optc/`).** (2) Cambios de cota
  incrementales en el revisado (`DFB_SX_FLIPINC=1`): mismas celdas, ×0.94–0.98
  → **por omisión** en los tres programas. (1) Tolerancia por motor
  (`DFB_SX_TOL_REVISADO=1e-6`): compo ×0.91 pero Bratu-0210 afín ×3.0 con la
  misma celda (cada pasada contrae menos y el punto fijo con HC4 hace 20
  pasadas en vez de 4) → no se adopta. En su lugar **tolerancia por
  relajación** en `ibexsolve_dfb`: 1e-6 con `--lin=compo`, 1e-9 con el afín.
  En opt con compo, 1e-6 cambia el árbol en los dos sentidos (concon ×0.60,
  launch ×0.64, hydro ×0.76; hs099 ×1.74, haldmads ×1.17; media ×0.94): se
  deja 1e-9, como la caja en opt.
- **Sistemas con las omisiones finales (`d60z`, brazos `z_sc_*`).** SC/caja:
  afín celdas 0.98–0.99, cpu 0.94 0.93 0.94, 362/362; compo celdas 0.96–1.01,
  cpu 1.07 1.12 1.14, 364/367 (antes 1.10–1.20 y 360/367). Peores con compo:
  DiscreteBoundary-0200 ×2.46, Katsurabis-14 ×2.16, CCR2-9 ×2.15,
  BroydenBanded-100/140 ×1.56–1.64, transistor-icse ×1.55.
- **Remedio 3: arranque tibio elegido midiendo (`DFBH_WARM=medir`,
  2026-10-21, `wm/`).** Alterna dual y primal+devex en los primeros 10 nodos,
  mide el costo por LP (promedio móvil) y usa el más barato, volviendo a probar
  el otro cada 25 nodos. Contra el dual (2 semillas, 60 s): compo media ×0.89
  (CCR2-9 ×0.52, CCR2-12/14 ×0.73–0.78, transistor-icse ×0.79, Katsurabis-14
  ×0.92; BroydenBanded ×1.11–1.13 en vez del ×2.25–2.5 del primal fijo;
  DiscreteBoundary-0200 ×1.46, pocos nodos); afín media ×1.03 (Geneig-icse
  ×1.14 con celdas ×1.16: el primal le cambia el árbol aunque cueste menos por
  LP). El primal fijo: compo ×1.01, afín ×1.07.
- `ibexsolve_dfb` con `--lin=compo` fija ahora `DFBH_WARM=medir` por omisión
  (con el afín sigue el dual). Verificación en curso: brazo `sc_compo` de
  sistemas (`d60w`, 484 corridas).
- **Verificación del arranque tibio medido con compo (`d60w`, `w_sc_compo`).**
  SC/caja en sistemas compo: celdas 0.98–1.00, cpu 1.03 1.05 1.06, 364/367
  (antes 1.07–1.14); contra SC compo sin medir ×0.92–0.96.
- **Contractor híbrido (`DFBH_HIBRIDO=T`, 2026-10-21).** Con más de T filas en
  la relajación (la del nodo anterior), la pasada la resuelve SPxSolver directo
  (la caja `DFBH_CAJA_SPX`), con el mismo orden de cota más cercana, la misma
  relajación fusionada y el mismo certificado; si no, SC. Corrección necesaria:
  la caja tiene que recibir la relajación fusionada aunque el motor se decida
  después de linealizar. Katsurabis-14 compo, seed 1: SC 14.3 s, híbrido
  (T=60) 7.1 s, la clase SoPlex 7.8 s, mismas celdas (171–177).
- Tamaño de la relajación por instancia (`filas_m2/`, `tamano.py`): contra la
  clase SoPlex, SC gana en todos los tamaños con el afín (opt ≤ 60 filas,
  sistemas hasta > 120); con compo empata desde ~30 filas en opt y ~15 en
  sistemas, y pierde ×1.54 sobre 120 en sistemas. La caja SPxSolver ahora usa
  también el tope de SC con `DFBH_CAJA_TOPE=igual`. Barrido `d60s2` en curso:
  SPxSolver directo, los dos programas, afín y compo (2888 corridas).
- **SC contra SPxSolver directo con las correcciones (`d60s2`, `s_spx_*`).**
  SC/SPx: opt afín cpu 0.91–0.92, 875/868; opt compo 1.00–1.03, 875/878;
  sistemas afín 0.97–0.98, 362/362; sistemas compo 1.12–1.17, 364/372.
  SPx/clase SoPlex: ×0.89–0.96. Por tamaño, SC gana solo hasta ~15 filas
  (×0.88–0.96); 15–30 empata; con compo pierde desde 15–30 (×1.1–1.5).
- **Híbrido con T=15 y T=30 (`hib/`, 16 instancias × 2 semillas).** En casi
  todas toma el mejor de los dos motores (Katsurabis-14 compo ×0.49 de SC,
  BroydenBanded-120 ×0.50). Fallas: T=15 en Geneig-icse afín ×2.18 de SC (SPx
  pierde cotas: SC era ×0.46 de SPx) y hhfair1 ×1.36; T=30 en hhfair1 ×1.41,
  concon compo ×1.35 de SPx, Katsurabis-14 afín ×1.34 de SPx.
- **Híbrido T=30 en los dos bancos (`d60h`, `h_hib_*`).** CPU (> 0.1, > 1,
  > 5 s), celdas ≈ 1.00 en todo:

  | | / clase SoPlex | / SPxSolver | / SC | terminadas (híbrido / SC / SPx / clase) |
  |---|---|---|---|---|
  | opt afín | 0.85 0.86 0.89 | 0.93 0.93 0.96 | 1.02 1.03 1.01 | 874 / 875 / 868 / 865 |
  | opt compo | 0.90 0.90 0.89 | 0.98 0.98 0.96 | 0.97 0.96 0.93 | 878 / 875 / 878 / 871 |
  | sistemas afín | 0.94 0.93 0.94 | 0.98 0.98 0.97 | 1.00 1.00 1.00 | 362 / 362 / 362 / 362 |
  | sistemas compo | 0.94 0.91 0.90 | 1.01 1.01 1.01 | 0.91 0.87 0.85 | 372 / 364 / 372 / 367 |

  El híbrido es igual o mejor que cada motor solo en todas las
  configuraciones (salvo opt afín, ×1.02 de SC). Peores contra la clase:
  DiscreteBoundary-0200 compo ×2.66, Katsurabis-14 afín ×1.24, ex7_2_3 ×1.29.
- **Banco de pasadas con SC corregido (`banco2/`, `ibexreplay_h2`, afín).**
  SC/caja por pasada (sistemas / opt), misma contracción y mismos vacíos:
  clase SoPlex sin presolve 0.60 / 0.53 (antes 0.62 / 0.56); SPxSolver 0.80 /
  0.71 (igual); Gurobi 0.26 / 0.20; HiGHS 0.18 / 0.11 (eigmaxc con HiGHS
  falló). SC más barato que SPxSolver en 26/36 sistemas y 30/35 opt.
- **Tope de pivotes por LP con compo (2026-10-21, `cap/`, `DFBH_LPCAP`).**
  Con el afín el truncamiento no pagaba (ESTADO §6), pero con compo sí: SC
  sistemas compo, media ×0.76 (tope 30) y ×0.71 (tope 60); BroydenBanded
  ×0.40–0.46 con celdas ×1.8–2.7, DiscreteBoundary-0200 ×0.42 con celdas ×0.58,
  Katsurabis-14 ×0.67; neutro en Geneig-icse, transistor-icse, CCR2-14;
  Bratu-0120 ×1.12 (2.6 s). Barrido en curso con tope 60 (`d60c`, 2888:
  `cap_art`, `cap_compo`).
