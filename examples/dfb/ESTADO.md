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
| **Motor híbrido**: simplex revisado sin tableau en relajaciones de más de 40 filas (§9) | `ibexsolve` cpu ×0.75–0.93 con las mismas celdas; neutro en `ibexopt` |
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
| `DFB_SX_PIVOTEO=primal\|dual\|mixta`, `DFB_SX_ENTRAK=1`, `DFB_SX_MONO=1` | re-optimización al cambiar de objetivo: `dual` es el dual puro (reubica también `k`), `mixta` la anterior; entrada explícita de `k`; sonda de monotonía (§9). **Omisión `dual`** |
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
| `DFB_SX_DEVEXP=1` | precio devex en los pasos primales del motor revisado (por omisión Dantzig; sin banco aún) (§4) |
| `DFB_SX_REVPRIMAL=1` | pasos primales en el motor revisado (omisión sólo el 2026-10-10 antes de la fusión; hoy apagados, la omisión es dual puro) (§4) |
| `DFB_SX_PIVOTEO=por_motor\|dual\|primal\|mixta` | **omisión `dual` (puro, en los dos motores) desde 2026-10-10**; `por_motor` (dual en revisados, primal en tableau) fue omisión 2026-10-08…10 (§1, §4, §9) |
| `DFBH_MOTOR=caja`, `DFBH_CAJA_TOL=t` | mismo contractor con SoPlex como caja negra y el mismo certificado gamma; tolerancia de SoPlex (omisión 1e-6) (§4) |
| `DFBH_CUENTAS=1` | `ibexsolve_dfb`: celdas, llamadas a DFB, LPs y pivotes al final |
| `IBEX_PH_STATS=1` | `CtcPolytopeHull`: estados de SoPlex por LP (óptimo, vacío probado / sin prueba, desconocido) |
| `DFBH_TESTIGO=archivo`, `DFBH_INFSTATS=1` | verifica un punto solución conocido contra la relajación y las cotas; cuenta filas donde el ±1e100 de SoPlex habría cortado (§9) |
| `DFBH_LATIDO=1`, `DFBH_LPLARGO=1` | latido del contractor cada 200 llamadas; LPs de más de 5000 vueltas y progreso de la pasada |
| `DFB_SX_MOTOR=tableau\|revisado\|auto`, `DFB_SX_MOTOR_M=m` | motor del simplex: tableau explícito, revisado sin tableau, o revisado por encima de `m` filas (**umbral 40** desde 2026-10-10; 80 entre el 08 y el 10) (§9). **Omisión `auto`** |
| `DFBH_LUSTATS=1` | cuenta factorizaciones, updates, BTRAN y FTRAN de la LU (§9) |
| `-b lsmeardfb`, `DFBH_SENS=redcost\|abs\|jac`, `DFBH_BISEC_CMP=1`, `IBEX_LSMEAR_DUMP=1` | bisector `LSmear` con sensibilidades de DFB guardadas en la celda; comparación con `LSmear` sobre el mismo nodo; volcado de puntajes (§9). Fuera de la omisión |
| `DFB_SX_DEVEXVIEJO=1`, `DFB_SX_DSE=1`, `DFB_SX_DSEABS=1`, `DFB_SX_HARRISABS=1` | restaura el defecto de devex (§9); dual steepest edge exacto y su variante con violación absoluta; Harris de dos pasadas. Apagados |
| `DFB_SX_BFRT=0`, `DFB_SX_BFRTINF=1` | ratio test de paso largo, **encendido por omisión** (`=0` lo apaga); la segunda restaura el defecto de la infactibilidad espuria (§9) |
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
   relajación tiene más de `DFB_SX_MOTOR_M` filas (umbral 40), tableau si no.
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
