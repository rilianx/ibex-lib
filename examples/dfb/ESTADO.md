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

- **`dual`** (omisión desde 2026-09-30): reubica **todas** las no básicas por
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

**Línea de base vigente (2026-09-30)**, con `eps_x` = 1e-7 de verdad (§8) y la
omisión actual de DFB. 138 instancias, 4 semillas, 8 cores, los dos brazos en el
mismo barrido (`results/linea_base_epsx1e-7_4semillas.csv`):

| corte | celdas | cpu | +rápido / −lento |
|---|---|---|---|
| > 0.1 s (n=245) | ×1.018 | **×0.617** | **217 / 10** |
| > 1 s (n=91) | ×0.973 | **×0.591** | 83 / 5 |
| > 5 s (n=45) | ×0.969 | **×0.575** | **45 / 0** |

cpu por celda ×0.584 (> 1 s), ×0.516 en las 22 que llegan al tope en los dos.
0 óptimos incompatibles. Resueltas en 4/4: 131 y 131; DFB resuelve además
`ex8_4_4bis` y `ship-1`, producción `bearing`. Bajar `eps_x` de 1e-6 a 1e-7, por
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
| Herencia de base entre nodos (`DFB_SX_NODOS=1`) | peor que arrancar en frío sobre la relajación nueva |
| Test de `γ` en las rodajas de ACID (`DFBH_ENACID=1`), re-medido tras el filtro de estabilidad | con `γ` propios acierta ≈ 0 %; con los de un ancestro mata rodajas que el lineal del nodo mata igual (§9) |
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

**La brecha que queda** está en las relajaciones grandes: `launch` y `ex5_3_2`
(×1.18–1.20) y `ex14_2_7` (×1.16). El dual reubica y reparte la infactibilidad
por muchas de las ~50 filas; el orden de cotas no es el problema (el predictor
por columnas acierta el costo del dual).

### Las opciones

0. **`bearing` con `eps_x` = 1e-7** (§4): verificar si el buscador de `loup`
   acierta antes de contraer y falla después.
1. **Usar en el resto de `ibexopt` lo que DFB ya calcula.** El punto primal `x*`
   como candidato de `loup` (hoy `LoupFinderXTaylor` resuelve otro LP con SoPlex
   por nodo; y el `loup` es lo que falla en `bearing`) y los multiplicadores del
   objetivo para `LSmear`, que hoy los pide a SoPlex. Techo ~12 % de cpu y
   `ibexopt` sin solver LP externo. Cambia la heurística de bisección: medir
   celdas.
2. **`ibexsolve`**: `DefaultSolver` compone `CtcPolytopeHull` + XTaylor en el
   mismo lugar. Sin cota superior, toda la poda del lineal es prueba de vacío.
   Conservar el linearizador de producción para aislar el cambio. Instancias
   naturales: las Brown/Brent de `data_tests` (que hacen caer a `ibexopt` por no
   tener objetivo).
3. **El rediseño para relajaciones grandes** (BTRAN/FTRAN sin tableau
   explícito): la brecha que queda y la condición para que `compo` (celdas
   ×0.83) sea rentable.
4. **Dual estricto**: eliminar los pasos primales que quedan sin que el Harris
   absoluto cueste pivotes.

No reabrir salvo evidencia nueva: cortes anytime por presupuesto o brecha,
caché de `γ`, HC4 intercalado, congelar la relajación, el test de `γ` en las
rodajas de ACID (re-medido el 2026-09-30: con los `γ` propios acierta ≈ 0 %;
con los de un ancestro mata rodajas que el lineal del nodo mata igual) y la cota
lagrangiana con monotonía (contrae de verdad, a precio de linealización). Todas
con la misma explicación (§6).
