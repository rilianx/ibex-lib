# DFB en Ibex — estado del trabajo

Documento de traspaso. Sólo lo que quedó en pie, con las cifras que lo
sostienen. El registro cronológico, con todas las idas y vueltas, está en
`MEDICIONES_DUELO.md`.

---

## 1. Qué es

`CtcDFBHull` reemplaza a `CtcPolytopeHull` en el hueco exacto del contractor
lineal: misma interfaz, misma composición, `--filtering=acidhc4 --lr=dfbhull`.

La relajación del nodo se escribe `Ā z = 0` con `z = (x, b)`, `na = nb_var + m`,
y la parte `b` de `Ā` es exactamente `−I`. De ahí salen dos propiedades:

- **Arranque dual-factible gratis**: con las `m` columnas de `b` como base queda
  `B = −I`, no hay que factorizar, y con `c = ±e_k` resulta `c_B = 0`. Sin fase 1.
- **Certificado**: de los multiplicadores `λ` se arma `γ = λᵀĀ` **en
  intervalos**, y como `Āz = 0` para toda solución, vale `0 ∈ γ·z`. Despejando,
  `z_k ∈ (−Σ_{i≠k} γ_i z_i) / γ_k`.

Tres consecuencias que definen el método:

1. **Se puede pivotear en punto flotante.** La cota se recertifica en intervalos
   desde `λ`, así que un `λ` malo da una cota débil y **nunca una falsa**.
2. **`γ` no depende de la caja**: un certificado vale para cualquier sub-caja de
   aquella en que se linealizó.
3. **Prueba de vacío por Farkas**: si `0 ∉ γ·z`, la caja no contiene ninguna
   solución de la relajación. Del mismo tipo que el test de infactibilidad de
   Neumaier–Shcherbina, verificado en intervalos.

### Qué es dual y qué es primal: tres ejes independientes

La palabra «dual» aparece en tres lugares y significa algo distinto en cada uno.

| eje | qué decide | hoy |
|---|---|---|
| **el acotamiento** | de dónde sale la cota certificada | **dual siempre**, por dualidad débil: cualquier `λ` da una cota válida |
| **la selección** (`DFBH_ORDEN`) | qué cota se resuelve después | **por costo dual** (omisión), §3 |
| **el pivoteo** (`DFB_SX_PIVOTEO`) | cómo se **re-optimiza** al cambiar de objetivo | **`primal`** (omisión) o `dual`. El valor describe las cotas 1 en adelante: **la cota 0 de cada nodo es dual en las dos** |

### Las dos estrategias de pivoteo

Lo común, que es la mayor parte:

1. **Por nodo** se linealiza, se hace `load` —que **reinicia la base**, porque la
   herencia entre nodos está apagada— y se recorre una pasada de hasta `2n`
   cotas con **un solo** simplex compartido.
2. **La cota 0 es idéntica en las dos**: base en frío `B = −I`, dual-factible sin
   fase 1 y primal-infactible, o sea pasos duales.
3. **La certificación es idéntica**, y **el bucle también**: paso dual mientras
   haya infactibilidad primal, paso primal mientras haya dual. Lo que cambia
   entre estrategias es **el estado con que se entra**, y por eso cuál de los dos
   pasos domina.

Difieren sólo en las cotas `1 … 2n−1`, donde cambia el objetivo y **la caja no se
mueve**. Por eso el interruptor nombra la re-optimización: la omisión es, en
rigor, «dual en la cota 0, primal de la 1 en adelante».

- **`primal` (omisión).** No se reubica nada. La base conserva la factibilidad
  **primal** —el poliedro no cambió— y pierde la dual, porque `d = c − Āᵀy` se
  recalcula con signos arbitrarios. Son **pasos primales**: el arranque tibio de
  un simplex primal, lo mismo que hace SoPlex.
- **`dual`.** Cada no básica se mueve a la cota que le da el signo correcto
  —sólo las violadas; las degeneradas se quedan— lo que restaura la factibilidad
  **dual** gratis y rompe la primal, porque cada columna reubicada corre `x_B` en
  `ā_{·c}·diam(z_c)`. Son **pasos duales**. La columna `k` del objetivo se
  excluye, por lo que sigue.

**Reubicar `k` nunca puede servir.** Con `c = ±e_k` y `k` no básica, mandarla a
la cota que el objetivo prefiere deja la base dual-factible con `k` **fuera** de
la base: `c_B = 0`, `y = 0`, `γ = 0`. El LP afirma que el óptimo es la cota de
entrada de `z_k` —ninguna contracción— y encima no deja certificado. Y no es un
caso de borde: cuando `k` es no básica, `y = 0` implica `d = c`, así que `k` es
la única columna con costo reducido no nulo. Medido: en el 82 % de las
resoluciones se movía exactamente una columna, siempre `k`. Excluyéndola, las
cotas sin certificado bajan de 86.0 % a 80.0 %, contra 77.5 % del primal. No era
la tolerancia de factibilidad: bajarla de 1e-9 a 1e-13 mueve las dos tasas dos
puntos y deja la brecha intacta. `DFB_SX_REUBK=1` restaura el defecto, sólo para
reproducir mediciones viejas.

### Qué se sigue

- **La poda es la misma, y no puede no serlo.** El vacío se decide en la cota 0,
  que es común: si una cota llega a `OPTIMAL` el poliedro no es vacío y ningún
  `λ` posterior puede probarlo. Medido: vacíos ×0.98–×1.01, medianas 4.2 %
  contra 4.8 % de los nodos.
- **La cota intermedia sólo sirve en la dual.** Con pasos duales la base es
  dual-factible en cada iteración, así que la cota certificada **sube
  monótonamente** desde el primer pivote y detenerse antes deja algo útil. Con
  pasos primales las bases intermedias dan cotas válidas sin garantía de calidad
  hasta el óptimo. **Es la propiedad que distingue a DFB de un simplex, y hoy
  sólo la tiene la estrategia dual.** No está medida (§9).
- **La diferencia que queda es costo.** Banco de 137 instancias, 4 semillas,
  dual contra primal:

  | corte | celdas | cpu | +rápido / −lento |
  |---|---|---|---|
  | > 0.1 s | ×1.000 | ×1.050 | 25 / 91 |
  | > 1 s | ×0.997 | ×1.050 | 7 / 31 |
  | > 5 s | ×0.990 | ×1.029 | 6 / 12 |

  Las mismas 130 instancias resueltas en 4/4 semillas, gana `ship-1`, no pierde
  ninguna, 0 óptimos incompatibles. Son pivotes, ×1.113 por LP, y de un tramo
  acotado: cuando `k` es no básica no se reubica nada, así que en el 72–82 % de
  los LPs el dual **es** el primal. Todo el exceso sale del 18–28 % restante,
  que compra celdas ×1.000, o sea nada medible.

Lo contrario también vale: cuando lo que cambia es la **caja** y no el objetivo
—propagación, ACID— el que hereda es el dual. Es el escenario de la propagación
propia archivada (§8) y la razón de conservarla.

---

## 2. Las optimizaciones que quedaron

Activas por omisión, sin variables de entorno.

| | efecto medido |
|---|---|
| **Prueba de vacío por Farkas** | 87 → 128 instancias resueltas |
| **Filtro de estabilidad en el test de razón primal** | celdas ×1.139 → ×0.935 contra producción (> 1 s) |
| **Orden de cotas por infactibilidades duales** | pivotes/LP ×0.78; cpu ×0.90 / ×0.92 / ×0.96 en los cortes > 0.1 s / > 1 s / > 5 s |
| **Pivoteo primal por omisión** | la estrategia dual cuesta cpu ×1.025–1.046 con celdas iguales |
| **Un solo simplex secuencial** | 6.8–12.3 % más rápido que `2n` ejemplares sembrados |
| **`x_B` sólo si está sucio en `primal_step`** | cpu ×0.99 con resultados **idénticos bit a bit** |
| **Dos costos fijos del linearizador afín** | cpu ×0.942 (> 1 s), 58 instancias más rápidas y **ninguna más lenta**, resultados **idénticos bit a bit** |

### El filtro de estabilidad — leer antes de tocar `primal_step`

El test de razón primal elegía la fila de salida aceptando **cualquier**
coeficiente no nulo, incluso `1e-18`; `pivot` lo rechazaba por `|p| < PIV_TOL`;
y ese rechazo llegaba al bucle **indistinguible de «soy óptimo»**:

```cpp
if (q < 0) return false;              // no hay entrante: soy ÓPTIMO
if (!pivot(r_out, q)) return false;   // el pivote FALLÓ  ← mismo valor
...
st = OPTIMAL; break;                  // la caída es optimismo
```

Medido: el 48 % de los `OPTIMAL` eran con cero pivotes, y de ésos el **83 % eran
dual-infactibles**, con violaciones de hasta `1.7e+06`. El paso **dual** ya tenía
el filtro (Harris); el primal nunca lo recibió.

### Dos costos fijos en el linearizador afín

Los encontró un perfilado, no una idea. Ninguno cambia una sola cuenta, así que
se validan con el criterio fuerte: **resultados idénticos bit a bit**.

- **`pow(2, -50)` estaba dentro del bucle de variables**, o sea `nb_ctr × nb_var`
  llamadas por linealización —736 en `ex5_3_2`—, para una constante. Y `rad` y
  `mid` de la caja se recalculaban una vez por restricción cuando sólo dependen
  de la caja. Izados y constante: la extracción de fila pasa de 263 a 51 µs por
  llamada en `ex5_3_2`, de 53 a 17 en `house`.
- **La aritmética afín pedía un bloque al heap por nodo del árbol de
  expresión**: entre 400 y 2200 asignaciones por linealización, de 7 a 24
  doubles cada una, 20 millones en 15 s en `ex14_2_7`. Ahora se reciclan en un
  pool por tamaño, con una cabecera de un double que guarda el tamaño para poder
  liberarlas sin conocer `n`. El modelo de propiedad lo permite: el constructor
  de copia y la asignación **copian** el contenido, no comparten el puntero.

Juntos dan cpu ×0.951 / ×0.942 / ×0.941 en los cortes > 0.1 s / > 1 s / > 5 s,
con 98, 58 y 29 instancias más rápidas y **ninguna más lenta**, y 0 pares con
celdas distintas en 521 comparaciones. `DFBH_LINSLOW=1` y `DFBH_AFNOPOOL=1`
restauran el estado anterior.

**Cuidado al tocar esto**: es código vendorizado de terceros en el camino
caliente. Las asignaciones de `ibex_AffineVar.cpp` tienen que ir al mismo
asignador que las de `ibex_Affine2_fAF2.cpp`, porque las libera el mismo
destructor; ruteadas a medias, el solver vuelca el core.

### `compute_basics` se llamaba dos veces por iteración

El bucle compuesto llama `dual_step` antes que `primal_step`, y aquél ya
recalcula `x_B` al entrar. Si devuelve falso es porque no había infactibilidad
primal que reparar, o sea que no tocó nada. `primal_step` lo recalculaba igual:
el mismo `O(m·nnz)` dos veces. Medido antes del arreglo: 2.06 y 2.25 llamadas
por pivote en `ex5_3_2` y `house`, cuando corresponde una.

Con la guarda `if (xb_sucio)`, las llamadas bajan 25–34 %, `compute_basics` pasa
del 20 % al 15 % del simplex y el banco da cpu ×0.99 en los tres cortes, **con
las celdas idénticas en los 521 pares resueltos**. No es la actualización
incremental del §44, que sigue apagada por deriva numérica: acá no cambia
ninguna cuenta, sencillamente no se ejecuta dos veces. `DFB_SX_NOXBSKIP=1`
restaura el comportamiento anterior.

### Arquitectura de estados

Hubo una versión con `2n` ejemplares del simplex sembrados desde la base óptima
de la cota 0. Usaba ×0.75 pivotes, pero **el ahorro no se traduce en tiempo**:
sembrar `2n` bases cuesta casi exactamente lo que ahorra. Medido sobre 3
relajaciones y 4 semillas: 6.8–12.3 % más lento, celdas iguales o peores, y las
12 combinaciones dan el mismo signo. Hoy hay **un solo `DFBSimplex`
secuencial**: la base de la cota anterior es el arranque tibio de la siguiente.

---

## 3. Selección de cota

### La regla por omisión: salteo de Achterberg + orden por costo dual

- **Salteo de Achterberg**: se da por terminada la cota cuyo `x*_j` ya está en el
  borde. Es exacto: si `x*_j = lb_j`, entonces `min x_j = lb_j`.
- **Orden por infactibilidades duales** (`DFBH_ORDEN=costo`). Para `x_j` básica
  en la fila `r`, los costos reducidos del objetivo `±e_j` sobre las no básicas
  son `d_c = −ck·ā_rc`; se cuentan las de signo violado. Es el número de columnas
  que el simplex primal tendría que hacer entrar, o sea una estimación directa
  del costo de esa cota, `O(nnz(fila r))` sobre el tableau explícito y sin
  pivotear (`DFBSimplex::dual_infeasibilities`). Se elige el conteo menor;
  desempate por la más lejana.

Contra el default anterior (Achterberg + «lejana»), 137 instancias, 4 semillas:

| corte | celdas | cpu | +rápido / −lento | LPs/nodo | piv/LP | cpu/LP |
|---|---|---|---|---|---|---|
| > 0.1 s | ×0.974 | **×0.903** | **128 / 21** | ×1.028 | **×0.777** | ×0.922 |
| > 1 s | ×0.971 | **×0.924** | 41 / 9 | | | |
| > 5 s | ×0.979 | ×0.960 | 17 / 7 | | | |

130 instancias en 4/4 semillas, gana `ship-1`, 0 óptimos incompatibles. Pivotes
por LP, mediana por instancia: 3.35 → 2.39.

**El mecanismo.** Con pivoteo primal la localidad del vértice pesa: el orden
«cerca» divide por dos los pivotes por LP contra «lejos» (10.9 → 5.4 en
`ex5_3_2`, 4.7 → 2.7 en `alkylbis`, 4.1 → 2.5 en `house`). Que aun así sea neutro
en cpu es porque los paga por otro lado: el vértice lejano deja más variables en
sus cotas y compra salteos; el cercano compra LPs baratos y resuelve más LPs por
nodo. La distancia `|x*_j − borde|` es un sustituto del costo; el conteo de
signos violados **es** el costo, y con él se ahorra el 22 % de los pivotes
conservando casi todos los salteos.

El ahorro no llega entero al tiempo (piv ×0.78, cpu/LP ×0.92): con 2.4 pivotes
por LP el pivoteo es la mitad del costo del simplex y `compute_basics`, que
`primal_step` rehace entero en cada pivote, es el 20–25 % (§9).

### Variantes medidas

Todas contra Achterberg + «lejana».

| variante | cpu (> 0.1 s / > 1 s / > 5 s) |
|---|---|
| **orden por costo dual** (`DFBH_ORDEN=costo`, **omisión**) | **×0.903 / ×0.924 / ×0.960** |
| hull + orden por costo | ×0.918 / ×0.943 / ×0.965 |
| hull solo (`DFBH_HULL=1`) | ×0.943 / ×0.956 / ×0.964 |
| orden más cercana (`DFBH_ORDEN=cerca`) | ×1.015 (banco) |
| orden natural (`DFBH_ORDEN=nat`) | ×1.038 (banco) |
| sin salteo de Achterberg (`DFBH_NOSKIP=1`) | ×1.133 (banco) |

El salteo es lo que más pesa. El orden importa cuando se ordena por **costo** y
no por distancia: «cerca» y «natural» son neutros porque lo que ahorran en
pivotes lo gastan en LPs.

### El hull de puntos primales, medido y no adoptado

Todo minimizador `x*` es factible, así que si `H` es la caja mínima que contiene
todos los `x*` vistos desde la última linealización, la contracción posible de la
cota `(j, min)` es a lo sumo `H.lb_j − X.lb_j`: una cota superior **exacta** de
la ganancia, `O(n)` por LP. Achterberg y el orden de la más lejana son el caso
particular con un solo punto. Se saltea si la ganancia relativa cae bajo `τ`
(`DFBH_HULLTAU`, omisión `1e-6`) y se ordena por esa misma magnitud; con la
ganancia absoluta las instancias de escalas dispares se desordenan (`alkylbis`,
70 → 700 celdas). El resultado **no depende de `τ`**: cuatro órdenes de magnitud
apenas lo mueven.

No quedó como default porque **no suma con el orden por costo**: el hull compra
celdas pagando LPs por nodo (×1.03–1.05) y el orden compra pivotes, y la
combinación queda entre los dos. El orden solo gana a las dos combinaciones en
todos los cortes.

---

## 4. Resultados contra producción

Producción = `--filtering=acidhc4` + `CtcPolytopeHull` + linearizador.

### El linearizador: afín es el nuevo default

Se incorporó el plugin `ibex-affine` (vendorizado en `affine/`, sin CMake ni
`make install`). `LinearizerAffine2` entra en la ranura `art` y
`LinearizerCompo(afín, xtaylor)` en `compo`. `DFBH_LR` elige el linearizador con
independencia del contractor; `DFBH_LR=xn` restaura xtaylor.

**Cambiar sólo el linearizador vale más que todo lo hecho sobre el contractor**:
con `CtcPolytopeHull` intacto da cpu ×0.294 en las instancias de más de un
segundo. Y la razón no es la esperada: la relajación afín no sólo es más
apretada, es **más chica** —de 34 filas a 23, pivotes por LP de 10.9 a 6.2,
costo por LP de 168 a 56 µs—. Contrae más y cuesta menos.

### DFB contra producción, misma relajación

137 instancias, 4 semillas, mismo binario, pares resueltos por ambos; cortes por
el tiempo de producción.

| | corte | celdas | cpu | +rápido / −lento |
|---|---|---|---|---|
| **DFB por omisión** | > 0.1 s (n=235) | ×0.974 | **×0.684** | **206 / 12** |
| | > 1 s (n=88) | ×0.942 | **×0.657** | 80 / 7 |
| | > 5 s (n=45) | ×0.924 | **×0.621** | **45 / 0** |
| con el orden anterior | > 1 s | ×0.960 | ×0.709 | 74 / 9 |

Resueltas en 4/4 semillas: producción 129, DFB 130, no pierde ninguna. 0 óptimos
incompatibles. La ganancia es costo por nodo, y con el orden por costo también
algo de poda en las pesadas (celdas ×0.92 en > 5 s), que antes era
indistinguible. Con `DFB_SX_PIVOTEO=dual` la ventaja contra producción baja a
cpu ×0.71 en > 1 s, que es el ×1.045 del §1 aplicado sobre esta tabla.

### `compo`

| | celdas | cpu |
|---|---|---|
| PolyHull + compo / PolyHull + afín | ×0.858 | ×1.110 |
| DFB + compo / PolyHull + afín | ×0.831 | ×0.922 |

Poda bastante más y no lo paga. DFB lo deja en empate donde PolyHull pierde
11 %, coherente con el régimen del §5.

### Sin punto fijo

`DFBH_NOFIX=1` saltea el `CtcFixPoint` del contractor lineal, en las dos ramas.

| | celdas | cpu |
|---|---|---|
| producción sin fix / con fix | ×1.389 | ×1.008 |
| DFB sin fix / con fix | ×1.380 | ×1.077 |
| DFB / prod, **ambos** con fix | ×0.985 | ×0.729 |
| DFB / prod, **ambos** sin fix | ×0.979 | ×0.780 |

El punto fijo conviene a los dos y **a DFB más**, consistente con que su ciclo de
LP es más barato. La ventaja de DFB no es un artefacto de una composición
sintonizada para `PolytopeHull`.

---

## 5. El régimen de DFB

La ventaja se sostiene con relajaciones de 23 y 34 filas y **desaparece con
54–81**. Sobre los mismos LPs, con `m` creciendo ×2.53:

| | SoPlex | DFB |
|---|---|---|
| µs por iteración | ×1.25 | ×2.36 |
| iteraciones por LP | ×1.30 | ×2.52 |

Dos mecanismos de peso parecido. El costo por pivote crece **proporcional a `m`**
porque `pivot` actualiza las `m` filas del tableau explícito, mientras SoPlex
trabaja desde la factorización LU; no es llenado, los no nulos por fila apenas
pasan de 40.8 a 46.4. Y la ventaja del arranque tibio en iteraciones cae de
2.60× a 1.34×.

**La certificación NO es el cuello de botella**: 2–3 % del tiempo y plana en `m`.

Esto define el terreno: DFB rinde con relajaciones **chicas**. La afín lo deja en
su mejor régimen; `compo` lo saca de él.

### De dónde viene la ventaja

No de pivotear menos: DFB hace **más** pivotes que SoPlex tibio (7.4 contra 3.1
en `ex7_3_4`). Viene de que el **ciclo completo del LP** —pivoteo más
certificación— es 3–4 veces más barato que resolver con SoPlex y certificar con
Neumaier–Shcherbina: 16 µs/LP contra 59 en `ex7_2_3`, 29 contra 116 en
`ex14_2_7`. Y DFB emite además 17–26 % menos LPs por celda, por el salteo y por
el corte temprano de la pasada.

---

## 6. Lo que NO funciona — no reintentar sin leer esto

Veinticuatro direcciones probadas y descartadas **con el banco completo**.

| idea | resultado |
|---|---|
| HC4 intercalado en la pasada, 4 disparadores | cpu ×1.020–1.034, celdas planas, hasta 299 instancias más lentas contra 86 |
| Cotas colaterales (`0 ∉ γ_j·z` para j≠k) | no contraen ni ahorran pivotes; como alimento del HC4 tampoco |
| Rayos alternos de `B⁻¹` (`DFBH_OTRORAYO=1`) | recuperan el **100 %** de los certificados perdidos y **cero** efecto en celdas; bajo pivoteo dual, idénticos, porque ahí no se pierden certificados |
| Cosecha de las filas de `B⁻¹` de las otras básicas (`DFBH_FILAS=1`) | contrae de verdad —11–13 % de las filas, 10–22 % del diámetro— y baja los LPs por nodo ×0.835, pero el `γ` por básica sube el costo por LP ×1.254: celdas ×1.02, cpu ×1.06–1.09, 16 más rápidas contra 98 |
| Cachear esos `γ` y reevaluarlos (`DFBH_FILASCADA=k`) | el ahorro se degrada con el período (18.6 → 22.5 → 24.3 → 27.3 LPs/nodo): `γ` no depende de la caja, así que un certificado viejo sigue **válido** pero no **fuerte** |
| Saltear el LP si la cota de cero pivotes ya basta (`DFBH_LIBRE=1`) | acierta en 30–84 % de los LPs que contraen y ahorra 1–3 % de los pivotes: iguala al óptimo justo cuando el LP iba a costar 0 o 1 pivote |
| Generación perezosa de filas | 34–56 % peor |
| Relajación adaptativa por nodo | aditiva sin descuento |
| Reuso de linealización (`DFBH_RELIN=1`) | catastrófico: 472 → 22 512 celdas |
| **Congelar `A` por nodo y actualizar sólo las cotas de `b`** (`DFBH_CONG=1`) | celdas ×1.20 en > 1 s y ×1.23 en > 5 s, cpu ×1.14, contra re-linealizar; el ahorro es sólo ×0.95 en pivotes por LP. Medido nodo a nodo: contracción ×0.85 y **pruebas de vacío ×0.66** |
| Afinar el umbral del punto fijo (`DFBH_RATIO`) | compra poda y la paga cara, en los dos contractores |
| Ciclo interno propio con re-linealización | gana ×0.97 en las pesadas, pierde 4–5 % en el banco |
| Apretar la parte `b` antes de certificar | neutro (×0.984 celdas, ×0.994 cpu) |
| Perturbar los costos nulos al reubicar (`DFB_SX_PERTURB`) | no acorta el camino: 23.0 contra 23.7 pivotes/LP en `ex5_3_2` |
| Ratio test de paso largo (`DFB_SX_BFRT=1`), remedido con la configuración actual | sigue catastrófico: `alkylbis` de 62 a 19 982 celdas y al tope de 60 s |
| Cambiar la regla de la **primera cota** del nodo (`DFBH_PRIMERA`) | la más ancha ×1.115 y la más angosta ×1.059 en cpu sobre 10 instancias; al azar parecía ×0.97 y el banco lo devuelve a ×0.99–1.00 en las dos estrategias de pivoteo. Sobre el mismo nodo ninguna regla contrae sistemáticamente más |
| Congelar `A` por nodo y sólo actualizar las cotas de `b` (`DFBH_CONG=1`) | celdas ×1.20 y cpu ×1.14 en > 1 s contra re-linealizar, con pivotes ×0.95. Nodo a nodo: contracción ×0.85 y pruebas de vacío ×0.66 |
| Re-linealizar sólo las restricciones cuyas variables se movieron (`DFBH_CONGVAR=τ`) | parecía cpu ×0.70, y era **incorrecto**: reusar `rango'_i` sin verificar contención da cotas inválidas, 38–53 óptimos incompatibles. Corregido, queda en cpu ×1.15 |
| Devex en la fila que sale (`DFB_SX_DEVEX=1`) | mixto: baja pivotes donde los LPs son caros (8.68 → 7.06 en `ex5_3_2`) y los sube donde son baratos. Sin umbral que lo salve |
| Presupuestos / truncamiento anytime (11 variantes, sobre el camino primal) | ninguna paga |
| **Corte anytime por brecha contra el techo del hull, sobre el pivoteo dual** (`DFBH_ANYTAU=τ`) | corta en 24–69 % de los LPs y ahorra pivotes de verdad (piv/LP ×0.84–0.85), pero celdas ×1.39–1.60 y cpu ×1.36–1.49, con 5 pares más rápidos contra 71 en > 1 s y 1 contra 46 en > 5 s |

### El patrón que las explica

**Adelantar trabajo dentro del nodo no mueve el árbol.** Los rayos alternos
recuperaban el 100 % de los certificados sin cambiar una celda; las colaterales
dan la misma cota unas cotas antes; el HC4 intercalado adelanta una propagación
que el punto fijo hace igual; la cosecha de filas saca un tercio de los LPs y
paga exactamente eso. Lo único que movió la aguja fue cambiar **qué** se calcula
—la relajación afín, el filtro de estabilidad, la regla de selección— no
**cuándo**.

Hay un corolario que aparece tres veces: **lo que se ahorra son los LPs
baratos**. La cota de cero pivotes iguala al óptimo justo cuando la base ya era
óptima; la cosecha elimina las cotas que Achterberg iba a saltear igual. La
calidad de la cota intermedia y el costo del LP están correlacionados, así que
cortar o saltear nunca cae del lado bueno.

### Qué aporta el punto fijo, descompuesto

El `CtcFixPoint` que envuelve al contractor lineal re-linealiza en cada vuelta.
Se puede separar cuánto de su valor viene de **re-linealizar** y cuánto de
simplemente **volver a acotar** sobre la caja ya apretada. Tres brazos con
pivoteo dual, 137 instancias, 4 semillas, celdas en el corte de > 1 s:

| | celdas |
|---|---|
| punto fijo re-linealizando | 1.00 |
| punto fijo con `A` congelada por nodo, sólo se actualizan las cotas | 1.20 |
| sin punto fijo | 1.53 |

**El 62 % del beneficio del punto fijo se consigue sin re-linealizar**, sólo
recalculando las cotas de `b` con una linealización afín fresca corregida por
la fila diferencia. El 38 % restante lo agrega cambiar las normales.

No sirve como default —congelar cuesta ×1.14 de cpu— pero acota el valor de
cada mitad, y dice que si alguna vez re-linealizar se volviera caro, congelar y
actualizar cotas es mucho mejor que resignar el punto fijo: ×0.78 celdas contra
no tenerlo, al mismo cpu.

La pieza que lo hace posible es `LinearizerAffine2::linearize_id`, que devuelve
fila, rango e identidad `(restricción, lado)` por cada restricción, incluidas
las redundantes que `linearize` descarta. Sin identidad estable no se puede
emparejar la fila vieja con la restricción nueva. Medido: 0 filas sin par en
seis instancias.

**El alcance del congelamiento es el nodo, y detectarlo importó mucho.** El
contractor no sabe dónde empieza un nodo, pero el `CtcFixPoint` vuelve a llamar
mientras `old_box.rel_distance(box) > 0.2`, así que una llamada que contrae
menos que ese umbral es la última de su nodo. Con un contador fijo de llamadas
en vez de esa señal, el mismo brazo medía ×1.4 en lugar de ×1.20.

### Sobre el anytime, con la garantía puesta

Las once variantes de truncamiento por presupuesto se midieron sobre el camino
primal, donde la cota intermedia es válida pero de calidad arbitraria, así que
el fracaso se explicaba por la falta de garantía. Con `DFB_SX_PIVOTEO=dual` la
cota **sí** sube monótonamente, y se probó el criterio que corresponde: cortar
por **brecha** contra el techo exacto que da el hull, no por presupuesto. O sea
el mismo criterio con que el hull saltea cotas enteras, aplicado por pivote.

Falla igual, y el modo de falla es el informativo. El corte hace lo que promete:
dispara en el 24–69 % de los LPs y baja los pivotes por LP un 15 %, con el costo
por nodo casi intacto (LPs/nodo ×1.07). **Todo el daño está en las celdas**,
×1.39 a ×1.60, y el árbol arrastra HC4, ACID, bisección y cota superior, que no
tienen nada que ver con el LP. En esta búsqueda la contracción vale mucho más
que los pivotes, y una pérdida pequeña por nodo se amplifica en el árbol.

No hay régimen que lo salve: los pivotes por LP del dual son 4.01 de media en
las siete instancias que ganan y 3.31 en las 42 que pierden, con contraejemplos
fuertes en las dos direcciones (`bearing` tiene 5.89 y pierde ×2.5, `like` tiene
2.05 y gana).

### Sobre el vacío

Es propiedad del **nodo**, no de la cota: si una cota llega a `OPTIMAL` el
poliedro no es vacío y **ningún** `λ` puede probar vacío después. Medido: 0 casos
de óptimo seguido de vacío. El simplex detecta la infactibilidad siempre y ya en
la primera cota; lo que varía es si el certificado en intervalos sobrevive.

DFB mata entre el 1 % y el 17 % de los nodos que ve, con apalancamiento enorme:
quitarle el test de vacío multiplica el árbol por 2.5 a 10.

---

## 7. Advertencias de método

### La actualización de cotas está verificada correcta

La duda razonable sobre la matriz congelada era si el mal resultado venía de un
error en el cálculo `b_i = (rango'_i + (A_i − a'_i)·X) ∩ A_i·X`. **No viene de
ahí.** El control es `DFBH_CONGSIEMPRE=1`: re-congelar en *cada* llamada, de modo
que `A_i = a'_i`, el término de corrección es exactamente cero y el camino
congelado tiene que reproducir el fresco. Con los conjuntos de filas igualados
(`DFBH_LINMISMO=1`) lo reproduce **bit a bit** en las 7 instancias × 4 semillas,
salvo en las que terminan por tiempo, donde las celdas dependen del reloj. La
degradación del congelamiento es antigüedad de las normales, no un defecto.

El control hizo falta porque `linearize` y `linearize_id` **no emitían las mismas
filas**, por dos diferencias independientes que casi se cancelaban en el conteo
de `m` y por eso pasaron desapercibidas:

| | `linearize` | `linearize_id` original |
|---|---|---|
| `LEQ`/`GEQ` redundante sobre la caja | la salta | la emite |
| `LT`/`GT` | sólo prueba de infactibilidad, sin fila | emite fila |
| `EQ` | **dos** filas unilateras | **una** fila bilátera |

En el sistema extendido hay exactamente una `EQ` por llamada —el objetivo
`y = f(x)`— y entre 0 y 2 redundantes, así que `m` coincidía dentro de una fila.
Hoy cada diferencia tiene su interruptor: `DFBH_LINRED` descarta redundantes y
`LT`/`GT`, `DFBH_LINEQ2` parte la `EQ` en dos, y `DFBH_LINMISMO` activa las dos.
Todos apagados por omisión; producción usa `linearize` y no los ve.

Medido sobre 4 semillas, ninguna de las dos formas de escribir el mismo poliedro
importa en el agregado: celdas ×1.011 colapsando la `EQ`, ×0.989 conservando las
redundantes, cpu ×1.004 en ambos casos.

### El acotamiento afín de `c` está verificado, y tenía un defecto de redondeo

El control de `CONGSIEMPRE` prueba el caso `dif = 0`, o sea que no dice nada del
término de corrección, que es justamente lo propio de la idea. Para eso está
`DFBH_CHKC=k`, que muestrea `k` puntos de la caja por fila y separa la
afirmación en sus dos partes:

- **identidad**, sin condición: `A_i.x − a'_i.x ∈ dif` para todo `x` de la caja;
- **implicación**, que es la que sostiene la cota:
  `a'_i.x ∈ rango'_i ⟹ A_i.x ∈ corr`.

La segunda se prueba contra `corr`, antes de intersectar con `b_acum`, porque
`b_acum` viene de cajas anteriores. El conjunto de puntos que la sonda condiciona
es un **superconjunto** del factible: todo `x` que satisface la restricción
original satisface `a'_i.x ∈ rango'_i`, así que pasar la prueba es más fuerte que
ser sano. Resultado sobre 6 instancias: **0 violaciones en 3.4 M de muestras**,
2.2 M de ellas condicionadas.

Al muestrear hay un cuidado: `A.x − a'.x` **no** se puede evaluar como
`Interval(A.x) − Interval(a'.x)`. Los dos comparten `x` y la cancelación infla el
resultado —la primera versión de la sonda reportó excesos de hasta `1.16e+18`,
todos falsos—. Hay que acumular término a término en precisión extendida.

La auditoría sí encontró un defecto real, ya corregido: `dif` se acumulaba como
`Interval(Ai[j] - ult_rows[i][j]) * box[j]`, con la **resta en punto flotante** y
el resultado ya redondeado envuelto como intervalo puntual. Perdía hasta un ulp
de la diferencia, lo que en un contractor riguroso no corresponde. Ahora va
`(Interval(Ai[j]) - Interval(ult_rows[i][j]))`. Los resultados no cambiaron en
las 5 instancias de diagnóstico.

Dos medidas más de la misma sonda, que confirman que el mecanismo opera como se
diseñó: **el 100 % de las filas recibe refresco de su `c` en cada llamada** —
ninguna se queda con el rango de una caja anterior— y el refresco aprieta de
verdad: el ancho de `corr` contra el de la cota trivial `A_i·X` queda entre
0.48 y 0.81 según la instancia.

O sea que la degradación del congelamiento no viene del cálculo de `c`. Viene de
que las normales congeladas `A_i` son peores direcciones que las frescas, y un
`c` más apretado no lo compensa.

### Cuánto aporta refrescar `c`, contra dejar el `c` original

La ablación propia de la idea es matriz congelada con el `c` **original** —el
rango capturado al congelar, intersectado cada vuelta con `A_i·X` y acumulado
mientras la caja se anide, que es lo que hace `DFBH_CONGBARATO=1`— contra el `c`
**refinado** por afín. Sobre 11 instancias × 4 semillas:

| | celdas | cpu |
|---|---|---|
| `c` original vs `c` refinado | ×1.269 | ×1.094 |

El refresco afín aporta, y se paga solo: dejar el `c` congelado cuesta 27 % más
celdas y 9 % más tiempo aunque ahorre toda la evaluación afín. El margen es mayor
en las instancias caras: `chembis` 165254 contra 87252, `ex6_2_10` 188472 contra
145584, `ex8_5_1bis` 2820 contra 1836.

**Cuidado con el estadístico.** La sonda `DFBH_CORIG=1` compara las dos cotas en
la misma caja y dice que son idénticas en el 70–98 % de las filas, con ancho medio
entre ×0.97 y ×1.01. Eso hace predecir un empate, y es falso. El ancho medio sobre
todas las filas está dominado por la mayoría que no cambia; lo que decide la
contracción es el 1–29 % donde el refresco sí aprieta, que es cuando esa cota
resulta ser la que ata. Para juzgar un cambio en la relajación, el promedio sobre
todas las filas no sirve.

**Y cuidado con el banco de diagnóstico.** Sobre las 6 instancias fáciles el
congelamiento con `c` refinado empataba con relinealizar (celdas ×1.006, cpu
×0.956, o sea parecía ser más rápido). Agregadas las 5 pesadas, el empate
desaparece: celdas ×1.143 y cpu ×1.045 **en contra** del congelamiento. La
conclusión que valía —relinealizar gana— sigue valiendo, y el episodio es otro
caso de la regla de incluir instancias caras antes de concluir.

### El ruido de trayectoria con una semilla llega a ±40 %

El mismo experimento da la advertencia más fuerte del documento. Los cuatro
brazos de arriba describen **el mismo poliedro** y sólo difieren en cómo está
escrito. Con la semilla 1, `ex14_2_7` mide 4128, 4808 y 2900 celdas según la
variante; con la semilla 2 el orden se da vuelta. `ex5_3_2` va de 154 a 196.

O sea que en esta familia una reescritura sin contenido matemático mueve las
celdas hasta ±40 % por instancia y semilla. Cualquier conclusión sacada de una
corrida sola es indistinguible de esa dispersión. La contracara: `ex7_3_5bis` da
4726 / 5668 / 5674 idéntico en las cuatro semillas, así que ahí la diferencia sí
es real. Antes de atribuir una diferencia a un mecanismo, comprobar en cuál de
los dos regímenes está la instancia.


Las correcciones de este trabajo tienen todas la misma causa: **generalizar desde
una muestra que no cubre la variación relevante**.

- **Controles chicos no filtran nada.** El HC4 intercalado daba ×0.57 en el
  control y ×1.03 en el banco. `bearing` invierte el orden entre brazos según la
  semilla y su base sola oscila ×3.0, así que no sirve como caso de diagnóstico. La cosecha de filas repitió el patrón con tres
  instancias y **cuatro semillas**: celdas ×0.89, ×0.98 y ×0.84 en el control
  contra ×1.02 en el banco, con dispersión por instancia de ×0.70 a ×1.50.
  Cuatro semillas no arreglan una muestra de tres instancias.
- **Los totales del banco** pesan una instancia de 28 000 celdas como 300 de 90.
  Usar media geométrica de razones, o ventaja acotada `(A−B)/max(A,B)`.
- **La caja raíz** no representa los nodos profundos, ni las instancias chicas a
  las grandes.
- **Al tocar algo, remedir lo que dependía de ello.** Dos decisiones se validaron
  una vez y se arrastraron mientras el código de alrededor cambiaba.
- **Reusar información entre cajas exige verificar contención.** Pasó tres
  veces: al reusar los `γ` dentro de ACID, al acumular las cotas de `b` con la
  matriz congelada, y al re-linealizar selectivamente. En la tercera se me pasó,
  y comparar diámetros en vez de contención produjo un resultado espectacular y
  falso, cpu ×0.70 con 38–53 óptimos incompatibles. Una caja puede encoger y a
  la vez desplazarse.
- **El perfilador propio distorsiona.** La instrumentación por fases con
  `clock_gettime` infla el tiempo total entre 22 % y 37 %, y ese sobrante se
  acumula en los ámbitos que envuelven a los demás. Eso creó un balde falso que
  parecía el segundo costo más grande del simplex. Antes de perseguir un rubro,
  comparar el tiempo con y sin instrumentación.
- **Una semilla no explica un mecanismo.** El pivoteo dual multiplicaba por diez
  el árbol de `ex7_2_3` en una corrida, y eso se explicó durante meses por la
  forma de sus rayos de Farkas. Era el loup: esa corrida no llegó al óptimo en
  60 s y sin cota superior el árbol no poda. Antes de explicar un efecto grande
  de una instancia, correr las otras semillas y mirar el loup final.

Para **diagnosticar** conviene lo contrario que para evaluar: el caso más chico
que exhiba el efecto y que **no** dependa de la semilla, elegido por el mayor
contraste entre configuraciones. Después hay que bajar hasta un nodo. La técnica
que sirve es correr las dos configuraciones **sobre la misma caja**, adoptando
una sola para seguir la búsqueda: comparar dos corridas completas no sirve,
porque divergen en el primer nodo distinto. Así se vio que el pivoteo dual
contraía igual en 4347 de 4400 nodos, lo que descartó la explicación agregada y
dejó un puñado de nodos profundos para volcar cota por cota.

---

## 8. Estado del código

Fuentes vivas en `examples/dfb/`:

```
ibex_CtcDFBHull.{h,cpp}     el contractor
ibex_DFBSimplex.{h,cpp}     el simplex
ibex_DFBBasis.{h,cpp}       la factorización
ibex_Optimizer05Config.cpp  el cableado de --lr y DFBH_LR
affine/                     el plugin ibex-affine vendorizado
```

Archivadas pero **no borrar** (hay intención de retomar una variante parcial):
`ibex_CtcDualFeasibleBounding.cpp`, `ibex_CtcDFBPropag.cpp`, modo `--lr=dfb`.
Arneses obsoletos en `obsoleto/`. Los datos de los barridos con 4 semillas están
en `results/*_4semillas.csv`, una fila por corrida, con su `results/analiza_*.py`
al lado.

### Interruptores vivos

| variable | |
|---|---|
| `DFBH_LR=xn\|art\|compo` | linearizador, independiente del contractor. Omisión `art` |
| `DFBH_ORDEN=costo\|lejos\|cerca\|nat` | orden de la próxima cota. Omisión `costo` |
| `DFB_SX_PIVOTEO=primal\|dual` | re-optimización al cambiar de objetivo (§1). Omisión `primal` = «dual en la cota 0, primal de la 1 en adelante»; `dual` es dual en todas |
| `DFBH_HULL=1`, `DFBH_HULLTAU=τ` | selección por hull (§3). No suma con el orden por costo |
| `DFBH_NOFIX=1` | sin `CtcFixPoint` en el contractor lineal |
| `DFB_SX_REUBK=1` | con `PIVOTEO=dual`, vuelve a reubicar la columna del objetivo. Es el defecto del §1 |
| `DFB_SX_NOXBSKIP=1` | recalcula `x_B` en `primal_step` aunque esté limpio, o sea el comportamiento previo al arreglo del §2 |
| `DFBH_PRIMERA=0\|ancha\|angosta\|azar` | qué cota se resuelve primero en cada nodo. Omisión `0`. Medido neutro (§6) |
| `DFBH_ABPRIM=R` | con `DFBH_ABDUAL`, compara dos reglas de primera cota sobre el mismo nodo en vez de dos estrategias de pivoteo |
| `DFBH_GOALPROBE=1` | observa en qué pivote la cota del objetivo cruza el loup, sin cortar |
| `DFB_SX_DEVEX=1`, `DFB_SX_BFRT=1`, `DFB_SX_PERTURB=eps` | variantes del simplex medidas en el §6 |
| `DFBH_FILAS=1`, `DFBH_FILASCADA=k`, `DFBH_LIBRE=1` | cosecha de certificados, medida en el §6 |
| `DFBH_ANYTAU=τ` | con `PIVOTEO=dual`, corte anytime por brecha contra el techo del hull. Medido, no paga (§6) |
| `DFBH_CONG=1`, `DFBH_CONGVAR=τ` | congela `A` por nodo y actualiza las cotas de `b`; con `τ`, re-linealiza sólo lo que se movió (§6). `DFBH_ABCONG=1` compara las dos contracciones sobre la misma caja |
| `DFBH_LINSLOW=1`, `DFBH_AFNOPOOL=1` | restauran el linearizador previo a las dos optimizaciones del §2 |
| `DFBH_PERFIL=1`, `DFBH_LINPERF=1`, `DFBH_AFALLOC=1` | reparto del contractor, interior de la linealización, y conteo de asignaciones de la afín |
| `DFBH_CRONO=1`, `DFB_SX_PERF=1`, `IBEX_PH_STATS=1` | reparto de costo por nodo, perfil del simplex, y estadísticas de `CtcPolytopeHull` |
| `DFB_SX_VERIF=1`, `DFB_SX_FEASTOL=t` | verifica cada `OPTIMAL` desde la factorización y cuenta los falsos; y fija la tolerancia de factibilidad primal, omisión `1e-9` |
| `DFBH_SINCOTA=1`, `DFBH_PREDICE=1` | tasa de cotas que terminan sin certificado, y conteo predicho por la heurística de orden contra los pivotes reales |
| `DFBH_ABDUAL=τ`, `DFBH_ABNODO=n` | corre las dos estrategias de pivoteo sobre el mismo nodo y vuelca las diferencias; con `ABNODO`, cota por cota |

Queda una veintena de interruptores de líneas descartadas (`DFBH_HC4*`,
`DFBH_COLAT`, `DFBH_LAZY`, `DFBH_ADAPT`, `DFBH_TECHO`, sondas). Cada uno tiene el
resultado que lo condena anotado al lado. Se pueden limpiar en una pasada
dedicada, pero varias sondas envuelven código que hay que conservar, así que no
sirve un borrado por rangos.

### Arreglos en el núcleo de Ibex

- `CtcPolytopeHull`: `n_soplex_iterations` y `n_soplex_calls` **no estaban
  inicializados** en ninguno de los dos constructores. Arreglado.
- Instrumentación de costo por LP en `LPSolver::minimize`, tras `IBEX_PH_STATS`.
  El reparto por rubros **no es de fiar** (acumuladores globales, suma LPs de
  otros usuarios); los agregados sí.

---

## 9. Lo que queda abierto

- **El costo fijo por LP.** Con el LP en 2.4 pivotes de mediana, el pivoteo es
  la mitad del costo del simplex y `compute_basics` el 15 %, ya sin la llamada
  duplicada (§2). Bajarlo más pide la actualización incremental, que está medida
  y **no paga**: acumula redondeo distinto que el recálculo y eso mueve las
  decisiones del árbol (cpu ×1.018, celdas distintas en 25 de 128). Por ahí no
  hay margen sin repensar la aritmética.
- **Qué hacer con el pivoteo dual.** Su propiedad propia ya se probó y no paga
  (§6), así que hoy `dual` es un empate a ×1.03–1.05 de cpu con celdas iguales,
  defendible como propuesta pero no por rendimiento. Lo que queda por explorar
  es usar la cota monótona para algo que **no** sea resignar contracción.
- **El pivoteo sobre el tableau explícito**: 49 % del tiempo del simplex con `m`
  chico y 66 % con `m` grande. Es el eje que limita a DFB fuera de su régimen. El
  camino sin tableau que existió (`DFB_SX_EXACT`) resultó **peor** al medirlo, así
  que pide rediseño real (BTRAN/FTRAN con actualización tipo Forrest–Tomlin).
- **La degradación del arranque tibio** con `m` grande. Pesa casi lo mismo que el
  pivoteo y no se arregla reimplementando: pide repensar qué se hereda.
- **`compo`** poda 33 % más y cuesta 3 % más con DFB. Abaratar el LP lo volvería
  la mejor opción.
- **Por qué el hull paga LPs por nodo** (×1.03–1.05), y si con `τ` más grande se
  suma al orden por costo en vez de competir con él.
- **Una variante parcial de la propagación propia de DFB.** Se archivó cuando la
  relajación tenía 34 filas y el defecto del pivote fallido estaba activo; ambas
  condiciones cambiaron.
