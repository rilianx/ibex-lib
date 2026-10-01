# Mediciones de la sesión del 2026-09-29/30

Registro textual de lo medido en esa sesión, extraído de `ESTADO.md` el
2026-10-01 al reducirlo. `ESTADO.md` conserva las conclusiones; acá están los
números completos, las hipótesis descartadas y los mecanismos. Las referencias a
secciones (§N) son las de `ESTADO.md` de entonces, respaldado en
`respaldo/ESTADO_2026-10-01_largo.md`.

---

## A. `bearing` con `eps_x` = 1e-7 (análisis completo)

**Mecanismo de `bearing`** (traza del optimizador, `IBEX_TRACE=1`, y
`IBEX_LOUPSTATS=1`, que cuenta intentos y éxitos del buscador de `loup` por
ancho de caja). La cota inferior de DFB es buena (1.9510660 contra 1.9510658 de
producción); lo que falta es el `loup`. `bearing` tiene 8 igualdades, que el
sistema normalizado relaja a `|h| ≤ 1e-8`, y el buscador **sólo tiene éxito en
cajas de ancho 1e-4 a 1e-3**: producción hace 1932 intentos ahí (2 éxitos) y
DFB 402 (1 éxito). En cajas más angostas que 1e-5, ninguno de los dos encuentra
nada, y ahí es donde DFB pasa el tiempo: 18 558 intentos en cajas de menos de
1e-6, contra 145 de producción. La contracción más fuerte lleva las cajas cerca
del óptimo de anchas a menos de 1e-6 sin detenerse en el rango donde el
buscador funciona, y después hay que bisecarlas hasta 1e-7 para descartarlas
por diminutas. Con `eps_x` = 1e-6 esas cajas se descartaban de inmediato. Es
una limitación del buscador de `loup` con igualdades en cajas angostas, que la
contracción de DFB deja a la vista; no un defecto del contractor.

Qué hay en esas cajas (`IBEX_TINYDUMP=1` vuelca las cajas diminutas con cota
inferior bajo `loup − 1e-6`; evaluadas con un verificador aparte): en las tres
que DFB descarta junto al óptimo, **las igualdades 3, 5 y 10 valen ≈ ±1e-8
sobre toda la caja** (la 3 en [9.9999e-9, 1.00001e-8]): no hay solución exacta
adentro, sólo del problema relajado, **en el borde de la holgura** `eqeps`.
Las otras igualdades contienen al 0 pero varían hasta 1e-4 en la caja, y en el
punto medio valen 1e-5 a 1e-7. O sea que la cota inferior 1.9510660 de esas
cajas se obtiene aprovechando la holgura de 1e-8 en tres igualdades, y un punto
que sirva de `loup` tendría que resolver las ocho a 1e-8 dentro de una caja de
1e-7. Un buscador por puntos «interiores» no puede.

**Corrección del análisis.** Bajar `eps_x` pide más resolución, así que se
descartan menos cajas, no más: el problema no es perder soluciones. Tampoco es
el orden de exploración: tras el `loup` de DFB, `y ≤ loup − eps` poda todo lo
demás y el 100 % de las 19 030 celdas siguientes está en la zona del óptimo
(cota a menos de 2e-6 del `uplo` final; `IBEX_CELLTRACE=1` traza cada celda).
El contraste que decide: con `eps_x` = 1e-6 el `loup` que cierra
(1.95106600) sale de una caja de **5.95e-7**, un éxito en 470 intentos en cajas
< 1e-6; con 1e-7, 0 éxitos en 18 633 intentos en cajas de ese tamaño o
menores. Hipótesis, sin verificar: con 1e-6 esas cajas se prueban una vez
recién salidas de la contracción y se descartan; con 1e-7 se siguen bisecando y
contrayendo, y la contracción las empuja contra el borde de la holgura, donde
las igualdades 3, 5 y 10 quedan en ±1e-8 y no queda margen para un punto
factible. 0 óptimos incompatibles. Bajar `eps_x` a 1e-7, por sí
solo: producción celdas ×0.99 / ×0.96 / ×0.93 y resueltas 130 → 131; DFB
celdas ×0.995 / ×0.978 / ×0.960.

Sin `bearing` los cortes no cambian (sus pares ya quedaban fuera por llegar DFB
al tope); resueltas en 4/4: producción 130, DFB 131, y DFB no pierde ninguna.

> **PROBLEMA ABIERTO — `bearing` con `eps_x` = 1e-7.** DFB no encuentra un
> `loup` suficiente en la zona del óptimo. Análisis y hipótesis sin verificar
> más abajo en esta sección. Dejado ahí por ahora (2026-09-30).

Las tablas que siguen son anteriores (`eps_x` = 1e-6, otra omisión de DFB).

## B. Tablas contra producción superadas por la línea de base del §4

Producción = `--filtering=acidhc4` + `CtcPolytopeHull` + linearizador.

**El linearizador afín es el nuevo default** (`ibex-affine` vendorizado en
`affine/`; `DFBH_LR=xn` restaura xtaylor). Cambiar sólo el linearizador vale
más que todo lo hecho sobre el contractor: con `CtcPolytopeHull` intacto da cpu
×0.294 en > 1 s. La relajación afín es más apretada **y más chica**: 34 → 23
filas, pivotes por LP 10.9 → 6.2, costo por LP 168 → 56 µs.

**DFB contra producción, misma relajación**, 137 instancias, 4 semillas, cortes
por el tiempo de producción:

| corte | celdas | cpu | +rápido / −lento |
|---|---|---|---|
| > 0.1 s (n=235) | ×0.974 | **×0.684** | **206 / 12** |
| > 1 s (n=88) | ×0.942 | **×0.657** | 80 / 7 |
| > 5 s (n=45) | ×0.924 | **×0.621** | **45 / 0** |

Resueltas en 4/4 semillas: producción 129, DFB 130, no pierde ninguna. 0
óptimos incompatibles.

**Con el tableau con guardia** (omisión desde 2026-09-30), mismo banco (138),
4 semillas, producción corrida en el mismo día
(`results/lu_tableau_4semillas.csv`, `analiza_lu.py <dir> prod tab`):

| corte | celdas | cpu | +rápido / −lento |
|---|---|---|---|
| > 0.1 s (n=244) | ×0.972 | **×0.573** | **228 / 2** |
| > 1 s (n=90) | ×0.951 | **×0.558** | 86 / 2 |
| > 5 s (n=45) | ×0.932 | **×0.536** | **45 / 0** |

cpu por celda ×0.571 (> 1 s); en las 26 que llegan al tope en los dos, ×0.545
(`chem*` hacen el doble de celdas en los 60 s). Resueltas en 4/4: producción
130, DFB 131; DFB resuelve además `ex8_4_4bis` y `ship-1` en alguna semilla que
producción no, y no pierde ninguna. 0 óptimos incompatibles. Los cortes de
arriba son pares resueltos por los dos, así que el n no coincide con la tabla
anterior. Con `DFB_SX_PIVOTEO=dual` la ventaja baja a cpu ×0.71
en > 1 s.

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

## 9. Lo que queda: el perfil de 2026-09-29 y las opciones

### El hallazgo

Perfil con `perf` y pila de llamadas sobre `ex6_1_3`, configuración por omisión,
una semilla, 30 s. Una instancia, así que son órdenes de magnitud:

| componente | tiempo inclusivo |
|---|---|
| `CtcDFBHull::contract` | 49 % |
| `DFBSimplex::pivot` | 23 % |
| de ese pivote, `DFBBasis::change_basis` (Forrest–Tomlin) | 10 % |
| `LPSolver::minimize` con SoPlex, **fuera** de DFB | 5.6 % |
| `LoupFinderXTaylor::find` | 4.2 % |
| `DFBBasis::btran_row` + `ftran_vec` | 3 % |

**DFB paga dos representaciones en cada pivote.** El tableau explícito se
actualiza con `axpy`, y además la LU se actualiza con Forrest–Tomlin en cada
`change_basis`. Esa segunda actualización se lleva algo más del 40 % del costo
del pivote, y sus únicos consumidores en el camino por omisión son el `λ` final
(`y_desde_factorizacion`, un BTRAN por LP) y el refresco al cambiar de cota
(`y`, `d` y `x_B` con `A_times` + FTRAN, una vez por LP; `compute_basics` lee el
tableau, no la LU), que juntos cuestan una fracción de lo que cuesta
mantenerla. Con 2.4 pivotes por LP de mediana, es probable que salga más barato
refactorizar una sola vez al final del LP, cuando hace falta `λ`, y obtener
`x_B` desde las filas del tableau que ya se tienen. Techo: en torno al 10 % del
tiempo total en esa instancia. Un `λ` refactorizado no es bit a bit igual al
actualizado, así que se valida con banco y 4 semillas, no con identidad.

**DFB depende hoy de SoPlex por esta factorización, no por su simplex.**
`DFBBasis` envuelve `soplex::SLUFactor`: `load` para factorizar, `solveLeft` y
`solveRight` como BTRAN y FTRAN, y `change` como actualización tipo
Forrest–Tomlin, con refactorización cada 50 actualizaciones
(`DFB_REFACTOR_PERIOD`). No se usa `SoPlex::solve` ni `LPSolver`. La razón de
sacar `λ` de la LU y no de la parte `λ` de la fila objetivo del tableau es que
ésta acumula el redondeo de todos los pivotes; afecta la calidad del
certificado, no su validez. Es medible: `λ` desde el tableau contra `λ` desde
la LU, en cotas sin certificado y en celdas. Si la diferencia es nula, la
dependencia sobra.

### La vía 1, medida (2026-09-30): LU diferida

`DFB_SX_LU=diferida` (fue omisión el 2026-09-30 hasta que la reemplazó el tableau; `actualizada` restaura el update por pivote): `pivot` ya no hace
`change_basis`, marca la LU sucia, y `asegurar_fact` la refactoriza desde
`basic[]` cuando alguien la lee —el `λ` final, el refresco de la cota siguiente,
las sondas—. Con base compartida, el `λ` final y el refresco leen la misma
base: **0.77 factorizaciones por LP** en vez de un update por pivote. Una
factorización fallida cae al respaldo del tableau, y a diferencia del update se
recupera en la cota siguiente. Test unitario idéntico en los dos modos.

Banco `easy+medium+hard`, 138 instancias, 4 semillas, 8 cores
(`results/lu_diferida_4semillas.csv`, `results/analiza_lu.py`):

| corte | celdas | cpu | +rápido / −lento | cpu por celda |
|---|---|---|---|---|
| > 0.1 s (n=211) | ×0.994 | **×0.950** | **106 / 28** | |
| > 1 s (n=79) | ×1.026 | ×0.975 | 35 / 18 | **×0.939** (n=107, con las del tope) |
| > 5 s (n=43) | ×1.059 | ×1.008 | 17 / 11 | ×0.934 |
| ambas al tope de 60 s (n=27) | | | | **×0.905** |

Resueltas en 4/4: 131 y 131. 0 óptimos incompatibles. `ship-1` la resuelve el
control en 1 de 4 semillas y el diferido en ninguna: es la instancia que
oscila ×5.2 por semilla, y decide el `loup`.

**Las celdas son caos de redondeo, no efecto.** Control nulo: el camino actual
con `DFB_REFACTOR_PERIOD=49` en vez de 50, que no cambia ningún mecanismo y
sólo perturba `λ` en los últimos bits, sobre las instancias de > 1 s:

| brazo | celdas > 1 s | celdas > 5 s | cpu por celda |
|---|---|---|---|
| período 49 (nulo) | ×1.051 | ×1.078 | ×1.007 |
| LU diferida | ×1.026 | ×1.059 | **×0.939** |

El nulo mueve el árbol **más** que el cambio. `himmel16`, insensible a la
semilla, da 4760 celdas con período 50, 5302 con 49, 8852 diferida, y entre
4760 y 8514 barriendo el período (20, 35, 100, 1000, siempre); `bearing` 7492
/ 17 574 / 5026. El default actual está parado en una configuración
afortunada de esas dos, que son las que dominan el corte > 5 s. Lo único
sistemático es el costo por celda, ×0.94 (×0.91 en `ex6_1_3`, la del perfil,
con las mismas celdas).

**Quedó por omisión** y después la reemplazó el tableau (abajo). Es el primer cambio del §2 que no es idéntico bit a bit,
así que el criterio no fue la identidad: con el nulo como referencia, la parte
de árbol queda dentro del ruido y la de costo fuera. Las cifras de §3–§5 son
anteriores y se midieron con `actualizada`.

### La vía 2, antes de implementarla: techo y deriva (2026-09-30)

**Techo.** `perf` con pila DWARF sobre `ex6_1_3`, LU diferida: lo que la LU
cuesta todavía dentro de DFB es `factorize` 4.1 % + `btran_row` 1.6 % +
`ftran_vec` y `A_times` ~1.1 % ≈ **6.8 % del tiempo total**. El `d` del
refresco no cuenta: el camino del tableau también lo calcula. Hoy `pivot`
inclusivo es 15 % y `axpy` 12.9 %.

**Deriva** (`DFB_SX_DERIVA=1`, no cambia la trayectoria: sigue usando la LU y
en los mismos puntos calcula lo que daría el tableau). 28 instancias de > 1 s
más `ex6_1_3`, semilla 1, 30 s:

| | resultado |
|---|---|
| LPs sondeados | 4.93 M |
| cota certificada con `λ` del tableau peor por > 1e-6 del ancho | 138 (0.003 %) |
| peor por > 1e-3 del ancho | 16 |
| mejor por > 1e-6 | 24 |
| vacíos probados, LU / tableau | 55 147 / 55 150 |
| refrescos con algún veredicto distinto (signo de `d` o factibilidad de `x_B`) | 1951 de 4.21 M (0.05 %) |

En 21 de las 28 la deriva es nula a efectos prácticos (`λ` con diferencia
relativa media 1e-12–1e-15, ningún veredicto distinto). **Se concentra en
cuatro**: `ex5_3_2` (5.2 % de los refrescos con veredicto distinto, `|Δd|`
hasta 2.7e11), `himmel16` (0.6 %), `bearing` (0.2 %) y `ship-1` (0.1 %), que
son también las más caóticas del banco. En `ex7_3_5bis` y `ex7_3_4` `λ` difiere
mucho (hasta 9e7 relativo) pero sin perder ninguna cota: pasa donde la cota no
contrae de todos modos.

Artefacto a evitar al implementar: con `k` no básica la LU da `y = 0` exacto y
el tableau un residuo, que deja `γ_k` diminuto y no nulo y una «cota» enorme e
inútil. Hay que forzar `y = 0` como hace `y_desde_factorizacion`.

Consecuencia: la vía 2 pura degradaría justo las instancias mal condicionadas.
La forma sensata es **tableau con guardia**: tomar `λ`, `d` y `x_B` del
tableau, verificar la consistencia de la fila objetivo (`d` contra
`c − Āᵀλ`, `O(nnz)`, lo mismo que ya cuesta el refresco) y refactorizar sólo
cuando falla. Los sistemas de `data_tests` (Brown, Brent) no sirven para esto:
`ibexopt` se cae con ellos en cualquier configuración porque no tienen objetivo.

### La vía 2, medida (2026-09-30): tableau con guardia

`DFB_SX_LU=tableau` (**omisión desde 2026-09-30**): `λ`, `d` y `x_B` salen del tableau. Una guardia de
`O(nnz)` verifica la fila objetivo (`d` contra `c − Āᵀλ`: residuo relativo
bajo `DFB_SX_GUARDIA`, omisión `1e-9`, y ningún veredicto de signo distinto en
las no básicas) y el residuo de `Āz = 0` para `x_B`; si falla, se refactoriza y
se sigue por la LU. Con `k` no básica, `y = 0` exacto. El arranque en frío no
factoriza. Test unitario idéntico.

La guardia salta donde la sonda de deriva decía: 0 fallos en 1.1 M de refrescos
de `ex8_4_4bis`, 5 en `ex6_1_3`, ~1 % en `himmel16` y 8–10 % en `ex5_3_2`.

Banco, 138 instancias, 4 semillas, 8 cores, contra la LU diferida, con los dos
brazos intercalados en el mismo barrido (`results/lu_tableau_4semillas.csv`,
`analiza_lu.py <dir> dif tab`):

| corte | celdas | cpu | +rápido / −lento |
|---|---|---|---|
| > 0.1 s (n=203) | ×0.991 | **×0.894** | **147 / 14** |
| > 1 s (n=78) | ×0.972 | **×0.895** | **62 / 4** |
| > 5 s (n=40) | ×0.956 | **×0.877** | **33 / 4** |

cpu por celda ×0.918 (> 1 s, n=106) y ×0.917 en las 26 que llegan al tope en
los dos. Resueltas en 4/4: 131 y 131; `ship-1` la resuelve el tableau en alguna
semilla y la diferida en ninguna. 0 óptimos incompatibles. Las celdas se mueven
con el mismo caos de redondeo que midió el control nulo, ahora a favor
(`himmel16` ×0.68, `ex14_2_7` ×0.80) y en contra (`bearing` ×1.30).
Reproducibilidad: el brazo diferido da las mismas celdas que el barrido del día
anterior en todas las corridas que terminan; las 28 distintas llegan al tope.

Contra la LU actualizada que era omisión hasta ayer, las dos vías juntas dan
del orden de cpu por celda ×0.86.

### La cota lagrangiana con monotonía (2026-09-30): mecanismo real, árbol indiferente

La idea: los multiplicadores del LP valen como multiplicadores de Lagrange de
las restricciones **originales**. Del certificado, `z_k = N/γ_k` con
`N = −Σ_{i≠k} γ_i x_i − Σ_j γ_{b_j} b_j` y `b_j = a_j·x`. En las filas activas
—lado de la restricción y signo del multiplicador tal que
`(−γ_{b_j}/γ_k)·g_c(x) ≥ 0` en todo factible— se sustituye
`b_j = g_c(x) − ρ_j(x)`, con `ρ_j = g_c − a_j·x` el resto no lineal, y se
descarta el término en `g_c`. Queda `Φ(x) = N′(x)/γ_k ≤ z_k`, función escalar
de `x`, y `min_caja Φ` es cota válida.

**Evaluada en afín sobre la misma caja, `Φ` es la cota del LP**, por
construcción: la forma afín de `ρ_j` es exactamente el intervalo que el LP usó
para `b_j`. Medido (`DFBH_LAGR=1`, 14 instancias, 4.9 M de LPs del objetivo):
`|b − a|` relativo ≤ 3.5e-8, casi siempre 1e-14. Lo único que el LP no puede
hacer es **monotonía**: si `∂Φ/∂x_i` tiene signo constante, el mínimo está en
una cara, y la forma afín de `ρ_j` sobre la cara tiene menos error. Se fija
alguna variable en el 80–100 % de los LPs del objetivo (1.5–6 por cota).

**La cota mejora; el árbol no lo nota.** Contracción rigurosa
(`DFBH_LAGRC`, intervalos en todo, óptimos compatibles en las 48 corridas), 8
instancias, 2 semillas, contra la omisión (`results/lagrangiana_sueltas.csv`):

| modo | dónde | celdas | cpu | mecanismo |
|---|---|---|---|---|
| 1: cota inferior del objetivo | mejora la del LP en el 50–91 % de las aplicaciones; en `ex8_4_4bis` mata el nodo en 30 826 de 143 012 LPs donde el LP no lo mataba | **iguales** en las deterministas (`ex6_1_3` ×0.99, `ex2_1_9` ×0.98, `house` ×1.00); caos en las caóticas | ×1.08–1.24 | los nodos que mata mueren igual un momento después: el punto fijo re-linealiza sobre la caja contraída. Es el §6: adelantar |
| 2: las `2n` cotas | mejora en el 60–97 % | **×0.81–0.96** en las deterministas (`ex8_4_4bis` 115 836 → 94 054 al tope, `ex2_1_9` ×0.956, `ex6_1_3` ×0.96, `alkylbis` ×0.94) | **×1.6–3.5** | contracción nueva sobre las `x`, que ninguna relajación lineal da; pero cada cota cuesta `m` gradientes por intervalos y una forma afín por fila activa, ~70 µs, o sea 2–4 LPs |

El piso de costo de la variante 2, con gradientes una vez por pasada y una
ronda, sigue siendo del orden de una forma afín por fila activa y por cota:
`~0.4·n` linealizaciones por nodo contra una hoy. Con `m = 23` eso pesa más que
el 4–19 % de celdas que compra. Queda apagada; el código se conserva porque el
mecanismo es correcto y podría valer donde los LPs sean caros (§5). Es la
tercera idea dual que se cierra con la misma explicación que las colaterales:
**en un contractor que resuelve las `2n` cotas al óptimo y itera a punto fijo,
la información dual del nodo ya está implicada por lo que se calcula**, y lo
que agrega tiene que pagarse a precio de linealización.

Esto toca una conclusión anterior: el camino sin tableau (`DFB_SX_EXACT`) se
midió peor con esta contabilidad, donde la LU ya se pagaba entera. Cualquier
remedición del rediseño para `m` grande va después de ordenar esto.

### Las opciones, en el orden en que conviene tomarlas

0. **Pendiente marcado: `bearing` con `eps_x` = 1e-7** (§4). Verificación
   propuesta: llamar al buscador de `loup` sobre cada caja antes y después de
   contraerla, y ver si acierta antes y falla después.

1. **Quitar la doble representación por pivote.** Dos vías: conservar
   `SLUFactor` pero factorizar una sola vez por LP —**hecha, medida arriba:
   cpu por celda ×0.94**, ya por omisión—; o sacar `λ` y `x_B` del
   tableau explícito y eliminar `DFBBasis` del camino por omisión, que es la
   única que deja a DFB sin SoPlex en el caso común —**hecha como tableau
   con guardia: cpu ×0.89 contra la diferida**, por omisión desde
   2026-09-30—.

2. **Sacar SoPlex del resto de `ibexopt`.** Con DFB en el contractor, SoPlex
   sigue resolviendo un LP por nodo en `LoupFinderXTaylor` y otro en `LSmear`
   (bisector por omisión `lsmearmg`). El techo de tiempo es chico, ~5 %, pero el
   valor es otro: un optimizador certificado sin solver LP externo, lo que
   exige además la segunda vía del punto 1. `LSmear` no
   necesita ni un LP nuevo: quiere los multiplicadores duales del objetivo, que
   DFB ya calculó en la cota del nodo. Cambia la heurística de bisección, así que
   hay que medir celdas.

3. **El ratio test de paso largo: el defecto era real y está arreglado
   (2026-09-30), medido abajo.** Con `DFB_SX_BFRT=1`, el **100 % de las
   resoluciones** de `ex7_2_6` y `ex9_2_5` terminaba en `INFEASIBLE` (en las dos
   estrategias), ningún LP contraía y el árbol explotaba (54 → 308 376 celdas
   al tope). `DFB_SX_VERIF` no lo veía porque sólo verifica los `OPTIMAL`. La
   causa: cuando ni cambiando de cota todas las candidatas alcanzaba para cubrir
   la violación `δ`, el test declaraba infactibilidad primal. Eso vale en
   aritmética exacta, pero no con tolerancias: las columnas con `|a|` bajo el
   umbral de pivote se excluyen de las candidatas y aun así mueven `x_B(r)`, y
   las violaciones del caso eran del orden de 1e-6. Ahora se pivotea en el
   último quiebre, cambiando de cota los anteriores, como la regla normal y las
   implementaciones robustas; la infactibilidad, si la hay, la prueba Farkas.
   `DFB_SX_BFRTINF=1` restaura el defecto.

   Con el arreglo (2 semillas, las dos estrategias, óptimos verificados 0 %
   infactibles): `ex9_2_5` 6/8 celdas como sin paso largo; `ex7_2_6` 54/64
   contra 54/68, piv/LP 1.52 → 1.24; `alkylbis` 54 → 60 en primal y 68/60 → 58
   en dual, piv/LP ×0.84–0.90. Test unitario 528/528 con paso largo.

   **Banco** (138 instancias, 4 semillas, 8 cores, tres brazos en el mismo
   barrido; `results/paso_largo_4semillas.csv`, `analiza_lu.py <dir> pri dbf`):

   | contra el primal (omisión) | celdas > 0.1 / > 1 / > 5 s | cpu > 0.1 / > 1 / > 5 s | +rápido / −lento (> 1 s) |
   |---|---|---|---|
   | dual | ×1.023 / ×1.018 / ×0.999 | ×1.062 / ×1.062 / ×1.031 | 6 / 28 |
   | **dual + paso largo** | ×1.003 / ×1.006 / **×0.945** | **×1.037 / ×1.042 / ×0.983** | 5 / 18 |

   Dual + paso largo contra dual: cpu ×0.976 / ×0.972 / ×0.953, 46 más rápidas
   y 25 más lentas en > 0.1 s. El paso largo recupera del orden del 40 % de la
   brecha del dual. Queda un 3–4 % en los cortes bajos, y en > 5 s el dual con
   paso largo ya es más rápido que el primal. Resueltas en 4/4: 131 en los tres
   brazos. 0 óptimos incompatibles. `ship-1` sigue decidiéndose por semilla, y
   `ex14_2_4` (×1.8 en celdas en los dos duales) es caos del `loup`: su
   óptimo es 0, y con el primal va de 48 a 1972 celdas entre semillas.

   **Pero ese «dual» no era dual.** Sonda `DFB_SX_MONO=1`: en la estrategia
   que se llamaba `dual`, entre el 38 y el 83 % de los LPs daban algún paso
   primal, y entre el 31 y el 88 % de los pasos eran primales. Causa: `k`
   estaba excluida de la reubicación, y en los LPs que empiezan con `k` no
   básica la base heredada es primal-factible, así que el bucle compuesto mete
   `k` con un paso primal y sigue en primal. Sólo el 0–24 % de los LPs era
   enteramente dual con cota monótona. El ×1.04 de arriba compara el primal con
   esa mezcla.

   **Estrategias rehechas (2026-09-30), `DFB_SX_PIVOTEO=`:**
   - `primal` (omisión): sin cambios. La cota 0 es dual (fase 1 gratis) y el
     resto primal; medido, el 78–92 % de los LPs dan pasos primales.
   - `dual`: el dual **puro**. Reubica todas las no básicas, `k` incluida, en
     frío y en tibio, así que la base es dual-factible en cada iteración. Si
     `k` termina no básica, el óptimo del LP es su cota y no había nada que
     contraer. Medido: 0–0.6 % de LPs con algún paso primal.
   - `dual` + `DFB_SX_ENTRAK=1`: si `k` no es básica, se la mete con un pivote
     en la fila de mayor `|a_rk|` y después se reubica todo. Dual puro con `k`
     básica desde el primer paso: el certificado es una cota útil desde el
     comienzo.
   - `mixta`: la que antes se llamaba `dual`, para reproducir mediciones.

   Instancias sueltas (7 pesadas, 2 semillas, las duales con paso largo;
   `results/estrategias_pivoteo_sueltas.csv`), tiempo contra primal:

   | | `ex6_1_3` | `ex8_4_4bis` | `ex2_1_9` | `ex7_2_3` | `ex5_3_2` | `bearing` | `himmel16` |
   |---|---|---|---|---|---|---|---|
   | piv/LP primal → dual → dual+entrada | 5.1 → 6.2 → 7.2 | 1.5 → 0.7 → 1.7 | 2.4 → 1.7 → 2.9 | 1.9 → 1.5 → 2.7 | 5.6 → 7.9 → 8.2 | 4.9 → 3.7 → 7.2 | 5.2 → 4.8 → 7.4 |
   | dual | ×1.06 | ×1.00 | ×1.00 | ≈ | ×1.23 | **×2.8, 18 240 celdas al tope** | **×1.7, celdas ×1.9** |
   | dual + entrada de `k` | ×1.10 | ×1.05 | ×1.07 | ×1.1–1.35 | ×1.27 | celdas 3690 (primal 8832/12 802) | celdas 4752 (primal 5300) |

   El dual puro es el más barato por LP y empata con el primal en cuatro de
   siete, pero **en `bearing` y `himmel16` el árbol crece ×2**. No es la
   semilla: `bearing` con el dual puro no termina en 60 s con ninguna de 8
   semillas ni con la LU diferida, y el primal la resuelve con las 8. La
   aparente mejora de la entrada de `k` en `bearing` sí era redondeo (con la
   LU diferida, una semilla llega al tope).

   **Mecanismo (sonda sobre el mismo nodo, `DFBH_ABDUAL=τ DFBH_ABPURO=1`,
   `DFBH_ABNODO=n` para trazar):** en 11 400 nodos de `bearing` el primal
   contrae más en 6363 y el dual puro en 500; reducción media de perímetro
   0.0616 contra 0.0522; 28 vacíos que sólo prueba el primal contra 1. Nodo
   665, cota `min x2`: el primal pivotea 10 veces, termina `INFEASIBLE` y su
   rayo prueba vacío en intervalos (`γ·z ∈ [−3.6e-5, −2.0e-8]`); el dual puro
   reubica `x2` en su cota, da 3 pasos y declara `OPTIMAL` con `x2` no
   básica, o sea `y = 0`, un certificado que no prueba nada. La caja mide
   ~5e-8 en `x0` y ~3e-7 en `x2`, y la tolerancia primal es
   `1e-9·(1+|cota|)` ≈ 7e-9: **del orden del 13 % del ancho**. En cajas
   angostas «factible dentro de la tolerancia» no significa nada, y cuando el
   LP termina con `k` no básica el certificado `y = 0` no puede probar ni
   contracción ni vacío. El primal lo sufre menos porque termina con `k`
   básica y `λ` no trivial. Explica también la regresión histórica que el §1
   atribuía a reubicar `k`.

   **Arreglo: tolerancia a escala de la caja (`DFB_SX_TOLANCHO=τ`).** La
   tolerancia primal de cada básica pasa a `min(1e-9·(1+|cota|), τ·ancho)`
   con un piso de redondeo, en el test de factibilidad del paso dual. Sobre el
   mismo nodo en `bearing`, con τ = 1e-3 el dual puro contrae como el primal
   (reducción media 0.0948 contra 0.0947, 0 vacíos exclusivos; sin τ, 0.0522
   contra 0.0616 y 1 contra 28). De punta a punta, `bearing` con el dual puro
   pasa de no encontrar nunca un `loup` a resolverse en 2934 / 15 668 celdas.

   Banco (138 instancias, 4 semillas, 8 cores, tres brazos intercalados;
   `results/tolancho_4semillas.csv`):

   | | celdas > 0.1 / > 1 / > 5 s | cpu > 0.1 / > 1 / > 5 s | +rápido / −lento (> 1 s) |
   |---|---|---|---|
   | primal con τ = 1e-3, contra primal sin τ | ×1.003 / ×1.008 / ×1.010 | ×1.004 / ×1.007 / ×1.004 | 5 / 7 |
   | **dual puro + paso largo + τ, contra primal + τ** | ×1.010 / ×1.027 / ×1.017 | **×1.037 / ×1.052 / ×1.038** | 6 / 25 |

   τ es neutro para el primal (487 de 552 pares con celdas idénticas) y es lo
   que hace funcionar al dual puro. El dual bien hecho resuelve las mismas 131
   instancias, 0 óptimos incompatibles, cpu por celda ×1.012 (×0.976 en las
   que llegan al tope), y queda a ×1.04–1.05 del primal. La pérdida se
   concentra en `ex14_2_7` (×1.31), `ex5_3_2` (×1.33, 7.9 contra 5.5 piv/LP),
   `himmel16` (×1.16) y `ex14_2_4` (caos del `loup`).

   **El ~1 % de pasos primales del dual puro.** Sonda (`DFB_SX_MONO=1`, que
   además clasifica cada paso primal): todos ocurren **después** de pasos
   duales, nunca al empezar el LP, nunca en la columna `k` y nunca con cota
   infinita; la violación de signo del costo reducido que los dispara va de
   1e-8 a más de 1e-2. O sea que algún paso dual rompe la factibilidad dual.
   Causa principal: la banda de Harris es **relativa al cociente**
   (`HARRIS_BAND = 1e-3`), así que el paso puede pasar hasta un 0.1 % del
   quiebre mínimo y las columnas que quedan atrás terminan con el signo violado
   en hasta 1e-3·|d|; el bucle compuesto lo repara con un paso primal. El
   Harris de libro, de dos pasadas con tolerancia absoluta sobre `d`
   (`DFB_SX_HARRISABS=1`), los elimina casi del todo —pero no en `bearing`
   (0.4–0.9 %)— y no gana con claridad: `ex5_3_2` 8.1 → 9.0 piv/LP (+15 %
   de tiempo), neutro o dentro del caos en el resto. Queda detrás del
   interruptor. Un paso primal esporádico no invalida nada: la cota se
   certifica igual; sólo corta la monotonía en ese LP.

   **La pérdida en `ex5_3_2` es estructural, no de orden.** Calibración de
   los predictores (`DFBH_PREDICE`, `DFBH_PREDICE2`): el conteo de columnas con
   signo violado, que es lo que usa el orden por omisión, predice bien el costo
   del dual (0.9 → 18.7 pivotes al crecer el conteo); el de filas no. Lo que
   cambia es el costo de cada caso: con una sola columna violada el primal
   gasta 5.45 pivotes y el dual 8.0, y las cotas tibias cuestan ~7.7 contra
   ~4.8. El primal mete esa columna a la base; el dual la reubica y deja
   infactibles muchas de las 49 filas, que repara una por una. El remedio
   estándar es el *pricing* de la fila que sale. Devex (`DFB_SX_DEVEX=1`), que
   se había medido «mixto» con la estrategia que en los hechos era primal, con
   el dual puro: `ex6_1_3` 6.17 → 4.16 piv/LP (primal 5.15), tiempo 24.5 →
   22.3 s contra 23.8 del primal; `ex5_3_2` 8.1 → 6.7 piv/LP, brecha contra el
   primal de +38 % a +14 %; `ex14_2_7` neutro; `himmel16` peor (caótica).

   Banco (138 instancias, 4 semillas, 8 cores, tres brazos intercalados;
   `results/devex_dual_4semillas.csv`), cpu > 0.1 / > 1 / > 5 s:

   | | con `bearing` | sin `bearing` |
   |---|---|---|
   | dual + devex contra dual | ×0.978 / ×0.950 / ×0.906 | **×0.991 / ×0.981 / ×0.987** (39/25 más rápidas/lentas en > 0.1 s) |
   | dual + devex contra primal | ×1.013 / ×0.992 / ×0.942 | **×1.027 / ×1.026 / ×1.025** |
   | dual contra primal | ×1.030 / ×1.045 / ×1.031 | |

   `bearing` sola hace ×0.37 en celdas y domina los cortes altos; sin ella,
   devex le quita al dual 1–2 % y lo deja a ~2.5 % del primal. Resueltas 131
   en los tres, 0 óptimos incompatibles. Donde más pierde el dual + devex:
   `himmel16` ×1.29, `launch` ×1.19, `ex5_3_2` ×1.19, `ex14_2_7` ×1.16.

   **Devex por omisión desde 2026-09-30** (`DFB_SX_DEVEX=0` lo apaga).

   **Con `eps_x` = 1e-7 la brecha se cerró, y DSE no suma** (banco, 138
   instancias, 4 semillas, 8 cores, tres brazos intercalados;
   `results/dse_epsx1e-7_4semillas.csv`). Cpu > 0.1 / > 1 / > 5 s:

   | | cpu | celdas | +rápido / −lento (> 0.1 s) |
   |---|---|---|---|
   | **dual + devex (omisión) contra primal** | **×1.014 / ×1.002 / ×0.959** | ×0.999 / ×0.997 / ×0.969 | 22 / 45 |
   | dual + DSE exacto (violación absoluta) contra dual + devex | ×0.988 / ×0.935 / ×1.043 | ×0.974 / ×0.915 / ×1.028 | 29 / 57 |

   Resueltas en 4/4: primal 131, dual 131, DSE 132. 0 óptimos incompatibles.
   La brecha del dual contra el primal, que con `eps_x` = 1e-6 era ~3 %, queda
   en ~1 % en los cortes bajos y el dual es más rápido en > 5 s; persiste en
   `launch` y `ex5_3_2` (×1.18–1.20) y `ex14_2_7` (×1.16). DSE exacto
   (`DFB_SX_DSE=1`; pesos `‖e_rᵀB⁻¹‖²` recalculados en `pivot` desde la parte
   λ del tableau, sin FTRAN; `DFB_SX_DSEABS=1` usa la violación absoluta)
   no gana: el ×0.935 en > 1 s lo hacen dos casos de caos del `loup`
   (`ex7_2_2`, que es la misma instancia que `ex8_1_8bis`, con una semilla en
   la que el dual con devex hace 2348 celdas en vez de ~70; y `bearing`); más
   instancias empeoran que mejoran y la cpu por celda sube ×1.03 por el cálculo
   de normas. Queda apagado. `bearing` con `eps_x` = 1e-7 también le cuesta al
   primal (una semilla al tope, otra 16 182 celdas): es caos del `loup` para
   todos, no sólo para DFB.

   **Por omisión desde 2026-09-30:** dual puro, paso largo y τ = 1e-3.
   `DFB_SX_PIVOTEO=primal`, `DFB_SX_BFRT=0` y `DFB_SX_TOLANCHO=0` restauran
   lo anterior. Las cifras de §2–§5 se midieron con el primal. Queda por
   atacar la pérdida en `ex5_3_2` y `ex14_2_7`.

   Omisión: sigue `primal`.

4. **El segundo hueco exacto: `ibexsolve`.** `DefaultSolver` compone
   `CtcPolytopeHull` + XTaylor en el mismo lugar. Sin cota superior, toda la poda
   del contractor lineal es prueba de vacío, donde DFB tiene su mayor ventaja.
   Conservar el linearizador de producción para aislar el cambio.

5. **La variante parcial de la propagación propia, dentro de ACID.** OJO: el
   test de `γ` sobre las rodajas **ya se construyó y se midió**
   (`CtcDFBGamma`, `DFBH_ENACID=1`, §28–§31 de `MEDICIONES_DUELO.md`): mata el
   11–17 % de las rodajas sin pivotear y es complementario con HC4, pero en el
   árbol pierde (celdas ×0.925 → ×0.976 contra xn). El costo del test es
   0.4–5 %; lo que se paga es que perturba la adaptación de ACID. Con el orden
   normal los `γ` venían del nodo anterior del recorrido y se rechazaban por
   contención en el 67–98 % de las llamadas; con el lineal primero se aplican
   siempre pero aciertan poco.

   Qué cambió desde entonces: el filtro de estabilidad (el 48 % de los
   `OPTIMAL` eran pivotes fallidos, así que muchos `γ` eran débiles), la
   relajación afín (23 filas, más apretada) y el costo del contractor (×0.57
   contra producción). Qué no cambió: el mecanismo del §31.

   Lo nuevo que no se probó: guardar los `γ` del **padre** en la celda
   (`Bxp`), de modo que todo hijo los tenga aplicables por contención. Es el
   régimen que en `alkylbis` acertaba el 40 %. Antes de construirlo, re-medir
   la tasa bajo la configuración de hoy; si no subió respecto del §28, no hay
   nada que ganar.

   **Re-medido (2026-09-30): no hay nada que ganar. Cerrada.** `DFBH_ENACID=1`
   con `DFBH_GDIST=1`, 10 instancias (las del §29–§30 más cinco pesadas),
   semilla 1:

   | | orden normal (`γ` de un nodo anterior que contiene al actual, casi siempre un ancestro) | lineal primero (`γ` del propio nodo) |
   |---|---|---|
   | aplicables | 2–35 % de las llamadas | 100 % |
   | aciertos sobre las aplicables | `alkylbis` 41 %, `ex14_2_1` 41 %, `ex8_4_4bis` 22 %, `bearing` 12 %, `ex6_1_3` 7 %, `himmel16` 4 %, `ex5_3_2` 4 %, `house` 2.5 %, `ex7_2_3` 2 %, `launch` 0 % | **≈ 0 % en las diez** (máx. 1.3 % en `bearing`); en el §29 era 13.8 % en `alkylbis` |

   Dos lecturas. (1) Con los `γ` del propio nodo el test ya no mata nada:
   con el filtro de estabilidad la contracción del lineal es completa, y el nodo
   queda consistente con sus certificados. (2) Con los `γ` de un ancestro sí
   mata, pero **lo que mata lo mataría igual el lineal del propio nodo**, que
   corre después de ACID: `alkylbis` mata 205 rodajas y termina en las mismas 64
   celdas; `ex8_4_4bis` mata 4331 y cambia las celdas −1 %. Donde las celdas no
   cambian, el tiempo sube 6–29 %. Guardar los `γ` del padre en la celda sólo
   multiplicaría las llamadas a un test cuyas muertes no mueven el árbol. Es el
   patrón del §6 en su forma más limpia. Las cubetas de `DFBH_GDIST` no
   distinguen padre de ancestro: miden la distancia de la **rodaja**, que HC4
   ya encogió en varias variables.

6. **El rediseño para relajaciones grandes** (BTRAN/FTRAN sin tableau). Único
   camino para que `compo` (celdas ×0.83) sea rentable. Va después del punto 1,
   porque el punto 1 cambia el reparto sobre el que se decide.

No reabrir salvo evidencia nueva: cortes anytime por presupuesto o brecha,
caché de `γ`, HC4 intercalado, congelar la relajación, el test de `γ` en las
rodajas de ACID (opción 5) y la cota lagrangiana con monotonía (arriba). Los cuatro tienen la
misma explicación (§6) y se midieron con el banco completo.

