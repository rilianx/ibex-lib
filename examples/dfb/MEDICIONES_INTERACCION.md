# La interacción DFB ↔ propagación: cómo funciona, y cuánto rinde

El diseño original quiere un **DFB incremental que interactúe con la propagación
de restricciones**. Este documento verifica si eso está ocurriendo, arregla lo
que estaba roto, y mide cuánto rinde.

Datos: [results/revive_ibexopt191.csv](results/revive_ibexopt191.csv).

---

## 1. Las dos vías de reactivación

En el bucle de `CtcDFBPropag::contract` un contractor DFB vuelve a la cola por
dos caminos distintos:

1. **Se reencola a sí mismo** mientras su `state` siga en `CONTRACTING`. Con
   `max_iters = 1` cada aplicación avanza **un pivote** y vuelve a la cola, así
   que dentro de una llamada cada contractor progresa incrementalmente sobre una
   caja que se va apretando. **Funciona**, en los dos modos: son las ~19
   aplicaciones por contractor que se miden en `Katsura-50`.
2. **La propagación HC4 lo reactiva** cuando aprieta una variable que influye en
   su cota, vía `CtcDFB::real_impact`. **Esta es la interacción DFB↔propagación**,
   el núcleo del diseño incremental.

Fuera de una llamada no hay incrementalidad: cada llamada relineariza y
reinicializa los `2n` contractores. Reusar la linealización se midió y no paga
(§7 de [MEDICIONES_INCREMENTAL.md](MEDICIONES_INCREMENTAL.md)).

## 2. La vía 2 estaba muerta en el camino flotante

`real_impact` leía `A[0]` —la matriz de **intervalos**— mientras la fila gamma
flotante vive en `Af[0]`, y `init_float` nunca llena `A`. Peor: una guarda
defensiva (`var >= A[0].size()`) lo volvía **silencioso** en vez de un fallo.

| instancia | reactivaciones en intervalos | en flotantes | llamadas a `real_impact` |
|---|---|---|---|
| ex14_2_4 | 875 | **0** | 88 860 |
| ex6_2_12 | 3380 | **0** | 389 120 |
| alkyl | 1354 | **0** | 30 630 |

Corregido: en modo flotante lee `Af[0][var]`. `DFB_NO_REVIVE=1` reproduce el
comportamiento roto, para poder aislar su efecto.

**Nota importante sobre el banco de una caja**: `bench_dfb` construye el
contractor con `only_dfb = true`, o sea **sin contractores HC4**. Todo ese banco
—los grupos A/B/C, las 153 instancias, las comparaciones contra PolytopeHull—
mide **DFB aislado** y **no ejercita esta interacción en absoluto** (0 llamadas a
`real_impact`, verificado con contadores). Es una comparación legítima contra
PolytopeHull, que también corre solo, pero no dice nada sobre el diseño
incremental.

## 3. Cuánto rinde la interacción: casi nada, y con mucha varianza

Sobre las 191 instancias (`easy`, `medium`, `hard`, `blowup`, `benchs-minlp`)
con límite de 30 s:

| | celdas | tiempo | resueltas |
|---|---|---|---|
| con reactivación vs sin | **×1.03** | **×1.03** | 116 vs 113 |

Mejora en 16 instancias, empeora en 21, **igual en 75**. Es decir: el mecanismo
que el diseño considera central es **neutro en promedio**, con dispersión alta.

| mejores | sin → con | | peores | sin → con |
|---|---|---|---|---|
| ex14_2_4 | 7012 → **5098** celdas | | house | 940 → **2410** |
| ex8_5_5-1 | 1688 → 1364 | | alkyl | 264 → **588** |
| ex7_3_5bis | 5776 → 4924 | | ex8_5_5 | 9584 → 13 140 |
| ex5_2_2_case1 | 388 → 298 | | launch | 218 → 344 |

**Se deja activado** porque es el comportamiento intencionado, corrige un error
real, y resuelve 3 instancias más. Pero **no es la palanca que cierra la
brecha**: es la cuarta hipótesis de la línea «costo por unidad de poda» que
resulta neutra o negativa, después del warm start entre nodos, la prioridad a la
variable objetivo y la aplicación periódica.

## 4. El banco pesado cambia los números

Comparando el mismo par de configuraciones en los dos bancos:

| | 132 instancias (easy+medium, 20 s) | **191 instancias (con hard, blowup, minlp, 30 s)** |
|---|---|---|
| producción, resueltas | 116 de 132 | 137 de 191 |
| `dfb`, resueltas | 102 | 116 |
| `dfb` vs producción | ×2.49 celdas, ×1.23 tiempo | **×2.83 celdas, ×1.43 tiempo** |
| `acid_dfb` vs producción | ×1.53 celdas, ×3.53 tiempo | ×1.47 celdas, **×3.75 tiempo** |

La brecha de `dfb` **se ensancha** con las instancias pesadas (×1.23 → ×1.43 en
tiempo, ×2.49 → ×2.83 en celdas). Confirma que el banco liviano la subestimaba,
y que las conclusiones de costo hay que sacarlas del banco ampliado.
