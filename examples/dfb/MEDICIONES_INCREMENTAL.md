# ¿Es DFB realmente incremental?

Primera medición de la propiedad que sostiene toda la arquitectura: **cuando la
propagación aprieta las cotas de la caja, ¿reoptimizar desde la base anterior
cuesta menos pivotes que arrancar de cero?**

Herramienta: [incr_dfb.cpp](incr_dfb.cpp) (`make incr_dfb`). Datos en
[results/incremental.csv](results/incremental.csv), lote en
[results/run_incr.sh](results/run_incr.sh).

---

## 1. Por qué hacía falta un arnés propio

Con el banco normal la pregunta es **vacía**: `CtcDFBPropag` construye los
contractores con `max_iters = 1`, así que los pivotes por llamada son ≤1 por
construcción y no hay nada que comparar. El arnés llama a `CtcDFB` directamente
con `max_iters` grande.

Por contractor se corren tres cosas:

1. **Frío** en la caja original, hasta el punto fijo.
2. **Tibio** en la caja apretada: se conserva la base de (1) y se sigue.
3. **Frío** en la caja apretada: contractor nuevo desde `refA`.

La caja apretada encoge cada variable un 25 % por un lado, como haría la
propagación al contraer; `refA` sigue siendo válida porque es una subcaja. La
pregunta es si (2) cuesta mucho menos que (3) **y llega a la misma cota**.

## 2. Resultado

148 instancias con datos, **7065 contractores medidos**.

| | pivotes |
|---|---|
| frío, caja original | 75 150 |
| **tibio, caja apretada** | **97 139** |
| frío, caja apretada | 97 168 |
| **razón agregada tibio/frío** | **1.000** |

| razón por instancia | |
|---|---|
| **mediana** | **0.032** |
| p25 / p75 | 0.000 / 0.222 |
| máximo | 1.982 |
| instancias con razón < 0.5 | **137 de 144** |

| calidad de la cota alcanzada | contractores |
|---|---|
| igual que en frío | 6692 (**95 %**) |
| tibio mejor | 11 |
| tibio **peor** | 362 (5 %) |

## 3. Lectura: sí es incremental, pero el agregado engaña

**En la instancia típica la incrementalidad funciona, y con holgura**:
reoptimizar cuesta el **3 %** de arrancar de cero (mediana), y en 137 de 144
instancias menos de la mitad. En un grupo entero de instancias cuesta
**literalmente cero pivotes** — la base anterior ya es óptima para la caja
apretada:

| instancia | tibio | frío |
|---|---|---|
| Brown-05 | **0** | 36 |
| Brown-06 | **0** | 59 |
| Brown-07 | **0** | 83 |
| Brown-10 / Brown-10-5 | **0** | 170 |
| Brown-10sp | **0** | 100 |
| Brent-8 / Brent-10 | **0** | 14 / 5 |

Eso es exactamente la propiedad dual que se venía invocando en teoría: apretar
cotas preserva la factibilidad dual, así que la base anterior sigue sirviendo.

**Pero la razón agregada es 1.000**, y no por casualidad: tres instancias
concentran la mayoría de los pivotes del banco y en ellas el warm start cuesta
**más** que arrancar de cero.

| instancia | tibio | frío | razón |
|---|---|---|---|
| Geneigbis | 16 124 | 8 136 | **1.982** |
| Fourbar-icse | 48 031 | 32 783 | **1.465** |
| Ex14-2-3 | 16 003 | 16 044 | 0.997 |

No están topeadas por el límite de seguridad (686 pivotes por contractor en
Fourbar-icse), así que es un fenómeno real: **la base heredada puede quedar
lejos en la geometría nueva** y llevar a un camino más largo que empezar de
nuevo. Citar solo el agregado (1.000) sería tan engañoso como citar solo la
mediana (0.032).

Y el 5 % de contractores que llegan a una **cota peor** tibios que fríos dice
que la incrementalidad no es gratis en calidad: las secuencias de pivotes
difieren y se detienen en óptimos duales locales distintos.

## 4. Consecuencia de diseño

El warm start entre nodos vale la pena, pero **necesita una salvaguarda**, que
además es práctica estándar en los solvers de LP: si la reoptimización tibia
supera un presupuesto —por ejemplo el número de pivotes que costó la corrida
anterior, o un múltiplo— **abandonar y reiniciar en frío**. Con eso se queda con
la mediana de 0.032 y se acota el daño de los casos como Geneigbis.

Sin esa salvaguarda, el warm start es neutro en el banco completo: ganaría en
137 instancias y perdería todo lo ganado en tres.

## 5. Lo que esta medición no cubre

- Solo mide **dentro de una linealización** (`refA` fija). El warm start *entre
  nodos* implica además que `refA` cambia al relinealizar, y ahí la base
  anterior no significa nada. Habrá que decidir si se relineariza en cada nodo
  (como PolytopeHull) o se reusa la linealización del padre mientras sirva.
- El aprieto del 25 % por un lado es arbitrario. Convendría repetirlo con
  aprietos más chicos (1 %, 5 %), que es lo que realmente hace la propagación
  entre llamadas consecutivas, y ver si la razón mejora todavía más.
- No mide tiempo, solo pivotes. Como el costo por pivote es ahora uniforme
  (flotantes, sin regeneraciones), los pivotes son un buen proxy, pero
  conviene confirmarlo.

---

## 6. Con aprietos realistas es todavía mejor

El 25 % del §2 es un aprieto agresivo. Repitiendo con lo que realmente hace la
propagación entre llamadas consecutivas
([results/incremental_sq01.csv](results/incremental_sq01.csv),
[results/incremental_sq05.csv](results/incremental_sq05.csv)):

| aprieto | razón agregada | mediana por instancia | razón < 0.5 | cota peor |
|---|---|---|---|---|
| **1 %** | 0.927 | **0.000** | 139 de 145 | **0 %** |
| **5 %** | 0.815 | **0.000** | 138 de 144 | **0 %** |
| 25 % | 1.000 | 0.032 | 137 de 144 | 5 % |

Con aprietos realistas **reoptimizar cuesta cero pivotes en la instancia
mediana**: la base anterior ya es óptima para la caja apretada. Y el 5 % de
cotas peores del §2 **era un artefacto del aprieto agresivo**: con 1 % y 5 %
desaparece por completo.

## 7. Warm start entre nodos: **medido y descartado**

La conclusión del §4 —que el warm start vale la pena con salvaguarda— resultó
equivocada al implementarlo, por una razón estructural que la métrica sola no
mostraba.

**Las bases tibias y la relinealización están acopladas.** Una base solo sigue
tibia si `refA` no cambia, o sea solo si se reusa la linealización. Reusar la
del padre **es sólido** —una relajación calculada sobre una caja mayor sigue
valiendo para cualquier subcaja— y se implementó verificando la contención
(`lin_box.is_superset(box)`), que es indispensable: en el optimizador las
llamadas consecutivas vienen de nodos distintos y dos hermanos no se contienen.

Pero una relajación calculada sobre una caja mayor es **más floja**, y eso se
paga en árbol ([results/reuse_ibexopt.csv](results/reuse_ibexopt.csv)):

| | celdas | tiempo total | peor en celdas |
|---|---|---|---|
| reuso mientras el perímetro no caiga a la mitad | **×1.167** | 89.2 s → 110.5 s | 45 de 98 |
| reuso mientras no caiga un 10 % | ×1.015 | 89.5 s → 92.3 s | 27 de 99 |

Y lo que se ahorraría ya no vale casi nada: **los pivotes de reoptimizar son
cero en la mediana** (§6), así que el warm start solo ahorra las `2n`
inicializaciones, y eso no compensa un 17 % más de celdas.

Queda implementado y **apagado** (`CtcDFBPropag::reuse_linearization`,
`DFB_REUSE_LIN=1`, con `DFB_RELIN_SHRINK` para el umbral). Serviría si alguna
vez la relinealización se volviera mucho más caro que hoy, o si se encontrara
cómo reusar la base **sin** reusar la relajación: conservar la estructura de
pivotes y refactorizarla contra la `refA` nueva es válido (cualquier `λ` da una
cota válida) y es lo que hace el simplex dual cuando cambia la matriz. Con los
pivotes ya en cero, no esperaría mucho.
