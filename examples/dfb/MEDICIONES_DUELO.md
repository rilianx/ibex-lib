# Duelo directo DFB vs SoPlex, sin propagación

Arnés: [duelo_dfb.cpp](duelo_dfb.cpp) (`make duelo_dfb`).
Banco: `benchs/optim/{easy,medium,hard}`, 139 instancias, **114** con
linealización utilizable (`m>1` y caja acotada). 2 046 cotas.

## 1. Por qué hacía falta

Todas las mediciones anteriores comparaban **DFB con su propagación** contra
**CtcPolytopeHull con la suya**. Eso mezcla tres efectos independientes:

1. el **precio por LP** de cada método;
2. el **volumen de LPs**: PolyHull resuelve el 47 % de las `2n` cotas porque
   `choose_next_variable` saltea las que no le sirven, y DFB visita el 145 %;
3. el **intercalado con HC4** dentro del punto fijo.

La conclusión que salía de ahí —«DFB cuesta ×64»— era la suma de los tres, y
resultó estar dominada por el segundo, que no dice nada sobre los métodos.

## 2. Qué se iguala

- **El mismo poliedro**: una sola linealización (`LinearizerXTaylor`, esquina
  `INF` fija, `HANSEN`), compartida por los dos. A SoPlex se le da en su forma
  natural —`nx` variables, `m` filas `lhs ≤ A_j·x ≤ rhs`— y a DFB la suya
  —`Ā z = 0` con `z=(x,b)`—, que son el mismo conjunto: `Ā[j][nx+j] = −1`
  significa `A_j·x = b_j`, y las cotas de `b_j` son `lhs/rhs`.
- **La misma caja para las `2n` cotas**: no se contrae sobre la marcha, así que
  los `2n` LPs son literalmente idénticos de los dos lados y las cotas se
  comparan **una a una**. (PolyHull en producción sí contrae progresivamente;
  eso es parte de su orquestación, no de su método.)
- **Las mismas `2n` cotas, en el mismo orden**, sin saltear ni reencolar.
- **El mismo criterio de parada**: los dos resuelven hasta el óptimo.

DFB es *anytime*: `DFBSimplex::solve` calcula `lambda` desde la factorización
pase lo que pase —también al agotarse el presupuesto— y `γ = λᵀĀ` recertifica
en intervalos, así que detenerlo en cualquier pivote da una cota **más débil
pero nunca inválida**. Un simplex no tiene eso: su cota aparece recién en el
óptimo. **Esa ventaja acá no se usa a propósito**, para que la comparación sea
pareja; se mide aparte en §5.

## 3. Solidez: verificación por testigos

Cada minimizador primal que devuelve SoPlex es un punto (casi) factible del
poliedro. Se guardan los 2 046 y se verifica que **ninguna** de las dos cajas
contraídas los deje afuera.

| | peor exceso relativo |
|---|---|
| residuo de factibilidad de los propios testigos | 2.24e-14 |
| exceso de la caja de PolytopeHull | 5.81e-14 |
| exceso de la caja de DFB | **2.24e-14** |

**0 testigos expulsados** por ninguno de los dos, en 0 instancias. Los excesos
están al nivel del residuo de los testigos, o sea que son el propio error del
testigo y no una poda indebida. Esto es lo que valida el §4: donde DFB poda más
que PolytopeHull, poda **más**, no **de más**.

## 4. Resultado

### Calidad de la cota, cota a cota

| | cotas | % |
|---|---|---|
| idénticas | 1 964 | **95.99** |
| DFB **mejor** | 66 | 3.23 |
| DFB **peor** | 16 | **0.78** |

**99 de 114 instancias dan las `2n` cotas exactamente iguales.** DFB pierde
alguna cota en 7 instancias y gana alguna en 11.

Que DFB gane el 3.2 % no es un error de medición: son cotas válidas (§3). Dos
causas, las dos del lado de PolytopeHull: SoPlex no logra `OptimalProved` en el
**8.2 %** de sus LPs (1 878/2 046, contra 2 032/2 046 de DFB) y esas cotas se
descartan enteras; y donde sí lo logra, la certificación Neumaier–Shcherbina
puede ser más floja que `γ = λᵀĀ` evaluado en intervalos.

(La tabla usa el `lambda` del óptimo. Guardar el mejor `lambda` del recorrido
—lo que hace producción— mejora la cota final en 2 instancias más: ver §5.)

### Costo

| | pivotes | por cota | tiempo |
|---|---|---|---|
| SoPlex | 6 485 | 3.17 | 0.297 s |
| DFB | **5 851** | **2.86** | **0.066 s** |
| razón | ×0.90 | | **×0.22** |

Por media geométrica sobre instancias (que pesa parejo a las chicas): pivotes
**×1.27**, tiempo **×0.32**.

La diferencia entre el total (×0.90) y la geométrica (×1.27) es que DFB usa algo
más de pivotes en las instancias chicas y bastantes menos en las grandes. En las
5 instancias con `n ≥ 20`:

| | por cota | tiempo | iguales | peor | mejor |
|---|---|---|---|---|---|
| SoPlex | 13.26 | 0.165 s | | | |
| DFB | **6.12** | **0.023 s** | 84.9 % | 4 | **35** |
| razón | ×0.46 | ×0.14 | | | |

`launch` (n=39, m=52): SoPlex 17.89 pivotes/cota, DFB 6.30. `ex2_1_10` (n=21):
SoPlex 36.74, DFB 8.36, tiempo ×0.04.

## 5. La propiedad anytime

`./duelo_dfb --anytime`. La curva sale de **una sola corrida tibia, la misma que
llega al óptimo**, leída en puntos anteriores: `DFBSimplex::traza_y` registra
`lambda` después de cada pivote, así que `traza_y[j]` es exactamente el
certificado que se obtendría deteniendo el bucle en el pivote `j`. Por
construcción

    pivotes(b) = Σ_i min(b, pivotes_óptimo_i)  ≤  pivotes(óptimo)

y la curva satura exactamente en el óptimo. **Medido: 0 violaciones de
`pivotes(b) ≤ pivotes(óptimo)` y 0 curvas no monótonas, en las 114 instancias.**

| presupuesto | pivotes | % del costo al óptimo | % de la contracción | instancias al 100 % |
|---|---|---|---|---|
| 0 | 0 | 0.0 % | 30.1 % | 23/114 |
| 1 | 1 282 | 21.9 % | 59.5 % | 43/114 |
| 2 | 2 215 | 37.8 % | 70.9 % | 55/114 |
| 4 | 3 352 | 57.2 % | 84.9 % | 71/114 |
| 8 | 4 535 | 77.4 % | **94.5 %** | 89/114 |
| 16 | 5 404 | 92.3 % | 98.8 % | 110/114 |
| 32 | 5 749 | 98.2 % | 99.7 % | 113/114 |
| 64 | 5 846 | 99.8 % | 100.0 % | 114/114 |
| **óptimo** | **5 857** | 100 % | 100 % | 114/114 |

La cota aparece temprano y el rendimiento decrece: con **0 pivotes** —la base
inicial de las `b`, que ya es dual-factible— se obtiene el 30 % de la
contracción; con 1 pivote por cota, el 22 % del costo da el 59.5 %; con 8, el
77 % del costo da el 94.5 %. El último cuarto del presupuesto compra el 5 %
final.

### Dos defectos que esta medición destapó, los dos arreglados

**1. `ITER_LIMIT` invalidaba la base.** `solve` hacía `base_lista = false` en
toda salida distinta de `OPTIMAL`, así que cortar por presupuesto obligaba a la
cota siguiente a arrancar en frío. Eso rompe justamente el caso de uso anytime:
detenerse antes salía **más caro** que llegar al óptimo (15 042 pivotes contra
5 851). Pero la base tras `ITER_LIMIT` es perfectamente válida —el tableau y la
factorización se mantuvieron consistentes en cada pivote—, solo que no es la
óptima. Ahora solo invalidan `INFEASIBLE`, `UNBOUNDED` y `SINGULAR`.

El cambio es inerte fuera del uso anytime: 0 instancias del §4 dan un resultado
distinto, y en el banco entero con la regla de simplex (`DFB_SIMPLEX_RULE=1`,
139 instancias) hay **0 `ITER_LIMIT` en 2 587 resoluciones**, porque el tope de
producción es `4m+100` y el óptimo está a unos pocos pivotes. Solo se activa
cuando uno pone un presupuesto chico a propósito, que es exactamente el caso que
la propiedad tiene que soportar.

**2. El `lambda` del óptimo no siempre es el mejor que se vio.** El bucle
compuesto intercala pasos **primales**, que restauran factibilidad primal a
costa de bajar el objetivo dual, así que la calidad puede empeorar al avanzar:
en `bearing` la contracción caía 10 puntos entre los presupuestos 2 y 4. Un
algoritmo anytime no puede empeorar al correr más, de modo que la curva se mide
guardando el **mejor certificado del prefijo**, que es lo que el contractor de
producción ya hace con `best_lam`. Con eso quedan 0 curvas no monótonas.

Ese segundo defecto no es solo de la curva: guardar el mejor `lambda` mejora
también la cota **final**, en 2 de 114 instancias (`bearing` 69.8 % → 89.8 % de
contracción, `ex7_3_6` 38.3 % → 39.1 %). O sea que la tabla del §4, que usa solo
el `lambda` del óptimo, **subestima levemente a DFB**.

### Qué vale la propiedad acá

Truncar no acelera en este régimen: el óptimo está a 2.86 pivotes de media, así
que un tope bajo corta poca cosa y hay que llegar al 77 % del costo para
conservar el 95 % de la poda. Lo que la propiedad da es otra cosa: **un tope de
presupuesto no puede producir una cota inválida**, así que se puede poner uno
sin analizar el caso, y en cualquier punto del bucle hay una cota utilizable.
Sería además una palanca de velocidad en un régimen donde el óptimo estuviera a
decenas de pivotes, que no es éste.

## 6. Qué cambia esto

El déficit de DFB en el árbol **no está en el método de LP**. Por LP, DFB iguala
o supera a SoPlex en las dos dimensiones: calidad (96 % idénticas, 3.2 % mejor,
0.8 % peor) y costo (×0.90 en pivotes, ×0.22 en tiempo; ×0.46 y ×0.14 en las
grandes).

Lo que queda por cerrar es **la orquestación**, y son las dos cosas que este
arnés apagó a propósito:

1. **El volumen de LPs.** PolyHull resuelve el 47 % de las `2n` cotas y DFB
   visita el 145 %: un factor ×3 que no tiene nada que ver con el simplex. El
   mecanismo de PolyHull es `choose_next_variable`, que elige la próxima cota
   mirando la solución primal y saltea las que ya están en su cota.
2. **La contracción progresiva.** PolyHull actualiza las cotas del LP después de
   cada resolución (`set_bounds(i, box[i])`), así que las cotas posteriores se
   resuelven sobre un poliedro ya más chico.

Las dos son transferibles a DFB tal cual: no tocan el simplex, solo el bucle que
lo llama.

## 7. Pseudocódigo del DFB anytime, y dónde se intercalan las cotas

La respuesta corta a «¿se aplican cotas intercaladas?» es **sí, entre las `2n`
cotas, pero no dentro de una**, y en el duelo del §4 están apagadas a propósito.
Hay tres niveles y conviene no confundirlos.

### Nivel 1 — dentro de una resolución: NO se intercala

```
resolver_cota(k, maximizar, z, presupuesto):
    # base inicial gratis: las m columnas de b, con B = -I.
    # c = ±e_k => c_B = 0 => y = 0 => d = ±e_k >= 0: dual-factible sin fase 1.
    si no hay base reusable: base <- columnas de b ;  no_basicas <- cota inferior

    refrescar d, y, x_B desde la factorizacion          # y = c_k * fila_r de B^-1
    mejor_lambda <- y                                    # el certificado de 0 pivotes

    repetir hasta presupuesto:
        si x_B es factible: OPTIMO, cortar              # factible primal = optimo
        paso dual   (sale la basica que mas viola, entra la del menor |d_j/a_rj|)
        o paso primal (si algun d_j viola su signo)
        y <- c_k * (fila r de B^-1)                     # un BTRAN, exacto
        si f(y) > f(mejor_lambda): mejor_lambda <- y    # <-- ver abajo

    devolver mejor_lambda                # SIEMPRE, llegue o no al optimo
```

Lo que hace a esto *anytime* es la última línea: `lambda` sale de la
factorización y existe en todo momento, así que **cortar en cualquier iteración
devuelve un certificado**. La cota se obtiene después, en intervalos, y por eso
un `lambda` malo da una cota débil y nunca una falsa.

`mejor_lambda` en vez del último hace falta porque los pasos **primales** bajan
el objetivo dual: sin eso, correr más puede empeorar (§5).

No se aplica la cota a mitad de la resolución, y no serviría: la cota que sale
es sobre `z_k`, que es justo lo que se está minimizando, así que apretarla no
cambia el mínimo del LP.

### Nivel 2 — entre las `2n` cotas: SÍ se intercala (en producción)

```
contraer_caja(x):                       # CtcDFBPropag, un nodo
    A, z_ref <- linearizar(x)           # una vez por nodo
    simplex.cargar(A)
    cola <- los 2n contractores (k, lado)

    mientras la cola no este vacia:
        (k, lado) <- sacar de la cola

        work[0..n-1]  <- x               # <-- INTERCALADO: ve lo ya podado
        work[n..n+m-1] <- z_ref          # las cotas de b NO se refrescan

        lambda <- resolver_cota(k, lado, work, 4m+100)
        gamma  <- lambda^T A             # en INTERVALOS: aca se certifica
        gaussSeidel(work, k, gamma)      # aprieta z_k y de paso otras variables
        x <- work                        # <-- publica de inmediato

        si mejoro >= 1%: reencolar (k, lado)
```

Cada resolución arranca de la caja que dejaron las anteriores, y `gaussSeidel`
sobre `gamma·z = 0` puede apretar **varias** variables de una vez, no solo `z_k`.
Es el mismo mecanismo que `CtcPolytopeHull::contract` con su
`set_bounds(i, box[i])`, con dos diferencias: DFB además reencola mientras
rinda, y las cotas de la parte `b` se quedan en las de la linealización.

### Nivel 3 — el duelo del §4: NO se intercala, a propósito

El arnés fija la caja para las `2n` cotas. No es un descuido: es lo que hace que
los `2n` LPs sean **los mismos** de los dos lados y que las cotas se puedan
comparar una a una. Intercalar haría que la primera diferencia cambiara todos
los LPs siguientes y ya no se sabría si una cota distinta viene del método o de
haber partido de otra caja.

### Lo que el anytime habilita y PolytopeHull no puede hacer

Como toda parada es válida, se puede recorrer las `2n` cotas **por rondas de
presupuesto creciente**:

```
para presupuesto b en 1, 2, 4, 8, ...:
    para cada (k, lado):
        lambda <- resolver_cota(k, lado, x, b)    # tibio: conserva la base
        x <- aplicar(x, lambda)                    # cota debil pero valida
    si x dejo de moverse: cortar
```

Cada ronda barata aprieta la caja para la siguiente, así que las rondas caras se
corren sobre un poliedro ya más chico. Un simplex no admite este esquema: su
cota aparece recién en el óptimo, de modo que una ronda truncada no devuelve
nada y el presupuesto gastado se pierde entero.

Con los números del §5 la ronda de 1 pivote cuesta el 22 % del presupuesto total
y ya entrega el 59.5 % de la poda, así que hay margen para que el esquema pague.
**Sin medir todavía**: es la continuación natural del §6, junto con el volumen
de LPs.

## 8. El round-robin de presupuesto creciente: medido, y no paga

`./duelo_dfb --rondas` corre tres estrategias sobre la misma linealización, las
tres **con intercalado** (cada cota se aplica apenas sale) y las tres iterando
hasta punto fijo con el mismo criterio (una pasada que mejora el perímetro menos
del 0.01 % corta). Difieren solo en quién resuelve el LP y con cuánto
presupuesto.

| estrategia | pasadas | pivotes | tiempo | contracción media |
|---|---|---|---|---|
| SoPlex al óptimo | 189 | 9 666 | 0.437 s | 31.94 % |
| **DFB al óptimo** | 192 | **9 396** | **0.095 s** | **34.01 %** |
| DFB por rondas 1,2,4,8… | 792 | 32 468 | 0.320 s | 31.37 % |

Media geométrica por instancia contra SoPlex: DFB al óptimo ×1.20 en pivotes y
**×0.29 en tiempo**; por rondas ×2.85 y ×0.71.

Comparación de contracción, instancia por instancia:

| | gana el primero | gana el segundo | empatan |
|---|---|---|---|
| DFB al óptimo vs SoPlex | **7** | **0** | 107 |
| rondas vs SoPlex | 6 | 9 | 99 |
| rondas vs DFB al óptimo | 1 | 9 | 104 |

**El round-robin cuesta ×3.46 en pivotes y ×3.37 en tiempo sobre el DFB al
óptimo, y encima contrae 2.63 puntos menos de media.** Gana en 1 instancia de
114 (`ex7_3_6`: 43.5 % contra 34.8 %, a cambio de ×3.4 en pivotes) y pierde en 9.

### Por qué no paga

El esquema no *reemplaza* la pasada exacta: la *precede*. Las rondas baratas
suman su costo y después hay que hacer las exactas igual, porque el criterio de
corte exige que la última ronda haya llegado al óptimo en todas las cotas. Para
que pagara, las rondas baratas tendrían que dejar la caja tan apretada que las
exactas salieran mucho más baratas — y no puede pasar, porque la pasada exacta
ya costaba **2.86 pivotes por cota** (§4). No hay de dónde ahorrar.

Se ve en la proporcionalidad: 792 pasadas contra 192 es ×4.1, y el costo sube
×3.46 en pivotes y ×3.37 en tiempo. O sea que **el costo va con el número de
resoluciones, no con los pivotes de cada una**: lo que domina es la parte fija
por resolución —refrescar el objetivo, reubicar las no básicas, y sobre todo
certificar `γ = λᵀĀ`, que son `m × na` operaciones de intervalo pase lo que
pase—. Multiplicar pasadas multiplica esa parte fija.

Es la misma conclusión del §6 desde el otro lado: el problema de DFB es el
**volumen** de LPs, no su precio. Cualquier esquema que agregue resoluciones
pierde, incluso uno que las haga más baratas cada una.

La propiedad anytime sigue valiendo por lo que decía el §5 —ningún tope de
presupuesto puede dar una cota inválida— y **como palanca de velocidad queda
descartada *para esta pregunta***: cuando el presupuesto alcanza para llegar al
punto fijo, truncar no conviene. La pregunta que sí la favorece es la del §12,
a presupuesto total fijo. La dirección que queda es la contraria:
resolver *menos* LPs, como hace `choose_next_variable`.

### Control de solidez, y por qué el de testigos no sirve con intercalado

El test de testigos del §3 **pierde poder en cuanto se intercala**. Una
contracción legítima puede dejar afuera un testigo *flotante* apenas la caja se
aprieta por debajo de su residuo, y de ahí en adelante la cadena amplifica: en
`ex7_3_6` la variable 8 colapsa a `[0,0]` contra un testigo que vale −2.6e-10
—dentro de su propio residuo de 1.8e-10, o sea una contracción válida— y con
`x₈ = 0` fijo los LPs siguientes derivan `x₁₆ ∈ [1,1]`, que deja al testigo
afuera por 0.189. No es un fallo: es un residuo de 1e-10 amplificado por la
cadena.

El control que sí sirve es per-paso y riguroso. Neumaier–Shcherbina garantiza
que el óptimo **verdadero** del LP cae dentro del encierro que devuelve
`LPSolver::minimum()`, así que toda cota inferior válida tiene que quedar por
debajo de su extremo superior. La estrategia *sombra* hace que las dos resuelvan
sobre la **misma** caja —que avanza solo con la cota de SoPlex, para que no haya
amplificación— y compara cota a cota:

| | |
|---|---|
| cotas vinculantes comparadas | 594 |
| DFB más apretada | 7 (1.2 %) |
| DFB más floja | 11 (1.9 %) |
| **peor exceso sobre el encierro certificado** | **0** |
| instancias con alguna violación | **0** |

Y el mismo control aplicado al duelo del §4, sobre sus 2 046 cotas: **peor
exceso 0, 0 violaciones**. La solidez queda establecida sin depender de testigos
flotantes.

## 9. ¿Las bases son independientes? No, y hacerlas independientes es peor

Hasta acá el arnés usaba **un solo** `DFBSimplex` para las `2n` cotas, con
`warm=true`: cada resolución hereda la base de la cota **anterior**, que tenía
otro objetivo. No es estado por cota. `duelo_dfb` ahora compara las tres
políticas posibles:

- **compartida** (omisión): un simplex para las `2n` cotas. Hereda entre
  objetivos distintos. En un LP de variables acotadas eso restaura la
  factibilidad **dual** gratis —basta reubicar las no básicas por el signo del
  costo reducido— pero destruye la **primal**.
- **propia** (`DUELO_PROPIO=1`): `2n` simplex, uno por cota. Cada uno hereda de
  **su propia** resolución anterior, o sea el mismo objetivo con la caja apenas
  cambiada, que es la situación en que un warm start sirve de verdad. Pero solo
  tiene algo que heredar si la misma cota se resuelve más de una vez.
- **fría** (`DUELO_FRIO=1`): sin herencia, siempre desde la base de las `b`.

### Una sola pasada sobre las `2n` cotas, caja fija (§4)

| política | pivotes | por cota | tiempo | óptimos | cotas iguales |
|---|---|---|---|---|---|
| **compartida** | **5 851** | **2.86** | **0.066 s** | **2 032** | 1 964 |
| propia | 15 042 | 7.35 | 0.135 s | 2 026 | 1 953 |
| fría | 15 042 | 7.35 | 0.097 s | 2 026 | 1 953 |
| *SoPlex* | *6 485* | *3.17* | *0.297 s* | *1 878* | — |

**`propia` y `fría` dan exactamente los mismos 15 042 pivotes**, y eso explica
todo: con una base por cota, cada simplex se usa **una sola vez** en la pasada,
así que no hay nada que heredar y `propia` degenera en `fría`. La compartida, en
cambio, se reusa `2n` veces y ahorra **×2.57**.

O sea que la herencia entre objetivos distintos —la que uno esperaría inútil,
porque el objetivo cambia— es justamente la que paga. Es la misma propiedad que
le da a `CtcPolytopeHull` sus pocas iteraciones por cota.

### Punto fijo con intercalado, varias pasadas (§8)

| política | `op` pasadas | `op` pivotes | `op` tiempo | `rd` pivotes | `rd` tiempo |
|---|---|---|---|---|---|
| **compartida** | 192 | **9 396** | **0.095 s** | 33 764 | 0.324 s |
| propia | 189 | 14 390 | 0.138 s | **15 943** | **0.180 s** |
| fría | 189 | 24 743 | 0.148 s | 61 934 | 0.357 s |
| *SoPlex* | *189* | *9 666* | *0.438 s* | — | — |

Acá la base propia sí tiene de dónde heredar, y se nota: contra la fría ahorra
×1.72 en `op` y ×3.9 en `rd`. Pero **sigue perdiendo contra la compartida** en la
estrategia que importa (`op`: 14 390 contra 9 396), porque 39 de las 114
instancias convergen en **una sola pasada** y ahí `propia` vuelve a ser `fría`.
Sobre las 75 que hacen más de una pasada la brecha se achica pero no se da
vuelta: 11 784 contra 8 457.

Donde la base propia gana claramente es en el round-robin (`rd`: 15 943 contra
33 764, ×2.1), que es coherente: con presupuestos chicos la base compartida se
revuelve a cada cambio de objetivo sin llegar a ningún lado, mientras que la
propia conserva el avance de **esa** cota entre rondas. Pero aun así el
round-robin con base propia (15 943) sigue costando más que el DFB al óptimo con
base compartida (9 396), así que no rescata el esquema del §8.

La contracción es prácticamente la misma en las tres (34.01 % contra 33.96 %):
la política de base es una decisión de **costo**, no de poda.

### Conclusión

La base compartida no es una simplificación que haya que arreglar: es la
política correcta, y por un factor de 2 a 2.6. Una base por cota cuesta además
`2n` factorizaciones LU en memoria en lugar de una. `DFB_SX_PROPIO` en el
contractor de producción queda confirmado como lo que es: una variante a
medir, no un camino a adoptar.

## 10. Elegir la siguiente cota mirando la anterior: sí, y cierra el hueco del §6

`./duelo_dfb --orden`. La heurística ya existe y es la de **Achterberg** (cf.
Baharev), la que `CtcPolytopeHull` aplica en `choose_next_variable`. Con la
solución primal `x*` de la resolución anterior hace **dos cosas a la vez**:

- **(a) da por terminada** toda cota cuyo `x*_j` ya esté pegado al borde: si
  `x*_j = box[j].lb()` entonces `min x_j` sobre el poliedro **es**
  `box[j].lb()`, y ese LP no puede aportar nada. De acá sale que PolytopeHull
  resuelva una fracción de las `2n` cotas;
- **(b) elige la siguiente** como la de `|x*_j − borde|` **mínimo**, la más
  cercana al vértice actual — que es la localidad que debería abaratar la
  herencia de base.

Se transplanta a DFB sin tocar el simplex: solo hace falta su punto primal, que
sale de `x_B` y de la ubicación de las no básicas (`DFBSimplex::primal_solution`).
Las seis variantes corren con intercalado, hasta punto fijo, con el mismo
criterio y el **mismo tope de iteraciones** para los dos solvers.

| variante | LPs | % de `2n` | pivotes | piv/LP | tiempo | contracción |
|---|---|---|---|---|---|---|
| SoPlex, orden natural | 3 544 | 173 % | 9 666 | 2.73 | 0.447 s | 31.94 % |
| SoPlex, Achterberg | 2 162 | 106 % | 5 468 | 2.53 | 0.306 s | 32.65 % |
| SoPlex, salteo + lejos | 1 777 | 87 % | 7 182 | 4.04 | 0.291 s | 32.57 % |
| DFB, orden natural | 3 642 | 178 % | 9 396 | 2.58 | 0.098 s | 34.01 % |
| **DFB, Achterberg** | **1 962** | **96 %** | 7 435 | 3.79 | **0.076 s** | 33.88 % |
| DFB, salteo + lejos | 1 714 | 84 % | 7 077 | 4.13 | 0.074 s | 34.17 % |

**Funciona, y es el hueco del §6.** DFB pasa de resolver el 178 % de las `2n`
cotas al **96 %**: ×0.56 en LPs (media geométrica), ×0.82 en pivotes, ×0.74 en
tiempo, **con la misma contracción** (gana 1, pierde 3, empata 110). El ×3 de
volumen de LPs que el §6 señalaba como la deuda de DFB desaparece con una
heurística que no toca el simplex.

Y contra la combinación que usa producción hoy, SoPlex + Achterberg:

| DFB+Achterberg vs SoPlex+Achterberg | |
|---|---|
| LPs | ×0.91 |
| pivotes | ×1.36 |
| **tiempo** | **×0.25** |
| contracción | **33.88 %** vs 32.65 % (gana 8, pierde 1) |

### El control: ¿el (b) aporta algo, o todo es el (a)?

La heurística mezcla dos efectos, así que la variante **«salteo + lejos»** usa
la misma regla (a) pero elige la cota **más lejana** al vértice. Si la localidad
del warm start importara, «lejos» tendría que salir claramente peor.

| mismo salteo, solo cambia a cuál se salta | LPs | pivotes | piv/LP |
|---|---|---|---|
| SoPlex cerca | 2 162 | **5 468** | 2.53 |
| SoPlex lejos | 1 777 | 7 182 | 4.04 |
| DFB cerca | 1 962 | 7 435 | 3.79 |
| DFB lejos | 1 714 | **7 077** | 4.13 |

**Para SoPlex el (b) sí paga** (×0.87 en pivotes por media geométrica); **para
DFB no aporta nada** —«lejos» sale incluso marginalmente mejor en pivotes, en
LPs y en contracción (gana 4, pierde 0)—. O sea que toda la ganancia de DFB
viene del salteo (a).

Es coherente con el §9. El warm start de SoPlex descansa en la factibilidad
**primal**, que un objetivo lejano destruye, así que le conviene saltar cerca. El
de DFB descansa en la **dual**, que en un LP de variables acotadas se restaura
gratis reubicando las no básicas por el signo del costo reducido, sin importar
qué objetivo venga: la herencia ya funciona incondicionalmente y la distancia al
objetivo anterior no la modula.

Nota sobre `piv/LP`: sube con la heurística (2.58 → 3.79) sin que nada empeore.
Es un efecto de composición: lo que (a) saltea son justamente los LPs baratos
—los que ya estaban en su cota y se resolvían en cero o un pivote— así que el
promedio de los que quedan sube aunque el total baje.

### Advertencia de medición

En la primera corrida SoPlex+Achterberg daba 1 005 148 pivotes (464.92 por LP).
El 99 % venía de **una** instancia, `m4wd`, que agotaba el `max_iter = 1e6` que
yo le había dado a SoPlex mientras DFB corría con tope `20m+200`. Con el mismo
tope para los dos el número baja a 5 468. El presupuesto de iteraciones es parte
de lo que hay que igualar, como la factorización o el tableau disperso.

## 11. Estados independientes **sembrados**: la mejor configuración medida

`./duelo_dfb --semilla`. El §9 concluyó que una base por cota no rinde. Esa
conclusión estaba mal sacada: mi versión de «base propia» **arrancaba en frío**,
y como en una pasada cada ejemplar se usa una sola vez, no tenía nada que
heredar. Estaba midiendo arranque en frío disfrazado de estado independiente.

La variante que faltaba: resolver la **primera cota al óptimo** y **copiar esa
base** a los `2n` ejemplares, que después cada uno evoluciona la suya sin que las
demás se la revuelvan (`DFBSimplex::sembrar_desde`). No se copia la
factorización LU —`soplex::SLUFactor` guarda punteros a las columnas del objeto
de origen— sino la lista de columnas básicas, y cada ejemplar refactoriza sobre
las suyas.

Las ocho variantes, con intercalado, hasta punto fijo, mismo criterio y mismo
tope de iteraciones:

| variante | LPs | % `2n` | pivotes | piv/LP | tiempo | contracción |
|---|---|---|---|---|---|---|
| *SoPlex + Achterberg (producción)* | *2 162* | *106 %* | *5 468* | *2.53* | *0.306 s* | *32.65 %* |
| DFB compartida | 3 642 | 178 % | 9 387 | 2.58 | 0.098 s | 34.05 % |
| DFB propia, fría (§9) | 3 572 | 175 % | 14 383 | 4.03 | 0.139 s | 34.00 % |
| **DFB propia, sembrada** | 3 670 | 179 % | 6 924 | **1.89** | 0.095 s | **34.17 %** |
| DFB compartida + Achterberg (§10) | 1 947 | 95 % | 7 506 | 3.86 | 0.078 s | 33.95 % |
| **DFB propia sembrada + Achterberg** | **1 891** | **92 %** | **5 540** | 2.93 | **0.074 s** | 34.04 % |
| DFB propia sembrada, presupuesto 4 | 4 766 | 233 % | 5 845 | 1.23 | 0.090 s | 28.17 % |
| DFB propia sembrada + Ach, presup. 4 | 2 460 | 120 % | 4 570 | 1.86 | 0.068 s | 29.87 % |

**La siembra sola da vuelta el §9**: de 14 383 pivotes (propia fría) a **6 924**,
×0.48, y también por debajo de la compartida (9 387), ×0.74. Y baja `piv/LP` a
**1.89**, el mejor número de todo el estudio —por debajo de los 2.53 de SoPlex—,
con la contracción apenas mejor (gana 2, pierde 0 contra la compartida).

El mecanismo es claro y es el que el §9 dejaba inactivo: la semilla le da a cada
ejemplar una base **buena** de entrada, y tener estado propio evita que los otros
`2n−1` cambios de objetivo se la revuelvan. Las dos cosas hacen falta juntas:
estado propio sin semilla es peor que compartir (×1.50), y la semilla sin estado
propio es simplemente la compartida.

Combinada con el salteo de Achterberg del §10, contra lo que usa producción hoy:

| DFB propia sembrada + Ach **vs** SoPlex + Ach | total | media geom. |
|---|---|---|
| LPs | 1 891 vs 2 162 (×0.87) | ×0.90 |
| pivotes | 5 540 vs 5 468 (×1.01) | ×1.11 |
| **tiempo** | 0.074 s vs 0.306 s (**×0.24**) | **×0.35** |
| contracción | **34.04 %** vs 32.65 % | gana 8, pierde **0** |

O sea que DFB queda **a la par de SoPlex en pivotes**, resolviendo menos LPs, a
**un cuarto del tiempo**, y sin perder poda en ninguna instancia.

### La parte anytime sigue sin pagar

Las dos últimas filas truncan el presupuesto a 4 pivotes por cota. Bajan los
pivotes (5 540 → 4 570) pero suben los LPs (1 891 → 2 460, porque hacen falta más
pasadas) y **pierden poda**: 29.87 % contra 34.04 %, perdiendo en 20 instancias y
ganando en 1. Es la misma conclusión del §8, ahora con la mejor base de partida:
con el óptimo a 2.93 pivotes no hay nada que truncar.

La propiedad anytime queda donde la dejó el §5 —garantiza que un tope nunca
produzca una cota inválida, lo que permite ponerlo sin analizar el caso— y no
como palanca de velocidad.

## 12. Contracción contra presupuesto: acá la propiedad anytime sí paga

`./duelo_dfb --presupuesto`.

Las mediciones de los §8 y §11 comparaban **«truncado hasta punto fijo» contra
«exacto hasta punto fijo»**, y eso le permite al truncado gastar más presupuesto
total para compensar. Mide otra cosa, y por construcción da negativo. La pregunta
que define *anytime* es la otra: **gastado un presupuesto total `P`, ¿cuánta
contracción tengo si corto ahí?**

Se corre la mejor configuración del §11 (estados propios sembrados + Achterberg)
con presupuesto por LP `b ∈ {1, 2, 4, 8, óptimo}`, registrando (pivotes, tiempo,
perímetro) después de **cada** cota, y se lee la curva en puntos de corte
comunes, expresados como fracción de lo que gasta la variante exacta hasta su
punto fijo. SoPlex entra como referencia: también es interrumpible, pero solo
**entre** LPs, no dentro de uno.

### Presupuesto medido en pivotes

| presupuesto | b=1 | b=2 | b=4 | b=8 | óptimo | SoPlex |
|---|---|---|---|---|---|---|
| 5 % | **0.59 %** | 0.32 % | 0.17 % | 0.04 % | 0.04 % | 0.15 % |
| 10 % | **1.73 %** | 0.74 % | 0.43 % | 0.26 % | 0.18 % | 0.70 % |
| 20 % | **3.02 %** | 2.52 % | 1.72 % | 0.84 % | 0.64 % | 2.25 % |
| 30 % | **5.57 %** | 4.95 % | 3.56 % | 2.05 % | 1.92 % | 5.46 % |
| 50 % | 13.40 % | 12.19 % | 9.04 % | 6.28 % | 5.05 % | **18.35 %** |
| 75 % | 16.84 % | 19.04 % | 18.06 % | 15.45 % | 10.35 % | **28.78 %** |
| 100 % | 18.51 % | 25.53 % | 28.64 % | 30.09 % | **34.04 %** | 31.94 % |

### Presupuesto medido en tiempo (el que decide)

El eje de pivotes favorece a SoPlex, cuyos pivotes cuestan unas **4 veces más**
(56 µs contra 13 µs, §11). En tiempo:

| presupuesto | b=1 | b=2 | b=4 | b=8 | óptimo | SoPlex |
|---|---|---|---|---|---|---|
| 5 % | 0.09 % | **0.13 %** | 0.04 % | 0.04 % | 0.04 % | 0.00 % |
| 10 % | 0.22 % | **0.35 %** | 0.18 % | 0.13 % | 0.09 % | 0.00 % |
| 20 % | 0.58 % | **0.92 %** | 0.56 % | 0.33 % | 0.22 % | 0.05 % |
| 30 % | 0.95 % | **2.00 %** | 1.77 % | 0.79 % | 0.70 % | 0.14 % |
| 50 % | 1.89 % | 4.54 % | **4.62 %** | 3.81 % | 2.43 % | 0.96 % |
| 75 % | 3.60 % | 14.22 % | **17.93 %** | 12.42 % | 9.92 % | 3.03 % |
| 100 % | 10.18 % | 24.60 % | 28.40 % | 28.68 % | **34.04 %** | 6.28 % |
| 200 % | 18.51 % | 25.62 % | 28.70 % | 30.15 % | **34.04 %** | 23.39 % |

**Corregido en el §13**: truncar gana contra el orden de Achterberg, pero no es
la palanca. Sin ningún tope, cambiando solo el orden en que se recorren las
cotas, se consigue bastante más. Lo que sigue compara contra Achterberg.

**Truncar gana en todo el rango donde el presupuesto es la restricción.** A la
mitad del tiempo, `b=4` contrae 4.62 % contra 2.43 % del exacto: **×1.9**. A tres
cuartos, 17.93 % contra 9.92 %: **×1.8**. En el eje de pivotes la ventaja llega a
**×4.7** (20 %: 3.02 % contra 0.64 %). El cruce está justo en el 100 %: si el
presupuesto alcanza para terminar, el exacto gana (34.04 %).

Cara a cara sobre las 114 instancias, en el eje de tiempo, `b=2` supera
estrictamente al exacto en **36** instancias al 50 % de presupuesto, y `b=4` en
**29** al 75 %.

### El mejor truncamiento no es el más agresivo

En pivotes gana `b=1`; en **tiempo** ganan `b=2` y `b=4`. La razón es la del
§8: el costo por resolución tiene una parte fija grande —refrescar el objetivo y
sobre todo certificar `γ = λᵀĀ`, que son `m × na` operaciones de intervalo pase
lo que pase— así que `b=1` compra pocos pivotes pero paga la parte fija entera.
El óptimo práctico está en **2 a 4 pivotes por cota**.

### Y SoPlex no compite acá

En el eje de tiempo SoPlex queda por debajo de todas las variantes de DFB en
todos los presupuestos: 6.28 % cuando DFB exacto ya lleva 34.04 %, y necesita el
200 % para llegar a 23.39 %. Es la consecuencia directa de no ser interrumpible
dentro de un LP más el costo por pivote: un presupuesto que se agota a mitad de
un LP no deja nada de ese LP.

### Conclusión, y corrección de los §8 y §11

La afirmación correcta no es «truncar no paga», sino:

- **si el presupuesto alcanza para llegar al punto fijo, truncar no conviene**
  (§8 y §11, que es lo que esas secciones miden);
- **si no alcanza, truncar da de 1.8 a 4.7 veces más contracción** por el mismo
  gasto, y esa es exactamente la situación de un contractor dentro de un árbol
  de búsqueda, donde el presupuesto por nodo es una decisión y no un dato.

O sea que la propiedad anytime no es solo una red de seguridad: es la que
permite elegir el punto de la curva, y DFB es el único de los dos que tiene esa
curva.

### Nota de método

La salida de este modo son 108 líneas por instancia. Con 16 procesos escribiendo
al mismo pipe eso supera el tamaño de escritura atómica (`PIPE_BUF`, 4 KiB) y las
líneas se entrelazan y corrompen. Los modos anteriores emitían pocas líneas por
instancia y nunca lo mostraron. La corrida definitiva escribe **un archivo por
instancia**.

## 13. No hace falta truncar: la palanca es el ORDEN

Todas las variantes del §12 son intercaladas —cada cota se aplica apenas sale— y
el esquema «intercalado hasta terminar, sin tope de pivotes» ya estaba ahí: es la
serie `óptimo`. Lo que estaba mal era su **orden**.

Achterberg elige la cota cuyo `|x*_j − borde|` es **mínimo**, la más cercana al
vértice actual. Esa es, por construcción, la que **menos** va a contraer: sirve
para gastar poco cuando uno llega hasta el final —es lo que mide el §10— pero es
la peor regla posible si a uno lo van a interrumpir. Con presupuesto limitado hay
que poner adelante las cotas que más rinden.

Se comparan tres reglas, las tres sin tope y hasta el óptimo en cada LP:

- **cercana**: Achterberg, `|x*_j − borde|` mínimo;
- **lejana**: el mismo salteo, pero `|x*_j − borde|` **máximo**;
- **más ancha**: la cota de la variable de mayor diámetro actual.

| presupuesto (tiempo) | b=2 | b=4 | cercana | lejana | **más ancha** |
|---|---|---|---|---|---|
| 20 % | 0.54 % | 0.62 % | 0.22 % | 0.76 % | **0.76 %** |
| 30 % | 0.91 % | 1.77 % | 0.66 % | 1.63 % | **2.18 %** |
| 50 % | 2.09 % | 4.36 % | 2.65 % | 17.34 % | **15.08 %** |
| 75 % | 4.19 % | 16.69 % | 9.86 % | 30.46 % | **31.53 %** |
| 100 % | 12.14 % | 28.27 % | **34.04 %** | 34.00 % | 33.88 % |
| 200 % | 25.62 % | 28.70 % | 34.04 % | **34.17 %** | 34.05 % |

Y en pivotes, donde el efecto es todavía más marcado:

| presupuesto (pivotes) | b=2 | b=4 | cercana | lejana | **más ancha** |
|---|---|---|---|---|---|
| 20 % | 2.52 % | 1.72 % | 0.64 % | 5.24 % | **6.25 %** |
| 50 % | 12.19 % | 9.04 % | 5.05 % | 23.17 % | **23.82 %** |
| 75 % | 19.04 % | 18.06 % | 10.35 % | 29.84 % | **31.04 %** |

**Sin truncar nada se consigue mucho más que truncando.** A la mitad del tiempo,
`más ancha` contrae 15.08 % contra 4.36 % del mejor truncado (**×3.5**) y 2.65 %
de Achterberg (**×5.7**). A tres cuartos llega a 31.53 %, que es el **93 %** de
toda la poda disponible.

Cara a cara contra Achterberg en el eje de tiempo, `más ancha` gana en **66**
instancias y pierde en 8 al 75 % de presupuesto; al 50 %, gana 50 y pierde 7.
Como mejor serie de las seis, al 75 % lo es en **94 de 114** (contra 56 de `b=4`
y 47 de Achterberg).

**Y no se paga nada por ello.** En el 100 % las tres reglas terminan iguales
—34.04 %, 34.00 %, 33.88 %— y en el 200 % `lejana` queda incluso arriba (34.17 %).
El orden **adelanta** la poda, no la cambia.

### Conclusión: el §12 medía contra la regla equivocada

La conclusión correcta es más simple que la del §12 y no necesita topes:

- **Intercalado, cada LP hasta el óptimo, sin tope de pivotes**, con las cotas
  recorridas **de la variable más ancha hacia la más angosta**.
- La propiedad anytime sigue siendo lo que lo habilita: el proceso se puede
  cortar en cualquier momento y la caja que hay es válida. Pero lo que llena la
  curva temprano es el **orden**, no truncar los LPs.
- Achterberg queda como lo que es: la regla que **minimiza el costo total** si se
  va hasta el final (§10), y la peor si se corta antes. Son dos objetivos
  distintos y conviene elegir la regla según cuál rige.

## 14. La estrategia, y cómo queda contra producción

El §13 dejaba una duda: el orden «más ancha» gana bajo presupuesto, pero
Achterberg era el que minimizaba el costo **total** (§10). Medido en la
configuración completa, **no hay disyuntiva**: el orden que rinde antes es
también el que llega más barato.

### La estrategia

1. **`2n` estados independientes**, uno por cota, **sembrados** copiando la base
   de la primera cota resuelta al óptimo (§11).
2. **Intercalado**: cada cota se aplica apenas sale y la siguiente parte de la
   caja ya apretada. Punto fijo: se repite mientras una pasada rinda ≥ 0.01 %.
3. **Cada LP hasta el óptimo, sin tope de pivotes** (§13).
4. **Salteo de Achterberg**: se da por terminada toda cota cuyo `x*_j` ya esté
   pegado al borde (§10, parte (a)).
5. **Orden de la más lejana primero**: entre las que quedan, la de
   `|x*_j − borde|` **máximo** —o la variable de mayor diámetro, que empata—, no
   la de mínimo como hace Achterberg (§13).

### Contra producción

Las 114 instancias, misma linealización, mismo criterio de punto fijo, mismo tope
de iteraciones:

| | LPs | % `2n` | pivotes | piv/LP | tiempo | contracción |
|---|---|---|---|---|---|---|
| *producción: SoPlex + Achterberg* | *2 162* | *106 %* | *5 468* | *2.53* | *0.307 s* | *32.65 %* |
| DFB, base compartida, sin orden | 3 642 | 178 % | 9 387 | 2.58 | 0.099 s | 34.05 % |
| DFB, sembrada + Achterberg (§11) | 1 891 | 92 % | 5 540 | 2.93 | 0.075 s | 34.04 % |
| **DFB, sembrada + más lejana** | **1 693** | **83 %** | **5 145** | 3.04 | **0.073 s** | **34.17 %** |
| DFB, sembrada + más ancha | 1 688 | 83 % | 5 229 | 3.10 | 0.072 s | 34.05 % |

| DFB (más lejana) **vs** producción | total | media geom. |
|---|---|---|
| LPs | ×0.78 | ×0.81 |
| pivotes | ×0.94 | ×1.01 |
| **tiempo** | **×0.24** | ×0.35 |
| contracción | **34.17 %** vs 32.65 % | gana **8**, pierde **0** |

DFB resuelve **menos LPs** (83 % de las `2n` contra 106 %), gasta **los mismos
pivotes**, tarda **un cuarto del tiempo** y no pierde poda en ninguna instancia.

Y encima tiene la curva del §13, que producción no puede tener: a la **mitad**
del tiempo lleva **17.34 %** de contracción contra **0.97 %** de SoPlex, porque un
simplex interrumpido a mitad de un LP no deja nada de ese LP.

### Solidez

Verificado con el control riguroso del §8 —toda cota inferior válida tiene que
quedar por debajo del extremo superior del encierro certificado de
Neumaier–Shcherbina, que contiene el óptimo verdadero—: **0 violaciones** en las
2 046 cotas del duelo aislado y en las 594 cotas vinculantes de la estrategia
sombra. Y 0 testigos expulsados en el régimen sin intercalado (§3).

### Lo que falta, y es lo importante

**Todo esto es sin propagación**: una sola linealización, sobre la caja raíz de
cada instancia. Mide el contractor lineal contra el contractor lineal, que era el
punto, pero **no es la medición en el árbol**.

La última medición en árbol es anterior a todo lo de este documento y daba
producción 137 instancias resueltas contra 120–123 de DFB. Aquellas
configuraciones de DFB **no tenían nada de esto**: ni estados sembrados, ni
salteo, ni orden. Hay que rehacerla con la estrategia del §14 puesta en el hueco
exacto de `CtcPolytopeHull`, que es el experimento que decide.

## 15. Round robin de un pivote: converge, pero cuesta ×3

Tercera estrategia, en el extremo opuesto a la secuencial: **una vuelta = un
pivote en cada una de las `2n` cotas**, certificando y aplicando después de cada
pivote, ciclando hasta que ninguna cota tenga nada que pivotear. Es el
intercalado en su granularidad más fina: cada pivote ve la caja que apretaron los
`2n−1` anteriores, en vez de resolver una cota entera sobre una caja congelada.
Sólo es posible porque DFB es anytime.

### Costo total hasta converger

| estrategia | llamadas | pivotes | pivotes/llamada | tiempo | contracción |
|---|---|---|---|---|---|
| secuencial, Achterberg | 1 838 | 5 491 | 2.99 | 0.070 s | 34.04 % |
| **secuencial, más lejana** | **1 695** | **5 135** | 3.03 | **0.067 s** | **34.17 %** |
| secuencial, más ancha | 1 691 | 5 169 | 3.06 | 0.067 s | 34.05 % |
| **round robin, 1 pivote** | **34 266** | 7 679 | **0.22** | **0.336 s** | **34.17 %** |
| *SoPlex + Achterberg* | — | *5 291* | — | *0.303 s* | *32.68 %* |

**Converge a la misma caja** —34.17 %, igual que la mejor secuencial, y le gana
en 2 instancias sin perder en ninguna— pero cuesta **×1.4 en pivotes** y
**×3.1 en tiempo**.

### Por qué: el costo está en las llamadas, no en los pivotes

Las cifras lo dicen sin ambigüedad: el round robin hace **34 266 llamadas** a
`solve` contra 1 695 de la secuencial, o sea **×20**, y sólo **0.22 pivotes por
llamada**. Es decir que el **78 % de sus llamadas no da ningún pivote**: revisan
una cota que ya está en su óptimo, comprueban que no hay nada que hacer, y pagan
igual el costo fijo por resolución —refrescar el objetivo, recalcular el punto
primal, y sobre todo certificar `γ = λᵀĀ`, que son `m × na` operaciones de
intervalo pase lo que pase.

Es el mismo mecanismo que hundió al round robin de presupuesto creciente (§8) y a
las variantes truncadas (§12), ahora aislado con un número: **en DFB el costo no
lo ponen los pivotes sino las resoluciones**. La secuencial gasta 39.4 µs por
llamada y hace pocas; el round robin gasta 9.8 µs y hace veinte veces más.

### Curva de convergencia

Con el presupuesto medido en tiempo, y el 100 % siendo lo que le cuesta a la
mejor secuencial llegar a su punto fijo:

| presupuesto | b=4 | cercana | **lejana** | más ancha | round robin | SoPlex |
|---|---|---|---|---|---|---|
| 30 % | 0.97 % | 0.66 % | **2.02 %** | 1.69 % | 0.97 % | 0.09 % |
| 50 % | 4.78 % | 2.56 % | **15.63 %** | 15.08 % | 4.04 % | 0.66 % |
| 75 % | 15.11 % | 11.26 % | **29.46 %** | 29.04 % | 10.29 % | 3.35 % |
| 100 % | 27.05 % | 27.63 % | **34.17 %** | 34.05 % | 20.13 % | 6.04 % |
| 200 % | 28.70 % | 34.04 % | 34.17 % | 34.05 % | 29.22 % | 20.77 % |

El round robin queda por debajo de la secuencial en todo el rango de tiempo.
En el eje de **pivotes** sí es competitivo —10.93 % contra 20.78 % al 50 %, pero
34.17 % al 200 %, empatando— lo que confirma que su problema no es hacer trabajo
de más sino pagar el costo fijo de más.

### Lo que sí saldría de acá

El arreglo es evidente y ya tiene nombre: **no revisar una cota que no cambió**.
Una cota vuelve a tener algo que pivotear sólo si se movió alguna variable que
entra en su LP, así que en lugar de ciclar sobre las `2n` habría que mantener una
cola disparada por los cambios. Eso es exactamente propagación, y es lo que hace
`CtcDFBPropag` — o sea que el round robin de un pivote es la versión ingenua de
algo que el contractor ya tiene, y que este documento apagó a propósito para
aislar el método. Medirlo es el trabajo del árbol, junto con el §14.

## 16. Las dos mejoras al round robin: saltear y certificar fuera del bucle

El §15 dejó dos arreglos evidentes. Los dos funcionan; ninguno cambia el ranking.

### (a) No volver a tocar una cota que ya está en su óptimo

Una cota que llegó a su óptimo no puede mejorar mientras la caja no cambie. Se
lleva un contador de versión `ver` que se incrementa **sólo cuando la caja
cambia de verdad**, y cada cota registra la versión en la que se la encontró
óptima; mientras coincidan, se la saltea sin siquiera llamar a `solve`. Es
exacto: no se pierde ninguna contracción.

| round robin | llamadas | pivotes | piv/llamada | tiempo | contracción |
|---|---|---|---|---|---|
| ingenuo (§15) | 34 266 | 7 679 | 0.22 | 0.337 s | 34.17 % |
| **+ salteo** | **18 573** | 7 679 | 0.41 | **0.216 s** | 34.17 % |

**×0.58 en llamadas, ×0.69 en tiempo, los mismos pivotes y la contracción
idéntica en las 114 instancias** (gana 0, pierde 0). Gratis.

### (b) Certificar fuera del bucle de pivoteo

Certificar `γ = λᵀĀ` son `m × na` operaciones de intervalo y es el grueso del
costo fijo. El round robin lo pagaba **después de cada pivote**.

La versión correcta tiene dos partes, y el primer intento fallaba por omitir la
segunda:

1. durante el bucle se certifica **sólo cuando una cota termina**;
2. al agotarse el presupuesto se hace un **barrido final** certificando todas las
   que quedaron a medias, con el `λ` que tengan.

Sin el barrido, los pivotes ya dados y no certificados se tiran — por eso la
variante `rr_salteo_certfin` del §15 salía peor (12 474 pivotes contra 7 679).
El barrido es válido justamente porque DFB es anytime.

Con las dos partes, y contando el costo del barrido:

| cupo de pivotes | pivotes | % del secuencial | tiempo | % del secuencial | contracción |
|---|---|---|---|---|---|
| 50 % | 2 662 | 52 % | 0.052 s | 60 % | 14.88 % |
| 75 % | 3 925 | 76 % | 0.075 s | 87 % | 23.42 % |
| 100 % | 5 194 | 101 % | **0.104 s** | 119 % | 28.29 % |
| 200 % | 7 340 | 143 % | 0.167 s | 193 % | 34.05 % |
| **secuencial (§14)** | **5 146** | 100 % | **0.087 s** | 100 % | **34.16 %** |

El round robin pasa de **0.337 s a 0.104 s** a igualdad de pivotes: certificar
por pivote era cerca del 70 % de su costo.

### Por qué aun así no gana

Porque **en el esquema secuencial la certificación ya está fuera del bucle de
pivoteo**: se paga una sola vez por LP, al terminarlo. La optimización sólo
existía en el round robin. Corregida, el round robin llega a la misma caja
—pierde en 2 de 114, gana en 0— pero necesita **143 % de los pivotes y 193 % del
tiempo**, y su curva queda por debajo en todo el rango: al 60 % del tiempo da
14.88 %, donde el secuencial al 50 % ya da 15.63 % (§13).

La conclusión del §14 no cambia: **secuencial, cada LP hasta el óptimo, orden de
la más lejana primero**. Lo que el §16 agrega es que las dos optimizaciones son
correctas y baratas, y que la de saltear cotas convergidas conviene igual —en el
esquema secuencial es lo que ya hace el salteo de Achterberg.

## 17. Selección de la próxima cota en el round robin

`./duelo_dfb --seleccion`. Las tres variantes llevan el salteo y la
certificación fuera del bucle del §16, y difieren sólo en cómo eligen.

| variante | llamadas | re-evaluaciones | pivotes | tiempo | contracción |
|---|---|---|---|---|---|
| orden fijo `0..2n−1` | 20 438 | — | 12 474 | 0.231 s | 34.05 % |
| **orden por variable más ancha** | 16 044 | — | **8 078** | **0.188 s** | 34.08 % |
| + re-evaluación de `γ` | 16 044 | 24 849 | 8 078 | 0.194 s | 34.08 % |

### El orden sí mejora la selección

Recorrer cada vuelta de la variable más ancha hacia la más angosta —la regla del
§13— da **×0.65 en pivotes** (12 474 → 8 078) y ×0.81 en tiempo, con la
contracción igual o mejor (gana 1, pierde 0). La media geométrica es sólo ×0.96,
o sea que el grueso viene de las instancias difíciles: en `ex7_3_6` son **35
pivotes contra 4 472**.

Es el mismo resultado del §13 en otro esquema: el orden es la palanca.

### La re-evaluación de `γ`: idea correcta, lugar equivocado

`γ = λᵀĀ` **no depende de la caja**: es función sólo de `λ`, y `γ·z = 0` vale
para cualquier caja que contenga al conjunto factible. Entonces un `γ` ya
calculado puede re-evaluarse sobre la caja apretada y dar una cota mejor a costo
`O(na)` en vez de `O(m·na)`, sin pivotear ni recertificar. Sirve además como
test exacto de «esta cota tiene algo para dar» y como prioridad.

**Medido, no aporta nada**: 24 849 re-evaluaciones, **cero** contracción extra,
las mismas llamadas y los mismos pivotes, y ×1.01 en tiempo.

La razón es que en el round robin no hay hueco donde meterla. La cota se aplica
justo después de certificar, y cuando la caja cambió lo suficiente como para que
la re-evaluación diera algo, esa cota ya quedó sucia y se re-resuelve igual. La
re-evaluación nunca entrega más que la re-resolución que va a ocurrir de todos
modos: **duplica trabajo en vez de reemplazarlo**.

Un detalle de implementación que costó: la re-evaluación barata aprieta la caja
en cantidades infinitesimales, cada una de las cuales volvía a marcar sucias las
`2n` cotas. El bucle pasaba de 5 vueltas a 4 860. Hace falta contar como cambio
sólo las mejoras **relativas y significativas** (`1e-6`).

### Lo que queda por probar

- **La re-evaluación en el esquema secuencial**, donde sí reemplazaría trabajo
  caro: la pasada 2 re-resuelve las `2n` cotas desde cero; si primero barriera
  los `γ` guardados y sólo re-resolviera las que ese barrido agotó, se ahorrarían
  resoluciones completas.
- **Ensuciar más fino**: hoy cualquier cambio de la caja ensucia las `2n` cotas.
  El test exacto es marcar sucia la cota `j` sólo si cambió alguna `z_i` con
  `γ_i ≠ 0` en su último certificado — información que ya está guardada.
- **Presupuesto adaptativo por visita**: dar más de un pivote a las cotas que
  vienen rindiendo, en vez de uno fijo para todas.

## 18. Las dos ideas sobre el esquema secuencial: casi neutras

`./duelo_dfb --secmejor`. Las dos ideas que el §17 dejó pendientes, ahora sobre
el esquema del §14, que es donde podían **reemplazar** trabajo caro en vez de
duplicarlo.

- **sucio fino**: una cota se ensucia sólo si cambió alguna `z_i` con
  `γ_i ≠ 0` en su último certificado. Se lleva una versión por variable y, por
  cota, el sello de esas versiones al certificar. Es exacto y mucho más fino que
  «cualquier cambio ensucia las `2n`».
- **γ barato**: antes de cada pasada cara, un barrido que re-evalúa los `γ`
  guardados sobre la caja de ahora —`O(na)` por cota contra `O(m·na)` de
  certificar, y cero pivotes— iterado hasta que no rinda. Recién después se pagan
  las resoluciones.

| variante | llamadas | saltadas | re-evaluaciones | pivotes | contracción |
|---|---|---|---|---|---|
| `sec` (§14) | 1 574 | 0 | 0 | 5 146 | 34.16 % |
| `sec_sucio` | **1 455** | 175 | 0 | 5 246 | 34.16 % |
| `sec_gamma` | 1 581 | 0 | 861 | 5 151 | **34.17 %** |
| `sec_ambos` | **1 453** | 184 | 861 | 5 251 | **34.17 %** |

**El sucio fino evita el 11.2 % de las llamadas** y la contracción queda idéntica
en las 114 instancias, pero cuesta **2 % más pivotes**: saltear una cota en la
pasada 2 la deja para más tarde, con la caja distinta. El neto es
indistinguible.

**El barrido de `γ` no aporta**: 861 re-evaluaciones, las mismas llamadas, los
mismos pivotes, y una sola instancia donde mejora la contracción.

La razón de fondo es la misma para las dos: **el esquema secuencial casi no tiene
segunda pasada que optimizar**. 39 de las 114 instancias convergen en **una sola
pasada**, y el salteo de Achterberg ya elimina el grueso del trabajo redundante
en las demás. No queda margen donde estas dos ideas puedan morder.

### Una trampa de medición que casi produce un número falso

La primera lectura daba **×0.68 en tiempo** para las tres variantes nuevas —
incluida `sec_gamma`, que no cambia ni las llamadas ni los pivotes. Eso no podía
ser una mejora, y no lo era: es la **posición** en el bucle de variantes. La que
corre primera paga la puesta en marcha (asignaciones, primer toque de los buffers
de la factorización LU).

| orden de las variantes | `sec` | `sec_sucio` | `sec_gamma` | `sec_ambos` |
|---|---|---|---|---|
| normal | **0.085 s** | 0.064 s | 0.065 s | 0.064 s |
| invertido (`SECMEJOR_REV=1`) | 0.064 s | 0.064 s | 0.064 s | **0.085 s** |

El sobrecosto sigue a la posición, no a la variante. En un arnés que corre varias
variantes en el mismo proceso hay que **invertir el orden** y comprobar que el
efecto no lo siga; si lo sigue, el tiempo de ese arnés no sirve y hay que quedarse
con las métricas deterministas (llamadas, pivotes, contracción), que es lo que
hace la tabla de arriba.

## 19. En el árbol: DFB en el hueco exacto de `CtcPolytopeHull`

Nuevo contractor [`CtcDFBHull`](ibex_CtcDFBHull.h): la estrategia del §14 con la
misma interfaz y la misma estructura que `CtcPolytopeHull` —**una sola pasada**
por llamada, el punto fijo lo pone el `CtcFixPoint` de afuera— enchufado en su
hueco exacto vía `--lr=dfbhull`. La configuración de producción es
`--filtering=acidhc4 --lr=xn` (`compo` pide ibex-affine, que no está compilado).
139 instancias, 40 s de límite.

### Primer intento: catastrófico

| | resueltas | celdas (geom) | cpu (geom) |
|---|---|---|---|
| producción | **128** | 1.000 | 1.000 |
| `dfbhull` | 87 | **×1.584** | ×1.282 |

Celdas ×1.58 quiere decir que DFB **poda menos en el árbol**, justo al revés de
lo que decía el duelo aislado. `ex14_1_5`: 3 870 celdas contra 16.

Ni el orden ni el salteo lo explicaban: con `DFBH_ORDEN=cerca|nat` y
`DFBH_NOSKIP=1` la instancia queda entre 554 y 3 870 celdas, siempre dos órdenes
de magnitud peor.

### La causa: DFB tiraba todas las pruebas de vacío

Con `DFBH_AB=1` cada nodo corre **también** `CtcPolytopeHull` sobre una copia de
la caja:

```
[ab] nodos=200  contrac media PH=0.8625 DFB=0.0477 | vacios PH=163 DFB=0 | gana PH=16 DFB=13
```

**PolytopeHull prueba la caja vacía en 163 de cada 200 nodos; DFB en ninguno.**
Cuando las dos contraen están parejas —16 contra 13—, así que la calidad de las
cotas nunca fue el problema. El problema es que en un árbol la mayoría de los
nodos no se contraen: **se descartan**, y esa dimensión el duelo del §4 al §18 no
la medía, porque todos sus arneses trabajan sobre cajas factibles.

El certificado estaba a mano y sin usar. `γ = λᵀĀ` y `Ā z = 0` para toda
solución, así que **`0 ∈ γ·z` es obligatorio**; si esa suma no contiene al cero,
la caja no tiene ninguna solución de la relajación. Es un certificado de Farkas
verificado en intervalos, del mismo tipo que el test de infactibilidad de
Neumaier–Shcherbina que usa `CtcPolytopeHull`. Vale para **cualquier** `λ`, así
que se evalúa siempre, llegue o no el simplex al óptimo.

Son cuatro líneas, y cambian esto:

| instancia | producción | dfbhull antes | dfbhull después |
|---|---|---|---|
| `ex14_1_5` | 16 | 3 870 | **20** |
| `house` | 240 | 4 994 | **142** |
| `ex8_1_8` | 48 | 1 396 | **78** |
| `ex6_1_2` | 48 | 812 | **36** |

### Resultado

| | resueltas /139 | celdas (geom) | cpu (geom) |
|---|---|---|---|
| producción `acidhc4+xn` | **128** | 1.000 | 1.000 |
| **`acidhc4+dfbhull`** | **128** | **×1.006** | **×0.951** |

- **Las mismas 128 instancias**, exactamente: ninguna que resuelva una y la otra
  no, en ninguna de las dos direcciones.
- **Celdas empatadas**: 214 518 contra 213 126 en total, ×1.006 por media
  geométrica. Gana 37, pierde 43, empata 48.
- **CPU ×0.951**: 254.7 s contra 268.3 s.
- **0 óptimos incompatibles** en las 128.

Es la primera vez que DFB iguala a producción en el árbol. La ventaja de ×0.24 en
tiempo de LP del §14 se diluye a ×0.95 en el total, porque el contractor lineal
es sólo una parte del costo del nodo: ACID y HC4 se llevan el resto.

Lo mejor y lo peor: `ex7_3_2` ×0.23, `ex4_1_5` ×0.48 con menos celdas (114 contra
142); `ex5_4_3` ×5.12 y `ex5_3_2` ×5.08, las dos con más celdas, que es donde hay
que mirar después.

### La lección de método

El duelo sin propagación (§4–§18) midió con cuidado **una** dimensión —cuánto
contrae y cuánto cuesta sobre cajas factibles— y sobre ésa la conclusión era
correcta y se sostiene. Pero en un árbol de búsqueda la dimensión que domina es
otra: **cuántos nodos se descartan**. Ningún arnés de este documento la medía, y
por eso un defecto que costaba 41 instancias resueltas quedó invisible hasta
ponerlo en el árbol.

## 20. Corte temprano en el árbol: medido, y no paga

Con DFB ya a la par de producción (§19), la propiedad anytime por fin tiene dónde
aplicarse: cortar la pasada antes de recorrer las `2n` cotas. Dos criterios.

**Rendimiento decreciente** (`DFBH_TOL`, `DFBH_PAC`). El §13 mostró que con el
orden de la más lejana la ganancia está muy adelantada, así que lo natural no es
un tope fijo sino cortar cuando deja de rendir: una cota «rinde» si contrae `z_k`
más de `τ` relativo, y se corta tras `pac` cotas consecutivas que no rindieron.
La paciencia hace falta porque una cota sola puede dar cero y la siguiente mucho.

**Tope de pivotes por nodo** (`DFBH_MAXPIV`): la palanca anytime cruda, para
comparar.

Sobre las 117 instancias que resuelven todas las configuraciones:

| config | celdas | ×xn | cpu | ×xn | resueltas /139 |
|---|---|---|---|---|---|
| producción `xn` | 165 856 | 1.000 | 180.8 s | 1.000 | 128 |
| **`dfbhull` sin corte** | **165 140** | **0.996** | **145.4 s** | **0.804** | **128** |
| τ = 0.001 | 201 884 | 1.217 | 160.7 s | 0.889 | 128 |
| τ = 0.01 | 202 766 | 1.223 | 158.9 s | 0.879 | 129 |
| τ = 0.1 | 226 478 | 1.366 | 166.5 s | 0.921 | 126 |
| tope 50 pivotes | 168 978 | 1.019 | 151.0 s | 0.835 | 127 |
| tope 20 pivotes | 220 050 | 1.327 | 205.8 s | 1.138 | 118 |

**Todas las variantes con corte quedan peor que sin corte**, en celdas y en cpu.

### Por qué, y por qué ningún `τ` lo arregla

Es aritmético. El contractor lineal es **~20 % del costo del nodo**: el trabajo de
LP es ×0.24 contra SoPlex (§14) y el total del árbol ×0.80, o sea que ACID y HC4
se llevan el resto. Cortar el contractor a la mitad ahorra ~10 % del nodo;
perder 20 % de poda cuesta 20 % más nodos. El cambio es malo por construcción, y
mover `τ` sólo desplaza el punto sobre una curva que nunca cruza.

A esto se suma lo del §19: cortar antes también pierde oportunidades de **probar
vacío**, que es lo que de verdad poda el árbol.

Y no rescata siquiera las instancias donde DFB es más lento, que es donde debería
haber ayudado más:

| instancia | `xn` | sin corte | tope 50 | τ = 0.01 |
|---|---|---|---|---|
| `ex5_3_2` | 58 / 0.64 s | 70 / 3.17 s | 244 / 2.81 s | 108 / 2.59 s |
| `ex5_4_3` | 14 / 0.15 s | 32 / 0.77 s | 174 / 0.91 s | 38 / 0.59 s |
| `launch` | 68 / 2.05 s | 84 / 6.32 s | — | 1 080 / 28.1 s |

### Conclusión

La propiedad anytime queda donde la dejó el §5: **red de seguridad, no palanca de
rendimiento**. Garantiza que un tope de presupuesto nunca produzca una cota
inválida, lo que permite ponerlo sin analizar el caso; pero usarlo como criterio
de sintonía pierde, y pierde por un margen que no depende de la sintonía.

La configuración recomendada sigue siendo la del §19, **sin corte**.

### Dos pendientes

**Una discrepancia sin explicar**: con `DFBH_MAXPIV=20`, `ex8_5_6-1` devuelve un
óptimo de `0.2399` donde todas las demás configuraciones —producción incluida—
dan `0.29444`. Es la única discrepancia de todo el barrido. La contracción y el
test de vacío están certificados en intervalos y un tope de pivotes sólo puede
*reducir* la poda, así que no debería poder producirla; la hipótesis es que el
buscador de loup acepte otro punto bajo un orden de búsqueda distinto. **No está
verificada.**

**El costo por nodo en instancias grandes**: `ex5_3_2` (n=23, m=30) cuesta ×5 por
nodo. El sospechoso concreto es la siembra: `una_pasada` llama a
`DFBSimplex::load` **`2n` veces por nodo**, o sea `O(2n·m·na)` sólo para cargar la
relajación en los `2n` ejemplares —unas 73 000 operaciones por nodo en esa
instancia— antes de dar un solo pivote. Copiar la carga ya hecha en vez de
recomputarla en cada ejemplar es la optimización evidente y no está intentada.

## 21. Revisando lo descartado, con cajas apretadas: la causa real del ×5

Varias cosas se descartaron midiendo sobre la **caja raíz**, donde el LP termina
en ~3 pivotes. Dentro del árbol el régimen es otro, así que se revisaron.

### Lo que se revisó, y no

| revivido | resultado en el árbol |
|---|---|
| `DFB_SX_EXACT` (medido antes con «4× menos pivotes») | **no reproduce**: 31.6 pivotes/LP contra 29.9. Sin efecto. |
| Tope de pivotes **por LP** (no por nodo, §20) | `ex5_3_2` de 70 a **832 celdas**. Mucho peor. |
| Base compartida en vez de sembrada | 70 → 98 celdas, 3.14 → 4.58 s. Peor. |
| Orden de Achterberg en vez de la más lejana | 70 → 76 celdas. Levemente peor. |
| **Devex** en la elección de la fila que sale | 29.91 → **29.21** pivotes/LP. Un 2 %: no es la causa. |

### Dónde está el tiempo, medido

Cronómetro dentro del contractor (`DFBH_CRONO=1`):

| instancia | n | m | carga | **simplex** | certificación | otro |
|---|---|---|---|---|---|---|
| `ex5_3_2` | 23 | 49 | 9 % | **85 %** | 2 % | 7 % |
| `launch` | 39 | 75 | 18 % | **74 %** | 1 % | 8 % |
| `house` | 9 | 19 | 19 % | 57 % | 4 % | 22 % |

**Corrige dos hipótesis previas.** La siembra (`load` llamado `2n` veces por
nodo) es sólo el 8–19 %, no el cuello de botella que se había supuesto. Y
certificar `γ = λᵀĀ`, que en el round robin era el grueso del costo fijo (§15–§16),
acá es el **1–4 %**: ahí había 20 veces más llamadas, cada una con casi ningún
pivote, y por eso el reparto era el opuesto.

### A/B contra SoPlex sobre las mismas cajas del árbol

| instancia | SoPlex iters/LP | DFB pivotes/LP | µs por pivote |
|---|---|---|---|
| `ex5_3_2` | 6.12 | **27.73** | 22.1 vs 28.7 |
| `launch` | 11.40 | **28.37** | 25.1 vs 53.9 |
| `alkylbis` | 4.07 | **10.07** | 24.2 vs 19.8 |

El costo **por pivote** es comparable; lo que se dispara es la **cantidad**: ×2.5
a ×4.5. En la caja raíz era ×1.2 (3.03 contra 2.53). El problema aparece sólo en
el régimen apretado y es de elección de pivote, no de implementación.

### La causa: estancamiento por degeneración

Contando los pivotes cuya columna entrante tiene costo reducido ≈ 0 —o sea que el
paso dual **no mejora el objetivo**:

| instancia | pivotes/LP | pivotes degenerados |
|---|---|---|
| `ex5_3_2` | 29.91 | **78 %** |
| `launch` | 32.47 | **88 %** |
| `house` | 6.05 | 71 % |
| `ex7_3_2` | 1.85 | 70 % |

Entre el 78 y el 88 % de los pivotes son degenerados. Eso explica por qué devex
no sirvió: devex arregla **qué fila sale**, y acá el problema es que el paso, una
vez elegida la fila, tiene longitud cero. En la caja raíz el LP terminaba antes
de que la degeneración importara.

### Lo que corresponde hacer

El remedio estándar para exactamente esto es el **ratio test dual de paso largo
con cambio de cota** (*bound-flipping ratio test*, Fourer / Maros / Koberstein):
en vez de detenerse en el primer punto de quiebre, se recorren varios cambiando
de cota las no básicas acotadas que se cruzan, y se da **un** paso dual largo en
lugar de muchos de longitud cero. Es lo que implementa SoPlex, y el caso ideal
para aplicarlo es justamente el nuestro: **todas** las variables de la relajación
están acotadas.

No está implementado. Es la única mejora identificada que ataca el ×5 sin ceder
poda, que es el intercambio que el §20 mostró perdedor.

## 22. Resumen (SUPERADO — ver el §45)

Este resumen quedó desactualizado: sus números fueron corregidos por el §38
(el «×0.925 en celdas» venía de una sola instancia), el §41 (la ventaja no era
específica de `acidhc4`) y el §43 (la cpu mejoró a ×0.830). **El resumen válido
es el §45**, al final del documento.

## 23. Paso largo activado sólo en relajaciones grandes: negativo

El §21 dejó el ratio test de paso largo en una posición rara: mejor per-LP
(×0.77 pivotes, mismo certificado) y peor en el árbol (celdas ×1.35). La
sospecha era que la pérdida estuviera confinada a instancias chicas, así que se
midió por franja de `m`:

| franja de `m` | n | ventaja en celdas | ventaja en cpu |
|---|---|---|---|
| 0–11 | 56 | **−0.150** | **−0.094** |
| 12–19 | 23 | −0.018 | +0.069 |
| 20–29 | 5 | −0.134 | +0.023 |
| 30–49 | 1 | −0.028 | +0.248 |
| ≥ 50 | 1 | **+0.357** | **+0.306** |

La sospecha era correcta: las catástrofes están todas en `m` chico —`ex9_2_5`
(m=10) pasa de 6 a 10 668 celdas, `ex7_2_6` (m=3) de 62 a 6 160—. Así que se
agregó `DFB_SX_BFRT_M=k`, que lo activa sólo cuando `m >= k`.

| config | resueltas | celdas ×xn | cpu ×xn | ventaja cpu contra `dfbhull` |
|---|---|---|---|---|
| `dfbhull` sin paso largo | 128 | 0.998 | 0.895 | — |
| paso largo siempre | 107 | 1.414 | 1.072 | −0.028 |
| `m ≥ 12` | 121 | 1.140 | 0.993 | −0.020 |
| `m ≥ 20` | 129 | 1.020 | 0.901 | −0.018 |
| `m ≥ 30` | 128 | **0.998** | **0.870** | **+0.006** |

**No hay ganancia.** El único umbral que no empeora es `m ≥ 30`, y ahí la ventaja
contra `dfbhull` es +0.006 —un empate— porque a ese umbral el paso largo casi
nunca se activa: son un puñado de instancias en todo el banco. El umbral converge
a «apagado».

El 129 de `m ≥ 20` es ruido de borde: gana `bearing` y `ex14_2_7` terminando en
37.7 s y 38.2 s contra un límite de 40 s, y pierde `himmel16`.

Queda por omisión apagado (`DFB_SX_BFRT_M=k` para habilitarlo).

### Un dato que este barrido deja más firme

Sobre las 40 instancias que cuestan ≥ 0.1 s, **todas** las variantes de DFB
quedan por debajo de producción:

| | ventaja cpu contra `xn` | geom |
|---|---|---|
| `dfbhull` | −0.092 | ×1.166 |
| `m ≥ 30` | −0.085 | ×1.140 |
| `m ≥ 20` | −0.112 | ×1.235 |
| paso largo siempre | −0.156 | ×1.350 |

La paridad global del §22 se sostiene sobre las instancias baratas. En las caras,
producción sigue adelante, y la causa es la del §21: 30 pivotes por LP contra 6.

## 24. El LP duro, resuelto por los dos: el §21 estaba mal

`DFBH_DUMP=prefijo` vuelca las relajaciones de las resoluciones que pasan de
`DFBH_DUMP_MIN` pivotes, y [`probe_lp.cpp`](probe_lp.cpp) las resuelve con SoPlex
y con DFB. Doce LPs de nodos de `ex5_3_2`, de 48×71 a 50×73.

| LP | SoPlex | DFB | degenerados |
|---|---|---|---|
| lp_00 | 30 | 42 | 84 % |
| lp_01 | 39 | 50 | 88 % |
| lp_03 | 83 | 69 | 85 % |
| lp_04 | 23 | 40 | 95 % |
| lp_07 | 36 | 37 | 98 % |
| lp_11 | 33 | 32 | 94 % |
| **media** | **~35** | **~41** | 84–98 % |

**Sobre el mismo LP, SoPlex necesita casi lo mismo que DFB: ×1.2, no ×4.5.** El
«30 pivotes contra 6» del §21 comparaba **poblaciones distintas de LPs**, no el
mismo LP, y de ahí salieron cuatro remedios —devex, paso largo, refresco exacto,
umbral por `m`— aplicados a un mal inexistente. La degeneración del 84–98 % es
una propiedad de **estas relajaciones**, no de nuestro simplex: SoPlex pivotea
igual de degenerado y le cuesta lo mismo.

### Entonces, ¿de dónde sale el 6.12 de SoPlex en el árbol?

Se descartaron dos explicaciones, midiendo:

- **el orden**: `cerca` contra `lejos` da 28.09 contra 29.91 pivotes/LP. No es.
- **la base compartida**: 31.09 contra 29.91. Tampoco.

Lo que queda, y los números lo cierran:

| | en frío, los mismos LPs | en el árbol | rendimiento del arranque tibio |
|---|---|---|---|
| SoPlex | ~35 iteraciones | **6.12** | **×5.7** |
| DFB | ~41 pivotes | 27.7 | ×1.5 |

Los dos resuelven unos 30 LPs por nodo y en frío cuestan casi lo mismo. **La
diferencia entera es cuánto les rinde el arranque tibio dentro del nodo.**

### Lo que falta, y es una sola cosa

SoPlex, tras el primer LP del nodo, conserva una base **primal-factible** y
re-optimiza el objetivo nuevo en una o dos iteraciones. DFB, al cambiar de
objetivo, reubica las no básicas por el signo del costo reducido: eso restaura la
factibilidad **dual** gratis —es lo que hace barata la base inicial y no pide
fase 1— pero **destruye la primal**, así que cada cota vuelve a correr un simplex
dual largo desde una base primal-infactible.

Ya estaba escrito en un comentario de `ibex_CtcDualFeasibleBounding.cpp`, sin
haberse hecho:

> *El warm start que si sirve —el que le da a PolytopeHull sus 1.08 iteraciones
> por cota— conserva la base Y la ubicacion de las no basicas, o sea que se queda
> primal-factible, y corre un simplex PRIMAL para el objetivo nuevo. Eso pide
> implementar el primal, que es el paso siguiente.*

`DFBSimplex` tiene un `primal_step` que se usa como respaldo dentro del bucle
compuesto, pero no hay un camino que **conserve la ubicación de las no básicas al
cambiar de objetivo y re-optimice con el primal**. Ése es el trabajo, y es el
único identificado que explica el factor que separa a DFB de producción en las
instancias caras.

### Lección de método

El §21 diagnosticó comparando **promedios sobre poblaciones distintas** —los LPs
que cada estrategia elige resolver— y lo leyó como una diferencia por LP. Cuatro
experimentos salieron de ahí y los cuatro fallaron, cada uno de forma
desconcertante, porque atacaban algo que no existía. Volcar un caso concreto y
resolverlo con los dos costó menos que cualquiera de esos cuatro y dio la
respuesta.

## 25. El arreglo: no reubicar las no básicas al cambiar de objetivo

El §24 dejó identificada una sola cosa. Resultó ser **una línea**.

### La asimetría

Una base parte las variables en básicas y no básicas, cada no básica apoyada en
una de sus cotas. Hay dos factibilidades:

- **primal**: los valores de las básicas caen dentro de sus cotas. Depende de la
  base y de dónde se apoyan las no básicas, **no del objetivo**;
- **dual**: cada costo reducido tiene el signo correcto para su cota. **Sí**
  depende del objetivo.

Al pasar de la cota `k` a la `k+1` el objetivo cambia de `±e_k` a `±e_{k+1}` pero
el poliedro no se mueve. La base óptima que se acaba de obtener **sigue siendo
primal-factible** y **deja de ser dual-factible**.

SoPlex se queda con lo que se conservó y corre un **simplex primal**, que arranca
primal-factible y trabaja hacia la dual; con objetivos tan parecidos le alcanza
una o dos iteraciones. DFB hacía lo contrario: **reubicaba** cada no básica a la
cota que le da el signo correcto, restaurando la dual gratis pero **rompiendo la
primal**, y el simplex dual tenía que remontar esa infactibilidad desde cero.

Estábamos tirando lo que se conserva para recuperar lo que era barato de
arreglar. El origen del diseño es legítimo pero vale sólo para el **primer** LP
del nodo, donde no hay base previa y la base de las columnas `b` es dual-factible
sin fase 1.

### El cambio

`refrescar_desde_factorizacion(k, maximize, z, reubicar)` ya aceptaba el
parámetro, y el bucle compuesto ya tenía `primal_step`. La ruta tibia pasaba
`reubicar = true`. Ahora pasa `false` (`DFB_SX_REUBICAR=1` restaura lo anterior).

| instancia | pivotes/LP | | µs/LP | | degenerados | |
|---|---|---|---|---|---|---|
| | antes | **ahora** | antes | **ahora** | antes | **ahora** |
| `ex5_3_2` | 29.91 | **7.32** | 685 | **171** | 78 % | **23 %** |
| `launch` | 32.47 | **8.39** | 1 294 | **336** | 88 % | **14 %** |
| `alkylbis` | 11.08 | **3.52** | 143 | **48** | 78 % | 30 % |
| `house` | 6.05 | 5.73 | 46 | 47 | 71 % | 24 % |

La degeneración cae de 78–88 % a 14–30 %: **no era propiedad de las
relajaciones** —como afirmaba el §24— sino consecuencia de arrancar
primal-infactible. Tercera corrección de la cadena §21 → §24 → §25.

### En el árbol

| config | resueltas | celdas | ×xn | cpu | ×xn | ventaja cpu |
|---|---|---|---|---|---|---|
| producción `xn` | 128 | 213 126 | 1.000 | 268.3 s | 1.000 | — |
| reubicando (§19) | 128 | 214 518 | 1.007 | 254.7 s | 0.950 | +0.049 |
| **sin reubicar** | **128** | **197 064** | **0.925** | **253.3 s** | 0.944 | **+0.086** |

Y en las **67 instancias que cuestan ≥ 0.1 s**, que es donde DFB venía perdiendo:

| | ventaja cpu | geom |
|---|---|---|
| reubicando | −0.013 | ×1.033 |
| **sin reubicar** | **+0.036** | **×0.974** |

Se dio vuelta. Las mismas 128 instancias, **0 óptimos incompatibles**, y el test
unitario sigue en 528/528 cotas coincidentes con SoPlex con los invariantes
intactos.

**Por primera vez DFB explora menos nodos que producción** (×0.925) y es más
rápido también en las instancias caras.

### Contrapartida

En LPs fáciles el arranque dual era mejor: el test unitario pasa de 717 a 810
pivotes. Y hay instancias que empeoran —`launch` va de 84 a 552 celdas— aunque su
tiempo total baje por lo barato de cada LP. La ventaja por instancia en celdas es
levemente negativa (−0.025 contra la versión que reubica) mientras el total baja
un 8 %: gana mucho en las grandes y pierde poco en las chicas. Elegir la regla
según el tamaño de la relajación es lo siguiente a probar.

## 26. DFB dentro de ACID: la sonda

Idea: ACID poda rebanando —parte el dominio de una variable en trozos y pregunta
si el trozo es infactible—, y DFB podría contestar sin re-linealizar ni tocar la
base, porque (a) la linealización del nodo vale para **todo** sub-trozo y (b)
`γ = λᵀĀ` **no depende de la caja**, así que el test de vacío del §19 sobre un
trozo cuesta un producto punto de `na` intervalos.

La sonda (`DFBH_SONDA=n`, `DFBH_SONDA_NS=s`) se corre en **nodos profundos, no en
la raíz** —toda esta sesión mostró que lo medido en la caja raíz no se traslada—
y compara, sobre los trozos que ACID probaría: cuántos mata el `γ` ya calculado,
cuántos matan corridas parciales de `b` pivotes, y cuántos mata resolver al
óptimo.

### Lo que se midió

| | |
|---|---|
| `house` (n=9, m=19), 27 000 trozos de 300 nodos | `γ` guardado **0.4 %**, b=1 1.5 %, b=8 2.4 %, al óptimo **2.4 %** |
| `launch` (n=39, m=74), nodos 200–204 | re-resolviendo mata **31–66 %** de los trozos; el `γ` guardado agarra **0–26 %** de eso |
| `launch`, nodos 150–153 | re-resolviendo mata **0–1.3 %** |

**La varianza entre nodos es enorme** —0 % en el nodo 150 y 31–66 % en el 200 de
la misma instancia— así que el promedio es lo que importa, y no se pudo obtener:
la sonda es tan cara que en `launch` **no completa 5 nodos en 280 s**, contra los
**6.3 s** que tarda la corrida entera sin sonda. Un factor **×500**, y eso ya es
un resultado sobre el costo.

### Veredicto, por partes de la idea

- **No re-linealizar: correcto y gratis.** Es un hecho estructural: todo trozo es
  subconjunto de la caja para la que se linealizó.
- **No tocar la base (reusar los `γ`): no alcanza.** 0.4 % de los trozos en
  `house`, y en `launch` agarra a lo sumo un cuarto de lo que agarra re-resolver.
  El certificado de la caja completa es demasiado flojo para un trozo angosto.
- **Corrida parcial: ataca el factor equivocado.** Con `b=1` el costo ya es ~50
  pivotes **por trozo**, porque son `2n` resoluciones por trozo. Truncar los
  pivotes de cada resolución no toca el `2n`. Es la misma conclusión de los §8,
  §15 y §20: **el costo está en el número de resoluciones, no en los pivotes**.

### Lo único que sobrevive

**«Sólo algunas variables»**, pero llevado al otro factor: probar cada trozo con
**una sola cota** —la de la variable objetivo, que es la que poda contra el
loup— en vez de con las `2n`. Eso divide el costo por `2n` y deja el test en el
mismo orden que un HC4 sobre el trozo. No está medido.

La comparación que además falta, y que decide si vale algo: **cuántos de esos
trozos mata HC4**, que es quien lo hace hoy dentro de ACID. Si HC4 ya los mata,
DFB no agrega nada por caro que sea el test.

## 27. Corrección del §26: la re-evaluación sí sirve, y el objetivo solo no

El §26 concluyó que reusar los `γ` «no alcanza» a partir del 0.4 % de `house`.
Ampliando la muestra a instancias más grandes, esa conclusión estaba mal sacada:

| instancia | sólo la `γ` del objetivo | **las `2n` `γ`, 0 pivotes** | al óptimo, con pivotes |
|---|---|---|---|
| `house` (n=9, m=19), 28 800 trozos | 0.09 % | 0.4 % | 2.3 % |
| `ex14_2_1` (n=6, m=16), 30 000 trozos | 0.41 % | 0.7 % | 2.5 % |
| **`alkylbis`** (n=15, m=30), 10 500 trozos | 0.63 % | **11.9 %** | 25.6 % |

**Re-evaluar las `2n` `γ` guardadas mata el 11.9 % de los trozos de `alkylbis`
sin dar un solo pivote.** El costo son `2n` productos punto de `na` intervalos
por trozo —unas 1 350 operaciones para esa instancia—, contra los ~10 pivotes que
cuesta la variante `b=1`, cada uno mucho más caro. Es aproximadamente dos órdenes
de magnitud más barato por trozo.

Y no contradice el §17, donde la misma re-evaluación no daba nada: allí la caja
apenas se movía entre visitas, así que el `γ` viejo no tenía información nueva.
Un trozo mueve el dominio de una variable en un factor de 10: es el régimen donde
sí la tiene.

**Restringir a la variable objetivo pierde el 95 % de los aciertos** (11.9 % →
0.63 %). Los trozos se matan con los certificados de **otras** variables. Tiene
sentido: la `γ` del objetivo prueba vacío sólo cuando la cota inferior de `y`
supera el loup, que es raro; las demás lo prueban por inconsistencia geométrica
del trozo, que es frecuente.

### El diseño que sobrevive

Una pasada de DFB por nodo como hoy, guardando las `2n` `γ`, y dentro de ACID un
test de `2n` productos punto por trozo, **sin pivotear**. Falta lo que decide si
vale algo: **cuántos de esos trozos mata HC4**, que es quien lo hace hoy. Si los
mata igual, DFB no agrega nada.

### Y una advertencia sobre este documento

Es la tercera vez en la sesión que una conclusión sacada de instancias chicas se
cae al mirar instancias grandes: el §21 (pivotes por LP), el §24 (degeneración) y
ahora el §26. El patrón es siempre el mismo —medir donde es barato medir— y el
remedio también: no concluir sin una instancia grande en la muestra.

## 28. Contra HC4: DFB aporta poda que ACID hoy no obtiene

La medición que decidía si la idea del §27 vale algo: sobre los mismos trozos,
cuántos mata el test de DFB (las `2n` `γ` guardadas, **cero pivotes**) y cuántos
mata **HC4**, que es quien hace ese trabajo hoy dentro de ACID.

| instancia | DFB (`2n` γ, 0 piv) | HC4 | **sólo DFB** | sólo HC4 |
|---|---|---|---|---|
| **`alkylbis`** (n=15, m=30), 10 500 trozos | **11.9 %** | 5.5 % | **10.41 %** | 4.1 % |
| `house` (n=9, m=19), 28 800 trozos | 0.4 % | 0.3 % | 0.32 % | 0.3 % |
| `ex14_2_1` (n=6, m=16), 30 000 trozos | 0.7 % | **3.3 %** | 0.46 % | 3.1 % |

En `alkylbis` DFB mata **más del doble** de trozos que HC4, y el **10.4 % de
todos los trozos los mata sólo DFB**: poda que ACID no está obteniendo hoy, por
`2n` productos punto y ningún pivote.

Y en las tres instancias **son complementarios**: hay «sólo DFB» y «sólo HC4»
sustanciales en todas. No se reemplazan, se suman — que es el mejor caso posible
para agregarlo, porque se obtiene la unión.

En `ex14_2_1` HC4 gana (3.3 % contra 0.7 %), así que el aporte de DFB es
instancia-dependiente. La instancia donde más aporta es también la más grande de
la muestra, lo que es coherente con todo lo anterior, pero **la muestra son tres
instancias** y la más grande (`ex2_1_8`) no llegó a completar. Eso hay que
ampliarlo antes de tomarlo como general.

### El diseño a construir

Una pasada de DFB por nodo como hoy, guardando las `2n` `γ`; dentro de ACID, cada
trozo se testea con esas `γ` re-evaluadas sobre el dominio del trozo —`2n`
productos punto de `na` intervalos, sin pivotear ni re-linealizar— y se descarta
si alguna da `0 ∉ γ·z`. El test es una condición **suficiente** de vacío,
certificada en intervalos, así que no puede dar falsos positivos.

## 29. La integración en ACID: el mecanismo funciona, la integración no

Se construyó: `CtcDFBGamma`, un contractor que expone el test barato de
`CtcDFBHull` —las `2n` `γ` re-evaluadas sobre la caja, sin pivotear ni
re-linealizar— y entra como sub-contractor de ACID junto a HC4
(`DFBH_ENACID=1`). El test verifica que la caja esté contenida en aquella para la
que se linealizó, porque el objeto se comparte en todo el árbol y los `γ` pueden
venir de otra rama; sin esa comprobación sería incorrecto.

### Con el orden normal: los `γ` casi nunca aplican

ACID corre **antes** que el contractor lineal, así que ve los `γ` del nodo
anterior del recorrido, que en profundidad suele ser un hermano y no el padre:

| instancia | llamadas | rechazadas por contención | vacíos de las aplicables |
|---|---|---|---|
| `ex5_3_2` | 5 548 | **98.0 %** | 0.00 % |
| `house` | 2 015 | 85.4 % | 5.43 % |
| `alkylbis` | 1 857 | 67.1 % | **40.53 %** |

Donde aplican funcionan —40.5 % en `alkylbis`, y esa instancia mejora de 64 a 46
celdas— pero aplican tan poco que el resto es costo puro.

| config | celdas ×xn | cpu ×xn | ventaja cpu | en las ≥ 0.1 s |
|---|---|---|---|---|
| sin el test (§25) | **0.925** | **0.944** | **+0.086** | **+0.039** |
| con el test en ACID | 0.962 | 0.957 | +0.068 | +0.001 |

Más lento en 51 instancias y más rápido en 25.

### Con el contractor lineal primero: aplican siempre, pero aciertan poco

`DFBH_LRFIRST=1` invierte la composición para que ACID vea los `γ` de su propio
nodo. La contención pasa a rechazar **0 %**, y la tasa de acierto cae:

| instancia | vacíos de las aplicables | celdas | cpu |
|---|---|---|---|
| `alkylbis` | 13.78 % | 64 → **46** | 0.300 → 0.290 s |
| `launch` | 2.00 % | 552 → **492** | 15.5 → 18.1 s |
| `house` | 0.19 % | 182 → 206 | 0.484 → 0.620 s |
| `ex5_3_2` | 0.72 % | 94 → **262** | 1.45 → **6.52 s** |

El 13.78 % coincide con el 11.9 % que había medido la sonda (§28), así que la
sonda era correcta. Lo que cambió es el contexto: con el contractor lineal
corriendo primero, los `γ` se calculan sobre una caja **ya contraída**, y los
trozos que quedan son los difíciles de matar. El 40.5 % del caso anterior venía
justamente de los casos raros en que los `γ` provenían de una caja mucho más
ancha.

### Conclusión

**El mecanismo funciona y la integración no.** El test mata el 11–14 % de los
trozos sin dar un pivote, y es complementario con HC4 (§28) — eso se sostiene.
Pero en el árbol ninguna de las dos composiciones paga: con el orden normal los
`γ` casi nunca son aplicables, y con el orden invertido la composición misma
cuesta más de lo que el test ahorra.

Queda apagado (`DFBH_ENACID=1`, `DFBH_LRFIRST=1` para probarlo).

La configuración recomendada sigue siendo la del §25: `--filtering=acidhc4
--lr=dfbhull`, sin el test en ACID.

## 30. Apretar la parte `b` con el trozo: gratis, y sube mucho la tasa

El §29 usaba el trozo sólo para la parte `x` de `z = (x, b)`, y para la parte `b`
las cotas guardadas de la linealización, deducidas para la caja **ancha**. Pero
`b_j = A_j·x`: sobre el trozo, la evaluación por intervalos de esa fila da un
rango más angosto. Intersectarlo hace el test estrictamente más fuerte por
`m×nx` operaciones de intervalo y **ningún pivote**. Usar la cota vieja
desperdiciaba la información del trozo en `m` de las `na` componentes.

| instancia | tasa antes | **tasa apretando `b`** |
|---|---|---|
| `launch` | 2.00 % | **16.79 %** (×8) |
| `alkylbis` | 13.78 % | **17.12 %** |
| `ex5_3_2` | 0.72 % | 0.77 % |
| `house` | 0.19 % | 0.19 % |

Y de punta a punta, con la composición **normal**, gana en las grandes:

| instancia | sin el test | **con el test** |
|---|---|---|
| `alkylbis` | 64 celdas / 0.300 s | **40 / 0.234 s** |
| `launch` | 552 / 15.59 s | **382 / 12.94 s** |
| `ex5_3_2` | 94 / 1.45 s | 112 / 2.13 s |
| `house` | 182 / 0.487 s | 194 / 0.514 s |

### Pero el agregado sigue negativo

| config | celdas ×xn | cpu ×xn | ventaja cpu | en las ≥ 0.1 s |
|---|---|---|---|---|
| sin el test (§25) | **0.925** | **0.889** | **+0.089** | **+0.045** |
| con el test en ACID | 0.976 | 0.912 | +0.072 | +0.010 |

Más lento en 50 instancias, más rápido en 26, y una instancia menos resuelta.
Gana donde la tasa es alta y pierde donde es baja, y las de tasa baja son
mayoría.

### El autoapagado tampoco

`m` no separa los casos —`ex5_3_2` tiene m=49 y pierde—, pero la tasa sí, y se
puede medir en ejecución: tras un período de prueba, apagar el test si no rinde.
Probado (`DFBH_GWARM`), sale **peor**: `launch` tiene una tasa temprana baja
aunque la global sea 16.8 %, se autoapaga y termina en **832** celdas contra las
552 de base. Y apagarlo no restaura el comportamiento base, porque las primeras
llamadas ya desviaron la trayectoria. Queda desactivado.

### Estado

El mecanismo quedó bastante mejor de lo que estaba —apretar `b` es gratis y
multiplica por 8 la tasa en `launch`— y la integración gana claramente en las
instancias grandes. Lo que no se encontró es una regla que la active sólo ahí.
Queda apagada (`DFBH_ENACID=1`), y la configuración recomendada sigue siendo la
del §25.

## 31. No es caro: cuesta el 0.4–5 % y aun así empeora

Se midió el costo **propio** del test, separado de su efecto sobre la búsqueda:

| instancia | costo del test | tiempo del run | share |
|---|---|---|---|
| `ex5_3_2` | **0.009 s** | 2.13 s | **0.4 %** |
| `house` | 0.005 s | 0.514 s | 1 % |
| `alkylbis` | 0.010 s | 0.234 s | 4 % |
| `launch` | 0.615 s | 12.9 s | 4.8 % |

**El test no es caro.** En `ex5_3_2` cuesta 9 milisegundos y el run se encarece
0.68 s. Lo que se paga no es el test sino lo que ACID y el simplex hacen
**después**: más nodos (94 → 112, +19 %) y además nodos más caros (+47 % de
tiempo), porque al rebanar distinto ACID ve otras cajas y su propia adaptación
—cuántas variables rebana— cambia.

Es una clase de problema distinta de la que veníamos persiguiendo, y peor: **no
se arregla abaratando el test**. Se está perturbando una heurística que ya estaba
sintonizada, y la poda extra a nivel de trozo no compensa esa perturbación.

### Y una advertencia que corresponde a todo el documento

**ESTA SECCIÓN ESTABA MAL, y se corrige en el §36.** Se afirmó acá que
`LinearizerXTaylor` usa esquina aleatoria y que por lo tanto las corridas no son
deterministas. **No es cierto**: repitiendo `alkylbis` cinco veces se obtiene
64, 64, 64, 64, 64 celdas, y con otra configuración 40, 40, 40, 40, 40. La
variación que motivó la advertencia (1 580 llamadas contra 2 185) venía de un
cambio de código entre las dos corridas, no de aleatoriedad. Las comparaciones de
una sola instancia **sí** son reproducibles.

## 32. Configuración de fábrica e interruptores

**Todos los interruptores tienen por omisión el valor recomendado**, así que
`--filtering=acidhc4 --lr=dfbhull` sin variables de entorno es la configuración
de la que hablan las tablas. Los que siguen existen para reproducir las
mediciones de este documento, no para sintonizar.

### En `CtcDFBHull`

| variable | efecto | veredicto |
|---|---|---|
| `DFBH_ORDEN=cerca\|lejos\|nat` | regla de orden de las cotas (omisión: `lejos`) | §13, §21 |
| `DFBH_NOSKIP` | desactiva el salteo de Achterberg | §10: el salteo es esencial |
| `DFBH_TOL`, `DFBH_PAC` | corte temprano por rendimiento decreciente | §20: pierde |
| `DFBH_MAXPIV` | tope de pivotes por nodo | §20: pierde |
| `DFBH_LPCAP` | tope de pivotes por LP | §21: pierde |
| `DFBH_COMPARTIDA` | una sola base en vez de `2n` sembradas | §11, §21: pierde |
| `DFBH_GOAL`, `DFBH_SOLOGOAL` | cota del objetivo primero / sólo ésa | §21: marginal / ×2.6 peor |
| `DFBH_ENACID`, `DFBH_GWARM`, `DFBH_GMIN` | test barato dentro de ACID y su autoapagado | §29–§31: pierde |
| `DFBH_AB` | corre `CtcPolytopeHull` en paralelo y compara por nodo | diagnóstico (§19) |
| `DFBH_CRONO` | reparto del tiempo del contractor | diagnóstico (§21) |
| `DFBH_DUMP`, `DFBH_DUMP_MIN` | vuelca LPs duros para `probe_lp` | diagnóstico (§24) |
| `DFBH_SONDA`, `DFBH_SONDA_NS` | sonda de trozos de ACID | diagnóstico (§26–§28) |

### En `DFBSimplex`

| variable | efecto | veredicto |
|---|---|---|
| `DFB_SX_REUBICAR` | reubica las no básicas al cambiar de objetivo | **§25: la omisión (no reubicar) es el cambio más importante de la sesión** |
| `DFB_SX_BFRT`, `DFB_SX_BFRT_M` | ratio test de paso largo, siempre o si `m >= k` | §21, §23: pierde |
| `DFB_SX_DEVEX` | devex en la fila que sale | §21: 2 %, no es la causa |
| `DFB_SX_EXACT` | refresco exacto en cada iteración | §21: no reproduce su ventaja |
| `DFB_SX_SCALE` | escalado de la relajación | medido dañino en sesiones previas |
| `DFB_SX_NODOS` | conserva la base entre nodos | sin medir en esta sesión |
| `DFB_SX_CHECK` | verifica factibilidad dual al declarar óptimo | diagnóstico |

### Herramientas

| programa | para qué |
|---|---|
| `duelo_dfb` | duelo sin propagación: `--anytime`, `--rondas`, `--orden`, `--semilla`, `--presupuesto`, `--certfuera`, `--seleccion`, `--secmejor` |
| `probe_lp` | resuelve con SoPlex y con DFB los LPs volcados por `DFBH_DUMP` |
| `test_dfbsimplex` | validación contra SoPlex: 528/528 cotas, invariantes de dualidad |

## 33. `ex8_5_6-1`: no era un problema de solidez

La discrepancia que venía arrastrándose desde el §20 —`0.2399` contra `0.29444`
en toda configuración que cambie mucho la trayectoria— tiene una explicación
simple: **ninguna de las dos configuraciones resuelve esa instancia**. Las dos
terminan con

```
possibly unbounded objective (f*=-oo)
 f* in  [-inf, 0.294438229351]     (best bound)
```

O sea que el optimizador no prueba optimalidad: el número que se reporta es el
**mejor punto factible encontrado**, un incumbente heurístico. Que difiera entre
configuraciones que recorren el árbol distinto es lo esperable, no un síntoma.
Producción da `0.2944382294` y DFB `0.2944382304`; con tope de 20 pivotes,
`0.2399`.

Queda cerrado: la contracción y el test de vacío están certificados en intervalos
y un tope de pivotes sólo puede reducir la poda, así que nunca podía producir una
cota falsa. Lo que producía era otro camino de búsqueda y otro incumbente.

### Un defecto del arnés de análisis que esto destapó

Mi contador de «instancias resueltas» contaba **toda instancia que produjera una
última línea**, no las que el optimizador declara resueltas. Instancias como
`ex8_5_6-1`, donde el objetivo queda posiblemente no acotado, entraban en las
128. El número sirve igual para **comparar** las dos configuraciones —las dos se
cuentan con el mismo criterio y coinciden instancia por instancia— pero no debe
leerse como «128 problemas resueltos a optimalidad».

## 34. Medición definitiva, en serie

Para el registro, el barrido con **un solo proceso**, sin contención:

| | resueltas | celdas | ×xn | cpu | ×xn | ventaja cpu |
|---|---|---|---|---|---|---|
| producción `acidhc4+xn` | 128 | 213 126 | 1.000 | 268.8 s | 1.000 | — |
| **`acidhc4+dfbhull`** | **128** | **197 064** | **0.925** | **249.8 s** | **0.929** | **+0.090** |

En las 64 instancias que cuestan ≥ 0.1 s: ventaja **+0.046**, media geométrica
**×0.961**.

### La contención del barrido en paralelo era despreciable

Comparando instancia por instancia el mismo barrido corrido en serie y con
`-P 8`:

| | celdas serie/paralelo | cpu serie/paralelo |
|---|---|---|
| producción | ×1.000 | ×0.996 |
| DFB | ×1.000 | ×0.989 |

**Un 1 %.** Todas las tablas de este documento vienen de barridos en paralelo y
quedan validadas tal cual. La precaución de correr en serie no estaba
justificada, y costó una hora de reloj para confirmar algo ya sabido.

## 35. Salir temprano por éxito: medido antes de implementarlo, y no hay nada

La última aplicación pendiente de la propiedad anytime, y la que tenía el mejor
argumento: todos los cortes que se probaron eran por **rendirse** (tope de
pivotes, τ, presupuesto por nodo) y todos perdieron porque sacrificaban poda. El
corte por **éxito** no sacrifica nada: se para apenas el certificado ya alcanza
para decidir. Y apuntaba a lo que el §19 identificó como decisivo, el vacío.

Antes de implementarlo, la pregunta barata: **¿en qué pivote se vuelve detectable
el vacío?** (`DFBH_VACIO_PROBE=1`, que recorre la traza de `lambda` pivote a
pivote y busca el primero donde `0 ∉ γ·z`).

| instancia | LPs con vacío detectable | primer pivote / total |
|---|---|---|
| `house` | 50 de 4 000 (1.3 %) | 10.45 de 11.21 — **94 %** |
| `ex5_3_2` | 2 de 4 000 (0.1 %) | 42.50 de 42.50 — **100 %** |
| `launch` | 100 de 20 000 (0.5 %) | 33.15 de 34.49 — **97 %** |

**El certificado que prueba vacío es esencialmente el último.** Cortar ahí
ahorraría el 3–6 % de los pivotes de ese LP, y sólo en el 0.1–1.3 % de los LPs.

Tiene sentido en retrospectiva: el vacío se prueba cuando `γ` es lo bastante
bueno, y `γ` mejora monótonamente hacia el óptimo dual. Pedirle que pruebe vacío
antes del óptimo es pedirle que sea bueno antes de serlo.

La otra mitad del argumento —«se ahorra el resto del nodo»— tampoco queda en pie,
pero por buena razón: la cota que detecta el vacío es la **0.8ª de 18**, la **2ª
de 46**, la **8ª de 78**, o sea el **5–11 %** del recorrido del nodo. El nodo ya
muere casi al principio, y ese ahorro **ya se está obteniendo**, porque
`una_pasada` retorna de inmediato. No hay resto que ahorrar.

### Balance de la propiedad anytime en todo el documento

| aplicación | resultado |
|---|---|
| truncar pivotes por LP | §12, §20, §21: pierde |
| corte temprano por rendimiento (τ) | §20: pierde |
| tope de pivotes por nodo | §20: pierde |
| round robin de 1 pivote | §15–§17: converge igual, ×1.4 pivotes |
| parar tras la cota del objetivo | §20: celdas ×2.6 |
| **salir temprano por éxito** | **§35: el vacío se detecta en el 94–100 % del LP** |
| **garantía de validez al cortar** | **la única que vale, y es la que está en producción** |

La propiedad anytime resultó ser exactamente lo que el §5 decía y nada más: **una
red de seguridad**, que permite poner un tope de presupuesto sin analizar el caso
porque ninguna parada puede producir una cota inválida. Siete intentos de
convertirla en palanca de rendimiento, siete negativos.

Lo que sí resultó decisivo del mismo aparato matemático no fue la
incrementalidad sino la **independencia de la caja**: `γ = λᵀĀ` no depende de `z`,
y de ahí salió la prueba de vacío de Farkas del §19 —que llevó a DFB de 87 a 128
instancias— y el test barato del §30. Son la misma ecuación, usada por el otro
lado.

## 36. Presupuesto graduado: exactas al principio, parciales después

Última variante del intercalado, y la mejor motivada: el §35 midió que el nodo
muere en la cota 0.8ª de 18, 2ª de 46, 8ª de 78 —el **5–11 %** del recorrido—,
así que ese tramo tiene que ser **exacto**, porque de ahí sale la prueba de
vacío. Las demás sólo mejoran la caja para bisecar, y ahí podría alcanzar con ser
parcial.

Todo lo probado antes era **uniforme**: todo exacto (§25), todo truncado (§12,
§20) o sólo la cota del objetivo (§20). Esto es **graduado**
(`DFBH_EXACTAS=-p` como porcentaje de `2n`, `DFBH_RESTO=b`).

| config | celdas ×xn | cpu ×xn | ventaja cpu | en las ≥ 0.1 s |
|---|---|---|---|---|
| **sin graduar (§25)** | **0.925** | **0.889** | **+0.089** | **+0.045** |
| 50 % exactas, resto `b=8` | 0.939 | 0.952 | +0.078 | +0.008 |
| 25 % exactas, resto `b=8` | 0.946 | 0.969 | +0.079 | −0.007 |

Peor en celdas y en cpu, más lento en 24 instancias contra 18, y en las caras cae
a empate con producción.

Lo que sí es real: `alkylbis` pasa de 64 a **40 celdas** de forma reproducible.
La graduación hace algo en algunas instancias. Es el mismo patrón de siempre —
pocos ganadores, muchos perdedores, agregado en contra.

### Y una corrección al §31

El §31 advertía que `LinearizerXTaylor` usa esquina aleatoria y que por lo tanto
las comparaciones de una sola instancia no son confiables. **Es falso.**
Repitiendo `alkylbis` cinco veces: 64, 64, 64, 64, 64 celdas; y con la
configuración graduada: 40, 40, 40, 40, 40. Las corridas son **deterministas**.
La variación que motivó aquella advertencia venía de un cambio de código entre
las dos corridas, no de aleatoriedad.

Corresponde entonces lo contrario de lo que decía el §31: las comparaciones de
una sola instancia de este documento **sí** son reproducibles, y sirven para
diagnosticar. Lo que no se puede es generalizarlas al banco, que es otra cosa.

## 37. 3BCID: confirma la hipótesis de la adaptación, y destapa una interacción

El §31 dejó una hipótesis: el test barato pierde dentro de ACID no por su costo
—9 ms contra 0.68 s de degradación— sino porque **perturba una heurística
adaptativa**, ya que ACID decide cuántas variables rebanar según cuánto le rinde.
3BCID rebana con política fija: si la hipótesis es cierta, ahí la poda extra
debería aparecer como ganancia neta.

**Se confirma a medias.** Dentro de 3BCID el test **reduce** celdas (+0.032 en
ventaja acotada contra no usarlo) donde dentro de ACID las **aumentaba**. La
adaptación era efectivamente lo que lo absorbía. Pero el tiempo igual empeora
(−0.039; más lento en 67 instancias contra 22): la poda extra es real y aun así
no paga el test.

### La interacción que apareció

| configuración | celdas | cpu |
|---|---|---|
| `acidhc4` + producción | 141 788 | 200.8 s |
| **`acidhc4` + DFB** | **125 356** | 202.8 s |
| `3bcidhc4` + producción | **111 464** | 203.0 s |
| `3bcidhc4` + DFB | 115 512 | 232.5 s |
| `3bcidhc4` + DFB + test γ | 114 674 | 242.3 s |

(126 instancias que resuelven todas las configuraciones.)

Bajo **ACID**, DFB le gana a producción: 125 356 celdas contra 141 788. Bajo
**3BCID**, producción le gana a DFB: 111 464 contra 115 512, y 203 s contra
232 s.

**La ventaja de DFB es específica de `acidhc4`.** Eso no invalida el resultado
—`acidhc4` es la configuración de producción y es en la que se hizo toda la
comparación— pero sí acota su alcance, y merece decirse: no es «DFB es mejor que
PolytopeHull», es «DFB es mejor que PolytopeHull dentro de acidhc4».

De paso, un dato sobre Ibex y no sobre DFB: **3BCID poda bastante mejor que ACID**
(111 464 contra 141 788 celdas con el mismo contractor lineal) al mismo tiempo de
cpu. No es lo que se esperaría de una heurística adaptativa que existe
precisamente para abaratar el 3BCID.

## 38. Corrección importante: el §37 leyó totales como comportamiento

Antes de meter DFB dentro de 3BCID convenía entender por qué el §37 decía que
DFB es peor que PolytopeHull ahí. **No lo es**, y la respuesta obliga a corregir
también el titular de este documento.

Midiendo **por instancia** en vez de por totales, los dos filtrados se comportan
igual:

| filtrado | totales celdas | por instancia (acotada) | geom celdas | geom cpu |
|---|---|---|---|---|
| `acidhc4` | ×0.925 | **−0.037** | **×1.049** | ×0.974 |
| `3bcidhc4` | ×1.036 | **−0.039** | **×1.045** | ×0.986 |

Prácticamente idénticos. La ventaja de DFB **no** era específica de `acidhc4`:
los totales diferían por **una sola instancia**. `ex7_2_3` aporta **−18 392
celdas** bajo `acidhc4` y **+710** bajo 3BCID.

Esa instancia es legítima, no un artefacto: `acidhc4+xn` hace 27 920 celdas en
21.4 s y `acidhc4+dfbhull` 9 528 en 8.9 s, con optimización exitosa y el mismo
óptimo. Es una ganancia real de ×2.9 — pero **una sola de 128**.

### Lo que hay que corregir del resumen

| | totales | por instancia |
|---|---|---|
| celdas | **×0.925** | geom **×1.049**, acotada −0.037 |
| cpu | **×0.944** | geom **×0.974**, acotada **+0.086** |

Las dos lecturas son ciertas y dicen cosas distintas:

- **En cpu DFB gana en las dos**: totales ×0.944 y por instancia ×0.974 / +0.086.
  Ese resultado es sólido.
- **En celdas no.** El ×0.925 de los totales viene de `ex7_2_3`; por instancia DFB
  usa **4.9 % más** celdas. El enunciado «por primera vez DFB explora menos nodos
  que producción» (§25, §22) **vale para el total del banco, no para la instancia
  típica**, y así hay que decirlo.

### Y una advertencia de método

Toda esta sesión reportó celdas y cpu como **totales** del banco. Un total es una
suma sin normalizar: una instancia de 28 000 celdas pesa lo mismo que 300
instancias de 90. La métrica acotada `(A−B)/max` que propuso Ignacio no tiene ese
problema, y fue la que destapó esto. Donde las dos discrepan hay que mirar qué
instancia manda, porque suele ser una sola.

## 39. Gastar el ahorro en podar más: las dos formas fallan

El §38 dejó el cuadro claro: DFB poda **4.9 % menos** por nodo y cada nodo le
sale **15 % más barato**, con saldo de **10.7 % menos tiempo**. Lo natural
entonces no es ahorrar más —ocho intentos, ocho negativos— sino **gastar el
ahorro en podar más**. Dos formas.

### (1) Apretar el punto fijo

`CtcFixPoint(cxn_compo, relax_ratio)` reitera mientras la mejora supere el
umbral (0.2 por omisión). Bajarlo compra más poda con el ahorro. `DFBH_RATIO=r`:

| instancia | `xn` | r=0.2 | r=0.1 | r=0.05 |
|---|---|---|---|---|
| `alkylbis` | 50 | 64 | **40** | 120 |
| `ex5_3_2` | 58 | **94** | 134 | 130 |
| `house` | 240 | 182 | 154 | **144** |
| `ex7_2_3` | 27 920 | **9 528** | 16 496 | 27 932 |

**No monótono y en general peor.** `ex7_2_3` pasa de 9 528 a 27 932 celdas al
apretar el umbral, que es lo contrario de lo que debería hacer más contracción
por nodo. Efectos de trayectoria otra vez.

### (2) Reusar la linealización dentro del punto fijo

`CtcFixPoint` llama al contractor lineal varias veces por nodo y entre esas
llamadas la caja **sólo se encoge**, así que la relajación anterior sigue siendo
válida. Reusarla evita re-linealizar, evita recargar los `2n` simplex, y deja los
LPs con el **mismo objetivo sobre una caja apenas más chica** — el caso donde el
§25 mostró que el arranque tibio rinde más. Es la incrementalidad en el régimen
ideal (`DFBH_RELIN=1`).

| instancia | base | reusando |
|---|---|---|
| `house` | 182 | **16 714** |
| `ex5_3_2` | 94 | 506 |
| `alkylbis`, `ex7_2_3`, `launch` | — | fuera de tiempo |

**Catastrófico**, y la razón importa más que el resultado.

### La razón, y lo que unifica

**La relajación no es un poliedro fijo: es un desarrollo de Taylor cuya calidad
depende fuertemente del ancho de la caja.** Una relajación calculada sobre una
caja más ancha es mucho más floja, y esa pérdida arrasa con cualquier ahorro.

Eso unifica varios negativos que veníamos tratando como independientes:

- **§17**: re-evaluar los `γ` guardados dio **cero** contracción extra.
- **§26–§30**: el test dentro de ACID acierta poco (0.2–17 %) y sólo en las
  instancias grandes.
- **§39**: reusar la linealización es catastrófico.

Los tres reusan información derivada de una relajación calculada en una caja más
ancha, y los tres son débiles **por la misma razón**. No son tres resultados,
son uno.

Y da una regla de diseño: en este contractor, **re-linealizar es barato
comparado con lo que cuesta no hacerlo**. Toda idea que empiece con «reusemos la
relajación anterior» arranca con una desventaja que hay que cuantificar antes de
seguir.

## 40. Re-linealizar pero heredar la base: el defecto estaba en otro lado

El §39 apagó la re-linealización **y** reusó las bases a la vez, así que su
desastre no distingue cuál de las dos fue. La variante que falta es la buena:
**re-linealizar** (poliedro fuerte) pero **arrancar de la base anterior** (warm
start). Existía ya, como `DFB_SX_NODOS=1`.

### Estaba roto

Medido, daba **resultados idénticos** a no usarlo. La causa: `load` reseteaba
`basic`, `rowof` y `atupper` con `assign` **antes** del bloque que decía
conservarlos, así que `cols_buf[r-1] = basic[r]` leía `-1`, `set_basis` fallaba y
la herencia no ocurría nunca. La opción existía desde hacía tiempo y no hacía
nada.

### Arreglado, es peor

Copiando la base antes de los `assign` y restaurándola después:

| instancia | base | heredando |
|---|---|---|
| `alkylbis` | 64 | **144** |
| `house` | 182 | 226 |
| `ex7_2_3` | 9 528 | 10 378 |
| `launch` | 552 | 718 |
| `ex5_3_2` | 94 / 1.42 s | 94 / 1.77 s |

Dos razones, y la primera es un defecto de implementación:

1. **El estado queda inconsistente.** `load` restaura la base pero deja el
   tableau vacío, porque `tidx`/`tval` se resetean arriba. La base dice que
   ciertas columnas son básicas y el tableau no lo refleja. Hacerlo bien exige
   **reconstruir el tableau** para la base heredada —`m` FTRAN—, que es trabajo
   real y no gratis.
2. **Ya hay un arranque mejor.** La siembra del §11 resuelve la primera cota en
   frío sobre la relajación **nueva** y copia esa base a los `2n` ejemplares: es
   óptima para *esta* relajación. La heredada es óptima para la anterior y para
   otro objetivo. Se estaría cambiando algo bueno por algo peor.

Queda apagado. La corrección del bug se conserva para que `DFB_SX_NODOS=1` haga
lo que dice, pero con la advertencia de que el tableau no se reconstruye.

La validación sigue intacta: 528/528 cotas coincidentes con SoPlex.

## 41. Corrección del §37, y las cuatro configuraciones por instancia

El §37 sacó dos conclusiones de **totales**, y las dos caen al mirar por
instancia. Es el mismo error que el §38 ya había corregido en otro lado.

| config | celdas tot | cpu tot | celdas geom | cpu geom | cpu acotada |
|---|---|---|---|---|---|
| `acidhc4` + PH | 141 788 | 200.8 s | 1.000 | 1.000 | +0.000 |
| **`acidhc4` + DFB** | 125 356 | 202.8 s | ×1.049 | **×0.901** | **+0.089** |
| `3bcidhc4` + PH | **111 464** | 203.0 s | ×0.991 | ×1.027 | −0.020 |
| `3bcidhc4` + DFB | 115 512 | 232.5 s | ×1.036 | ×0.945 | +0.066 |

(126 instancias que resuelven las cuatro.)

**«3BCID poda bastante mejor que ACID» es falso.** Por totales lo parece
(111 464 contra 141 788 celdas) pero por instancia poda **1 %** mejor
(×0.991) y es más **lenta** (×1.027). El 111 464 es la suma dominada por pocas
instancias, igual que el ×0.925 del §38.

**«DFB es peor bajo 3BCID» también es falso.** Dentro de cada filtrado, DFB
contra PolytopeHull:

| filtrado | celdas | cpu | tiempo por nodo |
|---|---|---|---|
| `acidhc4` | ×1.049 | ×0.901 | **×0.858** |
| `3bcidhc4` | ×1.045 | ×0.920 | **×0.880** |

Casi idénticos: DFB poda ~4.7 % menos y cada nodo le sale ~13 % más barato,
**bajo los dos filtrados**. No hay nada anómalo que entender en 3BCID, y no es
un hogar mejor para DFB.

**La mejor configuración es `acidhc4 + DFB`**, que es la que ya estaba.

### Un patrón de esta sesión que conviene dejar escrito

Tres veces —§37, §38, §41— una conclusión sacada de **totales del banco** se cayó
al recalcularla **por instancia**. Un total es una suma sin normalizar donde una
instancia de 28 000 celdas pesa como 300 de 90, así que mide «el trabajo total
del banco», no «cómo se comporta el método». Para comparar métodos hay que usar
la media geométrica de razones o la ventaja acotada `(A−B)/max`. Los totales sólo
sirven si lo que se quiere saber es cuánto cuesta correr el banco entero.

## 42. Revisión del documento con la métrica por instancia

El §41 dejó claro que los totales del banco engañan. Se revisaron todas las
afirmaciones importantes recalculándolas por instancia.

### Lo que se sostiene

| sección | afirmación | total | geom | ¿coincide? |
|---|---|---|---|---|
| §9/§11 | siembra vs base compartida (pivotes) | ×0.738 | ×0.754 | sí |
| §9/§11 | base fría vs compartida (pivotes) | ×1.532 | ×1.500 | sí |
| §10/§13 | Achterberg vs orden natural (LPs) | ×0.539 | ×0.564 | sí |
| §15/§16 | round robin vs secuencial (tiempo) | ×3.244 | ×2.107 | sí |
| §4, §14 | **tiempo** de DFB vs SoPlex | ×0.22 / ×0.24 | ×0.29 / ×0.34 | sí |

Y los §4 y §14 **ya reportaban las dos métricas** y explicaban la diferencia, así
que ahí no había nada que corregir.

### Donde totales y geométrica discrepan

Sólo en **pivotes**, y en las dos secciones ya estaba dicho:

| sección | pivotes total | pivotes geom |
|---|---|---|
| §4 | ×0.900 | **×1.265** |
| §14 | ×0.938 | **×1.012** |

DFB usa menos pivotes **en total** y algo más **en la instancia típica**: gasta
menos en las grandes y un poco más en las chicas. La conclusión que importa —el
tiempo— no depende de eso y es robusta en las tres métricas.

### Ninguna variante descartada se resucita

Se recalcularon las doce variantes rechazadas contra la configuración actual,
por instancia. **Todas pierden también así** (ventaja acotada entre −0.015 y
−0.057). Los descartes estaban bien hechos.

Pero aparece un patrón que no se había visto:

| variante descartada | celdas (geom) | cpu (geom) |
|---|---|---|
| objetivo primero (§21) | **×0.946** | ×1.042 |
| paso largo `m ≥ 30` (§23) | **×0.962** | ×1.052 |
| test en ACID + `b` (§30) | **×0.977** | ×1.034 |
| `3bcid` + DFB (§37) | **×0.984** | ×1.046 |

**Cuatro de las descartadas podan mejor por instancia** —2 a 5 % menos celdas— y
pierden igual, porque cuestan 3 a 5 % más de tiempo. No fallaban por no podar,
fallaban por el precio. En estas direcciones el sistema está limitado por
**costo**, no por poda, y eso acota dónde puede estar la próxima mejora: en algo
que baje el costo sin tocar la poda, no al revés.

## 43. Perfilado del simplex, y la primera optimización pura de costo

El §42 dejó una sola dirección abierta: **bajar el costo sin tocar la poda**. El
simplex es el 57–85 % del tiempo del contractor (§21) pero nunca se había mirado
adentro. `DFB_SX_PERF=1` mide cada fase en tiempo **exclusivo** —las fases están
anidadas, así que cada cronómetro descuenta lo de sus hijos; sin eso los
porcentajes pasan de 100 y no se puede leer nada.

| instancia | `pivot` | `compute_basics` | ratio test | primal | refresh | y_fact | resto | `load` |
|---|---|---|---|---|---|---|---|---|
| `alkylbis` | 38 % | **20 %** | 7 % | 11 % | 6 % | 6 % | 16 % | **32 %** |
| `ex5_3_2` | 53 % | **21 %** | 6 % | 8 % | 3 % | 3 % | 10 % | 24 % |
| `house` | 33 % | **16 %** | 9 % | 12 % | 6 % | 8 % | 20 % | 19 % |
| `launch` | 51 % | **26 %** | 5 % | 8 % | 3 % | 3 % | 8 % | **45 %** |
| `ex7_2_3` | 28 % | **16 %** | 9 % | 14 % | 7 % | 7 % | 23 % | 16 % |

Dos blancos, los dos de pura implementación:

1. **`load` se lleva del 16 % al 45 %**, y en `launch` casi la mitad. Se llama
   `2n` veces por nodo **con la misma matriz**, recomputando cada vez las
   escalas, la copia dispersa por filas y la estructura de `DFBBasis`.
2. **`compute_basics` se lleva 16–26 %** recalculando `x_B` desde cero en cada
   iteración del paso dual, `O(m·nnz)`, cuando un simplex lo actualiza
   incrementalmente tras el pivote en `O(m)`.

Una corrección: el §21 midió la carga como «8–19 % del contractor» y la descartó
por chica. Es el mismo número contra otro denominador — contra el tiempo del
**simplex** es 16–45 %. Se descartó un blanco grande por compararlo contra el
total equivocado.

### La carga compartida

`DFBSimplex::cargar_desde` y `DFBBasis::copiar_matriz_de` copian lo ya construido
en vez de recomputarlo: se carga `banco[0]` y los otros `2n−1` copian. No se
copia la factorización `lu` —`SLUFactor` guarda punteros a las columnas del
origen— sino sólo la estructura; la factorización la establece `set_basis`.

| config | celdas | geom | cpu | geom | acotada |
|---|---|---|---|---|---|
| producción `xn` | 213 126 | 1.000 | 268.3 s | 1.000 | — |
| DFB antes | 197 064 | ×1.049 | 253.3 s | ×0.898 | +0.091 |
| **DFB con carga compartida** | **197 064** | ×1.049 | **229.2 s** | **×0.830** | **+0.151** |

**Celdas idénticas en las 128 instancias**, como corresponde a un cambio que no
toca la poda. Más rápido en **73** instancias y más lento en **1**. Cero óptimos
incompatibles, y el test unitario sigue en 528/528 contra SoPlex.

La cuota de `load` baja del 16–45 % al 6–9 %.

Es la mejora más segura de la sesión: no cambia ninguna decisión del algoritmo,
sólo deja de repetir trabajo idéntico `2n` veces.

## 44. El segundo blanco: el bug encontrado, y por qué igual no paga

Con `load` bajado al 6–9 % (§43), el perfil queda `pivot` 38–53 % y
**`compute_basics` 20–26 %**. Ese 20–26 % recalcula `x_B` desde cero en cada
iteración del paso dual —`O(m·nnz)`— cuando puede actualizarse incrementalmente:
si `α` es la columna entrante transformada **antes** del pivote —que `pivot()` ya
recorre— y `Δ` el paso, entonces `x_B(i) ← x_B(i) − α_i·Δ`, que es `O(m)`.

### El bug, y es real

La primera versión daba **8 cotas distintas de 528** con los **mismos 810
pivotes**. Los casos concretos mostraron la firma: para la misma variable, las
respuestas de `up=0` y `up=1` salían **intercambiadas exactamente**. Es la misma
firma de un defecto que `ibex_CtcDualFeasibleBounding.cpp` documentaba desde
antes —*«el camino tibio devuelve para maximizar la respuesta de minimizar»*— y
que nunca se había explicado.

La pista decisiva fue cruzar los interruptores:

| | discrepancias |
|---|---|
| por omisión (sin reubicar, sin incremental) | 0 |
| `REUBICAR=1` (comportamiento previo al §25) | 0 |
| `XBINC=1` (incremental, sin reubicar) | **8** |
| `REUBICAR=1` + `XBINC=1` | 0 |

El incremental **sólo** fallaba combinado con «no reubicar», que es el default
del §25. Y esos dos caminos se distinguen en cuál motor re-optimiza: sin reubicar
el punto queda primal-factible y trabaja el **simplex primal**; reubicando queda
dual-factible y trabaja el dual.

El defecto estaba entonces en el camino primal, y ahí hay un caso que no pivotea:

```c
/* bound flip: la entrante puede recorrer a lo sumo su propio rango */
if (rango < lim) {
    atupper[q] = atupper[q] ? 0 : 1;
    return true;          // <-- sin llamar a pivot()
}
```

Un *bound flip* mueve `x_B` **sin pivotear**, así que `pivot()` —que es quien
marca el estado sucio— nunca se llama y `x_B` queda obsoleto. Con el camino que
recalcula eso no se notaba; con el incremental, la resolución devolvía la
respuesta de la anterior.

**Arreglado** (`xb_sucio = true` en esa rama): 528/528 en las tres
configuraciones. Es un arreglo que se conserva aunque el incremental quede
apagado, porque hace correcto un invariante que estaba mal.

### Y aun así no paga

| config | celdas | geom | cpu | geom | acotada |
|---|---|---|---|---|---|
| producción `xn` | 213 126 | 1.000 | 268.3 s | 1.000 | — |
| **DFB (recalculando)** | **197 064** | ×1.049 | **229.2 s** | **×0.830** | **+0.151** |
| DFB (incremental) | 197 402 | ×1.069 | 229.5 s | ×0.845 | +0.141 |

Contra la versión que recalcula: cpu geom **×1.018**, más lento en 14 instancias
y más rápido en 6, con **celdas distintas en 25 de 128**.

La causa es **deriva numérica**: la acumulación incremental redondea distinto que
el recálculo —matemáticamente equivalente, no bit a bit— y eso mueve las
decisiones del árbol. Es justo por lo que los simplex de libro recalculan o
refactorizan periódicamente. Y el ahorro resultó chico: `compute_basics` sólo
bajó del 20–26 % al 16–22 %, porque los bound flips y los pivotes primales la
ensucian seguido.

Queda apagada (`DFB_SX_XBINC=1`), ahora **con la causa entendida** y con el bug
del bound flip arreglado en el camino por omisión.

## 45. Resumen

*Este es el resumen válido. El §22 quedó superado.*

### Qué se hizo

Reemplazar `CtcPolytopeHull` —que resuelve un LP con SoPlex y certifica el
óptimo con Neumaier–Shcherbina— por **DFB**, que acota por dualidad: pivotea en
punto flotante, y de los multiplicadores `λ` construye `γ = λᵀĀ` para certificar
la cota **en aritmética de intervalos**. Un `λ` malo da una cota débil, nunca
falsa.

El contractor es [`CtcDFBHull`](ibex_CtcDFBHull.h), en el hueco **exacto** de
`CtcPolytopeHull` y con la misma composición. Se activa con `--lr=dfbhull`.

### Resultado

139 instancias de `benchs/optim` (easy, medium, hard), `--filtering=acidhc4`,
40 s:

**Sobre 4 semillas** (1, 7, 17, 42), mismo binario, límite 60 s, 484 corridas:

| semilla | instancias | celdas geom | **cpu geom** | cpu acotada |
|---|---|---|---|---|
| 1 | 129 | ×1.099 | ×0.887 | +0.106 |
| 7 | 128 | ×1.012 | ×0.823 | +0.149 |
| 17 | 129 | ×1.054 | ×0.835 | +0.149 |
| 42 | 130 | ×1.050 | ×0.839 | +0.143 |
| **TODAS** | **484** | **×1.053** | **×0.846** | **+0.137** |

- **DFB más rápido en 323 corridas de 484**, más lento en 86, empate en 75.
- Las mismas instancias resueltas por ambas, salvo 1 a 2 por semilla en las que
  DFB resuelve y producción no (nunca al revés).
- **0 óptimos incompatibles** salvo uno, en una instancia donde ninguna de las
  dos prueba optimalidad (§33).
- Validado contra SoPlex: **528/528 cotas**, invariante `γ = λᵀĀ` 200/200,
  dualidad fuerte 55/55.

**Cómo leerlo.** DFB es más rápido porque **cada nodo le sale ~20 % más barato**,
no porque pode más: por instancia usa 5.3 % **más** celdas (geom ×1.053) y 15 %
menos tiempo (geom ×0.846). El cociente da tiempo por nodo ×0.803.

**El agregado es robusto a la semilla** —la cpu geométrica cae en la banda
0.823–0.887— pero **una instancia sola no lo es**: `ship-1` oscila ×5.2 entre
semillas. Ver el §46, que corrige una conclusión de este resumen sacada de una
sola corrida.

### Las cuatro decisiones que lo hicieron posible

| decisión | efecto |
|---|---|
| **Prueba de vacío por Farkas** (`0 ∉ γ·z` ⟹ caja vacía) | **87 → 128 instancias resueltas** (§19) |
| **No reubicar las no básicas** al cambiar de objetivo | pivotes/LP de 29.9 a 7.3; la degeneración cae de 78–88 % a 14–30 % (§25) |
| **Estados propios sembrados** desde la primera cota | pivotes ×0.75 contra base compartida, ×0.48 contra fría (§11) |
| **Carga compartida** entre los `2n` ejemplares | cpu ×0.925, celdas idénticas (§43) |

Más el **salteo de Achterberg** (§10) y el **orden de la más lejana** (§13).

### Lo que no funcionó, y por qué

Nueve intentos de convertir la propiedad *anytime* en rendimiento, nueve
negativos: truncar por LP, corte por rendimiento, tope por nodo, round robin de
un pivote, parar tras la cota del objetivo, test barato dentro de ACID, salida
temprana por éxito, presupuesto graduado, herencia de base entre nodos.

Tres causas, y cada una explica varios:

1. **El contractor lineal es ~20 % del costo del nodo.** Ahorrar ahí no puede dar
   más de 20 %, y cada ahorro costó más poda de la que valía.
2. **La relajación es un Taylor que se degrada con el ancho de la caja.** Todo lo
   que reusa información de una caja más ancha —`γ` guardados, la linealización
   anterior— es débil por eso (§39).
3. **El sistema está en un óptimo local estrecho.** Casi cualquier empujón a ACID,
   al punto fijo o al orden lo empeora, aunque el empujón sea «contraer más»
   (§31: un test que cuesta 9 ms encarece el run 0.68 s).

La propiedad anytime quedó como lo que es: **una red de seguridad** —ningún tope
de presupuesto puede producir una cota inválida— y no una palanca de
rendimiento. Lo decisivo del mismo aparato no fue la incrementalidad sino la
**independencia de la caja** de `γ`, de donde salió la prueba de vacío.

### Lecciones de método

- **Los totales del banco engañan.** Tres veces (§37, §38, §41) una conclusión
  sacada de sumas se cayó al recalcularla por instancia; una de ellas era el
  titular. Un total pesa una instancia de 28 000 celdas como 300 de 90. Para
  comparar métodos: media geométrica de razones, o la ventaja acotada
  `(A−B)/max`.
- **Medir en el régimen que importa.** La caja raíz no predice los nodos
  profundos: el §21 diagnosticó «DFB usa 4.5× más pivotes» comparando promedios
  sobre poblaciones distintas de LPs, y de ahí salieron cuatro remedios contra un
  mal inexistente (§24).
- **Perfilar antes de optimizar.** La carga se descartó en el §21 por «8–19 % del
  contractor»; contra el tiempo del **simplex** era 16–45 %, y arreglarla dio la
  mejora más grande y segura de la sesión (§43).
- **Volcar un caso concreto cuesta menos que probar remedios.** Cuatro
  experimentos fallidos se habrían evitado resolviendo un LP duro con los dos
  solvers, que es lo que finalmente dio la respuesta (§24).

### Lo que queda abierto

1. **`pivot` es el 38–53 % del simplex.** Eliminarlo pediría trabajar sólo desde
   la factorización, sin tableau. Es rediseño, no optimización.
2. **La ventaja en celdas es negativa por instancia** (×1.049). Cerrar ese 4.9 %
   sin gastar tiempo daría la mejora más grande disponible.
3. **`compute_basics` (20–26 %) no es recuperable por la vía incremental**, y el
   §44 explica por qué: deriva numérica. Cualquier otro intento tendría que
   evitar la acumulación.
4. **`ex6_1_3`** es la única instancia dura donde DFB pierde de forma
   consistente entre semillas (cpu ×1.525). Es un caso concreto y reproducible
   para entender cuándo el acotamiento dual rinde menos que resolver el LP.

## 46. La semilla aleatoria invalida las conclusiones de una sola corrida

El §45 concluyó, con el barrido de 600 s, que **DFB pierde en las instancias más
duras**: sobre las 5 que pasan de 40 s daba cpu ×1.217 y celdas ×1.212, con
`ship-1` ×2.8 peor. **Es falso**, y la causa es que todas las mediciones de este
documento usan **una sola semilla**.

`LinearizerXTaylor` usa esquinas aleatorias (`RANDOM_OPP`), así que la semilla
cambia el poliedro de cada nodo y con él toda la búsqueda. Se fija con
`--seed=<n>`; la de omisión es **42**.

### `ship-1` con ocho semillas

| semilla | producción | DFB | DFB/prod |
|---|---|---|---|
| 1 | 120.0 s | **47.7 s** | ×0.40 |
| 7 | 98.2 s | **65.7 s** | ×0.67 |
| 17 | 93.2 s | **54.7 s** | ×0.59 |
| **42** | 87.4 s | **248.1 s** | **×2.84** |
| 99 | 69.0 s | 73.6 s | ×1.07 |
| 123 | 93.4 s | 102.3 s | ×1.10 |
| 777 | 50.8 s | 61.6 s | ×1.21 |
| 2024 | 105.3 s | 90.1 s | ×0.86 |

DFB va de **47.7 s a 248.1 s** según la semilla —un factor **×5.2**— y la 42, la
de omisión, resulta ser **su peor caso de las ocho**. La media geométrica es
**×0.925**: DFB es un 7 % *más rápido* en `ship-1`, no ×2.8 más lento.

### Las cinco duras, con ocho semillas

| instancia | cpu con s=42 | **cpu geom (8 sem)** | celdas geom | variación de DFB |
|---|---|---|---|---|
| `ship-1` | ×2.84 | **×0.925** | ×1.039 | ×5.2 |
| `ex8_4_4bis` | ×0.61 | ×0.618 | ×1.001 | ×1.0 |
| `ex6_1_3` | ×1.32 | ×1.525 | ×1.043 | ×1.2 |
| `ex14_2_7` | ×1.05 | ×1.048 | ×1.170 | ×1.6 |
| `bearing` | ×1.10 | **×0.661** | **×0.487** | ×3.7 |

**En conjunto: cpu ×0.904 y celdas ×0.909.** O sea que DFB **gana** en las
instancias duras, exactamente al revés de lo que el §45 concluyó con una sola
semilla. Sólo `ex6_1_3` pierde de forma consistente.

### Lo que esto obliga a revisar

**Una corrida no es una medición.** Y lo engañoso es que con la misma semilla las
corridas **sí** son deterministas y reproducibles —el §36 lo verificó repitiendo
`alkylbis` cinco veces con resultados idénticos— así que uno confirma la
reproducibilidad y concluye que el número es sólido. No lo es: es sólido para esa
semilla.

Esto afecta sobre todo a las conclusiones apoyadas en pocas instancias, y el §38
mostró que los totales del banco suelen estar dominados por una o dos. Las
comparaciones sobre 128 instancias promedian parte del ruido de semilla, pero no
todo, y no cuando una instancia manda.

Es la cuarta corrección de este documento por el mismo tipo de error: **generalizar
desde una muestra que no cubre la variación relevante**. Antes fueron la caja raíz
(§24), las instancias chicas (§27), los totales del banco (§38, §41). Ahora la
semilla.
