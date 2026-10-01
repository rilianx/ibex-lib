# El DFB perezoso: repartir el presupuesto de pivotes

Idea a probar (§8.0 B del [plan](PLAN_MEJORA_DFB.md)): gracias a la
certificación, **en DFB cada pivote ya entrega una cota válida**, así que se
puede gastar un presupuesto chico por cota y seguir, cambiando calidad por costo
de forma continua. `CtcPolytopeHull` no puede: su cota aparece cuando el LP
llega a `OptimalProved` y una corrida interrumpida no entrega nada por la API de
Ibex.

## 0. Lo primero que apareció: ya estaba a medias implementado

Los `2n` contractores se construían con `max_iters = 1`
([ibex_CtcDFBPropag.cpp:27](ibex_CtcDFBPropag.cpp#L27)) y, al reencolarse, la
tupla lleva `pq_order++`, así que el contractor vuelve al **final** de su grupo
en la cola. O sea que el comportamiento vigente **ya era** un round-robin de un
pivote por cota con HC4 intercalado. Lo que no existía era el mando para
variarlo, ni un tope total.

Se agregaron dos mandos independientes, y conviene no confundirlos porque miden
cosas distintas:

| mando | qué controla |
|---|---|
| `DFB_MAXITERS=k` | pivotes **por visita** de un contractor (`-1` = hasta el punto fijo de esa cota). Default 1. |
| `DFB_PIVOT_BUDGET=k`, `DFB_PIVOT_BUDGET_N=c` | tope **total** de pivotes DFB por caja (absoluto o `c·n`). Al agotarse se dejan de atender contractores DFB y la propagación sigue solo con HC4. Default sin tope. |

Lo ya contraído se conserva al agotar el presupuesto: cada pivote dio una cota
certificada, no hace falta llegar al punto fijo para quedarse con lo obtenido.

También se agregó `BENCH_WITH_HC4=1` a `bench_dfb`, porque el banco de una caja
construía el contractor con `only_dfb = true` y por lo tanto **medía DFB en
aislamiento**: el reparto solo tiene sentido si hay algo con que intercalar.
Verificado antes de medir: con los valores por omisión los dos mandos nuevos no
cambian nada (0 diferencias estructurales en las 153 instancias).

## 1. La curva, en el banco de una caja con HC4

137 instancias donde todas las variantes terminan `ok`, la linealización sirve
(`A_rows > 1`) y nadie vacía la caja. `t_dfb` y `contrac%` incluyen el trabajo
de HC4, porque la propagación es conjunta.

| variante | contrac % | `t_dfb` total | pivotes | ms/pivote |
|---|---|---|---|---|
| tope total 0 (solo HC4) | 24.85 | 2.717 s | 0 | — |
| tope total 0.5 n | 25.66 | 2.961 s | 1 820 | 1.63 |
| tope total 1 n | 26.42 | 3.198 s | 3 596 | 0.89 |
| tope total 2 n | 27.69 | 3.702 s | 6 692 | 0.55 |
| tope total 4 n | 28.27 | 4.848 s | 11 342 | 0.43 |
| 1 pivote/visita (default) | 36.92 | 5.789 s | 21 315 | 0.27 |
| 2 pivotes/visita | 37.69 | 5.713 s | 22 494 | 0.25 |
| 3 pivotes/visita | 37.90 | 5.744 s | 22 972 | 0.25 |
| 5 pivotes/visita | 38.10 | 5.716 s | 23 580 | 0.24 |
| punto fijo por visita | **38.10** | **5.742 s** | 24 653 | 0.23 |

### 1.1 La contracción es convexa en el presupuesto, y eso hunde la idea

Puntos de contracción ganados por segundo **adicional** sobre HC4 solo:

| presupuesto | Δcontrac | Δtiempo | puntos/s |
|---|---|---|---|
| 0.5 n | +0.81 | +0.244 s | 3.32 |
| 1 n | +1.57 | +0.481 s | 3.26 |
| 2 n | +2.84 | +0.985 s | 2.88 |
| 4 n | +3.42 | +2.131 s | 1.60 |
| sin tope, 1 piv/visita | +12.07 | +3.072 s | 3.93 |
| sin tope, punto fijo | **+13.25** | **+3.025 s** | **4.38** |

El presupuesto parcial no es un intercambio favorable: **es el peor punto de la
curva**. Con 11 342 pivotes (4 n) se ganan 3.42 puntos, y al duplicarlos a
24 653 se ganan 13.25. La contracción **crece más rápido que el gasto**, así que
recortar el presupuesto no compra eficiencia, la destruye.

La razón es el reparto: con un pivote por visita y round-robin sobre las `2n`
cotas, un tope de `c·n` pivotes le da a cada cota **c/2 pivotes**, o sea casi
nada. Y el primer pivote desde la base inicial es el que menos aporta —la cota
buena aparece después de varios—. El presupuesto se gasta entero en la parte
plana de todas las cotas a la vez.

### 1.2 Y el reparto fino tampoco: `max_iters = 1` es peor que el punto fijo

Este es el hallazgo que no se esperaba. Comparando el default vigente contra
llevar cada cota a su punto fijo en una sola visita, **con el mismo presupuesto
(sin tope)**:

| | contrac % | tiempo | pivotes | ms/pivote |
|---|---|---|---|---|
| 1 pivote/visita | 36.92 | 5.789 s | 21 315 | 0.272 |
| punto fijo/visita | **38.10** | **5.742 s** | 24 653 | **0.233** |

El reparto fino usa **menos pivotes** (21 315 contra 24 653) y sin embargo tarda
**lo mismo o un poco más**, y contrae **1.2 puntos menos**. O sea: está dominado
en las dos dimensiones. El costo por visita —reentrada, restaurar `b`, la
mantención de la cola, `real_impact`— se come lo que ahorra en pivotes.

La curva de `max_iters` es monótona y satura rápido: 1 → 36.92, 2 → 37.69,
3 → 37.90, 5 → 38.10, punto fijo → 38.10. **A partir de 5 no cambia nada**, así
que el punto fijo por cota no cuesta más que un tope de 5.

## 2. En el optimizador: el presupuesto parcial confirma su fracaso, y el reparto fino resulta ser una mala elección de default

191 instancias, límite de 30 s, `--filtering=dfb --lr=no`.

| variante | resuelve | timeout | otros |
|---|---|---|---|
| `it1` (1 pivote/visita, default) | 116 | 42 | 32 |
| `it2` | **117** | 40 | 33 |
| `it5` | 116 | 41 | 33 |
| `itinf` (punto fijo/visita) | 74 | 16 | **100** |
| `b1n` (tope total 1 n) | 77 | 83 | 30 |
| `b2n` (tope total 2 n) | 80 | 77 | 33 |
| `b4n` (tope total 4 n) | 88 | 70 | 32 |

### 2.1 El tope total: peor en el árbol, como en la caja

Sobre las 53 instancias que las siete resuelven, relativo a `it1`:

| variante | celdas (geom) | cpu (geom) |
|---|---|---|
| `b1n` | ×1.951 | ×1.574 |
| `b2n` | ×1.423 | ×1.276 |
| `b4n` | ×0.998 | ×0.968 |

Y resuelven 77, 80 y 88 contra 116. El tope total **recorta la poda mucho más de
lo que ahorra**, y el árbol paga la diferencia con creces. Cuanto más grande el
tope, mejor: la curva no tiene un óptimo interior, apunta a «sin tope». Queda
descartado.

### 2.2 `itinf` no es aceptable: 72 procesos hubo que matarlos

`DFB_MAXITERS=-1` produce **72 `proc_timeout`** —el proceso pasó de `TMO+30` s y
hubo que matarlo— contra 1 en todas las demás variantes. Es la trampa de la
no terminación otra vez (§7 del plan): sin tope por visita, `contract_float`
gira hasta que alguna condición de corte se cumpla, y si el pivoteo en
flotantes cicla por degeneración **no vuelve**. El límite de tiempo de
`ibexopt` se comprueba entre nodos, así que no puede actuar.

O sea que el 38.10 % del banco de una caja que `itinf` mostraba en §1.2 **no es
un punto usable**: se obtuvo en instancias donde una sola contracción alcanzó a
terminar. Un tope por visita no es solo una política de reparto, es lo que
garantiza que el contractor retorne.

### 2.3 Lo que sí quedó: subir el tope de 1 a 5

Con las tres variantes seguras el conjunto común es mucho mayor —112
instancias— y la tendencia es monótona:

| variante | celdas totales | ×  | cpu total | × | celdas (geom) | cpu (geom) |
|---|---|---|---|---|---|---|
| `it1` (default) | 305 786 | 1.00 | 262.89 s | 1.00 | 1.000 | 1.000 |
| `it2` | 283 204 | 0.93 | 227.25 s | 0.86 | 0.964 | 0.950 |
| `it5` | 281 736 | **0.92** | **226.29 s** | **0.86** | **0.929** | **0.934** |

Mismo número de instancias resueltas (116 / 117 / 116) y **0 óptimos
incompatibles**. Un **7 % menos de celdas y de tiempo** en media geométrica, con
tendencia monótona en las dos métricas, es una mejora chica pero creíble: no es
un caso aislado sino un desplazamiento de la distribución.

Encaja con el §1.2: el reparto de un pivote por visita paga costo de reentrada
sin comprar nada, y 5 pivotes por visita ya captura toda la contracción que
captura el punto fijo (38.10 % en las dos) **sin renunciar a la garantía de
terminación**.

## 3. Validación y cambio del default

Protocolo del §6 del [plan](PLAN_MEJORA_DFB.md) con `DFB_MAXITERS=5`:

| paso | resultado |
|---|---|
| `./tests_dfb_all`, 7 casos | pasan; los dos `NO SOLUTION TEST` se detectan vacíos |
| campaña de testigos, las 32 instancias de referencia | **32/32 conservados** |
| vaciados corroborables | 1 vaciado (`Prolog`), ya corroborado `infeasible` por `ibexsolve` |
| óptimos incompatibles en `ibexopt` | **0** en las 6 variantes de la curva |
| `valgrind` sobre `ibexopt --filtering=dfb` | 0 errores en 0 contextos |

**El default pasó de 1 a 5** en
[ibex_CtcDFBPropag.cpp](ibex_CtcDFBPropag.cpp). Banco de una caja con el default
nuevo (`results/baseline_153_it5.csv`): 148 ok, **grupo A 73** (era 72),
B 25, C 49, contracción media **30.29 %** (era 29.94 %), `t_dfb` 6.733 s (era
6.246 s).

Una pérdida menor y hay que anotarla: con tope 5 el banco **deja de probar
vacía a `Prolog-icse`** (con tope 1 la probaba). No es un problema de solidez
—ese vaciado era correcto, `Prolog-icse` es `infeasible`— sino menos poda: la
detección de vacío por la cota depende de la secuencia de pivotes, y al cambiar
el reparto se pierde. En `ibexopt` no se pierden instancias resueltas (116 con
tope 1, 116 con tope 5, 117 con tope 2).

### 3.1 Un efecto de borde que conviene arreglar: `step_cap` cuenta visitas

En el banco aislado, subir el tope de 1 a 5 hace **85 % más pivotes** en total,
lo que a primera vista contradice que ambas variantes corran hasta el punto fijo.
La causa está localizada y medida:

| | pivotes con tope 1 | pivotes con tope 5 |
|---|---|---|
| 143 instancias que **no** tocan el `step_cap` | 27 267 | 27 742 (+1.7 %) |
| **5 instancias que sí lo tocan** | 7 117 | 35 785 (×5.0) |

`step_cap` (§7 del plan, `max_propag_steps`, por omisión `50·(ctcs+1)`) limita
**visitas de contractores**, no pivotes. Con tope 1 eso equivale a un tope de
pivotes; con tope 5, al mismo número de visitas le corresponden 5× más pivotes.
Las cinco instancias afectadas —`Fourbar-icse`, `Redeco8/10/11`, `Ex14-2-3`—
llegan al 90 % del cap con tope 1, o sea que **no estaban en su punto fijo: las
truncaba el cap**.

Dos consecuencias:

1. Parte de la mejora en el árbol puede venir de que en esas instancias duras
   ahora se permite más trabajo, no solo de un reparto mejor. En las otras 143
   el pivoteo es el mismo (+1.7 %), así que el −7 % de celdas **no** se explica
   por ahí, pero la separación merece medirse.
2. **`step_cap` no es un límite de trabajo robusto**: su presupuesto efectivo
   depende de `max_iters`. Si su intención es acotar el trabajo por caja
   —y lo es: se agregó para cerrar una no terminación—, debería contar pivotes.
   Cambiarlo altera los resultados medidos, así que queda anotado como paso a
   seguir y no se tocó acá.
