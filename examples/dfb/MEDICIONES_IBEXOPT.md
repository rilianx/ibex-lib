# Medición dentro de `ibexopt`: DFB como reemplazo de PolytopeHull

Nivel 3 del §3.3 de [PLAN_MEJORA_DFB.md](PLAN_MEJORA_DFB.md): medir el
contractor **dentro del optimizador**, donde se llama en cada nodo y el tiempo
es la restricción real. Todo lo medido antes
([CORRECCION_M0.md](CORRECCION_M0.md), [MEDICIONES_M2_M3.md](MEDICIONES_M2_M3.md))
era un contractor sobre una caja.

- Banco: 132 instancias de `benchs/optim/easy` y `benchs/optim/medium`.
- Límite de 20 s por corrida, 12 procesos en paralelo.
- Runner: [results/run_opt.py](results/run_opt.py); datos en
  [results/ibexopt_barrido.csv](results/ibexopt_barrido.csv).

Configuraciones, donde `--lr=xn` agrega `CtcPolytopeHull` en punto fijo y
`--lr=no` lo omite:

| nombre | opciones | qué es |
|---|---|---|
| `acidhc4+ph` | `--filtering=acidhc4 --lr=xn` | **producción**: ACID+HC4 con PolytopeHull |
| `hc4+ph` | `--filtering=hc4 --lr=xn` | HC4 con PolytopeHull |
| `acidhc4` | `--filtering=acidhc4 --lr=no` | sin relajación lineal |
| `dfb` | `--filtering=dfb --lr=no` | **DFB en lugar de PolytopeHull** |
| `acid_dfb` | `--filtering=acid_dfb --lr=no` | DFB dentro de ACID, sin PolytopeHull |

---

## 1. Lo que este paso encontró: cinco defectos

Ninguno era visible midiendo el contractor sobre una caja aislada. Antes de
arreglarlos, `ibexopt --filtering=dfb` era inservible: 62 instancias resueltas,
15 abortos y la mitad del banco en timeout.

### 1.1 `bad_alloc`: no se chequeaba el `-1` del linealizador

El contrato de `Linearizer::linearize` es «devuelve el número de restricciones,
o **−1 si el sistema lineal es infactible**». `CtcDFBPropag::linearize` usaba
ese valor como número de filas: `A.resize(-1, nb_var-1)` intenta
`new Interval[-1]` y aborta. `CtcPolytopeHull` sí lo chequea y trata el −1 como
caja vacía. Corregido en `CtcDFBPropag` y `CtcDFBManager`, junto con el
`clear_constraints()` y el guard de caja no acotada que usa PolytopeHull.
Además se recupera una prueba de vacío válida.

### 1.2 La propagación HC4 no se ejecutaba sin linealización útil

`contract` tenía un `if (refA.nb_rows() <= 1) return;` antes de armar la cola,
así que con caja no acotada, sin restricciones o linealización trivial **no
corría ningún contractor**, aunque los HC4 no dependen de la matriz. En un
optimizador eso pasa en todo nodo donde la variable objetivo aún no tiene cotas,
o sea desde la raíz: `ibexopt --filtering=dfb` terminaba con «possibly unbounded
objective» en instancias que `hc4` resuelve. Ahora los HC4 corren siempre y solo
los DFB dependen de `refA` (`dfb_usable`).

### 1.3 `CtcDFB::contract` cambiaba el tamaño de la caja del llamador

Hacía `x_new.resize(x_ref.size())` para agregar las variables `b` y
`resize(nb_var)` al salir. Eso viola el contrato de `Ctc::contract`: `resize()`
libera el arreglo interno y `ibex::Optimizer::contract_and_bound` conserva
referencias a las componentes entre llamadas. Resultado: **lectura de memoria
liberada y segfault determinista a −O3**.

Vale registrar el método, porque dos caminos fueron falsos: con ASan **no se
reproducía** (ni a −O1 ni a −O3), porque sus redzones cambian el layout; y los
asserts de Ibex están activos, lo que descartaba un `operator[]` fuera de rango
en código propio. Lo resolvió **valgrind**, que señaló `operator delete[]` desde
`IntervalVector::resize` en `CtcDFB::contract`, sobre un bloque reservado por
`Cell::bisect`, y la lectura inválida posterior en `Optimizer::contract_and_bound`.
Ahora se trabaja sobre un buffer propio (`work`) y se copia de vuelta componente
a componente.

### 1.4 `real_impact` indexaba fuera de rango

Con la matriz marcador 1×1 (linealización inútil), `A[0][var]` está fuera de
rango. Lo introdujo el arreglo 1.2 al permitir que la propagación siguiera sin
linealización; se guarda con `dfb_usable` y con una comprobación defensiva en
`real_impact`.

### 1.5 La propagación no terminaba

Dos causas distintas:

1. En el reencolado, `old_box[v] = box[v]` estaba **dentro** del bucle de las
   restricciones de salida de `v`, así que no se ejecutaba cuando `v` no es
   salida de ninguna restricción — por ejemplo la variable objetivo del sistema
   extendido. La condición `ratiodelta >= ratio` quedaba verdadera para siempre
   y la cola se realimentaba sin fin.
2. Aun corregido eso, con cajas de cotas 1e100 `ratiodelta` sigue reportando
   cambios grandes con progreso absoluto despreciable. Se agregó un tope de
   aplicaciones de contractores por llamada (`CtcDFBPropag::max_propag_steps`,
   por defecto 50 × número de contractores). Cortar antes solo contrae menos.

Lo grave de la no terminación es que **el límite de tiempo de `ibexopt` no puede
actuar**, porque se comprueba entre nodos: `exinfinity2`, `exinfinity3`,
`ex8_5_3`, `ex8_5_4`, `ex8_5_5` quedaban colgadas indefinidamente y ahora
terminan con «optimization successful».

### 1.6 Trazas de depuración en el camino caliente

`CtcDFBPropag::linearize` imprimía la matriz completa (`cout << rows`) y
`CtcDFBManager` nueve líneas por contracción. Los ejemplos redirigen `cout` a
`/dev/null`; **`ibexopt` no**. El primer barrido midió I/O, no algoritmo. Quedan
comentadas, no borradas.

## 2. Resultados

### 2.1 Robustez

| config | resueltas | timeout | no acotada | infactible | aborta |
|---|---|---|---|---|---|
| `acidhc4+ph` (producción) | **116** | 10 | 5 | 1 | 0 |
| `hc4+ph` | 115 | 11 | 5 | 1 | 0 |
| `acidhc4` (sin relajación) | 60 | 64 | 5 | 1 | 2 |
| `dfb` | **103** | 22 | 5 | 1 | 0 |
| `acid_dfb` | 86 | 40 | 5 | 1 | 0 |

Antes de los arreglos: `dfb` resolvía 62 con 15 abortos. Ahora 103 sin abortos.

### 2.2 Costo, comparación pareada (mediana geométrica de la razón)

Sobre las 58 instancias que resuelven las cinco configuraciones:

| config | celdas × referencia | tiempo × referencia |
|---|---|---|
| `acidhc4+ph` | 1.00 | 1.00 |
| `hc4+ph` | 1.33 | 1.31 |
| `acidhc4` | 2.58 | 1.35 |
| `dfb` | 2.12 | 1.71 |
| `acid_dfb` | **1.32** | 6.05 |

### 2.3 DFB contra «sin relajación lineal»

Es la comparación que aísla lo que aporta DFB, con las 59 instancias que
resuelven ambas:

| | celdas | tiempo |
|---|---|---|
| `dfb` vs `acidhc4` | **×0.85** | ×1.37 |

**DFB reduce el árbol un 15 %** respecto de ACID+HC4 sin relajación, y sobre
todo lo hace mucho más robusto: 103 instancias resueltas contra 60. Hay casos
espectaculares: en `ex14_1_5`, `dfb` resuelve con **32 celdas** donde `acidhc4`
necesita **9146**.

### 2.4 Solidez

**0 instancias con óptimo incompatible** con la referencia, ni en `dfb` ni en
`acid_dfb` (criterio: intervalos `f*` disjuntos con tolerancia relativa 1e-4).

## 3. Veredicto sobre el objetivo

**DFB todavía no puede reemplazar a `CtcPolytopeHull`.** La configuración de
producción sigue ganando: `acid_dfb` necesita ×1.32 las celdas y ×6.05 el
tiempo, y `dfb` ×2.12 las celdas y ×1.71 el tiempo. Y resuelve 103 instancias
contra 116.

Pero el cuadro cambió de «inservible» a «competitivo en árbol, caro en tiempo»:

- **En árbol está cerca**: `acid_dfb` con ×1.32 celdas, e igual o menos celdas
  que la referencia en 30 de 86 instancias.
- **El problema es el tiempo**, ×6 en `acid_dfb`. Eso es consistente con el
  perfil del §2 de [MEDICIONES_M2_M3.md](MEDICIONES_M2_M3.md): el costo por
  pivote de DFB es alto (tableau denso de intervalos) y el pricing es el 54 %.
- **DFB sí supera claramente a no tener relajación lineal**, que es el escenario
  donde hoy tiene sentido usarlo.

### Consecuencia inmediata: la lista de candidatos de M2 pasa a tener sentido

Quedó apagada porque cambiaba contracción por tiempo (§5.2 de las mediciones) y
en el banco de una caja la contracción era la restricción. **Aquí el tiempo es
la restricción**, con ×6.05 de exceso, así que un −14 % de tiempo por −0.3 pp de
contracción probablemente sea buen negocio. Ese es el experimento que sigue:
repetir este barrido con `partial_pricing = true; partial_pricing_refresh = 1`.

## 4. El banco de una caja no se degradó

Tras todos los arreglos, sobre las 153 instancias
([results/baseline_153_final.csv](results/baseline_153_final.csv)):
2 cajas vaciadas (ambas correctas), grupo A **57**, B 40, C 48, 6 no terminan.
Idéntico a antes de esta ronda: los arreglos del optimizador no costaron nada en
el banco de una caja.

---

## 5. La lista de candidatos dentro del optimizador — **medida, no se activa**

Era el experimento que el §6.1 del plan proponía: en el banco de una caja la
lista de candidatos cambiaba contracción por tiempo y se apagó, y aquí el tiempo
es la restricción (×6 de exceso), así que el intercambio debía darse vuelta.

Datos en [results/ibexopt_partial.csv](results/ibexopt_partial.csv), con
`partial_pricing = true` y `partial_pricing_refresh = 1`:

| | celdas | tiempo | tiempo total |
|---|---|---|---|
| `dfb/partial` vs `dfb` (99 instancias) | ×0.982 | **×0.941** | 161.6 s → 146.8 s |
| `acid_dfb/partial` vs `acid_dfb` (83) | ×0.975 | **×0.948** | 183.1 s → 179.1 s |

Contra la configuración de producción, sobre las 80 instancias que resuelven las
cinco configuraciones:

| config | celdas × ref | tiempo × ref |
|---|---|---|
| `dfb` | 2.40 | 1.80 |
| `dfb/partial` | 2.38 | **1.70** |
| `acid_dfb` | 1.55 | 6.47 |
| `acid_dfb/partial` | 1.51 | **6.18** |

**Veredicto: el intercambio sí se da vuelta, pero la ganancia es chica.** Ya no
cuesta contracción —al contrario, las celdas bajan un 2 %— y el tiempo mejora un
5–6 % de forma consistente, sin ningún óptimo incompatible. Pero **no arregla el
×6**: lo baja a ×6.18. Y con `dfb` cuesta tres instancias resueltas (103 → 100).

No alcanza para activarla por defecto. Lo importante es la conclusión de
diagnóstico: **el ×6 no es un problema de pricing a nivel de optimizador**. El
pricing es el 54 % del tiempo *del contractor*, y optimizarlo da unidades
porcentuales; el exceso frente a PolytopeHull es estructural —costo por pivote
sobre un tableau denso de intervalos y `2n` matrices por caja— y eso es M1.
