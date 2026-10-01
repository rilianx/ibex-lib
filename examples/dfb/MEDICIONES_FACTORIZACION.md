# Factorización de la base: equivalencia verificada, y el lever resulta ser otro

Objetivo: darle a DFB la misma maquinaria que tiene `CtcPolytopeHull` —
factorización LU con actualizaciones en vez de un tableau denso— para poder
compararlos **en igualdad de condiciones**. La desventaja de costo medida
(×2.24 contra ×1.28) puede venir de que a DFB le falta una técnica estándar, no
de una debilidad del método, y mientras no se cierre esa asimetría la
comparación está confundida con la diferencia de ingeniería.

## 1. La equivalencia: verificada

Lo primero era establecer que el pivoteo de DFB **es** un cambio de base, para
poder reemplazar el tableau por la factorización sin cambiar el algoritmo. La
afirmación precisa: en todo momento

    Af = B⁻¹·[Ā | I] = [ B⁻¹Ā | B⁻¹ ]

donde `Ā` es la matriz flotante `m×na` (con la columna `k` negada si se contrae
la cota superior) y `B` la matriz de columnas básicas, una por fila. De ahí:

- la **parte λ** de una fila es `e_jᵀB⁻¹`, o sea **un BTRAN**;
- la **parte A** de esa fila es `λᵀĀ`;
- la fila objetivo no es una excepción: su columna básica es `k`, y por eso
  tiene 1 en `k` y 0 en las demás básicas, que es la condición de «costo
  reducido nulo sobre las básicas»;
- el paso `Af[0] += α·Af[j]` seguido de la eliminación de la columna entrante es
  exactamente la matriz elemental del cambio de base, con `α = −t₀ᵢ/t_jᵢ`.

[test_dfbbasis.cpp](test_dfbbasis.cpp) lo comprueba: construye matrices
aleatorias, aplica la eliminación del camino denso y la misma secuencia de
pivotes sobre [ibex_DFBBasis.h](ibex_DFBBasis.h), y compara **todas** las filas,
parte A y parte λ.

    casos comparados: 51016   fallos: 0   peor diferencia: 5.68e-13
    factorizaciones=60 updates=242 btran=2022

Las 242 actualizaciones tipo Forrest–Tomlin contra 60 factorizaciones muestran
que `SLUFactor::change()` mantiene la equivalencia, no solo `load()`.

**Conclusión**: la representación factorizada es un reemplazo válido del
tableau. Lo que sigue es si conviene.

## 2. Dónde está el costo, medido antes de escribir el motor

El plan asumía que la etapa 1 —calcular las filas a demanda con BTRAN en vez de
materializarlas— era el primer paso natural. Al derivar el pricing factorizado
apareció que no puede funcionar, y la medición lo confirma.

### 2.1 El pricing necesita las `m` filas transformadas, y eso es intrínseco a la regla

El bucle de `largest_impact_f` elige, para cada columna `i`, una esquina de
`work[i]`: mira `γᵢ = Af[0][i]` cuando `|γᵢ| ≥ 1e-5`, y mira `Af[j][i]` cuando
no. **Solo el primer caso es lineal en la fila `j`**, y por lo tanto solo ese
caso se puede calcular para las `m` filas a la vez con un único FTRAN
(`L = B⁻¹Āv` para un `v` fijo). El segundo exige la fila `j` transformada.

Medido con [probe_price.cpp](probe_price.cpp) sobre 68 instancias:

| | mediana | solo n ≥ 40 |
|---|---|---|
| fracción de columnas con `|γᵢ| ≥ 1e-5` | 0.090 | **0.016** |

O sea que en las instancias grandes la esquina la decide la fila `j` en el
**98.4 %** de las columnas. La descomposición lineal cubre el 1.6 %: no sirve.

La razón es estructural. Con `nnz(λ)` chico —y lo es, mediana 0.13— `γ = λᵀĀ`
hereda la dispersión de `Ā`, así que `γᵢ = 0` en casi todas las columnas. La
regla de DFB es entonces una **evaluación por intervalos de la fila**, del tipo

    acc_incr_j = Σᵢ ub(a_jᵢ·work[i]),   acc_decr_j = Σᵢ lb(a_jᵢ·work[i])

que separada en parte lineal y parte de radio queda
`L_j ± Σᵢ |a_jᵢ|·rad(work[i])`. El segundo término es una **norma L1 ponderada
de la fila**, o sea una cantidad tipo *steepest edge*: no se calcula sin la fila
y no tiene actualización exacta barata (la L2 sí, que es de donde salen devex y
steepest edge).

**Consecuencia**: calcular las filas a demanda costaría `m` BTRAN por pivote.
La etapa 1 del plan, tal como estaba escrita, es necesariamente más lenta. No
hace falta medirla: el motivo es la regla de pricing, no la implementación.

### 2.2 El lever real: el tableau denso hace 40–70× de trabajo sobre ceros

Las otras dos mediciones del mismo probe dicen dónde está el desperdicio.

| | mediana | solo n ≥ 40 |
|---|---|---|
| fracción de entradas **no nulas** que el pricing recorre | 0.144 | **0.023** |
| fracción de filas que un pivote **modifica** | 0.500 | **0.069** |

Y las filas modificadas son un número **constante**, no una fracción:

| instancia | n | m | filas que un pivote modifica |
|---|---|---|---|
| BroydenBanded-200 | 200 | 400 | 16.8 |
| BroydenBanded-140 | 140 | 280 | 16.7 |
| BroydenBanded-120 | 120 | 240 | 16.6 |
| BroydenTri-0100 | 100 | 200 | 7.5 |
| BroydenTri-0090 | 90 | 180 | 7.5 |
| CountercurrentReactors-44 | 44 | 176 | 18.9 |

`make_column_identity_f` resta `f·fila_pivote` de toda fila con `f ≠ 0`, o sea
de `nnz(columna entrante)` filas. Ese número **no crece con `m`**: 16.8 filas
con `m = 400` y 16.7 con `m = 280`.

Puesto en operaciones, sobre BroydenBanded-200 (`m = 400`, `na = 711`):

| por pivote | trabajo denso | trabajo útil |
|---|---|---|
| eliminación | `m·(na+m)` = 604 400 | 16.8 filas ≈ 25 000 |
| pricing | `(m−1)·na` = 283 689 | 2.3 % ≈ 6 500 |

**El tableau denso gasta entre 24× y 43× de más en las instancias grandes**, y
el sobrecosto es todo aritmética sobre ceros. En las instancias chicas el margen
es mucho menor (50 % de filas modificadas, 14 % de no nulos), lo que explica que
el perfil global mostrara el pricing en 37 % y la eliminación en 3.4 %: el
promedio está dominado por instancias donde no hay nada que ganar.

## 3. Lo que esto cambia en el plan

La prioridad estaba al revés. El plan decía «explotar la dispersión: solo
instancias grandes, y ahí el problema es memoria. Después de la factorización,
si acaso». Con lo medido:

1. **La dispersión del tableau transformado es el lever**, y vale 24–43× en las
   instancias grandes. No es un detalle posterior a la factorización.
2. **La factorización no abarata el pricing** —no puede, por la regla— pero sí
   es lo que impide que la dispersión se degrade: una fila recalculada como
   `λᵀĀ` no arrastra el *fill-in* acumulado de las eliminaciones sucesivas, y
   permite reconstruir una fila cuando su densidad crece. Sigue valiendo además
   por las otras dos razones del plan: la memoria (las 6 instancias que no
   terminan) y el cambio de objetivo barato, que es el habilitador de la cadena.
3. **La arquitectura que sale de esto** es un híbrido, y en este orden:
   - filas del tableau **dispersas**, pricing y eliminación que recorren solo
     los no nulos, y eliminación que toca solo las `nnz(columna entrante)` filas
     afectadas;
   - la factorización como fuente de `λ` y como reinicio de *fill-in*:
     recalcular una fila desde `λᵀĀ` cuando se densifica;
   - `change()` para el cambio de objetivo, que recién ahí vuelve compatibles la
     cadena de cotas y el intercalado fino con HC4.

El motor de base factorizada ya está escrito y verificado
([ibex_DFBBasis.h](ibex_DFBBasis.h), [ibex_DFBBasis.cpp](ibex_DFBBasis.cpp)),
con `load`, BTRAN, FTRAN y `change()` tipo Forrest–Tomlin, y con
`DFB_REFACTOR_ALWAYS` / `DFB_REFACTOR_PERIOD` para controlar cada cuánto se
refactoriza. Lo que falta es el tableau disperso que lo use.

## 4. El tableau disperso: implementado, equivalente y ×2.2

[ibex_DFBSparse.h](ibex_DFBSparse.h) reemplaza la matriz densa `Af` por filas
dispersas con la fila 0 densa, y `CtcDFB` lo usa con `DFB_SPARSE=1`. Decisiones
de representación:

- **Espacio de columnas único** `0..na+m-1` (parte A y parte λ juntas), porque
  toda combinación de filas tiene que tocar las dos para que λ siga siendo
  veraz.
- **Filas 1..m-1 dispersas**, índices ordenados; la eliminación recorre solo las
  filas con entrada no nula en la columna entrante.
- **Fila 0 densa**: es la fila objetivo, el pricing y el ratio test la leen en
  orden aleatorio por columna, y es la que acumula *fill-in*. Cuesta `na+m`
  doubles, no `m·(na+m)`.
- Se descarta el **cero exacto**; los valores chicos pero no nulos se
  conservan, porque el camino denso también los conserva y el objetivo es que
  las dos representaciones den la misma secuencia de pivotes.

### 4.1 La trampa que costó encontrar: `-x` no es `0 - x`

El test unitario [test_dfbsparse.cpp](test_dfbsparse.cpp) pasaba con **tolerancia
0** sobre 132 277 comparaciones, y sin embargo en el banco real **20 de 148
instancias** daban otra secuencia de pivotes.

Se agregó un modo autoverificante (`DFB_SPARSE_CHECK=1`) que mantiene las dos
representaciones a la vez y compara, en cada pivote, el tableau completo y las
decisiones de pricing y ratio test. Localizó la primera divergencia en una
entrada concreta: fila destino con `a = 0`, misma `f`, misma fila pivote
escalada (`48.921889824212422` en las dos), y aun así

    denso    0.99840591477984575
    disperso 0.99840591477984586

La causa: **Ibex corre con redondeo dirigido** (`-frounding-math`, y gaol deja
el modo hacia arriba), y con redondeo dirigido el redondeo **no es simétrico
respecto del signo**:

    round_up(-(f·b))  ≠  -round_up(f·b)

El camino denso calcula `Af[i][c] -= f*Af[row][c]` con `Af[i][c] == 0`, o sea
`0.0 - f·b`. Yo había escrito `-f·b`, que es lo mismo matemáticamente y **1 ulp
distinto** en aritmética dirigida. Ese ulp se amplifica pivote a pivote hasta
cambiar los desempates del ratio test.

La lección, anotada en las trampas del plan: **al reimplementar una operación
numérica hay que reproducir la expresión, no solo su valor matemático.** Y el
corolario metodológico: un test unitario con tolerancia 0 puede pasar y no
cubrir el caso; lo que lo encontró fue comparar las dos implementaciones
corriendo en paralelo sobre los datos reales.

### 4.2 Resultado

Corregido eso, sobre las 148 instancias comparables del banco de una caja:
**0 resultados distintos** —mismo perímetro, mismos pivotes, mismos vaciados—.

| | `t_dfb` denso | disperso | | DFB/PH denso | disperso |
|---|---|---|---|---|---|
| todas (148) | 6.686 s | **3.022 s** | ×2.21 | ×6.96 | **×3.15** |
| n ≥ 40 (23) | 4.842 s | **1.697 s** | ×2.85 | ×11.64 | **×4.08** |
| n ≥ 80 (10) | 4.301 s | **1.377 s** | ×3.12 | ×15.63 | **×5.00** |

Y **`Eiger-1000` pasa de `std::bad_alloc` a terminar**: la memoria era el
problema que el §5.6 del plan reportaba como «instancias que no terminan», y la
representación dispersa lo resuelve sin cambiar el algoritmo.

La ganancia crece con el tamaño, como predecía la medición del §2.2 (×2.2
global, ×3.1 en `n ≥ 80`), y queda todavía un factor ~3–5 contra PolytopeHull en
el banco de una caja. Ese banco corre DFB **hasta su punto fijo y sin HC4**, así
que no es la comparación del árbol; para eso está la medición en `ibexopt`.

## 5. En el árbol: la brecha de costo se cierra, y lo que queda es tamaño de árbol

191 instancias, 30 s, las cinco configuraciones en la misma corrida. Sobre las
97 que **todas** resuelven, relativo a producción (`acidhc4 --lr=xn` con la N-S
optimizada):

| config | resuelve | celdas | × | cpu | × | celdas (geom) | cpu (geom) |
|---|---|---|---|---|---|---|---|
| `acidhc4+ph` (producción) | **137** | 31 222 | 1.00 | 59.89 s | 1.00 | 1.000 | 1.000 |
| `dfb` denso | 116 | 97 254 | 3.11 | 95.27 s | 1.59 | 2.400 | 1.217 |
| **`dfb` disperso** | **118** | 97 254 | 3.11 | **92.00 s** | 1.54 | 2.400 | **1.171** |
| `acid_dfb` denso | 104 | 71 868 | 2.30 | 218.28 s | 3.64 | 1.434 | 3.876 |
| **`acid_dfb` disperso** | **105** | 71 868 | 2.30 | **193.00 s** | 3.22 | 1.434 | **3.385** |

**Celdas exactamente iguales** entre denso y disperso (97 254 y 71 868 en las
dos): la búsqueda es la misma, así que la diferencia de tiempo es puro ahorro.
0 óptimos incompatibles.

### 5.1 Por qué en el árbol el ×2.2 se vuelve un 3–12 %

`dfb` gana **3.4 %** de cpu y `acid_dfb` **11.6 %**, contra ×2.21 del contractor
en aislamiento. La razón está medida, y es el tamaño de las instancias:

| banco | `ext_n` mediana | filas que un pivote modifica | no nulos que recorre el pricing |
|---|---|---|---|
| `data_tests` con n ≥ 40 | — | **0.069** | **0.023** |
| **banco de `ibexopt`** (68 inst.) | **14** | **0.500** | **0.475** |

En el banco del optimizador la mediana de `ext_n` es **14** y solo 6 de 68
instancias llegan a `n ≥ 40`. A ese tamaño las matrices linealizadas son casi
densas —el 47.5 % de las entradas que el pricing recorre son no nulas—, así que
**no hay dispersión que explotar**. El ×2.2 es real y está donde se predijo (las
instancias grandes), pero el banco de `ibexopt` no tiene ninguna.

`acid_dfb` gana el triple que `dfb` porque ACID llama al contractor muchas más
veces, o sea que el contractor pesa más en su tiempo total.

### 5.2 La respuesta a «comparar en igualdad de condiciones»

Igualada la representación del tableau, el resultado es:

- **La brecha de costo por nodo está esencialmente cerrada**: `dfb` está a
  **×1.17 de cpu geométrica** contra producción, y su desventaja ya no es de
  implementación.
- **Lo que lo mantiene atrás es el árbol**: ×2.40 celdas geométricas sin ACID y
  ×1.43 con ACID. Con ACID el árbol casi se iguala, pero el tiempo se va a
  ×3.39 porque ACID multiplica las llamadas al contractor.
- **Resuelve 118 contra 137.**

O sea que la pregunta cambia de lugar. Ya no es «DFB cuesta más por cota», que
era en buena medida un problema de ingeniería y se corrigió; es **«DFB poda
menos por nodo que PolytopeHull, y la configuración que iguala la poda (ACID)
cuesta demasiado»**. Eso apunta a dos cosas concretas, ninguna de las cuales es
la factorización:

1. **Calidad de la cota por pivote** —el pricing y el criterio de terminación—,
   que es lo que decide cuánto poda cada llamada.
2. **El costo de `acid_dfb`**, que hoy es ×3.4 y donde el contractor sí es la
   mayor parte del tiempo (91.9 % en `alkyl`, §8.1 del plan). Ahí la
   factorización **sí** ayudaría, porque ACID cambia cotas repetidamente y el
   cambio de objetivo barato es exactamente lo que falta.

## 6. La sustitución directa, que era la comparación que faltaba

Las configuraciones medidas hasta acá (`--filtering=dfb` y
`--filtering=acid_dfb`) cambiaban **a la vez el componente y la composición**:
ponen DFB dentro de una propagación propia (`CtcDFBPropag`, con HC4
intercalado) y, en el segundo caso, bajo ACID. Ninguna aísla la sustitución que
el objetivo pide, que es **poner DFB en el lugar exacto de `CtcPolytopeHull`**.

Se agregó `--lr=dfb` en
[ibex_Optimizer05Config.cpp](ibex_Optimizer05Config.cpp): DFB ocupa el hueco de
`CtcPolytopeHull`, con el mismo envoltorio `CtcFixPoint(CtcCompo(·, hc44xn))` y
el mismo `relax_ratio`. Detalle que importa: se usa `CtcDFBPropag` con
`only_dfb = true`, **sin los HC4 internos**, porque PolyHull es un contractor
lineal a secas y el HC4 del bucle lo aporta el envoltorio; meter uno con HC4
adentro contaría HC4 dos veces y volvería a mezclar dos cambios.

Así, `--filtering=acidhc4 --lr=dfb` es la configuración de producción con **lo
único cambiado siendo el contractor lineal**.

### 6.1 DFB contra PolyHull en el mismo hueco

191 instancias, 30 s. Media geométrica de los cocientes **por instancia**, sobre
las 121 donde las dos resuelven:

| filtrado que acompaña | celdas | cpu | DFB más rápido en |
|---|---|---|---|
| `acidhc4` (producción) | ×1.384 | **×1.069** | 77 / 121 |
| `hc4` solo | ×1.268 | **×0.788** | 77 / 121 |

Instancias resueltas, de 190:

| config | resuelve | celdas (geom, vs producción) | cpu (geom) |
|---|---|---|---|
| `acidhc4` + PolyHull (producción) | **137** | 1.000 | 1.000 |
| **`acidhc4` + DFB** | **123** | **1.377** | **1.022** |
| `hc4` + PolyHull | 135 | 1.211 | 1.184 |
| **`hc4` + DFB** | **123** | 1.498 | **0.822** |
| `dfb` (composición propia) | 118 | 2.376 | 1.158 |
| `acid_dfb` | 104 | 1.440 | 3.487 |

0 óptimos incompatibles en las seis.

### 6.2 Qué dice

1. **La sustitución directa es, con diferencia, la mejor configuración de DFB
   medida**: 123 instancias contra 118 de `dfb` y 104 de `acid_dfb`, y una media
   geométrica de celdas de 1.38 contra 2.38 y 1.44. Las composiciones
   alternativas que se habían probado durante todo el trabajo **le estaban
   costando a DFB casi un factor 2 de árbol**, y eso no era una propiedad del
   método.
2. **El costo por nodo está a la par o mejor.** Con `hc4` DFB es **21 % más
   rápido** que PolyHull en el mismo hueco, y es más rápido en 77 de 121
   instancias. Con ACID la ventaja se pierde (×1.069), lo que es coherente con
   que ACID multiplique las llamadas al contractor y amplifique la diferencia de
   costo por llamada.
3. **Lo que queda es poda**: ×1.27–1.38 más celdas, y 12–14 instancias menos
   resueltas. Es el mismo diagnóstico del §5.2, pero ahora con la comparación
   limpia y con la brecha bastante más chica de lo que parecía.

**Consecuencia para el plan**: la configuración objetivo debería ser
`--filtering=acidhc4 --lr=dfb`, no `--filtering=dfb` ni
`--filtering=acid_dfb`. Y la pregunta abierta «¿`dfb` o `acid_dfb`?» queda mal
planteada: ninguna de las dos era la sustitución.
