# Revisión bibliográfica — mejoras de IbexOpt con Ipopt y cotas lagrangianas certificadas

Rama `ml-branching`, noviembre 2026. Basada en búsquedas web (resúmenes, páginas de
editorial, preprints y fuentes secundarias). **No leí el texto completo de ninguno de
los artículos citados**; la sección 4 lista los que hay que leer completos para
confirmar diferencias y novedad.

## 1. Qué aportamos (para ubicar cada componente)

| # | componente | archivo / flag |
|---|---|---|
| A | Ipopt sin umbral de ancho, punto de partida finito, sin arrancar desde puntos "en el infinito", corrección en cascada de la caja de certificación, protección contra cuelgues de Ipopt | `--ipopt-diam=-1 --ipopt-expansion --ipopt-expansion-first --ipopt-start-max --ipopt-call-limit` |
| B | Arranques guiados de Ipopt: presupuesto proporcional al tiempo de CPU, punto maximin dentro de la caja del nodo (lejos de soluciones previas), paciencia | `--ipopt-guided --ipopt-guided-patience` |
| C | Cota inferior por Lagrangiano convexo certificado: multiplicadores del Ipopt del upper bounding, Hessiano del Lagrangiano en intervalos, certificado PSD, plano tangente minimizado exacto sobre la caja, corrección por ε_h de las igualdades | `--cvx-close` |
| D | Desplazamiento diagonal tipo αBB del Lagrangiano (cota válida en toda caja), backoff | `--cvx-test 2 --cvx-backoff` |
| E | Bugs corregidos: relajación afín (`Affine2::operator=`), igualdades con constante de intervalo (hs118), tiempo de la raíz no contado, NaN en derivadas (MUMPS) | — |
| F | Métricas anytime (integrales primal y primal-dual) y medición de margen de la selección de nodos | `--bound-trace` |

## 2. Trabajos relacionados por componente

### A. Upper bounding con solver local en B&B de intervalos

- **Araya, Trombettoni, Neveu, Chabert (2014)**, *Upper bounding in inner regions for
  global optimization under inequality constraints*, J. Global Optim. 60(2):145–164
  ([RePEc](https://econpapers.repec.org/RePEc:spr:jglopt:v:60:y:2014:i:2:p:145-164),
  [HAL](https://hal-enpc.archives-ouvertes.fr/hal-01061701)); antecedente AAAI 2011,
  *Inner regions and interval linearizations for global optimization*
  ([AAAI](https://ojs.aaai.org/index.php/AAAI/article/view/7817)). Base de IbexOpt: regiones
  interiores (poliedro o caja totalmente factible) para el upper bound.
- **Kearfott (2014)**, *On rigorous upper bounds to a global optimum*, J. Global Optim.
  59(2):459–476 ([RePEc](https://econpapers.repec.org/RePEc:spr:jglopt:v:59:y:2014:i:2:p:459-476));
  **Kearfott (2005)**, *Improved and simplified validation of feasible points*
  ([preprint](https://www.reliable-computing.org/archive/preprints/2005_simplified_feasible_point_verification.pdf));
  documentación de GlobSol ([link](https://interval.louisiana.edu/GLOBSOL/whatisop/node12.html)).
  Punto aproximado → caja pequeña alrededor → prueba de existencia (Newton intervalar).
  Reconoce que con muchas igualdades la verificación es débil.
- **Relación con lo nuestro:** nuestra corrección en cascada (1e-6 y luego 1e-5) es una
  variante de la "caja pequeña alrededor del punto" de Kearfott/IbexOpt; la novedad es
  operativa (cuándo y desde dónde llamar a Ipopt, cajas no acotadas, cuelgues). Lo más
  cercano en espíritu en IbexOpt es Reyes & Araya (2024), ver G.

### B. Arranques múltiples y diversificación del solver local

- **Ugray, Lasdon, Plummer, Glover, Kelly, Martí (2007)**, *Scatter search and local NLP
  solvers: a multistart framework for global optimization*, INFORMS J. Comput.
  ([INFORMS](https://pubsonline.informs.org/doi/fpi/10.1287/ijoc.1060.0175),
  [preprint](https://www.uv.es/rmarti/paper/docs/global2.pdf)). OQNLP: puntos de
  scatter search, filtros que evitan arrancar cerca del borde o de mínimos ya
  encontrados (idea afín a nuestro maximin). Heurístico, fuera de un B&B riguroso.
- MLSL (multi-level single linkage) y sus variantes: misma familia (citado en Ugray et al.).
- **Relación con lo nuestro:** la diversificación por distancia a soluciones previas es
  conocida. Lo propio es hacerlo *dentro* de un B&B riguroso, eligiendo la caja del
  nodo best-first (región aún no descartada) y con presupuesto y paciencia.
  No encontré trabajos que hagan exactamente esto; búsqueda pendiente
  ("spatial branch and bound" + "local search" + "multistart").

### C. Cotas inferiores lagrangianas y certificación de convexidad

- **Van Voorhis (2002)**, *A global optimization algorithm using Lagrangian underestimates
  and the interval Newton method*, J. Global Optim. 24(3):349–370
  ([Springer](https://link.springer.com/article/10.1023/A:1020383700229)). **El trabajo
  más cercano.** B&B con cotas inferiores obtenidas minimizando un subestimador convexo
  del Lagrangiano restringido, más Newton intervalar cerca del óptimo; orientado a QCQP.
  Diferencias probables con lo nuestro (a confirmar leyendo el texto): ellos resuelven
  un problema convexo por nodo; nosotros usamos un único plano tangente minimizado en
  forma cerrada, con los multiplicadores que Ipopt ya entregó para el upper bound, en
  aritmética de intervalos y para NLP general.
- **Neumaier & Shcherbina (2004)**, *Safe bounds in linear and mixed-integer linear
  programming*, Math. Program. ([preprint](https://www.reliable-computing.org/archive/preprints/safe-bounds-review-revision-2.pdf),
  [Optimization Online](https://optimization-online.org/?p=9243)). Cota dual rigurosa
  desde multiplicadores aproximados con redondeo dirigido; requiere cotas finitas de
  las variables. Hay una revisión posterior que corrige la fórmula original
  (Guilbeau, Hossain, Karhbet, Kearfott et al.). Es el precedente directo de "dualidad
  débil + intervalos con multiplicadores aproximados", pero para LP.
- **Elloumi, Lambert, Neveu, Trombettoni (2025)**, *Global solution of quadratic problems
  using interval methods and convex relaxations*, J. Global Optim. 91(2):331–353
  ([IP Paris](https://researchportal.ip-paris.fr/en/publications/global-solution-of-quadratic-problems-using-interval-methods-and-/),
  [RePEc](https://ideas.repec.org/a/spr/jglopt/v91y2025i2d10.1007_s10898-024-01370-8.html)).
  QIBEX: relajación cuadrática convexa en cada nodo de un B&B de intervalos, con
  **certificación rigurosa de convexidad y de la validez de la cota**. Muy cercano en
  el aspecto "certificar convexidad dentro de un B&B de intervalos", pero para
  problemas cuadráticos y con otra relajación.
- **Hladík (2018)**, *Positive semidefiniteness and positive definiteness of a linear
  parametric interval matrix* ([arXiv 1704.05782](https://arxiv.org/abs/1704.05782v2));
  [diapositivas SWIM 2014](https://kam.mff.cuni.cz/~hladik/doc/2014-conf-SWIM-PsdParIntMat-talk.pdf);
  cotas de autovalores de matrices de intervalos ([publicaciones](https://kam.mff.cuni.cz/~hladik/publ/b2hd-Hla2013a.html)).
  Condiciones suficientes de PSD (p. ej. λ_min(A_c) ≥ ρ(A_Δ)), que es nuestro test
  "centro-radio". Probar PSD fuerte de una matriz de intervalos es NP-difícil.
- **Fourer, Maheshwari, Neumaier, Orban, Schichl (2010)**, *Convexity and concavity
  detection in computational graphs*, INFORMS J. Comput. 22(1):26–43
  ([preprint](https://arnold-neumaier.at/ms/convexity.pdf)). Detección simbólica de
  convexidad (alternativa/complemento al test con Hessiano en intervalos).
- **Hansen & Walster**, *Global Optimization Using Interval Analysis* (libro) y la
  literatura de Karlsruhe: **test de no convexidad** (si f no es convexa en una caja,
  esta no contiene un minimizador interior), tests de monotonía y punto medio. Es el uso
  inverso del Hessiano en intervalos (descartar en vez de acotar).
- **Neumaier (2004)**, *Complete search in continuous global optimization and constraint
  satisfaction*, Acta Numerica ([preprint](https://optimization-online.org/wp-content/uploads/2003/10/765.pdf)).
  Panorama; menciona chequeo de convexidad con intervalos dentro de B&B.
- Debate reciente sobre cotas válidas desde multiplicadores de Ipopt en NLP
  (dualidad débil vale con cualquier multiplicador): Pyomo/mpi-sppy
  ([PR #838](https://github.com/Pyomo/mpi-sppy/pull/838)). No es literatura revisada
  por pares, pero muestra que el problema está vigente.

### D. αBB y el desplazamiento diagonal

- **Androulakis, Maranas, Floudas (1995)**, *αBB: a global optimization method for
  general constrained nonconvex problems*, J. Global Optim.
  ([Springer](https://link.springer.com/article/10.1007/BF01099647)).
- **Adjiman, Dallwig, Floudas, Neumaier (1998)**, *A global optimization method, αBB, for
  general twice-differentiable constrained NLPs — I. Theoretical advances*, Comput. Chem.
  Eng. 22(9):1137–1158 ([PDF](https://www.ou.edu/class/che-design/che5480-11/A%20global%20optimization%20method,%20aBBTheory(Adjiman%20et%20al)-98.pdf));
  **Parte II** (Adjiman, Androulakis, Floudas), 1159–1179
  ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S009813549800218X)).
  Métodos rigurosos para α desde el Hessiano en intervalos (Gershgorin escalado, etc.).
- Refinamientos del cálculo de α: [arXiv 1307.2786](https://arxiv.org/pdf/1307.2786)
  (Gershgorin escalado eficiente), [arXiv 1307.2791](https://arxiv.org/pdf/1307.2791)
  (efecto de las evaluaciones del Hessiano), [J. Global Optim. 2018](https://link.springer.com/article/10.1007/s10898-018-0718-y),
  [Yugoslav J. Oper. Res. 2024](https://yujor.erc.fon.bg.ac.rs/index.php/yujor/article/view/776),
  perturbaciones no diagonales ([JOTA](https://unpaywall.org/10.1007%2Fs10957-012-0033-6)).
- Cotas para QCQP vía duales de relajaciones SDP que inducen subestimadores convexos:
  [arXiv 1902.08861](https://arxiv.org/pdf/1902.08861).
- **Relación con lo nuestro:** el desplazamiento d_i = max(0, Σ|h_ij| − h_ii) es el α de
  Gershgorin simple aplicado **al Lagrangiano** (no a cada función por separado) y con
  el término cuadrático centrado en el punto de Ipopt en vez de (x−l)(u−x). αBB
  clásico subestima cada término no convexo y resuelve un NLP convexo por nodo;
  nosotros hacemos un único plano tangente. No encontré una formulación "αBB del
  Lagrangiano con multiplicadores de un solver local + tangente" publicada, pero la
  búsqueda no es concluyente. Mejora inmediata disponible: Gershgorin escalado.

### E. Reducción de rangos basada en dualidad

- **Ryoo & Sahinidis (1995, 1996)** (Comput. Chem. Eng. 19(5):551–566; J. Global Optim.
  1996), base del *branch-and-reduce* de BARON: con multiplicadores de la relajación se
  ajustan cotas de variables ([diapositivas GAMS](https://www.gams.com/archives/ads/poster_baron.pdf),
  [manual BARON](https://minlp.com/baron-user-manual)).
- **Zhang, Sahinidis, Nohra, Rong (2020)**, reducción de dominio basada en optimalidad
  para NLP/MINLP con desigualdades, J. Global Optim.
  ([RePEc](https://ideas.repec.org/a/spr/jglopt/v77y2020i3d10.1007_s10898-020-00886-z.html)).
- **Relación con lo nuestro:** usamos los multiplicadores para acotar el objetivo, no para
  ajustar las variables. Una extensión natural: con la misma cota tangente, reducir el
  rango de cada variable (ya tenemos el gradiente del Lagrangiano). Pendiente.

### F. Relajación afín (contexto del bug)

- **Ninin, Messine, Hansen (2015)**, *A reliable affine relaxation method for global
  optimization*, 4OR 13:247–277 ([preprint](https://optimization-online.org/wp-content/uploads/2012/10/3650.pdf)).
- **Araya, Messine, Ninin, Trombettoni (2024)**, *Hybridizing two linear relaxation
  techniques in interval-based solvers*, J. Global Optim. ([arXiv 2308.05244](https://bohrium.dp.tech/paper/arxiv/2308.05244)).
  Combina la relajación afín (IBBA) con XTaylor (IbexOpt): es la configuración
  `--relax both`. **Conviene revisar si el código usado en ese trabajo tiene el mismo
  `Affine2::operator=`**: el bug que encontramos produce cortes inválidos cuando la
  forma afín se vuelve inactiva (división por un intervalo que contiene 0); en nuestro
  benchmark cambió el resultado de wall (declarada infactible) y de haldmads con loup
  temprano.

### G. Selección de nodos, variables y anytime

- **Neveu, Trombettoni, Araya (2016)**, *Node selection strategies in interval branch and
  bound algorithms*, J. Global Optim. 64:289–304 ([dblp](https://dblp.uni-trier.de/rec/journals/jgo/NeveuTA16.html)).
- **Reyes & Araya (2024)**, *Node selection through upper bounding local search methods in
  branch & bound solvers for NCOPs*, J. Global Optim.
  ([Springer](https://link.springer.com/article/10.1007/s10898-024-01403-2)). Usa
  soluciones factibles de un método iterativo barato para seleccionar nodos: cercano a
  nuestros arranques guiados.
- **Araya & Neveu (2017)**, *lsmear: a variable selection strategy for interval branch and
  bound solvers*, J. Global Optim. 71(3):483–500.
- **Berthold (2013)**, *Measuring the impact of primal heuristics*, Oper. Res. Lett.
  41(6):611–614 ([ZIB](https://opus4.kobv.de/opus4-zib/frontdoor/index/index/docId/1788)).
  Integral primal; base de nuestras métricas anytime.
- **Vanaret, Gotteland, Durand, Alliot (2020)**, *Certified global minima for a benchmark of
  difficult optimization problems* ([arXiv 2003.09867](https://arxiv.org/pdf/2003.09867)).
  Charibde: evolutivo + intervalos compartiendo el upper bound (otro modo de hibridar
  búsqueda local y B&B riguroso).

## 3. Evaluación de novedad (provisional)

- **C + D (cota lagrangiana certificada, con multiplicadores reutilizados del upper
  bounding, plano tangente cerrado, desplazamiento tipo αBB del Lagrangiano, rigor
  completo incluida la relajación ε_h, dentro de un B&B de intervalos):** las piezas son
  conocidas (dualidad débil con multiplicadores aproximados, Neumaier–Shcherbina; αBB;
  certificación de convexidad, QIBEX; subestimadores del Lagrangiano, Van Voorhis). **No
  encontré la combinación publicada.** Es el aporte con más potencial (234 de 298
  instancias cerradas, abel/fccu/polak3 nuevas, robot −53% de tiempo). Riesgo
  principal: Van Voorhis (2002); hay que leerlo completo.
- **B (arranques guiados en B&B riguroso):** probablemente incremental; comparar con
  Reyes & Araya (2024) y la literatura de multistart.
- **A, E, F:** ingeniería y corrección de bugs; valiosos para el software y como
  advertencia (bug afín), no como contribución metodológica.

## 4. Papers que necesito y a los que no pude acceder (texto completo)

Prioridad alta (deciden la novedad):
1. **Van Voorhis, T. (2002).** A global optimization algorithm using Lagrangian
   underestimates and the interval Newton method. *J. Global Optim.* 24(3):349–370.
   doi:10.1023/A:1020383700229. — ¿Cómo construye el subestimador del Lagrangiano? ¿Usa
   multiplicadores de un solver local? ¿Es riguroso en aritmética de intervalos?
2. **Elloumi, Lambert, Neveu, Trombettoni (2025).** Global solution of quadratic problems
   using interval methods and convex relaxations. *J. Global Optim.* 91(2):331–353.
   doi:10.1007/s10898-024-01370-8. — Cómo certifican convexidad y validez de la cota
   (comparar con nuestro test).
3. **Adjiman, Dallwig, Floudas, Neumaier (1998).** αBB Part I. *Comput. Chem. Eng.*
   22(9):1137–1158 (hay un PDF de curso enlazado arriba, no verificado). — Métodos de α
   (Gershgorin escalado, Hertz–Rohn) para mejorar nuestro desplazamiento.
4. **Kearfott (2014).** On rigorous upper bounds to a global optimum. *J. Global Optim.*
   59(2):459–476. — Por si incluye cotas inferiores con multiplicadores.
5. **Neumaier & Shcherbina (2004).** Safe bounds in linear and mixed-integer linear
   programming. *Math. Program.* 99:283–296 (preprint enlazado) + la revisión que corrige
   la fórmula (Guilbeau, Hossain, Karhbet, Kearfott et al.; no la ubiqué).

Prioridad media:
6. **Ryoo & Sahinidis (1996).** A branch-and-reduce approach to global optimization.
   *J. Global Optim.* 8:107–138. — Reducción de rangos por dualidad (extensión E).
7. **Zhang, Sahinidis, Nohra, Rong (2020).** *J. Global Optim.* 77(3). doi:10.1007/s10898-020-00886-z.
8. **Hladík (2018)** capítulo completo (arXiv 1704.05782 disponible) y **Hladík (2013)**
   sobre cotas de autovalores.
9. **Araya, Messine, Ninin, Trombettoni (2024).** Hybridizing two linear relaxation
   techniques. *J. Global Optim.* — Para ver si el bug afín pudo afectar sus resultados.
10. **Reyes & Araya (2024).** Node selection through upper bounding local search methods.
    *J. Global Optim.* doi:10.1007/s10898-024-01403-2.
11. **Hansen & Walster (2004).** *Global Optimization Using Interval Analysis*, 2ª ed.,
    Marcel Dekker (libro) — formulación exacta del test de no convexidad y de unicidad.
12. **Fourer et al. (2010)** (preprint disponible) — si conviene sumar detección simbólica.

Búsquedas aún pendientes: subestimadores del Lagrangiano en solvers deterministas
(ANTIGONE, Couenne, SCIP), "rigorous" + "Lagrangian bound" en COCONUT/GloptLab (Domes,
Neumaier), y multistart dentro de spatial B&B.


## 5. Lectura de los textos completos disponibles en arXiv (nov 2026)

Desde este entorno solo arxiv.org es accesible (reliable-computing.org,
optimization-online.org, arnold-neumaier.at, hal.science, core.ac.uk, springer y
acm devuelven 403 por la política de red). Descargados y leídos (abstract, métodos,
conclusiones):

- **Hladík, *On the efficient Gerschgorin inclusion usage in the global optimization
  αBB method*** (arXiv 1307.2786). Para el α de αBB con Gershgorin escalado basta la
  matriz puntual h̲_ii en la diagonal y −max(|h̲_ij|, |h̄_ij|) fuera de ella (= nuestra
  "matriz de comparación"). La escala estándar es d = ancho de la caja (= nuestro
  `--cvx-test 3`). Propone dos heurísticas de mejora local de d; la segunda cumple tres
  condiciones necesarias de optimalidad, corre en tiempo polinomial y en promedio
  ~1 iteración. En su ejemplo baja α₂ de 12 a 3. **Mejora directa para nuestro
  desplazamiento**: reemplazar la escala por la heurística 2 de Hladík.
- **Hladík, *The effect of Hessian evaluations in the global optimization αBB method***
  (arXiv 1307.2791). El cálculo simbólico del Hessiano, con simplificaciones, reduce
  mucho el sobreestimado del Hessiano en intervalos y con ello α. IBEX ya diferencia
  simbólicamente (`Function::DIFF`); ver si la simplificación de IBEX es suficiente.
- **Hladík, *Positive semidefiniteness and positive definiteness of a linear parametric
  interval matrix*** (arXiv 1704.05782). Condiciones suficientes/necesarias baratas
  para PSD de matrices de intervalos (incluye la de centro-radio que implementamos).
  Motiva explícitamente el uso en B&B de intervalos para chequear convexidad.
- **Zamani, *New bounds for nonconvex QCQP*** (arXiv 1902.08861). Cotas para QCQP con
  multiplicadores afines (cuadráticas) y cuadráticos (cúbicas); las relajaciones SDP
  son duales de cotas cuadráticas. Relacionado a nivel de idea (multiplicadores para
  cotas), específico de QCQP y sin aritmética de intervalos.
- **Vanaret et al., *Certified global minima ...*** (arXiv 2003.09867). Charibde
  (evolutivo + intervalos). Contexto de hibridación; no aborda cotas lagrangianas.
- arXiv 2308.05244 **no** es el artículo de hibridación afín/XTaylor (es *A hybrid
  approach to joint spectral radius computation*); el identificador vino mal de una
  fuente secundaria. El artículo de Araya, Messine, Ninin, Trombettoni (2024) sigue
  pendiente.

No están en arXiv (búsqueda por título y autor): Van Voorhis (2002), QIBEX (Elloumi et
al. 2025), Neumaier–Shcherbina (2004), αBB I/II (1998), Kearfott (2014), Ryoo–Sahinidis.


## 6. Lectura completa de los PDF entregados (Van Voorhis 2002, QIBEX 2025, revisión de Neumaier–Shcherbina)

### Van Voorhis (2002), *A global optimization algorithm using Lagrangian underestimates and the interval Newton method*, JOGO 24:349–370

- **Problema:** solo QCQP con cotas finitas (f_i cuadráticas, Hessianos constantes).
- **Cota inferior:** minimiza un subestimador convexo del Lagrangiano restringido
  φ(x,v) = f⁰ + Σ v_i (fⁱ − b_i). El subestimador es término a término (ec. 3): los
  bilineales se escriben como (x_j ± x_k)² − (x_j² + x_k²) y los −p_j x_j² se
  linealizan con (u_j − x_j)(x_j − l_j) ≥ 0. Es una idea del tipo αBB, pero para
  cuadráticas y sin Hessiano en intervalos. Por nodo resuelve un QP convexo sobre la caja.
- **Multiplicadores:** iteraciones de subgradiente / direcciones de deflexión dentro de
  cada nodo (`cnt` iteraciones), heredados del padre.
- **Upper bound:** Newton sobre las KKT desde (x̂, v̂); Newton intervalar sobre las
  condiciones de Fritz-John para descartar cajas o probar unicidad.
- **Rigor:** no trata el redondeo de punto flotante (los intervalos solo aparecen en el
  Newton intervalar).
- **Conclusiones, cita textual relevante:** "it requires several iterations to establish
  good estimates of the Lagrange multipliers... Perhaps Lagrangian underestimates could
  be strategically employed at later stages of the branch-and-bound process, when good
  multiplier estimates have been obtained via RLT or some other relaxation strategy."
  **Esto es justamente lo que hacemos**: los multiplicadores vienen gratis de Ipopt (el
  upper bounding), sin iterar en el nodo, y la cota es un plano tangente (sin resolver un
  QP), con backoff para usarla solo donde rinde.

### Elloumi, Lambert, Neveu, Trombettoni (2025), QIBEX, JOGO 91:331–353

- B&B de intervalos sobre IbexOpt para MIQCQP con una relajación cuadrática convexa
  (tan fuerte como SDP) obtenida de un SDP resuelto una vez.
- **Rigor (QIBEX-R):** (1) certifica la convexidad del objetivo de la relajación
  corrigiendo S₀* + |δ|I con un algoritmo que encierra los autovalores evaluando el
  polinomio característico en intervalos; (2) certifica la cota inferior con el
  **Teorema 1** (Jansson 2004; Messine & Trombettoni 2019): para un problema **convexo**
  y cualquier (x̃, μ̃ ≥ 0, λ̃), f(x*) ≥ L(x̃,μ̃,λ̃) + (z − x̃)ᵀ∇ₓL(x̃,μ̃,λ̃), con z_i en la
  cota que minimiza cada término. **Es la misma fórmula de nuestro plano tangente.**
- Costo: certificar la cota multiplica el tiempo por 4–8 (implementación en AMPL
  interpretado).

**Consecuencia para la novedad:** la cota "Lagrangiano en x̃ + gradiente minimizado sobre
la caja" para problemas convexos es conocida (Jansson 2004) y ya se usa en IBEX
(QIBEX, Messine–Trombettoni 2019) sobre una relajación convexa construida aparte. Lo que
no aparece en ninguno de los tres:
1. aplicarla **directamente al problema no convexo original**, certificando por caja que
   su Lagrangiano es convexo (Hessiano en intervalos + test PSD) o convexificándolo con
   un desplazamiento diagonal tipo αBB del Lagrangiano;
2. con los **multiplicadores del solver local que ya hace el upper bounding** (Ipopt),
   sin SDP, sin subgradientes y sin resolver un problema por nodo;
3. para **NLP general** (no solo cuadráticos), con el control de costo (backoff) y la
   corrección por la relajación ε_h de las igualdades.

### Guilbeau, Hossain, Karhbet, Kearfott, Sanusi, Zhao, *A review of computation of mathematically rigorous bounds on optima of linear programs*

- Revisa el método de Neumaier–Shcherbina (2004) para cotas duales rigurosas en LP desde
  multiplicadores aproximados, **corrige dos errores tipográficos de fórmulas** del
  original y muestra su impacto. Solo LP. Remite a referencias nuevas que importan
  (abajo).

## 7. Papers que faltan (actualizado)

Prioridad alta:
1. **Jansson, C. (2004).** A rigorous lower bound for the optimal value of convex
   optimization problems. *J. Global Optim.* 28(1):121–137. — Origen del Teorema 1.
2. **Messine, F., Trombettoni, G. (2019).** Reliable bounds for convex relaxation in
   interval global optimization codes. *AIP Conf. Proc.* 2070(1):020050,
   doi:10.1063/1.5090017. — Su uso en IbexOpt/IBBA; lo más cercano a lo nuestro dentro
   de IBEX.
3. **Kearfott, R.B., Castille, J., Tyagi, G. (2013).** A general framework for convexity
   analysis in deterministic global optimization. *J. Global Optim.* 56(3):765–785. —
   Detección/certificación de convexidad en B&B riguroso.
4. **Adjiman, Dallwig, Floudas, Neumaier (1998)**, αBB I (Comput. Chem. Eng. 22(9)).

Prioridad media:
5. **Jansson, C. (2009).** On verified numerical computations in convex programming.
   *Japan J. Indust. Appl. Math.* 26(2–3):337–363.
6. **Kearfott (2014)**, On rigorous upper bounds to a global optimum, JOGO 59(2).
7. **Araya, Messine, Ninin, Trombettoni (2024)**, hibridación afín + XTaylor, JOGO (por el
   bug afín).
8. **Ryoo & Sahinidis (1996)**, branch-and-reduce, JOGO 8 (si se extiende a reducción de
   rangos).
9. **Li & Sun (2000)**, Local convexification of the Lagrangian function in nonconvex
   optimization, JOTA 104:109–120 (citado por Van Voorhis; convexificación del
   Lagrangiano).
10. **Dür & Horst (1997)**, Lagrange duality and partitioning techniques in nonconvex
    global optimization, JOTA 95:347–369 (convergencia del gap dual con particiones).

Ya no hace falta: Van Voorhis (2002), QIBEX (2025) y Neumaier–Shcherbina (cubierto por
la revisión de Guilbeau et al.).
