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

## 7. Lectura completa del segundo lote de PDF (Messine–Trombettoni 2019, Jansson 2004 LP, Kearfott et al. 2013, Adjiman et al. 1996)

### Messine, F., Trombettoni, G. (2019), *Reliable bounds for convex relaxation in interval global optimization codes*, AIP Conf. Proc. 2070, 020050

- **Teorema 1 (solo cajas):** para f convexa en X y cualquier x̃ ∈ X,
  f(x*) ≥ f(x̃) + (z − x̃)ᵀ∇f(x̃), con z_i el extremo de X_i según el signo de ∂f/∂x_i.
- **Teorema 2 (restricciones):** f y g convexas, igualdades **lineales** Ax = b,
  μ̃ ≥ 0, λ̃ libre: f(x*) ≥ L(x̃, μ̃, λ̃) + (z − x̃)ᵀ∇ₓL(x̃, μ̃, λ̃).
  La prueba usa Taylor de 2º orden del Lagrangiano y que Hₓ ⪰ 0 en X.
- **Corolario 1 (versión con intervalos):** si el signo de ∂L/∂x_i no se decide con
  intervalos, toma z_i = [x_i, x̄_i] completo. Es exactamente lo que hacemos.
- Reconocen que Jansson (JOGO 2004) da el mismo teorema por otra vía.
- Anuncian resultados numéricos futuros; el artículo no trae experimentos.
- **Qué nos aporta.**
  - En la prueba, la convexidad de f y g y la linealidad de las igualdades solo se usan
    para tener Hₓ L(·, μ̃, λ̃) ⪰ 0 en X **para los multiplicadores fijos**. Es una
    condición suficiente, no necesaria.
  - Nuestro aporte es justamente debilitar la hipótesis: basta **certificar en la caja**
    que ∇²ₓL(·, μ̃, λ̃) ⪰ 0 (Gershgorin, matriz de comparación, Cholesky intervalar).
    Así entran objetivo y restricciones no convexos e igualdades no lineales.
  - Si no se certifica, se desplaza con d_i al estilo αBB. La cota pasa a ser válida
    para el NLP original no convexo.
  - Ellos lo plantean para una **relajación convexa** resuelta con Ipopt. Nosotros
    usamos el Ipopt del **upper bounding** sobre el problema original, así que los
    multiplicadores salen gratis.
  - Detalle técnico: nuestras igualdades están relajadas a |h| ≤ ε_h, por eso restamos
    Σ|λ_j|ε_h. Ellos no lo necesitan porque sus igualdades son exactas y lineales.
- **Errata menor:** escriben λ ∈ ℝⁿ y μ ∈ (ℝⁿ)⁺ donde corresponde ℝ^q y (ℝ^p)⁺. Además,
  en el desarrollo final el término −λ̃ᵀb debe ir con λ̃ᵀA x̃ dentro del gradiente; la
  forma final es correcta.

### Jansson, C. (2004), *Rigorous lower and upper bounds in linear programming*, SIAM J. Optim. 14(3):914–935

- **Ojo: no es el paper que pedí.** Este es el de **LP** del SIAM. El que cita
  Messine–Trombettoni como origen del teorema convexo es *A rigorous lower bound for the
  optimal value of convex optimization problems*, JOGO 28(1):121–137 (2004). Sigue
  faltando.
- **Contenido.**
  - Cota inferior rigurosa en O(n²) desde un dual aproximado (Teorema 6.1). Es
    equivalente a Neumaier–Shcherbina y fue desarrollada en paralelo.
  - Cota superior rigurosa: deflaciona el problema y encierra un punto factible en el
    interior relativo (método de Hansen/Kearfott).
  - Soporta datos intervalares y cotas simples infinitas.
- **Para nosotros:**
  - Es la base de la corrección del LP que IBEX ya usa en el contractor de relajación
    lineal.
  - Su conclusión anuncia la generalización a programación convexa, que es el paper de
    JOGO 28.
  - Su esquema de deflación para la cota superior es análogo a nuestra corrección de
    la solución de Ipopt (búsqueda de una caja factible certificada alrededor del
    punto).

### Kearfott, Castille, Tyagi (2013), *A general framework for convexity analysis in deterministic global optimization*, JOGO 56:765–785

- **Análisis estructural sobre el grafo computacional (code list).** Identifica
  variables "no convexas", incluidas variables intermedias, cuyo subdividirlas basta
  para que el resto del problema sea convexo o quede bien acotado.
- **Objetivo:** reemplazar el B&B adaptativo por una malla no adaptativa sobre pocas
  variables, más predecible y paralelizable. Hace experimentos exploratorios sobre
  COCONUT Lib-1.
- **No es certificación numérica de la Hessiana por caja**, ni cota lagrangiana. Es
  simbólica y global: el problema completo o por subexpresión, no el Lagrangiano con
  multiplicadores fijos.
- **Relación con nuestro trabajo.**
  - Es complementario: su análisis podría decir a priori qué instancias tienen
    estructura convexa, y nuestra certificación es numérica, local a la caja y depende
    de λ̃.
  - Hay casos que su análisis no capta y el nuestro sí: el Lagrangiano puede ser
    convexo en una caja aunque f o g no lo sean (candidatos a verificar: polak3 y abel, cerradas
    con la cota).
  - También sugiere una línea futura: ramificar preferentemente sobre las variables
    responsables de d_i > 0, que son las "no convexas" en la caja.

### Adjiman, Androulakis, Maranas, Floudas (1996), *A global optimization method, αBB, for process design*, Comput. Chem. Eng. 20(Suppl.):S419–S424

- **Precursor de αBB I/II (1998).**
  - Descompone f en términos lineales, convexos, bilineales (envolvente de McCormick),
    cóncavos univariados y no convexos generales.
  - A los no convexos les suma Σα(x_i^L − x_i)(x_i^U − x_i).
- **Cómo obtiene α en esta versión:**
  - polinomio característico intervalar de la Hessiana más el teorema de Kharitonov
    (4 polinomios);
  - la raíz mínima da λ_min; el ejemplo muestra que es conservador (−4 frente a −2.5
    real).
- **No incluye Gershgorin escalado ni el α no uniforme por variable.** Esos (con la
  matriz de comparación, Hertz, etc.) están en Part I de 1998, que sigue faltando para
  citar nuestro d_i.
  - Nuestra fórmula d_i = max(0, Σ_{j≠i}|h_ij| − h̲_ii) es el "Gershgorin escalado con
    d = 1" de Part I.
  - Nuestra variante escalada corresponde a la de Part I con d = x̄ − x̲.
- **Diferencia clave con nosotros.**
  - αBB construye un **subestimador convexo de cada función** y resuelve un NLP convexo
    por nodo.
  - Nosotros desplazamos solo el **Lagrangiano agregado** y no resolvemos nada extra: la
    cota es la tangente en el punto de Ipopt más la penalidad cuadrática, con costo
    O(n²) por nodo.

### Evaluación de novedad (tras este lote)

- **El teorema base es conocido** (Jansson 2004 en JOGO; Messine–Trombettoni 2019):
  cota tangente del Lagrangiano para problemas convexos con un punto aproximado.
- **Lo nuevo:**
  1. **Hipótesis debilitada:** se certifica la convexidad del Lagrangiano con
     multiplicadores fijos **en cada caja**, no la del problema. Eso permite NLP no
     convexo con igualdades no lineales.
  2. **Extensión αBB del Lagrangiano agregado,** con el término cuadrático acotado en la
     caja. Es válida siempre, aunque no se certifique.
  3. **Los multiplicadores vienen del Ipopt de upper bounding,** sin relajación convexa
     ni subgradiente. Van Voorhis lo sugiere como trabajo futuro.
  4. **Corrección ε_h** para igualdades relajadas, y backoff adaptativo.
  5. **Evaluación empírica en un B&B de intervalos riguroso:** +12 instancias cerradas
     sobre 298 y −50 % de nodos.
- Messine–Trombettoni no reportan experimentos. Si existe una implementación de su
  teorema en IbexOpt (alguna rama de Trombettoni), habría que comparar.

## 8. Jansson (2004, JOGO) y αBB Part II (1998)

### Jansson, C. (2004), *A rigorous lower bound for the optimal value of convex optimization problems*, JOGO 28:121–137

**Problema (10).**
- min f(x) s.a. G(x) ≤ 0, Hx = h, x̲ ≤ x ≤ x̄.
- f y todas las G_i son convexas; las igualdades son **lineales**.
- Las cotas simples pueden ser infinitas.

**Lema 1.**
- Linealiza **cada función por separado** con la desigualdad de gradiente:
  f(x) ≥ f(x̃) + ∇f(x̃)ᵀ(x − x̃), y lo mismo para cada G_i.
- Obtiene un LP cuya solución dual factible (ỹ ≥ 0, z̃, u = d⁺, v = −d⁻), con defecto
  d = ∇f + ∇G ỹ − Hᵀz̃, da la cota
  f* ≥ f(x̃) + (G(x̃) − ∇G(x̃)ᵀx̃)ᵀỹ + hᵀz̃ + x̲ᵀd⁺ + x̄ᵀd⁻ − ∇f(x̃)ᵀx̃.

**Teorema 1 y Teorema 2.**
- El Teorema 1 es la versión con intervalos: datos inciertos p ∈ **p** y variables
  libres (J^±), donde exige d_j = 0 vía un sistema intervalar.
- El Teorema 2 muestra que la cota es **exacta en un punto KKT**.
- Sección 4: QP convexo, f* ≥ hᵀz̃ − ½x̃ᵀQx̃ + x̄ᵀd⁻.
- Lo plantea como **postprocesamiento** de la salida de cualquier solver no lineal,
  sin tocar su código.

**Equivalencia con nuestra cota.**
- Algebraicamente, el Lema 1 es nuestra cota tangente del Lagrangiano:
  L(x̃) + Σ_i d_i (z_i − x̃_i), con z_i el extremo de la caja según el signo de d_i = ∂L/∂x_i.
  Messine–Trombettoni lo reconocen ("same theorem, different proof").

**Diferencias con lo nuestro.**
1. **Convexidad.**
   - Jansson necesita **cada** f y G_i convexas, porque linealiza término a término, y
     las igualdades lineales.
   - Nosotros solo necesitamos que el **Lagrangiano agregado** L(·, ỹ, z̃) sea convexo en
     la caja, y lo certificamos.
   - Si no se certifica, agregamos −Σ d_i/2·max(y_i − x̃_i)². Con eso cubrimos f y g no
     convexas e igualdades no lineales.
2. **Origen de (x̃, ỹ, z̃).**
   - En Jansson vienen de resolver el problema convexo (o la relajación) con un solver
     local.
   - En nosotros vienen del Ipopt que ya corre para el upper bounding del problema
     original no convexo. El costo marginal es solo la Hessiana intervalar, O(n²) por
     llamada.
3. **Igualdades relajadas** |h| ≤ ε_h: restamos Σ|z̃_j|ε_h. Jansson no lo necesita.
4. **Cotas infinitas.**
   - Jansson las maneja con J^± y d_j = 0, lo que exige resolver un sistema
     intervalar.
   - Nosotros solo aplicamos la cota en cajas acotadas (los nodos del B&B), así que es
     más simple.
   - Si quisiéramos usarla en la raíz con dominios infinitos, su Teorema 1 es la
     receta.
5. **Parámetros inciertos.**
   - Jansson los modela.
   - En IBEX las constantes gruesas [a, b] cumplen ese papel: las tratamos como rangos
     en `NormalizedSystem`.

Jansson anuncia aplicarla en B&B ("Branch and bound algorithms can be made reliable by
using such rigorous lower bounds"), pero como cota de **relajaciones convexas**, no del
problema original.

### Adjiman, Androulakis, Floudas (1998), *αBB … II. Implementation and computational results*, Comput. Chem. Eng. 22(9):1159–1179

**Ojo: es Part II, no Part I.** Part I (pp. 1137–1158) trae la teoría y las fórmulas de
los métodos de α. Part II los compara empíricamente y usa la numeración de Part I:
- I.1 Gerschgorin;
- I.2 E-matrix;
- I.3 Mori–Kokame;
- I.4 lower bounding Hessian;
- I.5 Kharitonov;
- I.6 Hertz;
- II.1 Gerschgorin escalado (d_i = 1 o d_i = x_i^U − x_i^L);
- II.2 H-matrix;
- II.3 min. max. distance (SDP).

**Resultados relevantes.**
- **Costo del cálculo de α:** el **~2.5 % del CPU** en casi todos los métodos. Kharitonov
  sube a 11 % y la versión con SDP a 30 %. Coincide con nuestro cvx < 5 % del tiempo.
- **Mejores métodos:** E-matrix, Hertz y **Gerschgorin escalado con d = x^U − x^L**
  (ej. 3.1.3: 330 iteraciones contra 353 de Gerschgorin simple).
  - Las diferencias entre métodos son modestas, del orden del 10 %.
  - Esto es consistente con lo que observamos al probar nuestro Gerschgorin escalado:
    sin ganancia clara frente a d = 1.
  - Respalda dejar Gerschgorin sin escalar (test 2) por defecto.
- **Branching:**
  - La estrategia 4 ramifica según la "influencia global de cada variable en la calidad
    del subestimador".
  - Es el antecedente directo de la idea pendiente de ramificar sobre las variables con
    d_i > 0.
  - También proponen actualizar las cotas de las variables con un problema convexo por
    variable. Solo en una fracción de ellas, porque es caro: análogo a la reducción de
    rangos.
- **Comparación con αBB:**
  - αBB resuelve un **NLP convexo por nodo** con un α por término no convexo.
  - Nosotros no resolvemos nada adicional: tangente más penalidad cuadrática, en forma
    cerrada sobre la caja.

## 9. Papers que faltan (actualizado)

Prioridad alta (solo uno):
1. **Adjiman, Dallwig, Floudas, Neumaier (1998).** A global optimization method, αBB,
   for general twice-differentiable constrained NLPs — **I. Theoretical advances**.
   *Comput. Chem. Eng.* 22(9):**1137–1158**. — Para citar el Teorema del Gerschgorin
   escalado (nuestro d_i). Lo entregado fue Part II (pp. 1159–1179).

Prioridad media (opcionales):
2. **Jansson (2009)**, On verified numerical computations in convex programming, Japan
   J. Indust. Appl. Math. 26:337–363.
3. **Kearfott (2014)**, On rigorous upper bounds to a global optimum, JOGO 59(2).
4. **Araya, Messine, Ninin, Trombettoni (2024)**, hibridación afín + XTaylor, JOGO.
5. **Ryoo & Sahinidis (1996)**, branch-and-reduce, JOGO 8.
6. **Li & Sun (2000)**, Local convexification of the Lagrangian, JOTA 104:109–120.
7. **Dür & Horst (1997)**, Lagrange duality and partitioning, JOTA 95:347–369.

Ya leídos: Van Voorhis (2002), QIBEX (2025), revisión de Neumaier–Shcherbina,
Messine–Trombettoni (2019), Jansson (2004, SIAM, LP), **Jansson (2004, JOGO, convexo)**,
Kearfott–Castille–Tyagi (2013), Adjiman et al. (1996), **Adjiman et al. (1998, Part II)**.

## 10. Local solver calls in other global solvers (upper bounding)

- **BARON** (GAMS/minlp.com option docs): `NumLoc` = number of local searches in
  preprocessing (default -2: BARON's own rule; i.e. a multistart before/at the
  root), `DoLocal` = local search for upper bounding during the tree (default on,
  BARON's rule), `NLPSol` = dynamic choice among NLP solvers during the search.
  The rule itself is not documented.
- **ANTIGONE** (GAMS docs): feasible solutions by local NLP solves (CONOPT by
  default, SNOPT optional), 30 s per NLP solve by default; schedule not documented
  beyond that.
- **Couenne**: `local_optimization_heuristic` (default yes) and
  `log_num_local_optimization_per_level` (default 2): on average 2^k local NLP
  solves per tree level at randomly chosen nodes (all nodes up to depth k, then
  half, a quarter...): a depth-based decreasing schedule.
- **SCIP 8** (Bestuzheva et al., arXiv 2301.00587, sec. 2.8): subNLP heuristic
  (fix integers, presolve, Ipopt from the LP solution or points of other
  heuristics); whether to run at a node depends on the nodes since the last
  run, the iteration limit it would use and its past success (adaptive,
  success-based); iteration limit = 2x the average of previous runs. A
  multistart heuristic (sampling + constraint consensus + clustering) feeds
  starting points to subNLP for continuous problems.
- Consequence: an adaptive, success-based schedule for local solves is already
  standard (SCIP), and multistart at the root too (BARON, SCIP). Our call
  strategy is not a contribution by itself.
