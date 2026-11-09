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
