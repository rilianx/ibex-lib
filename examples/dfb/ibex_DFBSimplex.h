//============================================================================
//                                  I B E X
// File        : ibex_DFBSimplex.h
// Author      : (DFB) simplex dual de variables acotadas
//============================================================================

#ifndef __IBEX_DFB_SIMPLEX_H__
#define __IBEX_DFB_SIMPLEX_H__

#include "ibex_Matrix.h"
#include "ibex_IntervalVector.h"
#include "ibex_DFBBasis.h"

#include <vector>

namespace ibex {

/**
 * \brief Simplex dual de variables acotadas sobre la relajacion de DFB.
 *
 * ## Por que hace falta
 *
 * El pivoteo de DFB es un ascenso por coordenadas sobre
 * `f(lambda) = -sum_{i!=k} ub(gamma_i(lambda)*[z_i])`, que es concava y lineal
 * a trozos y cuyo maximo **es** la cota que calcula `CtcPolytopeHull` (dualidad
 * LP). El ascenso por coordenadas sobre una funcion no suave se detiene en
 * puntos no optimos, y esta medido que eso pasa en 33 de 148 instancias, sin que
 * lo arregle mas presupuesto (27x mas pivotes, los mismos grupos), ni el ratio
 * test, ni el umbral de parada. Ver MEDICIONES_PODA.md.
 *
 * El premio de cerrar esa brecha esta medido: DFB cuesta **x0.754 por nodo**
 * contra PolytopeHull en el mismo hueco, asi que igualar la poda lo dejaria
 * corriendo el mismo arbol a tres cuartos del tiempo.
 *
 * ## El LP y por que la base inicial es gratis
 *
 * Para la cota inferior de `z_k`:
 *
 *     min z_k    sujeto a   Abar z = 0,   l <= z <= u
 *
 * con `z = (x, b)`, `na = nx + m`. La linealizacion de DFB pone
 * `Abar[j][nx+j] = -1` exactamente, o sea que **la parte `b` de `Abar` es `-I`**.
 * Tomando como base las `m` columnas de `b`:
 *
 *  - `B = -I`, asi que `B^-1 = -I` y el tableau inicial es `-[Abar | I]`: no hay
 *    nada que factorizar ni invertir;
 *  - el objetivo es `c = e_k` con `k < nx`, o sea `c_B = 0`, de donde `y = 0` y
 *    los costos reducidos son `d = e_k >= 0`. Poniendo **todas** las no basicas
 *    en su cota inferior, el punto de partida es **dualmente factible**.
 *
 * O sea que no hace falta fase 1: se arranca dual-factible y se itera el simplex
 * dual hasta factibilidad primal, que es la condicion de optimalidad.
 *
 * ## Estructura del tableau
 *
 * `m+1` filas sobre un espacio de columnas unico `0..na+m-1` (primero la parte
 * A, despues la parte lambda, una por restriccion):
 *
 *  - **fila 0, densa**: la fila objetivo, `[d | -y]`, inicializada a `[e_k | 0]`.
 *    A diferencia del DFB actual, que sacrifica una fila de restriccion para
 *    llevar `gamma`, aca la fila objetivo es aparte y las `m` filas de
 *    restriccion quedan disponibles para la base.
 *  - **filas 1..m, dispersas**: `-[Abar_{r-1} | e_{r-1}]`, de modo que la
 *    columna basica inicial tiene coeficiente `+1` y el invariante
 *    «parte A = (parte lambda)^T * Abar» se cumple.
 *
 * Al terminar, `lambda() = y` y `gamma = y^T*Abar = e_k - d`, que es la relacion
 * valida con la que se certifica la cota en intervalos.
 */
/** \brief Vuelca el perfilado interno del simplex (DFB_SX_PERF=1). */
void perf_volcar();

class DFBSimplex {
public:
    enum Status {
        OPTIMAL,     //!< factible primal: la cota es el optimo del LP
        INFEASIBLE,  //!< el ratio test dual no encuentra candidato: LP infactible
        ITER_LIMIT,  //!< se agoto el presupuesto de pivotes
        UNBOUNDED,   //!< el LP no esta acotado: no hay cota que dar
        SINGULAR     //!< no se pudo pivotear (pivote degenerado)
    };

    DFBSimplex();

    /** \brief DFB_SX_NOSCALE=1 desactiva el escalado (para comparar). */
    static bool escalar;

    /** \brief Guarda la relajacion. Una vez por linealizacion. */
    void load(const Matrix& Abar, int nx);


    /**
     * \brief Resuelve `min z_k` (o `max`, con `maximize`) sobre la caja `z`.
     *
     * `z` tiene dimension `na` e incluye las cotas de las `b`.
     */
    /**
     * \param warm  si es true y ya hay una base cargada, se **reusa**: solo se
     *              recalcula la fila objetivo y se reubican las no basicas.
     *
     * En un LP de variables acotadas **cambiar de objetivo es gratis**: la
     * factibilidad dual es solo la condicion de signo `d_j >= 0` en la cota
     * inferior y `d_j <= 0` en la superior, y uno es libre de elegir en que
     * cota se apoya cada no basica. Asi que dada cualquier base, se recalcula
     * `d = c - Abar^T y` y se reubica cada no basica segun el signo: el punto
     * queda dual-factible sin un solo pivote. Es exactamente el mecanismo por
     * el que `CtcPolytopeHull` gasta 1.08 iteraciones por cota, y es lo que
     * hace que una sola base sirva para las 2n cotas.
     */
    Status solve(int k, bool maximize, const IntervalVector& z, int max_iter,
                 bool warm = true);

    /**
     * \brief Si no es NULL, `solve` deja aca el `lambda` despues de CADA
     * pivote, empezando por el estado inicial de 0 pivotes.
     *
     * Es la instrumentacion de la propiedad *anytime*: `traza_y[j]` es
     * exactamente el certificado que se obtendria deteniendo el bucle en el
     * pivote `j`, asi que la curva presupuesto -> cota sale de UNA sola
     * resolucion, la misma que llega al optimo, y por construccion el costo con
     * presupuesto `b` nunca supera al del optimo. Cuesta un BTRAN por pivote:
     * es para medir, no para produccion.
     */
    std::vector<std::vector<double> >* traza_y;

    /**
     * \brief Techo para el corte anytime por BRECHA. `NaN` lo desactiva.
     *
     * Sólo tiene sentido con `pivoteo_dual` y con `k` BASICA: ahi la base es
     * dual-factible en cada iteracion, asi que el valor basico de `z_k` es el
     * objetivo dual —una cota VALIDA— y sube monotonamente hacia el optimo del
     * LP. El contractor pone aca el techo que da el hull de puntos primales,
     * corrido por la tolerancia: como todo `x*` visto es factible, el optimo no
     * puede pasar de ese techo, asi que cuando la cota flotante lo alcanza lo
     * que queda por ganar ya esta por debajo de la tolerancia y seguir
     * pivoteando no puede aportar.
     *
     * Cortar no afecta la CORRECCION: la cota que se aplica se certifica en
     * intervalos desde el `lambda` que haya, como siempre.
     */
    double tope_anytime;
    static long n_vueltas;   //!< vueltas del bucle compuesto, sin reloj
    static long n_anytime;   //!< resoluciones cortadas por ese criterio
    static long n_hereda;    //!< cargas que heredaron las cotas de apoyo del nodo anterior
    static long n_base_rechazada; //!< bases heredadas descartadas por tableau inconsistente

    /** \brief OBSERVA en que pivote la cota flotante cruza `obs_tope`, sin
     *  cortar. `NaN` lo desactiva. Deja el resultado en `obs_pivote`, o -1 si
     *  no cruzo. Sirve para medir el techo de un corte antes de implementarlo. */
    double obs_tope;
    int    obs_pivote;

    /** \brief Multiplicadores `y` de la solucion, para certificar. */
    const std::vector<double>& lambda() const { return y; }

    /**
     * \brief Punto primal de la base actual, en coordenadas ORIGINALES.
     *
     * `out[c] = x_B(fila)` si la columna `c` es basica, y la cota en la que se
     * apoya si es no basica. Es lo que necesita la heuristica de Achterberg
     * para elegir la proxima cota: la que el vertice actual ya tiene mas cerca.
     *
     * \param z  la misma caja con la que se resolvio (dimension `na`).
     */
    void primal_solution(const IntervalVector& z, std::vector<double>& out);

    /**
     * \brief Infactibilidades duales de la base ACTUAL para el objetivo `+-e_k`.
     *
     * Si `k` es basica en la fila `r`, los costos reducidos del nuevo objetivo
     * son `d_c = -ck * abar_rc` sobre las no basicas, y se cuentan las que
     * tienen el signo violado por su ubicacion. Es el numero de columnas que
     * el simplex primal tendria que hacer entrar: una estimacion directa del
     * costo de resolver esa cota desde aqui, `O(nnz(fila r))`, leida del
     * tableau explicito sin pivotear. Sirve para ORDENAR la seleccion de la
     * proxima cota (`DFBH_ORDEN=dual`).
     */
    int dual_infeasibilities(int k, bool maximize) const;

    /**
     * \brief Infactibilidades PRIMALES que dejaria reubicar, para `+-e_k`.
     *
     * Es el modelo de costo del camino DUAL, en contraposicion a
     * `dual_infeasibilities`, que lo es del primal. Un paso primal mete una
     * columna a la base por pivote, asi que al primal lo predice el numero de
     * columnas con el signo violado. Un paso dual repara UNA FILA por pivote, y
     * esas columnas el dual no las pivotea: las REUBICA gratis. Lo que le
     * cuesta es la infactibilidad primal que esos saltos producen.
     *
     * Se calcula `dx_B = -sum_{c en F} T(.,c) * dx_c` sobre el tableau
     * explicito, con `F` el conjunto que se reubicaria, y se cuentan las filas
     * que quedan fuera de sus cotas. Cuesta `O(nnz(tableau))`, o sea un pivote.
     *
     * \param z  la caja sin escalar, dimension `na`.
     */
    int primal_infeasibilities(int k, bool maximize, const IntervalVector& z) const;

    /** \brief Si la columna `c` esta en la base actual. */
    bool es_basica(int c) const { return c >= 0 && c < (int)rowof.size() && rowof[c] >= 0; }

    /** \brief Fila 1..m en que la columna `c` es basica, o -1. */
    int fila_de(int c) const { return (c >= 0 && c < (int)rowof.size()) ? rowof[c] : -1; }

    /** \brief Fila `r` de `B^-1` como candidato a multiplicador, desescalada.
     *  El test `0 not-in gamma.z` es valido para CUALQUIER `lambda`, no solo
     *  para el rayo que entrega el simplex: cada fila de la inversa de la base
     *  es un certificado candidato y cuesta un BTRAN, ningun pivote. */
    bool rayo_de_fila(int r, std::vector<double>& out);


    /** \brief Costos reducidos, de donde sale `gamma = e_k - d`. */
    const std::vector<double>& reduced_costs() const { return d; }

    /** \brief Cota flotante obtenida (solo diagnostico; la que vale se
     *  certifica en intervalos a partir de `lambda()`). */
    double bound() const { return bnd; }

    int iterations() const { return iters; }

    /* Contadores globales de diagnostico. */
    static long n_solves;
    static long n_pivots;
    static long n_optimal;
    static long n_infeasible;
    static long n_iterlimit;
    static long n_warm;      //!< resoluciones que reusaron la base
    static long n_cold;      //!< resoluciones que arrancaron en frio
    static long n_fact_fallo; //!< veces que la factorizacion no sirvio
    static long n_refact_dif; //!< refactorizaciones hechas por la LU diferida
    static bool lu_diferida;  //!< `pivot` no actualiza la LU: `diferida` o `tableau`, ver `fact_sucia`
    /**
     * \brief `DFB_SX_LU=tableau` (omision): `lambda`, `d` y `x_B` salen del TABLEAU, y
     * la LU se usa solo cuando una guardia detecta deriva.
     *
     * La fila objetivo lleva `d` en la parte A y `-y` en la parte lambda, y
     * en un tableau sin deriva cumplen `d = c - Abar^T y`. La guardia rehace
     * ese producto —`O(nnz)`, lo mismo que el refresco desde la LU— y exige
     * residuo chico y ningun veredicto de signo distinto sobre las no basicas.
     * Para `x_B` exige residuo relativo chico en `Abar z = 0`. Si falla, se
     * refactoriza y se sigue por el camino de la LU. Medido con
     * `DFB_SX_DERIVA`: la deriva se concentra en pocas instancias mal
     * condicionadas, y la guardia existe para ellas.
     *
     * `DFB_SX_GUARDIA=tol` fija la tolerancia relativa, omision `1e-9`.
     */
    static bool lu_tableau;
    static long n_guardia_ok, n_guardia_fallo_d, n_guardia_fallo_xb, n_guardia_fallo_fin;
    /**
     * \brief Sonda de deriva del tableau (`DFB_SX_DERIVA=1`). No cambia la
     * trayectoria: se sigue usando la LU, y en los mismos puntos se calcula lo
     * que daria el tableau solo.
     *
     *  - al terminar cada LP, `y_tab` = `lambda` leido de la fila objetivo del
     *    tableau (`-ylam`), desescalado; el contractor lo certifica al lado del
     *    de la LU sobre la misma caja;
     *  - en cada refresco sin reubicar, `d` y `x_B` del tableau
     *    (`set_objective` + `compute_basics`) contra los de la LU: signos de
     *    costo reducido y factibilidad de las basicas que cambian de veredicto.
     */
    static bool sonda_deriva;
    std::vector<double> y_tab;
    static long dv_refrescos, dv_d_signo, dv_d_nobas, dv_xb_fact, dv_xb_bas, dv_refr_dist;
    static double dv_d_max, dv_xb_max;
    static long n_base_heredada; //!< cargas que conservaron la base del nodo anterior
    /** \brief Cuantas veces `z_k` queda basica o no basica en el optimo. Si
     *  queda NO basica, `c_B = 0` y por lo tanto `y = 0`: no hay cota que
     *  certificar. Es la hipotesis del §7.7 de MEDICIONES_PODA sobre por que el
     *  warm start pierde contraccion. */
    static long n_k_basica;
    static long n_k_no_basica;
    /** \brief Pivotes DEGENERADOS: la columna entrante tiene costo reducido
     *  ~0, asi que el paso dual no mejora el objetivo. Es la firma del
     *  estancamiento por degeneracion, y el sospechoso de que en el arbol haga
     *  falta 2.5-4.5 veces mas pivotes que SoPlex. */
    static long n_basics;   //!< llamadas a compute_basics
    static long n_piv_degen;
    /** \brief Cambios de cota del ratio test de paso largo. */
    static long n_flips;
    /** \brief Perturbacion de costos al reubicar (`DFB_SX_PERTURB=eps`). */
    static long n_pert_lps;          //!< LPs que arrancaron con costos perturbados
    static long n_pert_limpieza;     //!< pivotes gastados tras quitar la perturbacion
    static long n_piv_degen_pert;    //!< pivotes cuya entrante tenia |d| <= 3 eps: degenerados enmascarados
    /** \brief Verificacion del optimo (`DFB_SX_VERIF=1`): x_B y d exactos desde la factorizacion. */
    /**
     * \brief Estrategia de RE-OPTIMIZACION al cambiar de objetivo.
     *
     * OJO CON EL NOMBRE: el valor describe las cotas `1 .. 2n-1`, no la cota 0.
     * **La cota 0 de cada nodo es DUAL en las dos estrategias** y no depende de
     * este interruptor: `load` reinicia la base, se arranca en frio con
     * `B = -I` —las `m` columnas de `b`—, que con `c_B = 0` es dual-factible
     * sin fase 1 y primal-infactible, o sea pasos duales. Recien de la cota 1
     * en adelante hay una base previa que re-optimizar, y ahi el interruptor
     * decide como.
     *
     * `DFB_SX_PIVOTEO=primal` (omision) no reubica las no basicas: la base
     * conserva la factibilidad PRIMAL heredada —la caja no se movio— y pierde
     * la dual, asi que se re-optimiza con pasos primales. La estrategia
     * completa es entonces «dual en la cota 0, primal en el resto».
     *
     * `DFB_SX_PIVOTEO=dual` las reubica por el signo de su costo reducido, lo
     * que restaura la factibilidad DUAL y rompe la primal, asi que se
     * re-optimiza con pasos duales: dual en todas. Es la unica que mantiene una
     * cota certificada monotona durante todo el LP.
     *
     * Es escribible para poder comparar las dos sobre el MISMO nodo
     * (`DFBH_ABDUAL`).
     */
    static bool pivoteo_dual;
    /**
     * \brief `DFB_SX_PIVOTEO=dual`: el dual PURO. Reubica todas las no
     * basicas, `k` incluida, asi que la base es dual-factible en cada
     * iteracion y solo hay pasos duales. `pivoteo_dual` vale true tambien
     * para `mixta`, la variante anterior que excluye `k` de la reubicacion y
     * termina en primal en los LPs que empiezan con `k` no basica (medido con
     * `DFB_SX_MONO`: 38-83 % de los LPs con algun paso primal).
     */
    static bool dual_puro;
    static long n_entrak;   //!< pivotes de entrada explicita de `k` (`DFB_SX_ENTRAK=1`)
    static long n_verif;             //!< optimos verificados
    static long n_verif_pinf;        //!< de esos, con violacion primal relativa > 1e-6
    static long n_verif_dinf;        //!< de esos, con algun costo reducido de signo violado > 1e-7
    static void reset_counters();

private:
    int m, nx, na, nc;

    /* Relajacion ESCALADA, dispersa por filas. */
    std::vector<std::vector<std::pair<int,double> > > arow;

    /**
     * \brief Escalas de filas y columnas (equilibrado en norma infinito).
     *
     * Sin escalar, el simplex **no es robusto**: con el regimen de Brown-*
     * —coeficientes de hasta 1e72 sobre cajas de radio 1e9— falla contra SoPlex
     * en 217 de 400 casos (§8 de MEDICIONES_PODA.md). No es un problema de
     * solidez —la cota se recertifica en intervalos desde lambda, asi que un
     * lambda malo da una cota debil y nunca falsa— pero devuelve optimos
     * equivocados.
     *
     * Se resuelve sobre `A_s = R A C` con `z = C z_s`. Dos detalles:
     *
     *  - las columnas de `b` se escalan con `C[nx+j] = 1/R[j]`, de modo que la
     *    parte `b` de `A_s` sigue siendo **exactamente -I** y la base inicial
     *    sigue siendo trivial y dual-factible;
     *  - las escalas se redondean a **potencias de 2**, asi que multiplicar y
     *    dividir por ellas es exacto y el escalado no agrega error propio.
     *
     * La vuelta atras es `y = R y_s` (los multiplicadores originales, que es lo
     * que necesita la certificacion) y `z_k = C[k] * z_s,k`.
     */
    std::vector<double> rscale, cscale;
    IntervalVector      zs;
    void calcular_escalas(const Matrix& Abar);

    /* Tableau: fila 0 densa, filas 1..m dispersas. */
    std::vector<double> d;      //!< parte A de la fila objetivo (na)
    std::vector<double> ylam;   //!< parte lambda de la fila objetivo (m) = -y
    std::vector<std::vector<int> >    tidx;
    std::vector<std::vector<double> > tval;

    /* Base y ubicacion de las no basicas. */
    std::vector<int>  basic;     //!< columna basica de cada fila 1..m
    std::vector<int>  rowof;     //!< fila de cada columna basica, o -1
    std::vector<char> atupper;   //!< no basica en su cota superior

    /* Perturbacion de costos (ver `perturbar_costos`). */
    bool     perturbado;   //!< `d` lleva perturbacion en esta resolucion
    double   pert_eps;     //!< magnitud base de la perturbacion vigente
    unsigned pert_estado;  //!< generador congruencial propio: determinista
    void perturbar_costos();

    /**
     * \brief Pesos de referencia DEVEX para la eleccion de la fila que sale.
     *
     * La regla actual es Dantzig: sale la basica que mas viola sus cotas. En
     * cajas anchas da lo mismo, pero medido DENTRO del arbol —cajas apretadas y
     * degeneradas— el simplex necesita 2.5 a 4.5 veces mas pivotes que SoPlex
     * con un costo por pivote comparable, que es la firma conocida de Dantzig
     * contra devex/steepest-edge.
     *
     * Devex elige la fila que maximiza `violacion^2 / w_r`, o sea que normaliza
     * la violacion por la norma (aproximada) de la fila en el espacio de
     * referencia, y actualiza los pesos en cada pivote con la columna entrante
     * ya transformada, que `pivot()` recorre de todos modos.
     */
    std::vector<double> dw;

    /**
     * \brief Pesos EXACTOS de dual steepest edge (`DFB_SX_DSE=1`).
     *
     * `dse[r] = ||e_r^T B^-1||^2`, la norma de la fila `r` de la inversa de la
     * base. Con el tableau explicito esa fila es la parte lambda de la fila `r`
     * del tableau, asi que no hace falta la FTRAN de la actualizacion de
     * Forrest-Goldfarb: `pivot` recalcula la norma de cada fila que modifica.
     * La fila que sale maximiza `viol_r^2 / dse[r]`. Reemplaza a devex en la
     * seleccion cuando esta encendido.
     */
    std::vector<double> dse;
    double norma_lambda(int r) const;

    /** \brief Candidata del ratio test de paso largo. */
    struct Cand { int c; double a; double ratio; double rango; };
    static bool cand_menor(const Cand& a, const Cand& b);
    std::vector<Cand> bfrt_cand;
    std::vector<int>  bfrt_flip;

    /* Buffers. */
    std::vector<double> xN, xB, colq, vbuf, xhat;
    std::vector<int>    bidx;
    std::vector<double> bval;
    std::vector<double> y;
    double bnd;
    int iters;
    bool  base_lista;   //!< hay una base reusable cargada

    /** \brief Factorizacion LU de la base, para obtener `lambda` sin deriva.
     *  Se mantiene en paralelo al tableau: `change_basis` por pivote (tipo
     *  Forrest-Tomlin) y refactorizacion periodica. */
    DFBBasis fact;
    bool fact_ok;

    /**
     * \brief Politica de la LU (`DFB_SX_LU=actualizada|diferida`).
     *
     * `actualizada` sigue cada pivote con un update Forrest-Tomlin,
     * de modo que la LU esta siempre al dia. Pero sus unicos lectores en el
     * camino por omision son el `lambda` final del LP y el refresco al cambiar
     * de objetivo, que ocurren una vez por LP, mientras que el update se paga
     * en cada pivote: el tableau explicito ya lleva la base.
     *
     * `diferida` no toca la LU en `pivot`: la marca sucia, y `asegurar_fact`
     * refactoriza desde `basic[]` recien cuando alguien la lee. Con base
     * compartida el `lambda` final y el refresco de la cota siguiente leen la
     * misma base, asi que es una factorizacion por LP que pivoteo. El
     * interruptor, `lu_diferida`, es publico.
     */
    bool fact_sucia;
    bool asegurar_fact();
    std::vector<double> lam_buf;
    std::vector<int>    cols_buf;

    /** \param reubicar  si es false, se conserva la ubicacion de las no
     *  basicas. Conservarla mantiene la factibilidad PRIMAL (que es lo que
     *  hace util el warm start) a costa de perder la dual, que es lo que
     *  despues arreglan los pasos primales. */
    /**
     * \brief `y = B^-T c_B` calculado desde la FACTORIZACION, no del tableau.
     *
     * Es el arreglo del problema del §7.6 de MEDICIONES_PODA: manteniendo el
     * tableau por operaciones de fila acumuladas, `gamma_k` —que tiene que
     * valer exactamente +-1 porque la columna k es basica— deriva hasta un
     * 2.5 % entre cotas, y como la certificacion usa `lambda` y no la cota
     * flotante, la cota certificada queda debil (valida, nunca falsa: los 32/32
     * testigos lo confirman).
     *
     * Con `c = +-e_k` la cuenta es **un solo BTRAN**: si `k` es basica en la
     * fila `r`, `c_B = c_k*e_r` y por lo tanto `y = c_k * (fila r de B^-1)`; si
     * `k` es no basica, `c_B = 0` y `y = 0`.
     */
    bool y_desde_factorizacion(int k, bool maximize);

    /** \brief Camino `DFB_SX_LU=tableau`, ver `lu_tableau`. Devuelven false si
     *  la guardia detecta deriva, y el llamador cae al camino de la LU. */
    bool fila_objetivo_consistente(int k, bool maximize);
    bool basicas_consistentes();
    bool y_desde_tableau(int k, bool maximize);
    bool refrescar_desde_tableau(int k, bool maximize, const IntervalVector& z, bool reubicar);
    /** \brief El tableau con guardia si corresponde, y si no la LU. */
    bool refrescar(int k, bool maximize, const IntervalVector& z, bool reubicar);
    /** \brief Reubicacion minima de las no basicas por el signo de `d`, sin
     *  mover `k`. Devuelve cuantas movio, y en `n_k` si una fue `k`. */
    int  reubicar_minima(int k, int& n_k);
    std::vector<double> dchk;

    void registrar_traza(int k, bool maximize);

    void set_objective(int k, bool maximize, bool reubicar);

    /**
     * \brief Fila objetivo y valores basicos calculados desde la
     * FACTORIZACION, sin pasar por el tableau.
     *
     * El tableau se mantiene por operaciones de fila acumuladas y deriva. Eso
     * no importaria si solo decidiera pivotes, pero `d` es tambien el **test de
     * optimalidad** y `x_B` el de factibilidad: con los dos derivados el bucle
     * corta antes del optimo, y entonces `lambda` —que sale exacto de la
     * factorizacion— certifica una cota valida pero **debil**. Es lo que hacia
     * que el warm start entregara el costo y no la poda.
     *
     *     y   = c_k * (fila r de B^-1)        un BTRAN
     *     d   = c - Abar^T y                  un producto disperso
     *     x_B = -B^-1 (Abar x_N)              un producto + un FTRAN
     */
    /** \param reubicar  true solo al cambiar de objetivo. Reubicar en cada
     *  iteracion fue un error: es ambiguo cuando `d_i = 0` —el caso degenerado,
     *  que es omnipresente— y la ubicacion de una degenerada no cambia el
     *  objetivo pero SI el test de factibilidad primal, con lo que el criterio
     *  de optimalidad queda mal planteado (grupo B de 3 a 30). Un simplex dual
     *  de verdad conserva la ubicacion como ESTADO: la que sale va a la cota que
     *  violaba. Y cuando hay que reubicar se mueven solo las columnas cuyo signo
     *  esta violado, dejando las degeneradas donde estan, que es lo que preserva
     *  la factibilidad primal del warm start. */
    bool refrescar_desde_factorizacion(int k, bool maximize,
                                       const IntervalVector& z, bool reubicar);

    /** \brief Columna `q` del tableau transformado, `B^-1 a_q`. */
    void tab_column(int q, std::vector<double>& out) const;

    /** \brief Un paso primal: entra la no basica que viola el signo de su
     *  costo reducido, sale la basica que primero toca una cota. Devuelve
     *  false si no hay ninguna que viole (optimo dual). */
    bool primal_step(const IntervalVector& z, bool& unbounded);

    /** \brief Un paso dual: sale la basica que mas viola sus cotas, entra la
     *  del menor cociente |d_j/a_rj|. Devuelve false si no hay violacion
     *  (factible primal). */
    bool dual_step(const IntervalVector& z, bool& infeasible);
    /** \brief Igual que dual_step pero SIN recalcular x_B: usa el exacto que
     *  dejo refrescar_desde_factorizacion. */
    bool dual_step_exacto(const IntervalVector& z, bool& infeasible);
    bool dual_step_comun(const IntervalVector& z, bool& infeasible);

    double entry(int r, int col) const;
    void   axpy(int r, int rp, double f, int drop_col);
    bool   pivot(int r, int col);
    void   compute_basics(const IntervalVector& z);
    /** \brief Repara la ubicacion de las no basicas con cota infinita. */
    void   reparar_cotas_infinitas(const IntervalVector& z);

    /**
     * \brief `x_B` al dia sin recalcularlo, actualizandolo tras cada pivote.
     *
     * `compute_basics` rehace `x_B` desde cero —`O(m*nnz)`— en cada iteracion
     * del paso dual, y eso es el 20-26 % del tiempo del simplex (§43). Pero
     * `x_B(i) = -sum_N t_ij x_j`, asi que si la entrante `q` se mueve `delta`
     * desde su cota, cada basica cambia en `-alpha_i * delta`, con `alpha` la
     * columna `q` transformada ANTES del pivote — que `pivot()` ya recorre para
     * hacer las operaciones de fila. Actualizar cuesta `O(m)`.
     *
     * `xb_sucio` marca cuando hay que recalcular igual: al cambiar de objetivo,
     * al refrescar desde la factorizacion, o tras cambiar cotas de no basicas.
     */
    bool xb_sucio;
    bool pivote_fallido;  //!< un pivote no se pudo aplicar: no confundir con optimo
    std::vector<double> alpha_piv;   //!< columna entrante antes del pivote
};

} /* namespace ibex */

#endif /* __IBEX_DFB_SIMPLEX_H__ */
