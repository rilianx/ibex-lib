//============================================================================
//                                  I B E X
// File        : ibex_CtcDFBHull.h
// Author      : DFB en el hueco exacto de CtcPolytopeHull
//============================================================================

#ifndef __IBEX_CTC_DFB_HULL_H__
#define __IBEX_CTC_DFB_HULL_H__

#include "ibex_Ctc.h"
#include "ibex_Linearizer.h"
#include "ibex_LPSolver.h"
#include "ibex_DFBSimplex.h"
#include "ibex_CtcPolytopeHull.h"

#include <vector>

namespace ibex {

/**
 * \brief Contractor de casco poliedral con DFB, sustituto directo de
 *        `CtcPolytopeHull`.
 *
 * Misma interfaz, misma estructura de `contract` y **una sola pasada** sobre las
 * `2n` cotas, igual que `CtcPolytopeHull::optimizer`: el punto fijo lo pone el
 * `CtcFixPoint` de afuera, no este contractor. Lo unico distinto es quien
 * resuelve el LP y como se obtiene la cota.
 *
 * La estrategia es la del §14 de MEDICIONES_DUELO.md, que es la que gano el
 * duelo sin propagacion (114 instancias, misma linealizacion):
 *
 *  1. **`2n` estados de simplex independientes**, uno por cota, **sembrados**
 *     copiando la base de la primera cota resuelta al optimo. Una base por cota
 *     sin sembrar degenera en arranque en frio y cuesta x2.6; compartir una sola
 *     cuesta x1.4. Las dos cosas hacen falta juntas.
 *  2. **Intercalado**: cada cota se aplica apenas sale, y la siguiente parte de
 *     la caja ya apretada. Es lo mismo que hace `CtcPolytopeHull` con su
 *     `set_bounds(i, box[i])`.
 *  3. **Cada LP hasta el optimo**, sin tope de pivotes. Truncar no paga cuando
 *     el presupuesto alcanza para terminar.
 *  4. **Salteo de Achterberg**: se da por terminada toda cota cuyo `x*_j` ya
 *     este pegado al borde, porque ese LP no puede aportar nada. Es lo que baja
 *     de 178 % a ~83 % de las `2n` cotas.
 *  5. **Orden de la mas lejana primero**: entre las que quedan, la de
 *     `|x*_j - borde|` **maximo**, no la de minimo como hace Achterberg. La
 *     cercana minimiza el costo total si se llega al final; la lejana ademas
 *     adelanta la poda, y resulto mejor en las dos dimensiones.
 *
 * La cota se certifica en intervalos desde `lambda`: `gamma = lambda^T Abar` y
 * `gamma.z = 0` para todo `z` factible, de modo que un `lambda` malo da una cota
 * debil y **nunca falsa**. Si la cota certificada no intersecta a `box[k]`, eso
 * es una prueba rigurosa de vacio.
 */
class CtcDFBHull : public Ctc {
public:
    /** \param goal_var  indice de la variable objetivo del sistema extendido,
     *  o -1 si no la hay. Lo unico que poda un nodo es la cota inferior de esa
     *  variable contra el loup; las otras `2n-1` solo mejoran la caja para
     *  bisecar. Saberlo permite atar el corte a lo que la busqueda usa. */
    /** \param sonda_ctc  contractor de referencia (HC4) que usa SOLO la sonda
     *  del §28, para comparar cuantos trozos mata cada uno. No se usa en la
     *  contraccion normal. */
    CtcDFBHull(Linearizer& lr, int goal_var=-1,
               double eps=LPSolver::default_tolerance, Ctc* sonda_ctc=NULL);

    virtual ~CtcDFBHull();

    virtual void contract(IntervalVector& box);

    virtual void contract(IntervalVector& box, ContractContext& context);

    /**
     * \brief Test de vacio BARATO con los `gamma` ya calculados, sin pivotear ni
     * re-linealizar. Devuelve true si prueba que `box` no contiene ninguna
     * solucion de la relajacion.
     *
     * Dos hechos lo habilitan: `gamma = lambda^T Abar` **no depende de la caja**,
     * y la relajacion del nodo vale para **toda sub-caja**. Entonces
     * `0 in gamma.z` es obligatorio para todo punto factible, y si falla la caja
     * esta vacia. Cuesta `2n` productos punto de `na` intervalos.
     *
     * **La validez exige que `box` este contenida en la caja para la que se
     * linealizo**, y eso se VERIFICA: este objeto se comparte en todo el
     * recorrido del arbol, asi que los ultimos `gamma` pueden venir de otra rama
     * y ahi no valdrian. Sin esa comprobacion el test seria incorrecto.
     */
    bool gamma_prueba_vacio(const IntervalVector& box) const;

    int nb_var_publico() const { return nb_var; }

    /* Contadores de diagnostico, analogos a los de CtcPolytopeHull. */
    static long n_lps;        //!< resoluciones de LP (llamadas a solve)
    static long n_pivots;
    static long n_congela;   //!< veces que se re-congelo la matriz
    static long n_ctr_evaluadas, n_ctr_total;
    static long n_diag, n_diag_vieja, n_diag_nueva, n_diag_igual;
    static double acc_diag, acc_corr_v, acc_corr_n;
    static long n_barata;    //!< vueltas resueltas sin linealizar
    static long n_m_tot, n_m_llam;
    std::vector<Interval> b_orig;    //!< cota sin refresco afin, para la sonda DFBH_CORIG
    std::vector<Interval> b_orig_cong; //!< rango capturado al congelar, punto de partida de b_orig
    Interval gn_cong(int i, const std::vector<int>& idn, const std::vector<Interval>& gn);
    static long n_orig, n_orig_mejor, n_orig_peor, n_orig_igual, n_orig_vacio;
    static double acc_orig;
    static long n_fila_tot, n_fila_vieja, n_apr, n_apr_nada;
    static double acc_apr;
    static long n_chk, n_chk_ident, n_chk_cond, n_chk_viola;  //!< sonda del acotamiento afin de c
    static double chk_peor, chk_peor_ident;
    void chequear_c(int i, const Vector& Ai, const Vector& ap, const Interval& rango_p,
                    const IntervalVector& box, const Interval& dif, const Interval& corr);   //!< filas de la relajacion, para comparar los dos caminos
    static long n_acum_gana; //!< cotas de b que la acumulacion logro apretar
    static long n_sinpar;    //!< filas congeladas sin restriccion correspondiente     //!< pivotes del simplex
    static long n_nodes;      //!< llamadas a contract con linealizacion util
    static long n_bounds;     //!< cotas disponibles (2n por nodo)
    static long n_empty;      //!< vacios probados por la cota certificada
    static long n_cortes;     //!< pasadas cortadas por el criterio temprano
    static long n_gtest;      //!< llamadas al test barato
    static long n_gvacio;     //!< vacios probados por el test barato
    static long n_gfuera;
    static void reset_counters();

private:
    Linearizer& lr;
    int nb_var;
    int goal_var;
    Ctc* sonda_ctc;

    /** \brief Solo se usa para recibir la linealizacion; no se resuelve con el. */
    LPSolver mylinearsolver;

    /* Relajacion del nodo: `Abar z = 0` con `z = (x, b)`, `na = nb_var + m`. */
    IntervalMatrix A;
    IntervalVector z;
    Matrix Af;

    /** \brief Filas de `A` que se usan en la pasada. 0 = todas.
     *  El sistema truncado es el submatriz principal: las `m1` primeras filas y
     *  las `nb_var + m1` primeras columnas, porque las componentes `b` estan
     *  ordenadas por fila. Sigue siendo una relajacion valida, sólo mas floja. */
    int m_uso;

    /** \brief La ultima pasada activo el HC4 intercalado al menos una vez.
     *  Es el disparador del ciclo externo: si HC4 contrajo, la caja cambio lo
     *  suficiente como para que una relajacion nueva valga la pena. */
    /** \brief Certificados de fila cacheados, uno por variable.
     *
     *  `gamma = lambda^T Abar` NO depende de la caja, asi que un certificado
     *  construido con una base vale para cualquier sub-caja posterior y se
     *  puede reevaluar en `O(na)` en vez de reconstruirlo en
     *  `O(nnz(lambda)*na)`. Se reinician al linealizar. */
    std::vector<IntervalVector> gfila;
    std::vector<char>           gfila_ok;

    /** \brief Traza cota por cota dentro de `una_pasada` (`DFBH_ABNODO`). */
    /** \brief Regla de la primera cota del nodo (`DFBH_PRIMERA`). Miembro para
     *  poder alternarla en la sonda A/B sobre el mismo nodo. */
    /** \brief Matriz CONGELADA y la identidad de sus filas (`DFBH_CONG=k`).
     *
     *  Con `A` congelada la cota de `b_i` se recalcula sin re-linealizar:
     *  `b_i = A_i.x = (A_i - a'_i).x + a'_i.x`, y como toda solucion cumple
     *  `a'_i.x` en el rango que da la afin nueva, sale
     *  `b_i` en `rango'_i + (A_i - a'_i).X`, con la fila diferencia evaluada
     *  por intervalos, que para una forma lineal sobre una caja es exacta. */
    bool                reusa_A;   //!< esta llamada reusa la matriz congelada
    /** \brief Cota ACUMULADA de cada `b_i`, y la caja en que se acumulo.
     *
     *  Con `A` congelada, `b_i = A_i.x` es siempre la misma cantidad, asi que
     *  una cota valida en una llamada sigue siendo valida en la siguiente
     *  mientras la caja solo se encoja: se INTERSECTA en vez de recalcularse, y
     *  se estrecha monotonamente acumulando todo lo aprendido. Con la
     *  relajacion normal esto no existe, porque al re-linealizar `b_i` pasa a
     *  significar otra cosa. */
    std::vector<Interval> b_acum;
    IntervalVector        caja_acum;
    /** \brief Ultima linealizacion vista de cada fila congelada, y los
     *  diametros de la caja en que se evaluo.
     *
     *  Permite re-linealizar SOLO las restricciones cuyas variables se movieron
     *  mas que un umbral (`DFBH_CONGVAR=r`). Las demas conservan su `a'_i` y su
     *  `rango'_i`, que siguen siendo validos porque la caja solo se encoge, y
     *  cuya correccion incluso mejora sola al evaluarse sobre una caja menor. */
    std::vector<Vector>              ult_rows;
    std::vector<Interval>            ult_rango;
    std::vector<std::vector<double> > ult_diam;
    /** \brief Caja en que se evaluo por ultima vez cada fila. Reusar
     *  `rango'_i` sin verificar CONTENCION es incorrecto: el rango vale para
     *  los `x` factibles de aquella caja, y si la de ahora no esta contenida
     *  —por ejemplo al pasar a un hermano del arbol— la cota deja de ser
     *  valida. Comparar diametros no alcanza, porque una caja puede encoger y
     *  a la vez desplazarse. */
    std::vector<IntervalVector>      ult_caja;
    std::vector<Vector> cong_rows;
    std::vector<int>    cong_ids;
    /** \brief Si la llamada anterior NO contrajo, el punto fijo termino y la
     *  proxima llamada es de un nodo nuevo: hay que re-congelar. Es la unica
     *  señal de frontera de nodo que el contractor tiene, porque el punto fijo
     *  lo llama repetidas veces sin avisar. */
    bool                nodo_nuevo;
    int  regla_primera;
    bool traza_cotas;
    /** \brief Resumen de cotas sin certificado (`DFBH_SINCOTA`). */
    bool traza_resumen;
    bool hc4_actuo;

    /** \brief Un unico simplex para las `2n` cotas, secuencial: la base de la
     *  cota anterior es el arranque tibio de la siguiente. Es la arquitectura
     *  de `CtcPolytopeHull`, y medida sobre 3 relajaciones y 4 semillas gana
     *  entre 6.8 % y 12.3 % contra `2n` ejemplares sembrados. */
    DFBSimplex sx;
    std::vector<char>       hecho;
    std::vector<double>     xstar;

    /** \brief Caja minima que contiene TODOS los puntos primales vistos desde
     *  la ultima linealizacion.
     *
     *  Todo minimizador primal `x*` es un punto factible de `P cap X`, asi que
     *  si `H` los contiene a todos, la contraccion posible de la cota (j,min)
     *  es a lo sumo `H.lb_j - X.lb_j` (simetrico para max). Es una cota
     *  SUPERIOR de la ganancia, exacta y O(n) por LP, y generaliza tanto el
     *  salteo de Achterberg como el orden de la mas lejana: ambos miran solo el
     *  punto del ultimo LP. */
    IntervalVector hull_primal;
    bool hull_vacio;
    IntervalVector          gamma;

    /** \brief Construye `A` y `z` del nodo. Devuelve el numero de filas, 0 si no
     *  hay nada que linealizar y -1 si la relajacion es infactible. */
    int linearize(const IntervalVector& box, ContractContext& context);

    /** \brief Linealizacion con la matriz CONGELADA: se queda con las filas
     *  viejas y recalcula solo el rango de `b` (ver `cong_rows`). */
    int linearize_congelada(const IntervalVector& box, int cong);

    /** \brief El camino normal: re-linealiza y reconstruye `A`. */
    int linearize_fresca(const IntervalVector& box, ContractContext& context);

    /** \brief Con DFBH_AB=1: en cada nodo se corre tambien CtcPolytopeHull
     *  sobre una COPIA de la caja y se comparan las dos contracciones. Es la
     *  unica forma de ver por que DFB poda menos DENTRO del arbol, donde la
     *  medicion aislada decia lo contrario. */
    CtcPolytopeHull* ph_ab;

    /** \brief `gamma` de cada cota, y la caja para la que se linealizo. */
    std::vector<IntervalVector> gs;
    std::vector<char> gs_ok;
    IntervalVector caja_lin;

    /**
     * \brief Cota LAGRANGIANA con monotonia para `z_k` (`DFBH_LAGRC=1|2`).
     *
     * Del certificado `gamma` del LP `min/max z_k`: `z_k = N/gamma_k` con
     * `N = -sum_{i!=k} gamma_i x_i - sum_j gamma_{b_j} b_j`, `b_j = a_j.x`.
     * En las filas ACTIVAS —las que el LP acoto por el lado de la
     * restriccion y cuyo multiplicador tiene el signo que hace
     * `(-gamma_{b_j}/gamma_k) g_c(x) >= 0` (<= 0 para la cota superior) en
     * todo factible— se sustituye `b_j = g_c(x) - rho_j(x)` y se descarta el
     * termino en `g_c`. Queda `Phi(x) = N'(x)/gamma_k <= z_k` (>= para lado 1),
     * con `rho_j = g_c - a_j.x` el resto NO lineal de la linealizacion afin.
     *
     * Lo que el LP no puede hacer: si `dPhi/dx_i` tiene signo constante en la
     * caja, el extremo de `Phi` esta en una cara, y la forma afin de `rho_j`
     * sobre la cara tiene menos error que sobre la caja. Se fijan esas
     * variables (hasta `DFBH_LAGRR` rondas, omision 2) y se evalua `N'` en
     * afin sobre la caja reducida. Todo en intervalos: `gamma` es el
     * certificado en intervalos, `a_j` son los flotantes exactos de la fila y
     * `rho_j` sale del evaluador afin del linearizador.
     *
     * Modo 1: solo la cota inferior del objetivo. Modo 2: las `2n` cotas.
     * Devuelve la semirrecta, o ALL_REALS si no aplica.
     */
    Interval cota_lagrangiana(const IntervalVector& gamma, int k, int lado);
    static long n_lagr_aplic, n_lagr_mejora, n_lagr_vacio, n_lagr_fijas;
    bool lin_valida;
    /** \brief El test se apaga solo si no esta rindiendo. */
    mutable bool g_apagado;

    /** \brief Sonda: cuanto podria podar ACID reusando estos `gamma` sin
     *  resolver ningun LP. Ver §26 de MEDICIONES_DUELO.md. */

    /** \brief Una pasada sobre las `2n` cotas, con salteo y orden. */
    void una_pasada(IntervalVector& box);
    /** \brief Pasada que reusa la relajacion y las bases del punto fijo. */
    void una_pasada_reusada(IntervalVector& box);
    bool sin_recarga;
};

/**
 * \brief Envoltorio que expone el test barato de `CtcDFBHull` como contractor,
 * para meterlo dentro de ACID.
 *
 * ACID poda rebanando, y este test contesta «¿el trozo esta vacio?» con `2n`
 * productos punto y ningun pivote. Medido sobre los trozos que ACID prueba, mata
 * el 11.9 % en `alkylbis` contra el 5.5 % de HC4, y el 10.4 % de todos los
 * trozos los mata **solo** DFB (§28). Los dos son complementarios.
 */
class CtcDFBGamma : public Ctc {
public:
    CtcDFBGamma(CtcDFBHull& hull) : Ctc(hull.nb_var_publico()), hull(hull) {}
    virtual void contract(IntervalVector& box) {
        if (hull.gamma_prueba_vacio(box)) box.set_empty();
    }
private:
    CtcDFBHull& hull;
};

} /* namespace ibex */

#endif /* __IBEX_CTC_DFB_HULL_H__ */
