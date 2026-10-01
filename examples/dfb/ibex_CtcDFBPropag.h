#ifndef __IBEX_CTC_ALL_PROPAG_H__
#define __IBEX_CTC_ALL_PROPAG_H__

#include <queue>
#include <map>
#include <set>
#include <vector> // Para std::vector
#include "ibex_Ctc.h"
#include "ibex_Linearizer.h"
#include "ibex_CtcDualFeasibleBounding.h"
#include "ibex_DFBSimplex.h"
#include "ibex_DirectedHyperGraph.h"
#include "ibex_CtcFwdBwd.h"
#include "ibex_ExtendedSystem.h"
#include "ibex_LinearizerXTaylor.h"


using namespace std;

namespace ibex {

/**
 * \ingroup contractor
 *
 * \brief Propagation contractor.
 *
 * This class is an implementation of the classical interval variant of the AC3 constraint propagation
 * algorithm.
 *
 */
class CtcDFBPropag : public Ctc {
public:

    CtcDFBPropag(ExtendedSystem& sys, Linearizer& lr, double ratio=0.1, bool stand_alone=true, bool only_hc4=false, bool only_dfb=false);

    virtual ~CtcDFBPropag();

    void update_ref(const IntervalVector& box);

    void linearize(const IntervalVector& box, IntervalMatrix& A, IntervalVector& x);

    virtual void contract(IntervalVector& box);

    void init_dfb_contractors(IntervalMatrix& A, IntervalVector& x_ref);

    /** \brief Inicializa los 2n contractores DFB, en pares (ver M5a). */
    void init_all_dfb(IntervalMatrix& A, IntervalVector& x);

    double compute_rhs_ub(b_constraint& b_ctr, const IntervalVector& box);

    

    /**
     * \brief Contract a box.
     */
    virtual void contract(IntervalVector& box, ContractContext& context);


    /**
     * \brief The linearization technique
     */
    Linearizer& lr;

    /**
     * \brief  The linear solver that will be used
     */
    LPSolver mylineardummysolver;

    std::vector<CtcDFB*> dfb_ctc; // Vector de punteros a CtcDFB
    std::vector<Ctc*> hc4_ctc;    // Vector de punteros a Ctc

    double ratio;

    IntervalMatrix refA;

    /** \brief refA en doubles, compartida y de solo lectura por los 2n
     *  contractores del modo liviano (una sola copia por caja). */
    Matrix refAf;
    IntervalVector refbox;

    DirectedHyperGraph g; // constraint network (hypergraph)

    //ctc2id
    std::map<Ctc*, int> ctc2id;

    bool stand_alone;

    int count_dfb;
    int count_hc4;

    list< pair<int, double> > history; //dfb contractors

    ExtendedSystem& sys;

    map<int, list <pair<int, b_constraint*> > > adj_b; // x -> list(b_id,b_constraint)

    static bool b_contraction;

    /**
     * \brief true si la ultima linealizacion probo que el sistema lineal es
     * infactible (lr.linearize devolvio -1).
     *
     * Ese valor de retorno es parte del contrato de Linearizer: "el numero de
     * restricciones, o -1 si el sistema lineal es infactible". Antes no se
     * chequeaba y se usaba como numero de filas, con lo que resize(-1, ...)
     * intentaba reservar new Interval[-1] y abortaba con std::bad_alloc (asi
     * moria ibexopt --filtering=dfb). Ademas se perdia una prueba de vacio
     * valida: si la relajacion lineal es infactible, la caja no tiene
     * soluciones.
     */
    bool linearization_infeasible;

    /**
     * \brief Si es false, no se usa el -1 del linealizador para vaciar la caja
     * (solo se evita el resize(-1,...)). Diagnostico: DFB_NO_LIN_EMPTY=1.
     */
    static bool trust_linearizer_infeasible;

    /**
     * \brief Tope de aplicaciones de contractores por llamada a contract.
     *
     * El bucle de propagacion se realimenta segun ratiodelta, que con cajas
     * enormes (cotas de 1e100 en un optimizador) puede seguir dando cambios
     * "grandes" con progreso absoluto despreciable, y entonces contract no
     * retorna nunca. Ni el limite de tiempo de ibexopt puede actuar, porque se
     * comprueba entre nodos. Cortar antes solo contrae menos, nunca de mas.
     *
     * 0 = automatico: 50 * (numero de contractores + 1).
     */
    static int max_propag_steps;

    /**
     * \brief Warm start entre llamadas: reusar la linealizacion y las bases.
     *
     * Hoy cada llamada a contract() relinealiza y re-inicializa los 2n
     * contractores en frio. La linealizacion del padre sigue siendo VALIDA para
     * cualquier subcaja, asi que se puede reusar — pero solo si la caja actual
     * esta contenida en la que se linealizo. En un optimizador las llamadas
     * consecutivas vienen de nodos distintos y NO son anidadas (dos hermanos no
     * se contienen), asi que la contencion hay que verificarla, no suponerla.
     *
     * Si se reusa, cada contractor conserva su Af (y su lambda), que es el warm
     * start medido en MEDICIONES_INCREMENTAL.md: mediana de 0.032 pivotes
     * respecto de arrancar de cero.
     *
     * relin_shrink acota lo contrario: una relajacion calculada sobre una caja
     * mucho mas grande es demasiado floja, asi que si el perimetro cayo por
     * debajo de esa fraccion se relineariza igual.
     */
    /**
     * \brief Prioridad para la cota inferior de la variable objetivo.
     *
     * En un optimizador lo que poda es la cota inferior de la variable
     * objetivo del sistema extendido, y CtcPolytopeHull la ataca en forma
     * dirigida (optimizer() / set_contracted_vars). DFB, en cambio, reparte el
     * esfuerzo entre las 2n cotas por igual. Con esto ese contractor entra
     * primero en la cola y recibe mas iteraciones por invocacion.
     *
     * Hipotesis a medir: la brecha con la configuracion de produccion es de
     * arbol (x2.48 celdas) y no de tiempo por nodo (x1.31), asi que podria
     * venir de repartir el esfuerzo en vez de concentrarlo donde poda.
     */
    /**
     * \brief Periodo de aplicacion de la parte DFB: 1 = en cada llamada.
     *
     * La descomposicion del arbol (MEDICIONES_TECHOS.md) mostro que DFB corta
     * el arbol a la mitad respecto de ACID+HC4 y cubre el 74% de lo que logra
     * PolytopeHull, pero cuesta x2.24 de tiempo contra x1.28 de PolytopeHull:
     * el deficit es costo por unidad de poda, no poder de poda.
     *
     * Si la contraccion que consigue DFB en un nodo sirve para todo su subarbol,
     * aplicarlo cada K nodos deberia conservar buena parte de la poda a 1/K del
     * costo. Con periodo > 1 las llamadas intermedias no linealizan ni
     * inicializan los 2n contractores: solo corre la propagacion HC4.
     *
     * Es la misma idea que en ibex-ipopt, donde Ipopt rinde aplicado cada cierto
     * numero de nodos y no en cada uno.
     */
    /**
     * \brief Barrido por lotes: las cotas desde UNA base en cadena.
     *
     * Es el mecanismo de CtcPolytopeHull, validado en chain_dfb: -30% de
     * pivotes y una sola copia de matriz por caja en vez de 2n
     * (MEDICIONES_TECHOS.md §8). En vez de 2n contractores independientes en la
     * cola, se usa UN tableau: para cada variable k se renormaliza la columna k
     * en la fila 0 y se pivotea desde donde quedo k-1.
     *
     * gaussSeidel contrae las DOS cotas de x_k de una vez, asi que un barrido
     * por variable alcanza; lo que se pierde es la optimizacion del pivoteo
     * hacia el lado superior, que tenia su propio contractor.
     *
     * Se alterna con la propagacion HC4 en rondas gruesas, en vez del
     * intercalado fino de la cola. Ese intercalado rinde x1.03
     * (MEDICIONES_INTERACCION.md), asi que el intercambio deberia ser favorable.
     */
    static bool batch_sweep;
    static int  sweep_pivots;   //!< pivotes por cota en el barrido
    static int  sweep_rounds;   //!< alternancias barrido/HC4

    /** \brief Contractor unico del barrido (dueno del tableau en cadena). */
    CtcDFB* sweeper;

    /** \brief Un barrido: las n cotas desde una base en cadena. */
    void dfb_sweep(IntervalVector& box);

    /** \brief Punto fijo de la propagacion HC4, sin contractores DFB. */
    void hc4_fixpoint(IntervalVector& box, ContractContext& context);

    /** \brief Camino alternativo de contract() cuando batch_sweep esta activo. */
    void contract_sweep(IntervalVector& box, ContractContext& context);

    static long n_sweeps;

    static int dfb_period;
    long call_count;

    static bool goal_priority;
    static int  goal_iters;

    static bool   reuse_linearization;
    static double relin_shrink;

    /** \brief Caja sobre la que se calculo la linealizacion vigente. */
    IntervalVector lin_box;
    bool           has_lin;

    /* Contadores de diagnostico del warm start. */
    static long n_relin;        //!< relinealizaciones
    static long n_reuse;        //!< reusos de la linealizacion
    static long n_reuse_denied; //!< reusos rechazados por no haber contencion
    static long n_revive;       //!< reactivaciones de un DFB desde la rama HC4
    static long n_impact_calls; //!< llamadas a real_impact

    /** \brief Simplex dual compartido por los 2n contractores (DFB_SIMPLEX=1).
     *  Se carga una vez por linealizacion: copiarlo por contractor costaria 2n
     *  veces la matriz. */
    /** \brief Simplex por contractor (DFB_SX_PROPIO=1) en vez de uno compartido.
     *
     * Compartir uno solo hace inutil el warm start: cada resolucion hereda la
     * base de OTRO objetivo, que es el caso donde menos vale, y encima `d` se da
     * vuelta entero y hay que reubicar todas las no basicas. Con uno por
     * contractor, cada cota arranca tibia desde **su propia llamada anterior**:
     * mismo objetivo y la caja apenas cambiada, que es el regimen donde un
     * simplex dual reoptimiza en dos o tres pivotes. Es la incrementalidad que
     * el proyecto tenia como premisa.
     *
     * Cuesta 2n copias de la matriz; para el banco del optimizador (mediana
     * n=14) es trivial, y si rinde se comparte la matriz y se guarda solo base y
     * factorizacion por contractor. */
    DFBSimplex simplex;
    std::vector<DFBSimplex*> simplex_por_ctc;

    /** \brief Metodo perezoso: cuanto se gasta por cota y en total.
     *
     * Gracias a la certificacion, en DFB *cada pivote ya entrega una cota
     * valida*, asi que se puede repartir el esfuerzo en vez de llevar cada cota
     * a su punto fijo. Dos mandos independientes:
     *
     *  - \c lazy_pivots (DFB_MAXITERS): pivotes por visita de un contractor.
     *    Al terminar la visita, si el contractor sigue en CONTRACTING vuelve al
     *    *final* de su grupo en la cola, o sea que el reparto es un round-robin
     *    sobre las 2n cotas con HC4 intercalado. Por omision 1.
     *
     *  - \c pivot_budget (DFB_PIVOT_BUDGET): tope *total* de pivotes por caja.
     *    Al agotarse se dejan de atender contractores DFB y la propagacion
     *    termina solo con HC4. Negativo = sin tope (comportamiento historico:
     *    se corre hasta el punto fijo colectivo). Si es 0 se expresa como
     *    multiplo de nb_var con DFB_PIVOT_BUDGET_N.
     *
     * PolytopeHull no puede hacer esto: su cota aparece cuando el LP llega a
     * OptimalProved, y una corrida interrumpida no entrega nada por la API.
     */
    static int    lazy_pivots;
    static int    pivot_budget;
    static double pivot_budget_n;

    /** \brief Pivotes DFB gastados en la caja en curso (para el tope total). */
    long budget_spent;


};

} // namespace ibex
#endif // __IBEX_CTC_ALL_PROPAG_H__
