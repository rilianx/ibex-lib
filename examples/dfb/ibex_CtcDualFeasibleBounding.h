#ifndef __IBEX_CTC_DFB_H__
#define __IBEX_CTC_DFB_H__

#include "ibex_Ctc.h"
#include "ibex_DFBSparse.h"
#include "ibex_DFBSimplex.h"
#include <map>
#include <list>
#include <vector>

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
class CtcDFB : public Ctc {

    private:
    int max_iters;
    map<int, int> identity_rows;
    
public:

	/**
	 * \brief Create a DFB contractor for contracting bound(x) in A.x=0
	 */
	CtcDFB(int nb_var, int k, bool upper_contract=false, bool contract_all=false, int max_iters=2): Af(1,1), refA_ref(NULL), A(1,1),
    x_ref(1), upper_contract(upper_contract), contract_all(contract_all), k(k), 
    Ctc(nb_var), max_iters(max_iters), init_ok(false), iters(0),
    ratio_test_status(INCONCLUSIVE), work(1), cand_age(0), sx(NULL), sx_skip(0),
    gamma_k_target(1.0), gam(1), lam(1), refAf(NULL) {  };   

	/**
	 * \brief Contract a box.
	 */
	virtual void contract(IntervalVector& x_new) ;

    /**
     * \brief Resultado del ratio test (ver calculateAlpha).
     *
     * Distinguir estos tres casos es indispensable para la solidez del
     * contractor: solo NO_CANDIDATE autoriza a declarar la caja vacia.
     */
    enum RatioTest {
        BLOCKING_ROW,   //!< se encontro fila de bloqueo: el pivote se aplica
        NO_CANDIDATE,   //!< ninguna fila puede bloquear: paso dual no acotado
        INCONCLUSIVE    //!< indeterminado: no permite concluir nada
    };

    /** \brief Resultado de la ultima llamada a calculateAlpha. */
    RatioTest ratio_test_status;

    /** \brief false si init()/regenerateA() no pudo construir la base. */
    bool init_ok;

    /**
     * \brief Desempate tipo Harris entre candidatos del ratio test.
     *
     * Si es true, entre los candidatos cuyo cociente esta dentro de
     * harris_rel_band del minimo se elige el de mayor |A[j][i]|, que es mas
     * estable numericamente.
     *
     * ACTIVADO por defecto: medido sobre las 153 instancias, mejora la
     * contraccion en 8 y la empeora en 6, con un neto de +63 puntos
     * porcentuales, -3% de pivotes y +1% de tiempo. Ver MEDICIONES_M2_M3.md.
     * Se puede desactivar (DFB_HARRIS=0 no sirve; poner la variable a false)
     * para reproducir la comparacion.
     */
    static bool harris_tie_break;
    static double harris_rel_band;

    /**
     * \brief Si es true, la cota de gaussSeidel que no intersecta a x[k] se
     * usa como prueba de vacio. Ver la discusion en HALLAZGOS_ANALISIS.md:
     * la validez depende de que la fila gamma sea una relacion valida, lo que
     * el forzado A[0][i]=0 tras cada pivote no garantiza en intervalos.
     */
    static bool prove_empty_by_bound;

    /**
     * \brief M2: pricing fusionado, sin asignaciones por fila.
     *
     * calculateImpacts construia dos IntervalVector por fila y hacia un
     * producto punto aparte, o sea 2m asignaciones y dos pasadas por pivote.
     * Ademas reevaluaba m veces la eleccion de esquina, que depende solo de
     * gamma y no de la fila. La version fusionada la calcula una vez por
     * pivote y acumula el producto en la misma pasada, en el mismo orden de
     * indices, de modo que el resultado es identico.
     */
    static bool pricing_fused;

    /**
     * \brief M2: partial pricing con lista de candidatos.
     *
     * El pricing recorre las m filas en cada pivote y es el 54% del tiempo de
     * DFB. Con lista de candidatos se evalua solo un subconjunto de filas
     * (partial_pricing_size, por defecto ~sqrt(m)), refrescandolo cada
     * partial_pricing_refresh pivotes o cuando ninguna candidata mejora. Es la
     * tecnica estandar del simplex para el mismo problema.
     *
     * Cambia la fila elegida, asi que puede requerir mas pivotes: hay que
     * medir tiempo, no solo pivotes.
     */
    static bool partial_pricing;
    static int  partial_pricing_size;     //!< 0 = automatico (~sqrt(m))
    static int  partial_pricing_refresh;  //!< 0 = automatico (= size)

    /**
     * \brief Refresco adaptativo de la lista de candidatos.
     *
     * Un refresco cada K pivotes es ciego: gasta pasadas completas cuando la
     * lista funciona y las escatima cuando no. Con esto, la lista se invalida
     * en cuanto un pivote deja de mejorar la cota, que es la senal de que la
     * lista se quedo sin filas utiles. Es el equivalente a "si el partial
     * pricing se estanca, ampliar" del simplex.
     */
    static bool partial_pricing_adaptive;

    /* Contadores globales de diagnostico (reiniciables con reset_counters). */
    static long n_no_candidate;   //!< vaciados de caja declarados
    static long n_inconclusive;   //!< pivoteos detenidos por indeterminacion
    static long n_init_failed;    //!< inicializaciones fallidas
    static long n_empty_by_bound; //!< vaciados probados por la cota (gaussSeidel)
    static long n_regenerations;  //!< llamadas a regenerateA

    /* Tiempo acumulado por fase, en segundos. Permite saber si el cuello es el
     * pricing, el ratio test, la eliminacion o la regeneracion, que es lo que
     * decide la prioridad entre M1, M2 y M4. */
    static double t_pricing;      //!< largestImpact + calculateImpacts
    static double t_ratio_test;   //!< calculateAlpha
    static double t_pivot;        //!< A[0] += alpha*A[j] y makeColumnIdentity
    static double t_bound;        //!< gaussSeidel
    static double t_regen;        //!< regenerateA
    static double t_init;         //!< init() de los 2n contractores
    static double t_init_copy;    //!< dentro de init: la copia de refA
    static double t_init_elim;    //!< dentro de init: la eliminacion
    static long   n_inits;        //!< llamadas a init()
    static long   n_applied;      //!< contractores DFB efectivamente aplicados

    static void reset_counters() {
        n_no_candidate = 0; n_inconclusive = 0; n_init_failed = 0;
        n_empty_by_bound = 0; n_regenerations = 0;
        t_pricing = 0; t_ratio_test = 0; t_pivot = 0; t_bound = 0; t_regen = 0;
        t_init = 0; t_init_copy = 0; t_init_elim = 0; n_inits = 0; n_applied = 0;
        n_certifications = 0; t_certify = 0;
    }

    /**
     * \brief Virtual destructor to clean up resources.
     */
    virtual ~CtcDFB() {
        //cout << "[CtcDFB] Destroying instance with k=" << k << endl;
        identity_rows.clear();
        A.clear();
        //cout << "[CtcDFB] Instance destroyed successfully." << endl;
    }
    
    
    void init(IntervalMatrix& A, IntervalVector& x_ref);

    /**
     * \brief M5(a): inicializa el contractor de cota superior a partir del de
     * cota inferior de la misma variable, sin repetir la eliminacion.
     *
     * Los dos contractores de una variable k parten de la misma matriz de
     * referencia y difieren solo en que el de cota superior niega antes la
     * columna k. Desarrollando makeColumnIdentity con la columna k negada:
     *
     *   - fila 0:  se divide por -R[0][k] en vez de R[0][k], asi que queda
     *              exactamente la fila 0 del caso inferior con las entradas
     *              i != k negadas (la entrada k se fuerza a 1 en ambos);
     *   - filas jj >= 1: el factor de eliminacion tambien cambia de signo, y
     *              (-factor)*(-fila 0) = factor*fila 0, de modo que quedan
     *              IDENTICAS.
     *
     * En aritmetica de intervalos la igualdad es exacta: negar es simetrico
     * bajo redondeo dirigido, luego -(a/b) y (-a)/b dan el mismo intervalo.
     *
     * Resultado: una eliminacion por variable en vez de dos. Los resultados
     * deben ser identicos bit a bit; DFB_NO_PAIR_INIT=1 vuelve al camino
     * anterior para comprobarlo.
     */
    void init_from_lower(const CtcDFB& lower, IntervalVector& x_ref);

    static bool pair_init;

    /**
     * \brief Pivoteo en punto flotante con certificacion de la cota final.
     *
     * La cota es valida para CUALQUIER eleccion de multiplicadores lambda
     * (dualidad debil), asi que buscarlos es una heuristica sin requisitos de
     * rigor. Se pivotea en double y, al usar la cota, se recalcula
     * gamma = lambda^T . refA en intervalos y se evalua sobre la caja. Un
     * lambda malo da una cota debil, nunca una falsa.
     *
     * Para eso el tableau flotante se mantiene AUMENTADO con la identidad:
     * Af = [ A | I ], y la parte aumentada de la fila j es el lambda de esa
     * fila (invariante: lambda_j^T . refA_k == parte A de la fila j, con refA_k
     * = refA con la columna k negada cuando se contrae la cota superior). El
     * forzado de entradas a 0 exactos se aplica solo a la parte A: lambda
     * queda siempre veraz, y por eso la certificacion absorbe cualquier
     * suciedad numerica del pivoteo.
     *
     * Consecuencia: desaparece el crecimiento de ancho, y con el get_Aerror()
     * y regenerateA().
     *
     * ACTIVADO por defecto. Medido sobre las 153 instancias: tiempo del banco
     * 15.62 s -> 5.65 s (-64%), regeneraciones 731 -> 0, grupo A 57 -> 64, las
     * mismas 2 cajas vaciadas (correctas). En ibexopt: acid_dfb pasa de 86 a 96
     * instancias resueltas y de x6.47 a x3.44 el tiempo de la configuracion de
     * produccion, sin ningun optimo incompatible.
     * DFB_INTERVAL_PIVOT=1 vuelve al pivoteo en intervalos.
     */
    static bool float_pivoting;

    /**
     * \brief Modo liviano: sin tableau. Solo gamma y lambda.
     *
     * El tableau denso aumentado existe para mantener las filas transformadas
     * B^-1.refA, y es lo que vuelve CARO cambiar de objetivo: renormalizar la
     * columna k cuesta una eliminacion O(m.(n+2m)). Esa carestia es la que hace
     * incompatibles la cadena de cotas y el intercalado fino con HC4
     * (MEDICIONES_TECHOS.md §9), y lo que una factorizacion vendria a resolver.
     *
     * Este modo lo resuelve por otra via: **no transformar nada**. Las
     * direcciones de pivoteo son las filas ORIGINALES de refA —compartidas y de
     * solo lectura entre los 2n contractores— y el estado de cada cota es solo
     *
     *   gam : la combinacion vigente, gam = lam^T . refA        (n+m doubles)
     *   lam : los multiplicadores                                (m doubles)
     *
     * Consecuencias:
     *  - cambiar de objetivo es O(n+m): se elige otra fila de arranque;
     *  - la memoria por contractor pasa de O(m.(n+2m)) a O(n+m), unas 100 veces
     *    menos con n=m=100, lo que ataca las 6 instancias que no terminan;
     *  - gamma_k NO se fuerza a 1: la cota divide por gamma_k, que sirve con
     *    cualquier valor no nulo.
     *
     * Es un algoritmo distinto: ascenso por coordenadas en el dual usando las
     * filas originales, en vez de simplex sobre un tableau transformado. Puede
     * dar cotas mas debiles por pivote; la certificacion garantiza que sean
     * validas de todos modos.
     */
    static bool light_mode;

    Vector gam;                  //!< combinacion vigente (n+m)
    Vector lam;                  //!< multiplicadores (m)
    const Matrix* refAf;         //!< filas originales, compartidas y de solo lectura

    void init_light(const Matrix& R, IntervalMatrix& A, IntervalVector& x_ref);
    void set_bound_light(int k_new);
    void contract_light(IntervalVector& x_new);
    bool certify_gamma_light(IntervalVector& gamma_out);

    /**
     * \brief M2(b): normalizacion del pricing.
     *
     * 0 = ninguna (regla tipo Dantzig: se elige la fila de mayor mejora
     *     estimada, sin mirar cuanto se puede avanzar en esa direccion);
     * 1 = se divide por la norma 2 de la fila.
     *
     * El fundamento es el mismo del steepest edge: el paso alpha que permite el
     * ratio test es aproximadamente inversamente proporcional a la magnitud de
     * la fila, asi que delta/||A_j|| aproxima la mejora REAL (delta*alpha) en
     * vez de la mejora por unidad de paso. En el camino flotante sale casi
     * gratis porque el pricing ya recorre todas las columnas de cada fila.
     *
     * Ojo: con normalizacion delta queda escalado y no sirve como criterio de
     * parada; la normalizacion cambia solo el ORDEN de eleccion y la parada se
     * evalua sobre la mejora sin normalizar de la fila elegida.
     *
     * Por defecto 1 (aceptado). Medido sobre las 153 instancias: contraccion
     * media 24.03% -> 30.03%, grupo A 64 -> 72, suma de deltas +870 puntos
     * porcentuales (ganan +901, pierden -31), a cambio de +29% de pivotes y
     * +6% de tiempo. DFB_PRICING_NORM=0 vuelve a la regla tipo Dantzig.
     */
    static int pricing_norm;

    /** \brief DFB_NO_REVIVE=1: real_impact devuelve 0 y la propagacion HC4 no
     *  reactiva contractores DFB. Reproduce el comportamiento que tenia el
     *  camino flotante por el error de leer A[0] en vez de Af[0]. */
    static bool no_revive;

    /** \brief Tableau flotante aumentado [A | I], m x (nb_var+m + m). */
    Matrix Af;

    /** \brief Indices con lambda distinto de cero (coste de certificar). */
    std::vector<int> lam_nz;

    /** \brief Matriz de referencia, para certificar. No se copia. */
    const IntervalMatrix* refA_ref;

    static long   n_certifications;
    static double t_certify;

    /** \brief Dispersion del pricing (para decidir si la factorizacion puede
     *  abaratarlo). Se acumulan en largest_impact_f, una vez por pivote:
     *  - g_total: columnas de la parte A recorridas
     *  - g_nz   : de esas, cuantas tienen |gamma_i| >= 1e-5, o sea cuantas
     *             tienen la esquina determinada por la fila 0 y no por la fila j
     *  - a_nz   : entradas no nulas del tableau transformado recorridas
     *  Si g_nz/g_total es chico, la esquina depende de la fila en casi todas las
     *  columnas y el pricing necesita las m filas transformadas completas. */
    static long   n_price_gtotal;
    static long   n_price_gnz;
    static long   n_price_anz;
    static long   n_sx_fallback;  //!< veces que el simplex no llego al optimo
    /** \brief Densidad de lambda en la certificacion. Importa porque la cota
     *  certificada **no depende solo del optimo del LP**: entre los lambda
     *  optimos, los mas dispersos certifican mas apretado, porque
     *  `gamma = lambda^T refA` se evalua en INTERVALOS y cada termino agrega
     *  ancho. */
    static long   n_lam_nz;
    static long   n_lam_cert;
    static long   n_hib_barato;   //!< cotas resueltas solo con la regla de impacto
    static long   n_hib_simplex;  //!< cotas que ademas pagaron el simplex
    static long   n_hib_memo;     //!< cotas que se saltearon por memoria

    /** \brief El simplex ya se probo en esta caja para esta cota y no dio nada.
     *
     * Separa las dos razones por las que la regla de impacto no contrae:
     * **se trabo** (grupo B, y el simplex lo arregla) o **no hay nada que
     * sacar** (grupo C, 50 de 147 instancias, donde PolyHull tampoco contrae).
     * A priori no se distinguen, pero a posteriori si: si el simplex ya corrio
     * y no contrajo, no va a contraer mientras la caja no cambie.
     *
     * Y la informacion util es **entre nodos**: si para la variable `k` la
     * relajacion no da nada en un nodo, es probable que tampoco de en los
     * vecinos. Por eso NO se limpia en el init —dentro de una caja cada
     * contractor corre el simplex una sola vez, asi que limpiarlo lo volvia
     * inutil— sino que se programa un salteo de los proximos
     * `DFB_SX_MEMO_N` nodos. */
    int sx_skip;

    /** \brief Filas del tableau que un pivote modifica realmente.
     *  make_column_identity_f resta f*fila_pivote de toda fila con f != 0, o
     *  sea de nnz(columna entrante) filas. Si ese numero es chico, casi todas
     *  las filas del tableau siguen siendo las originales de refA y
     *  rematerializarlas es trabajo perdido. */
    static long   n_elim_rows;
    static long   n_elim_calls;

    void init_float(IntervalMatrix& A, IntervalVector& x_ref);

    /**
     * \brief Camino de tableau disperso (DFB_SPARSE=1).
     *
     * Mismo algoritmo, misma secuencia de pivotes y **mismos resultados bit a
     * bit** que el camino flotante denso: solo cambia la representacion del
     * tableau, que pasa de una matriz `m x (na+m)` a filas dispersas con la
     * fila 0 densa. La justificacion esta medida en
     * MEDICIONES_FACTORIZACION.md: en las instancias con `n >= 40` el pricing
     * recorre un 2.3 % de entradas no nulas y un pivote modifica el 6.9 % de
     * las filas, o sea que el tableau denso hace 24-43x de aritmetica sobre
     * ceros.
     *
     * La equivalencia bit a bit la comprueba test_dfbsparse.cpp (132 277
     * comparaciones, diferencia 0).
     */
    static bool sparse_tableau;

    /** \brief Tableau disperso, usado cuando sparse_tableau es true. */
    DFBSparseTableau St;

    /** \brief DFB_SPARSE_CHECK=1: mantiene ademas el tableau denso y compara,
     *  en cada pivote, el tableau completo y las decisiones de pricing y ratio
     *  test. Sirve para localizar cualquier divergencia entre las dos
     *  representaciones; es carisimo, solo para diagnostico. */
    /** \brief Umbral de mejora estimada por debajo del cual el pricing declara
     *  que no hay fila que mejore y DFB se detiene (`DFB_STOP_TOL`).
     *
     *  Estaba fijo en 1e-6. Importa porque en el grupo B —las instancias donde
     *  DFB contrae menos que PolytopeHull— DFB se detiene SIEMPRE por aca
     *  (`n_nocand = n_incon = 0`) y no por el ratio test, y subir el tope de
     *  pivotes de 5 a 200 no cambia ni un pivote. Bajar el umbral separa dos
     *  causas: si el grupo B se encoge, era el umbral; si no se mueve, DFB se
     *  traba de verdad, porque su ascenso esta restringido a `m` direcciones
     *  por paso y el ascenso por coordenadas sobre una funcion concava lineal a
     *  trozos se detiene en puntos no optimos. */
    static double stop_tol;

    static bool sparse_check;
    static long n_check_difs;
    static int  chk_last_j;
    static int  chk_last_i;
    static double check_worst;

    void compare_tableaus(const char* donde);
    bool check_ready() const {
        return sparse_check && Af.nb_rows() == St.nb_rows()
                            && Af.nb_cols() == St.nb_cols();
    }

    /**
     * \brief DFB_DUAL=1: seleccion de fila por violacion de cota, que es el
     * criterio del simplex dual, en vez del "impacto" por evaluacion de
     * intervalos de la fila.
     *
     * El pricing actual calcula, para la fila j,
     *
     *     acc = sum_i a_ji * esquina(i),   esquina segun signo de gamma_i
     *                                      o, si |gamma_i| < 1e-5, de a_ji
     *
     * y eso es *casi* la violacion de cota de la variable basica de la fila j.
     * Se aparta en dos puntos, y los dos son absolutos y no relativos al tamano
     * de la caja:
     *
     *  - si |a_ji| < 1e-5 la contribucion se **descarta**. Con cajas de ancho
     *    1e9 —y las hay— un coeficiente de 1e-6 aporta 1e3, asi que descartarlo
     *    puede esconder una violacion entera. Es un mecanismo de **falso
     *    negativo**, o sea de estancamiento.
     *  - las columnas con |gamma_i| < 1e-5 se evaluan con la esquina mas
     *    favorable *para cada fila*, o sea contra un punto no basico distinto
     *    en cada fila y distinto del que define la cota. Es un mecanismo de
     *    **falso positivo**, o sea de pivotes en falso.
     *
     * Este modo calcula en cambio el valor basico con la **unica** asignacion
     * consistente de las no basicas (la que dicta el signo de gamma, que es la
     * que alcanza la cota) y mide su violacion respecto de las cotas reales de
     * la variable basica; para las columnas artificiales, respecto de 0.
     * Termina cuando ninguna variable basica viola: eso si es optimalidad.
     */
    /**
     * \brief DFB_SIMPLEX=1: la cota sale de un simplex dual de variables
     * acotadas que **alcanza el optimo del LP**, en vez del ascenso por
     * coordenadas que se traba (ver MEDICIONES_PODA.md).
     *
     * El simplex esta validado contra SoPlex: 400/400 optimos coincidentes con
     * error relativo maximo 1.5e-15, invariante `gamma = y^T*Abar` en todos los
     * casos y dualidad fuerte donde hay cota que certificar, con 1.7 pivotes por
     * cota (test_dfbsimplex.cpp).
     *
     * La solidez se preserva igual que siempre: el pivoteo es flotante y lo
     * unico que se usa de el son los multiplicadores, con los que se recalcula
     * `gamma = lambda^T*refA` **en intervalos** y se evalua sobre la caja. Un
     * lambda malo da una cota debil, nunca una falsa.
     */
    static bool simplex_mode;

    /** \brief Simplex compartido por los 2n contractores: lo carga
     *  CtcDFBPropag una vez por linealizacion, porque copiar `Abar` en cada
     *  contractor costaria 2n veces la matriz. */
    DFBSimplex* sx;

    void contract_simplex(IntervalVector& x_new);

    /** \brief `gamma = lambda^T*refA` en intervalos, sin la negacion de la
     *  columna k del camino viejo: el simplex trabaja sobre `refA` tal cual. */
    bool certify_from_lambda(const std::vector<double>& lam,
                             IntervalVector& gamma_out);

    static bool dual_pricing;

    /** \brief Columna basica de cada fila, y si una columna es basica. */
    std::vector<int>  brow;
    std::vector<char> colbasic;

    /** \brief Buffers del pricing dual (asignacion no basica y violaciones). */
    std::vector<double> zN;

    /**
     * \brief Cota flotante que da la fila 0 actual, con signo tal que **mas
     * grande es mejor**.
     *
     * De `sum_i gamma_i z_i = 0` sale `z_k = (-sum_{i!=k} gamma_i z_i)/gamma_k`,
     * asi que la cota inferior de `z_k` es `-S_ub` si `gamma_k = 1` y `S_lb` si
     * `gamma_k = -1`, donde `S_ub` y `S_lb` son las sumas de los extremos de
     * `gamma_i*[z_i]`. Es un proxy en flotantes: sirve para **elegir** el mejor
     * lambda, y la cota que se usa se certifica igual en intervalos.
     */
    double bound_proxy() const;

    /** \brief Mejor lambda visto en la secuencia de pivotes, y su cota.
     *
     * La certificacion se hace una sola vez al final, asi que el algoritmo
     * tenia que ser monotono en la cota para no perder terreno. El pricing por
     * impacto lo era aproximadamente; el criterio del simplex dual no lo es
     * —se mueve por optimalidad del LP, no por mejora de la cota— y puede
     * terminar en un lambda peor que uno intermedio. Guardar el mejor vuelve el
     * metodo monotono por construccion y no le cuesta nada a ninguna regla.
     */
    std::vector<double> best_lam;
    double best_bound;
    bool   have_best;

    int  largest_impact_dual(int& j_out, double& delta_out, double& dir_out);
    void set_basic(int row, int col);
    void init_basis_bookkeeping(int na_cols, int m_rows, int k_col, int swapped);

    void init_sparse(IntervalMatrix& A, IntervalVector& x_ref);
    void init_from_lower_sparse(const CtcDFB& lower, IntervalVector& x_ref);
    void contract_sparse(IntervalVector& x_new);
    int  largest_impact_s(int& j_out, double& delta_out, double& dir_out);
    int  ratio_test_s(int j, double dir, double& alpha_out);
    bool certify_gamma_s(IntervalVector& gamma_out,
                         const std::vector<double>* lam_use = NULL);

    /**
     * \brief Deduce el tableau flotante de cota superior del de cota inferior.
     *
     * La identidad del par (ver init_from_lower) se extiende al tableau
     * AUMENTADO. Con la columna k negada, el pivote de la fila 0 pasa de
     * R[0][k] a -R[0][k], asi que la fila 0 completa —parte A y parte lambda—
     * queda negada, salvo la entrada k de la parte A que se fuerza a 1 en
     * ambos. Y las filas jj >= 1 quedan identicas, porque el factor de
     * eliminacion tambien cambia de signo:
     *   e_jj - (-R[jj][k])*(-(e_0/R[0][k])) = e_jj - R[jj][k]*e_0/R[0][k].
     *
     * Ahorra una eliminacion por variable. Los resultados deben ser identicos.
     */
    void init_from_lower_float(const CtcDFB& lower, IntervalVector& x_ref);

    /**
     * \brief Renormaliza el tableau vigente para acotar otra variable.
     *
     * Hace identidad la columna k_new en la fila 0 sobre el tableau que ya
     * viene pivoteado, sin copiar ni rearmar desde refA. Es el mecanismo con
     * el que CtcPolytopeHull resuelve las 2n cotas desde una sola base tibia:
     * cambiar de cota es cambiar el objetivo, no empezar de nuevo.
     *
     * El invariante lambda_j^T . refA = parte A de la fila j se conserva,
     * porque las operaciones de fila se aplican al tableau aumentado completo.
     */
    bool renormalize_to(int k_new);

    /**
     * \brief Cambia el objetivo al lado opuesto de la misma variable.
     *
     * Niega la fila 0 completa (parte A y parte lambda, con lo que el
     * invariante se conserva y gamma_k pasa de +1 a -1). gaussSeidel da una
     * envoltura de dos lados con cualquier signo de gamma_k, y el pricing
     * pasa a mejorar el otro lado. Es lo que en el camino de la cola hacian dos
     * contractores separados por variable.
     */
    void flip_side();

    /** \brief Valor al que se fuerza gamma_k tras cada pivote (+1 o -1). */
    double gamma_k_target;

    /**
     * \brief true si la convencion de signo para la cota superior vive en la
     * COLUMNA k de la matriz (camino de la cola), en cuyo caso hay que negar
     * x_k en la caja de trabajo.
     *
     * Con flip_side el signo vive en la FILA 0 (gamma_k = -1) y negar la caja
     * ademas es un error: se evaluaria un gamma valido contra una caja
     * transformada. La certificacion protege los multiplicadores, no la
     * evaluacion.
     */
    bool box_sign_flip() const { return upper_contract && gamma_k_target > 0; }
    void contract_float(IntervalVector& x_new);

    /** \brief gamma = lambda^T . refA en intervalos. false si no sirve. */
    bool certify_gamma(IntervalVector& gamma_out);

    int  largest_impact_f(int& j_out, double& delta_out, double& dir_out);
    int  ratio_test_f(int j, double dir, double& alpha_out);
    bool make_column_identity_f(int col, int row, bool interchange);

    std::pair<IntervalVector, IntervalVector> calculateImpacts(
        const IntervalMatrix& A, const IntervalVector& x_new, const IntervalVector& gamma);
    
    std::tuple<int, Interval, Interval> largestImpact(
        const IntervalMatrix& A, const IntervalVector& x_new, const IntervalVector& gamma);

    std::pair<Interval, int> calculateAlpha(
        const IntervalVector& Aj, const IntervalVector& gamma, const Interval& direction);

    void changeSigns(IntervalMatrix& A, IntervalVector& x_new);

    int makeColumnIdentity(IntervalMatrix& A, const int k, bool interchange, int j);

    void makeColumnsIdentity(std::list<int> columns);

    Interval gaussSeidel(IntervalVector& x, int k, IntervalVector& gamma);

    std::pair<Interval, int> getMaxValue(const IntervalVector& vector);

    void regenerateA(IntervalMatrix& Aref);

    double get_virtual_bound();

    double get_perc_impr(int iter);

    double get_Aerror();

    double real_impact(Interval& x_k, int var, double eps=0.01);

    IntervalMatrix A;
    IntervalVector x_ref;
    bool upper_contract;
    bool contract_all;
    int k;
    int iters;

    Interval virtual_x;
    list<double> perc_imprs; //porcentajes de mejora

    /* Buffers reutilizables del pricing fusionado (evitan asignar por fila). */
    std::vector<char>     col_case;   //!< 0: usar ub, 1: usar lb, 2: depende de A[j][i]
    std::vector<Interval> col_lb, col_ub;

    /**
     * \brief Buffer de trabajo con la dimension extendida (x seguido de b).
     *
     * contract() operaba redimensionando la caja del llamador (x_new.resize)
     * para agregarle las variables b y luego volviendola a su tamano. Eso
     * viola el contrato de Ctc::contract: quien llama puede conservar
     * referencias a las componentes de la caja entre llamadas (ibex::Optimizer
     * lo hace) y resize() libera el arreglo interno, con lo que esas
     * referencias quedan colgando. Valgrind lo reporta como lectura de memoria
     * liberada en Optimizer::contract_and_bound y a -O3 termina en segfault.
     */
    IntervalVector work;

    /** \brief Copia el resultado de work a la caja del llamador, sin resize. */
    void write_back(IntervalVector& x_new) const;

    /* Estado del partial pricing */
    std::vector<int> cand_rows;   //!< filas candidatas
    int cand_age;                 //!< pivotes desde el ultimo refresco

    /** \brief Impacto de la fila Aj en las dos direcciones (pricing fusionado). */
    void row_impact(const IntervalVector& Aj, int n, Interval& g_incr, Interval& g_decr) const;

    /** \brief Prepara col_case/col_lb/col_ub para el pivote actual. */
    void prepare_columns(const IntervalVector& x_new, const IntervalVector& gamma, int n);

    //enum State
    enum State {
        INITIAL,
        CONTRACTING,
        FINAL
    };
    State state;

};

} // namespace ibex
#endif // __IBEX_CTC_DFB_H__
