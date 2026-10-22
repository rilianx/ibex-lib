//============================================================================
//                                  I B E X
// File        : ibexsolve_dfb.cpp
// Author      : (DFB) ibexsolve con DFB en el hueco de CtcPolytopeHull
//============================================================================
//
// Misma composicion que DefaultSolver (ibex_DefaultSolver.cpp):
//   HC4, ACID(HC4), Newton si el sistema es cuadrado, y
//   CtcFixPoint(CtcCompo(<contractor lineal>(LinearizerXTaylor), HC4)),
// y el mismo buffer (pila: profundidad primero). Cambia el contractor lineal:
// CtcPolytopeHull (--lr=xn, produccion) o CtcDFBHull (--lr=dfbhull), y el
// linearizador (--lin=art por omision, --lin=xt el de DefaultSolver), y el
// bisector: -b lsmear (omision desde 2026-10-06), -b smearsumrel, el de
// DefaultSolver (SmearSumRelative si es cuadrado, RoundRobin si no), o
// -b lsmeardfb, la regla de LSmear con las sensibilidades de DFB (requiere
// --lr=dfbhull). --lin=compo escribe las filas del afin y de XTaylor.
//
// Uso: ibexsolve_dfb archivo.bch [--lr=xn|dfbhull] [--lin=art|xt|compo] [-b lsmear|smearsumrel|lsmeardfb] [-e eps_min] [-t timeout]
//                                [--seed=s] [--bfs]
// Ultima linea: celdas tiempo estado soluciones frontera desconocidas pendientes

#include "ibex.h"
#include "ibex_CtcDFBHull.h"
#include "affine/ibex_LinearizerAffine2.h"
#include "ibex_LSmearDFB.h"
#include "ibex_LinearizerCompo.h"
#include <cstring>
#include <cstdlib>

using namespace std;
using namespace ibex;

/* `DFBH_ACIDTRACE=N`: cada N llamadas a ACID imprime cuantas variables recorta
 * (`nbcidvar`, el valor que ACID aprende en sus fases de ajuste). */
class AcidSonda : public CtcAcid {
public:
    AcidSonda(const System& sys, Ctc& c) : CtcAcid(sys, c) {}
    long n = 0, vac = 0;
    void contract(IntervalVector& box, ContractContext& context) {
        static const long cada = getenv("DFBH_ACIDTRACE") ? atol(getenv("DFBH_ACIDTRACE")) : 0;
        CtcAcid::contract(box, context);
        ++n; if (box.is_empty()) ++vac;
        if (cada > 0 && n % cada == 0)
            fprintf(stderr, "[acid] llamadas=%ld nbcidvar=%d vacias=%ld\n", n, nbcidvar, vac);
    }
    void contract(IntervalVector& box) { ContractContext c(box); contract(box, c); }
};

int main(int argc, char** argv) {
    /* Orden de las cotas por omision (2026-10-20): el conteo del tableau donde
     * lo hay y la mas cercana con el motor revisado. Antes el conteo caia en
     * silencio a la MAS LEJANA con el revisado, el peor orden para el warm
     * start dual en relajaciones grandes. `DFBH_ORDEN` lo cambia. */
    if (!getenv("DFBH_ORDEN")) setenv("DFBH_ORDEN", "costo_o_cerca", 0);
    /* CORRECCIONES POR OMISION (2026-10-20; cada una se apaga desde fuera):
     *  - DFB_SX_RAYO=fila: el certificado de vacio con la fila infactible de B^-1;
     *  - DFB_SX_PISO=error: piso de la tolerancia primal al error de la fila;
     *  - DFBH_FILAS_PARALELAS=paralelas: fusion de filas proporcionales (compo);
     *  - tolerancias de SC en 1e-9 sin el tope por ancho (2026-10-21: con 1e-6,
     *    la de la caja, SC perdia contraccion con el afin; lo que arreglaba
     *    BroydenBanded era quitar el tope); tope de iteraciones de la caja igual al de SC. */
    if (!getenv("DFB_SX_RAYO"))           setenv("DFB_SX_RAYO", "fila", 0);
    if (!getenv("DFB_SX_PISO"))           setenv("DFB_SX_PISO", "error", 0);
    if (!getenv("DFBH_FILAS_PARALELAS"))  setenv("DFBH_FILAS_PARALELAS", "paralelas", 0);
    /* Tolerancia por relajacion (2026-10-21): 1e-9 con el afin (mas fina =
     * mas contraccion por pasada) y 1e-6 con compo (sus filas casi paralelas
     * cuestan pivotes sin ganar contraccion); se decide tras leer `--lin`. */
    const bool tol_externa = getenv("DFB_SX_FEASTOL") || getenv("DFB_SX_DUALTOL");
    if (!getenv("DFB_SX_FLIPINC"))        setenv("DFB_SX_FLIPINC", "1", 0);
    if (!getenv("DFB_SX_FEASTOL"))        setenv("DFB_SX_FEASTOL", "1e-9", 0);
    if (!getenv("DFB_SX_DUALTOL"))        setenv("DFB_SX_DUALTOL", "1e-9", 0);
    if (!getenv("DFB_SX_TOLANCHO"))       setenv("DFB_SX_TOLANCHO", "0", 0);
    if (!getenv("DFBH_CAJA_TOPE"))        setenv("DFBH_CAJA_TOPE", "igual", 0);
    /* Tope de 60 pivotes por LP en SC (2026-10-21): cualquier base certifica, asi
     * que cortar da una cota valida. Las cajas negras siguen con 20m+200: una caja
     * cortada por iteraciones no deja certificado y pierde la cota entera. */
    if (!getenv("DFBH_LPCAP") && !getenv("DFBH_LPCAP_M")) setenv("DFBH_LPCAP", "60", 0);

    if (argc < 2) { cerr << "uso: ibexsolve_dfb archivo.bch [--lr=xn|dfbhull] [-e eps] [-t s] [--seed=s] [--bfs]" << endl; return 1; }
    string lr = "xn", lin = "art", bis = "lsmear"; bool sols = false; double eps = 1e-3, timeout = -1, seed = 1.0; bool bfs = false;
    for (int i = 2; i < argc; ++i) {
        string a = argv[i];
        if (a.rfind("--lr=", 0) == 0) lr = a.substr(5);
        else if (a == "-b" && i+1 < argc) bis = argv[++i];
        else if (a == "-e" && i+1 < argc) eps = atof(argv[++i]);
        else if (a == "-t" && i+1 < argc) timeout = atof(argv[++i]);
        else if (a.rfind("--lin=", 0) == 0) lin = a.substr(6);
        else if (a.rfind("--seed=", 0) == 0) seed = atof(a.substr(7).c_str());
        else if (a == "--bfs") bfs = true;
        else if (a == "--sols") sols = true;
        else { cerr << "argumento desconocido: " << a << endl; return 1; }
    }
    if (lr != "xn" && lr != "dfbhull") { cerr << "--lr debe ser xn o dfbhull" << endl; return 1; }
    if (lin != "xt" && lin != "art" && lin != "compo") { cerr << "--lin debe ser art, xt o compo" << endl; return 1; }
    if (!tol_externa && lin == "compo") { setenv("DFB_SX_FEASTOL", "1e-6", 1); setenv("DFB_SX_DUALTOL", "1e-6", 1); }
    /* Con compo, el arranque tibio se elige midiendo (`DFBH_WARM=medir`); con el
     * afin, dual (el primal le cambia el arbol a Geneig-icse). */
    if (lin == "compo" && !getenv("DFBH_WARM") && !getenv("DFB_SX_PIVOTEO")) setenv("DFBH_WARM", "medir", 0);
    if (bis != "lsmear" && bis != "smearsumrel" && bis != "lsmeardfb") { cerr << "-b debe ser lsmear, smearsumrel o lsmeardfb" << endl; return 1; }
    if (bis == "lsmeardfb" && lr != "dfbhull") { cerr << "-b lsmeardfb requiere --lr=dfbhull" << endl; return 1; }

    System sys(argv[1]);
    RNG::srand(seed);

    // Sistema cuadrado de igualdades, como en DefaultSolver
    int nb_eq = 0;
    for (int i = 0; i < sys.nb_ctr; i++) if (sys.ctrs[i].op == EQ) nb_eq += sys.ctrs[i].f.image_dim();
    System* eqs = (sys.nb_var == nb_eq) ? new System(sys, System::EQ_ONLY) : NULL;

    CtcHC4 hc4(sys.ctrs, 0.01);
    CtcHC4 hc4diag(sys.ctrs, 0.01);
    if (getenv("DFBH_PIVDIAG") || getenv("DFBH_HC4RED") || getenv("DFBH_HC4VAR") || getenv("DFBH_HC4ACUM")) CtcDFBHull::hc4_diag = &hc4diag;
    CtcHC4 hc4acid(sys.ctrs, 0.1, true);
    AcidSonda acid(sys, hc4acid);
    CtcNewton* newton = eqs ? new CtcNewton(eqs->f_ctrs, 5e8, eps, 1.e-4) : NULL;
    /* --lin=art (omision desde 2026-10-05: el linearizador afin, una fila por
     * lado de cada restriccion; mas rapido con los dos contractores) o --lin=xt
     * (el de DefaultSolver: XTaylor con dos esquinas). Para reproducir al
     * ibexsolve oficial: --lr=xn --lin=xt. */
    Linearizer* linp;
    if (lin == "art") linp = new LinearizerAffine2(sys);
    else if (lin == "xt") linp = new LinearizerXTaylor(sys);
    else linp = new LinearizerCompo(*new LinearizerAffine2(sys), *new LinearizerXTaylor(sys));
    Ctc* clin = (lr == "dfbhull") ? (Ctc*) new CtcDFBHull(*linp) : (Ctc*) new CtcPolytopeHull(*linp);
    if (bis == "lsmeardfb") ((CtcDFBHull*) clin)->activar_sensibilidades();
    CtcHC4 hc4xn(sys.ctrs, 0.01);
    CtcCompo compo_xn(*clin, hc4xn);
    CtcFixPoint fix(compo_xn);

    Array<Ctc> lista(4); int k = 0;
    lista.set_ref(k++, hc4);
    lista.set_ref(k++, acid);
    if (newton) lista.set_ref(k++, *newton);
    lista.set_ref(k++, fix);
    lista.resize(k);
    CtcCompo ctc(lista);

    /* -b lsmear: sobre las igualdades si el sistema es cuadrado, sobre todo el
     * sistema si no. Sin objetivo, cada LP de LSmear minimiza o maximiza una
     * variable al azar. -b smearsumrel reproduce a DefaultSolver. */
    Bsc* bsc;
    if (bis == "lsmear") bsc = new LSmear(eqs ? *eqs : sys, eps);
    else if (bis == "lsmeardfb") bsc = new LSmearDFB(sys, eps);
    else bsc = eqs ? (Bsc*) new SmearSumRelative(*eqs, eps) : (Bsc*) new RoundRobin(eps);
    CellBuffer* buf = bfs ? (CellBuffer*) new CellList() : (CellBuffer*) new CellStack();

    /* `DFBH_SOLO=1`: aplica cada contractor UNA vez, por separado, a la caja
     * inicial e informa si la vacia (para cajas volcadas con DFBH_DUMPNODO). */
    if (getenv("DFBH_SOLO")) {
        const char* nom[4] = {"hc4", "acid", "newton", "fix(dfb,hc4)"};
        Ctc* cs[4] = {&hc4, &acid, newton, &fix};
        for (int q = 0; q < 4; ++q) {
            if (!cs[q]) continue;
            IntervalVector b(sys.box);
            cs[q]->contract(b);
            cerr << "[solo] " << nom[q] << ": " << (b.is_empty() ? "VACIA" : "no vacia") << endl;
        }
        IntervalVector b(sys.box); clin->contract(b);
        cerr << "[solo] dfb una vez: " << (b.is_empty() ? "VACIA" : "no vacia") << endl;
        return 0;
    }
    Solver s(sys, ctc, *bsc, *buf, Vector(sys.nb_var, eps), Vector(sys.nb_var, POS_INFINITY));
    if (timeout > 0) s.time_limit = timeout;
    /* `DFBH_SOLTRACE=1`: imprime cada solucion, incierta y pendiente al encontrarla. */
    if (getenv("DFBH_SOLTRACE")) s.trace = 1;
    Solver::Status st = s.solve(sys.box);
    const CovSolverData& d = s.get_data();
    s.report();
    if (lr == "dfbhull" && getenv("DFBH_CUENTAS"))
        cerr << "[cuentas] celdas=" << s.get_nb_cells() << " llamadas=" << CtcDFBHull::n_nodes
             << " LPs=" << CtcDFBHull::n_lps << " pivotes=" << CtcDFBHull::n_pivots << " fusionadas=" << CtcDFBHull::n_filas_fusionadas << endl;
    if (CtcDFBHull::modo_cortes > 0)
        cerr << "[cortes] generados=" << CtcDFBHull::n_cortes_gen << " filas agregadas=" << CtcDFBHull::n_cortes_filas
             << " en " << CtcDFBHull::n_cortes_llam << " llamadas" << endl;
    if (bis == "lsmeardfb")
        cerr << "[lsmeardfb] llamadas=" << LSmearDFB::n_llamadas << " sin_sens=" << LSmearDFB::n_sin_sens
             << " fuera=" << LSmearDFB::n_fuera << " eligio_dfb=" << LSmearDFB::n_dfb << endl;
    if (sols) cout << d << endl;
    cout << s.get_nb_cells() << " " << s.get_time() << " " << (int) st << " "
         << d.nb_solution() << " " << d.nb_boundary() << " " << d.nb_unknown() << " " << d.nb_pending() << endl;
    delete clin;   /* los destructores imprimen las sondas */
    return 0;
}
