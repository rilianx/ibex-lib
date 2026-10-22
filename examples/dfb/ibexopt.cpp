#include <cstdio>
//============================================================================
//                                  I B E X
//
//                               ************
//                                  IbexOpt
//                               ************
//
// Author      : Ignacio Araya, Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Last Update : Jul 09, 2017
//============================================================================

#include "ibex.h"
#include "../src/bin/parse_args.h"

#include <sstream>
#include <regex>
#include "ibex_Optimizer05Config.h"
#include "ibex_CtcDualFeasibleBounding.h"
#include "ibex_CtcDFBPropag.h"
#include "ibex_DFBSimplex.h"


using namespace std;
using namespace ibex;

void printArgsSummary(const std::vector<std::string>& args) {
    if (args.size() < 11) {
        std::cerr << "Número insuficiente de argumentos." << std::endl;
        return;
    }

    std::cout << "Resumen de Argumentos:" << std::endl;
    std::cout << "Filename: " << args[1] << std::endl;
    std::cout << "Filtering: " << args[2] << std::endl;
    std::cout << "Linear Relaxation: " << args[3] << std::endl;
    std::cout << "Bisection: " << args[4] << std::endl;
    std::cout << "Strategy: " << args[5] << std::endl;
    std::cout << "Beamsize: " << args[6] << std::endl;
    std::cout << "Precision (eps_x): " << args[7] << std::endl;
    std::cout << "Goal Precision (eps): " << args[8] << std::endl;
    std::cout << "Time Limit: " << args[9] << std::endl;
    std::cout << "Random Seed: " << args[10] << std::endl;
}

// Función para parsear el string a IntervalVector
IntervalVector parseIntervalVector(const std::string& box_str) {
    std::regex interval_regex(R"(\[\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*,\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*\])");
    std::sregex_iterator next(box_str.begin(), box_str.end(), interval_regex);
    std::sregex_iterator end;

    std::vector<std::pair<double, double>> intervals;
    while (next != end) {
        std::smatch match = *next;
        double lower = std::stod(match[1].str());
        double upper = std::stod(match[2].str());
        intervals.emplace_back(lower, upper);
        ++next;
    }

    int n = intervals.size();
    double bounds[n][2];
    for (int i = 0; i < n; ++i) {
        bounds[i][0] = intervals[i].first;
        bounds[i][1] = intervals[i].second;
    }

    return IntervalVector(n, bounds);
}

int main(int argc, char** argv) {
    /* Omision de `ibexopt` (2026-10-13, ESTADO §4): estrategia primal, con
     * pasos primales tambien en el motor revisado y revisado por encima de 15
     * filas. En los LPs tibios de la pasada el primal conserva el vertice
     * anterior; el dual puro lo descarta al reubicar y hace ~2x pivotes.
     * Elegido en el banco del paper (cpu x0.94) y validado en `all` (x0.91).
     * `ibexsolve_dfb` queda en dual puro con umbral 30. Si `DFB_SX_PIVOTEO`
     * esta definida, mandan las variables de entorno como antes. */
    if (!getenv("DFB_SX_PIVOTEO")) {
        DFBSimplex::pivoteo_dual = false;
        DFBSimplex::dual_puro = false;
        if (!getenv("DFB_SX_REVPRIMAL")) DFBSimplex::rev_primal = true;
        if (!getenv("DFB_SX_MOTOR_M")) DFBSimplex::motor_umbral = 15;
    }
    /* Orden de cotas del optimizador (2026-10-18): la mas cercana primero. Con
     * el arranque tibio primal la base anterior conserva su vertice, y la cota
     * mas cercana es la que menos pivotes primales necesita; el orden por costo
     * (conteo de infactibilidades duales) es el adecuado para el arranque dual
     * del solver. En los 240 problemas: cpu x0.98 / x0.97 / x0.95, mismas
     * celdas. `DFBH_ORDEN` definida manda como antes (`costo` es el anterior). */
    if (!getenv("DFBH_ORDEN")) setenv("DFBH_ORDEN", "cerca", 0);
    /* CORRECCIONES POR OMISION (2026-10-20; cada una se apaga desde fuera):
     *  - DFB_SX_RAYO=fila: el certificado de vacio con la fila infactible de B^-1;
     *  - DFB_SX_PISO=error: piso de la tolerancia primal al error de la fila;
     *  - DFBH_FILAS_PARALELAS=paralelas: fusion de filas proporcionales (compo);
     *  - tolerancias de SC iguales a las de la caja negra en optimizacion (1e-9),
     *    sin el tope por ancho; tope de iteraciones de la caja igual al de SC. */
    if (!getenv("DFB_SX_RAYO"))           setenv("DFB_SX_RAYO", "fila", 0);
    if (!getenv("DFB_SX_PISO"))           setenv("DFB_SX_PISO", "error", 0);
    if (!getenv("DFB_SX_FLIPINC"))        setenv("DFB_SX_FLIPINC", "1", 0);
    if (!getenv("DFBH_CAJA_TOL"))         setenv("DFBH_CAJA_TOL", "1e-9", 0);   /* la mejor de la caja en opt (hibrido) */
    if (!getenv("DFBH_FILAS_PARALELAS"))  setenv("DFBH_FILAS_PARALELAS", "paralelas", 0);
    if (!getenv("DFB_SX_FEASTOL"))        setenv("DFB_SX_FEASTOL", "1e-9", 0);
    if (!getenv("DFB_SX_DUALTOL"))        setenv("DFB_SX_DUALTOL", "1e-9", 0);
    if (!getenv("DFB_SX_TOLANCHO"))       setenv("DFB_SX_TOLANCHO", "0", 0);
    if (!getenv("DFBH_CAJA_TOPE"))        setenv("DFBH_CAJA_TOPE", "igual", 0);
    /* Tope de 60 pivotes por LP en SC (2026-10-21): cualquier base certifica, asi
     * que cortar da una cota valida. Las cajas negras siguen con 20m+200: una caja
     * cortada por iteraciones no deja certificado y pierde la cota entera. */
    if (!getenv("DFBH_LPCAP") && !getenv("DFBH_LPCAP_M")) setenv("DFBH_LPCAP", "60", 0);
    // Menú interactivo para configurar las opciones
    args::ArgumentParser parser("Optimizer04 configuration", "Configure and run Optimizer04");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});


    // Continuando con la configuración de tus argumentos existentes...
    args::Group filteringGroup(parser, "Filtering Options", args::Group::Validators::DontCare);
    args::ValueFlag<std::string> _filtering(filteringGroup, "filtering", "Filtering option (hc4|acidhc4*|3bcidhc4|dfb)", {"filtering"});

    args::Group hc4Group(filteringGroup, "HC4 Options", args::Group::Validators::DontCare);
    args::ValueFlag<double> _precHc4(hc4Group, "prec_hc4", "Precision for hc4 (default=1e-7)", {"hc4_p"});

    args::Group bcidhc4Group(filteringGroup, "3BCIDHC4 Options", args::Group::Validators::DontCare);
    args::ValueFlag<double> _prec3Bcidhc4(bcidhc4Group, "prec_3bcidhc4", "Precision for 3bcidhc4 (default=1e-7)", {"prec_3bcidhc4"});
    args::ValueFlag<int> _slices3Bcidhc4(bcidhc4Group, "slices_3bcidhc4", "Number of slices for 3bcidhc4", {"slices"});


    //args::ValueFlag<std::string> _filtering(parser, "filtering", "Filtering option (hc4|acidhc4*|3bcidhc4)", {'f', "filtering"});
    args::ValueFlag<std::string> _linearRelaxation(parser, "linear relaxation", "Linear relaxation option (xn*)", {"lr"});
    args::ValueFlag<std::string> _bisection(parser, "bisection", "Bisection option (roundrobin|largestfirst|largestfirstnoobj|smearsum|smearmax|smearsumrel|smearmaxrel|lsmear|lsmearmg*)", {'b'});
    args::ValueFlag<std::string> _strategy(parser, "strategy", "Search Strategy (bfs|dh|bs*)", {'s'});
    args::ValueFlag<int> _beamsize(parser, "beamsize", "Beamsize (default=1)", {'B'});
    args::ValueFlag<double> _prec(parser, "prec", "Precision (default=1e-7)", {"eps_x"});
    args::ValueFlag<double> _goalPrec(parser, "goal precision (default=1e-6)", "Goal precision", {"eps"});
    args::ValueFlag<double> _initial_loup(parser, "initial loup", "Initial loup", {"loup"});
    args::ValueFlag<std::string> _initial_box(parser, "initial box", "Initial Box", {"box"});

    args::ValueFlag<double> _timeLimit(parser, "time limit", "Time limit (default=1000)", {'t'});
    args::ValueFlag<int> _randomSeed(parser, "random seed", "Random seed (default=42)", {"seed"});

	args::Positional<std::string> filename(parser, "filename", "The name of the MINIBEX file.");


    try {
        parser.ParseCLI(argc, argv);
    } catch (const args::Help&) {
        std::cout << parser;
        return 0;
    } catch (const args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    } catch (const args::ValidationError& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // Aquí se recogen los valores de las opciones del menú
    // Si alguna opción es obligatoria, se debe verificar que haya sido proporcionada
    // Por ejemplo:
    if(!filename) {
        std::cerr << "Filename is required." << std::endl;
        return 1;
    }

    // Valores por defecto
    std::string defaultFiltering = "acidhc4";
    std::string defaultLinearRelaxation = "xn";
    std::string defaultBisection = "lsmearmg";
    std::string defaultStrategy = "bs";
    int defaultBeamsize = 1;
    double defaultPrec = 1.e-7;
    double defaultGoalPrec = 1.e-6;
    double defaultTimeLimit = 1000;
    int defaultRandomSeed = 42;


    /* Los double se pasan con precision completa. `std::to_string` usa `%f`
     * (6 decimales), y como gaol deja el redondeo hacia arriba desde el
     * arranque, 1e-7 salia como "0.000001": eps_x valia 1e-6, igual que eps,
     * en vez de ser menor. Cualquier valor por debajo de 1e-6 se perdia. */
    auto num = [](double v) { char b[64]; snprintf(b, sizeof b, "%.17g", v); return std::string(b); };

    // Construir los argumentos para el constructor de Optimizer04Config
    std::vector<std::string> args = {
        "optimizer04", filename.Get(),
        _filtering ? _filtering.Get() : defaultFiltering,
        _linearRelaxation ? _linearRelaxation.Get() : defaultLinearRelaxation,
        _bisection ? _bisection.Get() : defaultBisection,
        _strategy ? _strategy.Get() : defaultStrategy,
        _beamsize ? std::to_string(_beamsize.Get()) : std::to_string(defaultBeamsize),
        num(_prec ? _prec.Get() : defaultPrec),
        num(_goalPrec ? _goalPrec.Get() : defaultGoalPrec),
        num(_timeLimit ? _timeLimit.Get() : defaultTimeLimit),
        _randomSeed ? std::to_string(_randomSeed.Get()) : std::to_string(defaultRandomSeed)
    };

    //printArgsSummary(args);

    // Convertir std::vector<std::string> a char*[]
    std::vector<char*> argv_new;
    for (auto& arg : args) {
        argv_new.push_back(&arg[0]);
    }

    // Crear el objeto Optimizer04Config
    Optimizer05Config config(argv_new.size(), argv_new.data());

    // Build the default optimizer
    Optimizer o(config);

    // display solutions with up to 12 decimals
    cout.precision(12);


    // Search for the optimum
    double initial_loup = _initial_loup ? _initial_loup.Get() : POS_INFINITY;
    if (_initial_box) {
        IntervalVector initial_box = parseIntervalVector(_initial_box.Get());
        //cout << initial_box << endl;
        o.optimize(initial_box, initial_loup);
    } else
        o.optimize(config.sys->box, initial_loup);

    o.report();
    cout << o.get_nb_cells() << " " << o.get_time() << " " << o.get_loup() << " " << o.get_loup() << endl;

    /* DFB_PHASES=1: desglose del tiempo que se fue en el contractor DFB, para
     * saber el techo de cualquier optimizacion suya frente al tiempo total. */
    if (getenv("DFB_PHASES")) {
        const double T = o.get_time();
        const double c = CtcDFB::t_pricing + CtcDFB::t_ratio_test + CtcDFB::t_pivot
                       + CtcDFB::t_bound + CtcDFB::t_regen + CtcDFB::t_init
                       + CtcDFB::t_certify;
        printf("PHASES,%s,%g,%g,%g,%g,%g,%g,%g,%g,%ld,%ld,%ld,%ld\n",
               filename.Get().c_str(), T, c,
               CtcDFB::t_init, CtcDFB::t_pricing, CtcDFB::t_pivot,
               CtcDFB::t_certify, CtcDFB::t_bound, CtcDFB::t_ratio_test,
               CtcDFB::n_inits, CtcDFB::n_applied, CtcDFB::n_certifications,
               (long)CtcDFBPropag::n_relin);
        printf("REVIVE,%s,%ld,%ld,%ld,%ld\n", filename.Get().c_str(),
               CtcDFBPropag::n_revive, CtcDFBPropag::n_impact_calls,
               CtcDFB::n_applied, CtcDFB::n_inits);
    }


    return 0;
}
