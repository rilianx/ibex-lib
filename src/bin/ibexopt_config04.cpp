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
#include "parse_args.h"

#include <sstream>

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

int main(int argc, char** argv) {
    // Menú interactivo para configurar las opciones
    args::ArgumentParser parser("Optimizer04 configuration", "Configure and run Optimizer04");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});


    // Continuando con la configuración de tus argumentos existentes...
    args::Group filteringGroup(parser, "Filtering Options", args::Group::Validators::DontCare);
    args::ValueFlag<std::string> _filtering(filteringGroup, "filtering", "Filtering option (hc4|acidhc4*|3bcidhc4)", {"filtering"});
    //args::ValueFlag<std::string> _filtering(parser, "filtering", "Filtering option (hc4|acidhc4*|3bcidhc4)", {'f', "filtering"});
    args::ValueFlag<std::string> _linearRelaxation(parser, "linear relaxation", "Linear relaxation option (no|xn*)", {"lr"});
    args::ValueFlag<std::string> _bisection(parser, "bisection", "Bisection option (roundrobin|largestfirst|largestfirstnoobj|smearsum|smearmax|smearsumrel|smearmaxrel|lsmear|lsmearmg*)", {'b'});
    args::ValueFlag<double> _ratio_bisect(parser, "ratio_bisect", "Ratio for bisection (default=0.5)", {"ratio_bisect"});
    args::ValueFlag<std::string> _strategy(parser, "strategy", "Search Strategy (bfs|dh|bs*)", {'s'});
    args::ValueFlag<int> _beamsize(parser, "beamsize", "Beamsize (default=1)", {'B'});
    args::ValueFlag<double> _prec(parser, "prec", "Precision (default=1e-7)", {"eps_x"});
    args::ValueFlag<double> _goalPrec(parser, "goal precision (default=1e-6)", "Goal precision", {"eps"});
    args::ValueFlag<double> _timeLimit(parser, "time limit", "Time limit (default=1000)", {'t'});
    args::ValueFlag<int> _randomSeed(parser, "random seed", "Random seed (default=42)", {"seed"});
	args::Flag quiet(parser, "quiet", "Print no report on the standard output.",{'q',"quiet"});

    // Propagation ratios
    args::Group propagationGroup(parser, "Propagation Ratios", args::Group::Validators::DontCare);
    args::ValueFlag<double> _hc4_ratio_propag(parser, "hc4_ratio_propag", "Ratio for hc4 propagation (default=0.01)", {"hc4_ratio_propag"});
    args::ValueFlag<double> _hc44cid_propag(parser, "hc44cid_propag", "Ratio for hc44cid propagation (default=0.1)", {"hc44cid_propag"});
    args::ValueFlag<double> _hc44xn_propag(parser, "hc44xn_propag", "Ratio for hc44xn propagation (default=0.01)", {"hc44xn_propag"});
    args::ValueFlag<double> _hc43bcidhc4_propag(parser, "hc43bcidhc4_propag", "Ratio for hc43bcidhc4 propagation (default=0.1)", {"hc43bcidhc4_propag"});
    args::ValueFlag<double> _poly_hc4_propag(parser, "poly_hc4_propag", "Ratio for poly_hc4 propagation (default=0.1)", {"poly_hc4_propag"});
    args::ValueFlag<double> _all_propag(parser, "all_propag", "Ratio for all propagation (default=0.1)", {"all_propag"});


    //ACID & 3BCID
    args::Group acid3bcidGroup(parser, "ACID & 3BCID Options", args::Group::Validators::DontCare);
    args::ValueFlag<int> _s3b(parser, "s3b", "Number of slices for 3bcid (default=10)", {"s3b"});
    args::ValueFlag<int> _scid(parser, "scid", "Number of slices for cid (default=1)", {"scid"});
    args::ValueFlag<double> _var_min_width(parser, "var_min_width", "Minimum width for variables (default=1e-11)", {"var_min_width"});
    args::Flag _optim(parser, "optim", "Optimize the solution (default=true)", {"optim"});
    args::ValueFlag<double> _ct_ratio(parser, "ct_ratio", "CT ratio (default=0.002)", {"ct_ratio"});

    //group linear relaxation... add group
    args::Group polyRelaxationGroup(parser, "Polynomial Relaxation Options", args::Group::Validators::DontCare);
    args::ValueFlag<int> _poly_maxiter(parser, "poly_maxiter", "Maximum iterations for polynomial relaxation (default=100)", {"poly_maxiter"});
    args::ValueFlag<double> _poly_eps(parser, "poly_eps", "Epsilon for polynomial relaxation (default=1e-9)", {"poly_eps"});
    args::ValueFlag<double> _relax_ratio(parser, "relax_ratio", "Relaxation ratio (default=0.2)", {"relax_ratio"});

    //LoupFinder  
    args::Group loupfinderGroup(parser, "Loupfinder Options", args::Group::Validators::DontCare);
    args::Flag _in_hc4(parser, "in_hc4", "Use hc4 in loupfinder (default=true)", {"in_hc4"});



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

    Optimizer04Config::hc4_ratio_propag=_hc4_ratio_propag? _hc4_ratio_propag.Get() : 0.01;
    Optimizer04Config::hc44cid_propag= _hc44cid_propag ? _hc44cid_propag.Get() : 0.1;
    Optimizer04Config::hc44xn_propag= _hc44xn_propag ? _hc44xn_propag.Get() : 0.01;
    Optimizer04Config::hc43bcidhc4_propag= _hc43bcidhc4_propag ? _hc43bcidhc4_propag.Get() : 0.1;
    Optimizer04Config::poly_hc4_propag= _poly_hc4_propag ? _poly_hc4_propag.Get() : 0.1;
    Optimizer04Config::all_propag= _all_propag ? _all_propag.Get() : 0.1;

    Optimizer04Config::ratio_bisect= _ratio_bisect ? _ratio_bisect.Get() : 0.5;

    //3bcid & acid
    Optimizer04Config::s3b= _s3b ? _s3b.Get() : 10;
    Optimizer04Config::scid= _scid ? _scid.Get() : 1;
    Optimizer04Config::var_min_width= _var_min_width ? _var_min_width.Get() : 1e-11;
    //only acid
    Optimizer04Config::optim= _optim ? _optim.Get() : true;
    Optimizer04Config::ct_ratio= _ct_ratio ? _ct_ratio.Get() : 0.002;

    //lr=xn
    Optimizer04Config::poly_maxiter= _poly_maxiter ? _poly_maxiter.Get() : 100;
    Optimizer04Config::poly_eps= _poly_eps ? _poly_eps.Get() : 1e-9;
    Optimizer04Config::relax_ratio= _relax_ratio ? _relax_ratio.Get() : 0.2;

    //loupfinder
    Optimizer04Config::in_hc4= _in_hc4 ? _in_hc4.Get() : true;


    // Construir los argumentos para el constructor de Optimizer04Config
    std::vector<std::string> args = {
        "optimizer04", filename.Get(),
        _filtering ? _filtering.Get() : defaultFiltering,
        _linearRelaxation ? _linearRelaxation.Get() : defaultLinearRelaxation,
        _bisection ? _bisection.Get() : defaultBisection,
        _strategy ? _strategy.Get() : defaultStrategy,
        _beamsize ? std::to_string(_beamsize.Get()) : std::to_string(defaultBeamsize),
        _prec ? std::to_string(_prec.Get()) : std::to_string(defaultPrec),
        _goalPrec ? std::to_string(_goalPrec.Get()) : std::to_string(defaultGoalPrec),
        _timeLimit ? std::to_string(_timeLimit.Get()) : std::to_string(defaultTimeLimit),
        _randomSeed ? std::to_string(_randomSeed.Get()) : std::to_string(defaultRandomSeed)
    };

    printArgsSummary(args);

    // Convertir std::vector<std::string> a char*[]
    std::vector<char*> argv_new;
    for (auto& arg : args) {
        argv_new.push_back(&arg[0]);
    }

    // Crear el objeto Optimizer04Config
    Optimizer04Config config(argv_new.size(), argv_new.data());

    // Build the default optimizer
    Optimizer o(config);

    // display solutions with up to 12 decimals
    cout.precision(12);

    if (!quiet)
        cout << "\n\nrunning............" << endl << endl;

    // Search for the optimum
    o.optimize(config.sys->box);

    // Report some information (computation time, etc.)
    if (!quiet)
        o.report();


    return 0;
}
