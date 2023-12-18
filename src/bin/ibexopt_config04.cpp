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

    args::ValueFlag<std::string> _filtering(parser, "filtering", "Filtering option (hc4|acidhc4*|3bcidhc4)", {'f', "filtering"});
    args::ValueFlag<std::string> _linearRelaxation(parser, "linear relaxation", "Linear relaxation option (xn*)", {"lr"});
    args::ValueFlag<std::string> _bisection(parser, "bisection", "Bisection option (roundrobin|largestfirst|largestfirstnoobj|smearsum|smearmax|smearsumrel|smearmaxrel|lsmear|lsmearmg*)", {'b'});
    args::ValueFlag<std::string> _strategy(parser, "strategy", "Search Strategy (bfs|dh|bs*)", {'s'});
    args::ValueFlag<int> _beamsize(parser, "beamsize", "Beamsize (default=1)", {'B'});
    args::ValueFlag<double> _prec(parser, "prec", "Precision (default=1e-7)", {"eps_x"});
    args::ValueFlag<double> _goalPrec(parser, "goal precision (default=1e-6)", "Goal precision", {"eps"});
    args::ValueFlag<double> _timeLimit(parser, "time limit", "Time limit (default=1000)", {'t'});
    args::ValueFlag<int> _randomSeed(parser, "random seed", "Random seed (default=42)", {"seed"});
	args::Flag quiet(parser, "quiet", "Print no report on the standard output.",{'q',"quiet"});

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
