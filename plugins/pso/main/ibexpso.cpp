/*
 * optimsolver.cpp
 *
 *  Created on: 21 dic. 2017
 *      Author: iaraya
 */

#include "ibex.h"
#include "args.hxx"

using namespace std;
using namespace ibex;

int main(int argc, char** argv){

	try {

		args::ArgumentParser parser("********* PSO *********.", "Solve a Minibex file.");
		args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
		args::ValueFlag<std::string> _strategy(parser, "string", "the search strategy", {'s', "search"});
		args::ValueFlag<int> _seed(parser, "int", "seed", {"seed"});
		args::Flag trace(parser, "trace", "Activate trace. Updates of loup/uplo are printed while minimizing.", {"trace"});

		args::Positional<std::string> filename(parser, "filename", "The name of the MINIBEX file.");

		try
		{
			parser.ParseCLI(argc, argv);
		}
		catch (args::Help&)
		{
			std::cout << parser;
			return 0;
		}
		catch (args::ParseError& e)
		{
			std::cerr << e.what() << std::endl;
			std::cerr << parser;
			return 1;
		}
		catch (args::ValidationError& e)
		{
			std::cerr << e.what() << std::endl;
			std::cerr << parser;
			return 1;
		}

		string strategy= (_strategy)? _strategy.Get() : "feasdiv";
		if (trace) cout << "Search strategy: " << strategy << endl;
		int nseed= (_seed)? _seed.Get() : 1 ;
		if (trace) cout << "seed: " << nseed << endl;
		PSOSwarm::trace = trace;

		double eqeps= 1.e-8;
		double default_relax_ratio = 0.2;
		srand(nseed);

		System* orig_sys,*sys;
		LoupFinderDefault* loupfinder;
		std::size_t found = string(filename.Get().c_str()).find(".nl");
		if (found!=std::string::npos){
			   AmplInterface interface (argv[1]);
			   orig_sys= new System(interface);
		 }else
			   orig_sys = new System(filename.Get().c_str());

		// the extended system
		if(strategy=="solver") sys=orig_sys;
		else {
			sys=new ExtendedSystem(*orig_sys,eqeps);
			NormalizedSystem* norm_sys = new NormalizedSystem(*orig_sys,eqeps);
			loupfinder = new  LoupFinderDefault(*norm_sys,true);
		}

		Timer timer;
		timer.start();

		double c1 = 2;
		double c2 = 2;
		double constrictorMod = 1;
		int particles = 200;
		int iterations = 2000;
		PSOSwarm* swarm = new PSOSwarm(c1,c2,particles,iterations);

		Vector valueSwarm = swarm->executePSO(orig_sys, false, constrictorMod);
		cout << "\033[0mPSO result:" << valueSwarm << endl;

		timer.stop();
		// Display the cpu time used
		cout << "cpu time used=" << timer.get_time() << "s."<< endl;
	}


	catch(ibex::SyntaxError& e) {
	  cout << e << endl;
	}
}



