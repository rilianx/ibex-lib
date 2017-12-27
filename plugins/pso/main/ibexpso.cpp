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

		/*if (trace) cout << *orig_sys << endl;

		IntervalVector box = orig_sys->box; //dominios de las variables
		if (trace) cout << "dominio de la primera variable: " << box[0] << endl;

		Vector v = box.random();
		if (trace) cout << "vector aleatorio:" << v << endl;


		if (trace) cout << "eval(obj):" << orig_sys->goal->eval(v).ub() << endl;

		for(int i=0; i<orig_sys->ctrs.size();i++){
			Interval eval=orig_sys->ctrs[i].f.eval(v);
			if((orig_sys->ctrs[i].op == LEQ || orig_sys->ctrs[i].op == LT) && eval.ub()>0.0)
				cout << "error("<< i <<"):" << eval.ub() << endl;
			else if((orig_sys->ctrs[i].op == GEQ || orig_sys->ctrs[i].op == GT) && eval.lb()<0.0)
				cout << "error("<< i <<"):" << -eval.lb() << endl;
			else if(orig_sys->ctrs[i].op == EQ )
				cout << "error("<< i <<"):" << abs(eval).ub() << endl;
		}*/
		if(trace) cout << "cduartes" << endl;
		PSOSwarm* swarm = new PSOSwarm(2,2,100,1000);

		Vector valueSwarm = swarm->executePSO(orig_sys);
		if (trace) cout << "PSO result:" << valueSwarm << endl;
	}


	catch(ibex::SyntaxError& e) {
	  cout << e << endl;
	}
}



