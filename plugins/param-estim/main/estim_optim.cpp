#include "ibex.h"
#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <set>
#include "args.hxx"
using namespace std;
using namespace ibex;

/*
 * 
 * Detects plan in real 3d point clouds using q-intersection.
 * From original work of Bertrand Neveu, Martin de la Gorce and Gilles Trombettoni
 * "Improving a Constraint Programming Approach for Parameter Estimation"
 * 
 */

int main(int argc, char** argv) {
	try{
		args::ArgumentParser parser("********* QInter project (QInter Optimization) *********", "Solve a problem's file with  Q-Inter.");
		args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
		args::Positional<std::string> filename(parser, "filename", "The name of the problem file.");
		args::ValueFlag<double> _epseq(parser, "double", "epsilon tolerance", {'e',"epseq"});
		args::ValueFlag<int> _Q(parser, "int", "known number of inliers", {'q',"Qinter"});
		args::ValueFlag<double> _precbc(parser, "double", "precision for bc", {"precbc"});
		args::ValueFlag<double> _precd(parser, "double", "precision for c", {"precd"});
		args::ValueFlag<int> _flist(parser, "int", "f list", {"flist"});
		args::ValueFlag<int> _gaplimit(parser, "int", "gap limit", {'g',"gaplimit"});
		args::ValueFlag<int> _nbrand(parser, "int", "n brand", {'n', "nbrand"});
		args::ValueFlag<std::string> _bisector(parser, "string", "the bisection method", {'b', "bis"});
		args::ValueFlag<int> _seed(parser, "int", "seed integer", {'s', "seed"});
		args::Flag trace(parser, "trace", "Activate trace.", {"trace"});

		try{
			parser.ParseCLI(argc, argv);
		}catch (args::Help&){
			std::cout << parser;
			return 0;
		}catch (args::ParseError& e){
			std::cerr << e.what() << std::endl;
			std::cerr << parser;
			return 1;
		}catch (args::ValidationError& e){
			std::cerr << e.what() << std::endl;
			std::cerr << parser;
			return 1;
		}
	
		vector<double> *x = new vector<double>;
		vector<double> *y = new vector<double>;
		vector<double> *z = new vector<double>;

		/*
		* Parameter section
		* */
		double cputime = 0;
		double totaltime = 0;
		string input_file_name = (filename)? filename.Get() : "randomize";
		double epseq = (_epseq)? _epseq.Get(): 1.e-3;
		int Q = (_Q)? _Q.Get(): 0;
		double precbc = (_precbc)? _precbc.Get(): 1.e-4;
		double precd = (_precd)? _precd.Get(): 1.e-4;
		int flist = (_flist)? _flist.Get(): 1;
		int gaplimit = (_gaplimit)? _gaplimit.Get(): 0;
		int nbrand = (_nbrand)? _nbrand.Get(): 0;
		string bisect = (_bisector)? _bisector.Get(): "rr";
		int seed = (_seed)? _seed.Get(): 1;
		cout << input_file_name << endl;
		ifstream input(input_file_name.c_str());

		/**
		* Load data from file 
		**/
		while (!input.eof()){
			double in;
			input >> in;
			x->push_back(in);
			// cout << in << " " ;
			input >> in;
			y->push_back(in);
			// cout << in << " " ;
			input >> in;
			z->push_back(in);
			// cout << in << " " << endl;
		}
		x->pop_back();y->pop_back();z->pop_back();
		srand(seed);

		cout << "Configuration:" << endl;
		cout << "\tfilename: " << input_file_name << endl;
		cout << "\tepseq: " << epseq << endl;
		cout << "\tQ: " << Q << endl;
		cout << "\tprecbc: " << precbc << endl;
		cout << "\tprecd: " << precd << endl;
		cout << "\tflist: " << flist << endl;
		cout << "\tgap limit: " << gaplimit << endl;
		cout << "\tnbrand: " << nbrand << endl;
		cout << "\tbisector method: " << bisect << endl;
		cout << "\tseed: " << seed << endl;

		if(trace){
			cout << "press a key to continue..." << endl;
			getchar();
		}

		if(trace) cout << "nb points " << x->size() << endl;
		
		int n = 3;
		Variable u(3);
		Variable w(3);
		int p = x->size();
		int nz = 0;
		for (int i=0;i<p;i++)
			if ((*z)[i] >=2.79 && (*z)[i] <3)
				nz++;
		if(trace) cout << "nz " << nz << endl;
		int K = 1;
		int np = p;
		double ** linfun;
		linfun = new double*[p];
		for (int i=0; i<p; i++)
			linfun[i] = new double[n];

		Array<Ctc> m_ctc(p);
		Array<Ctc> m_ctc1(np);
		Function* m_func[p] ;	
		Function* m_f0;
		Ctc* ctc0;
		Variable v(3);

		clock_t start, start0, end;
		int nb_cells = 0;
		start = clock();
		start0 = clock();
		int Q0 = Q;
		int Qvalid = Q;

		if(trace) cout << np << " " << Q0 << " " << Qvalid << " " << epseq << endl;

		int bestsolpointnumber=0;
		int Qoct=Q;
		Vector bestsol (3);
		/**
		 * Main optimization loop
		 **/
		for (int oct=0; oct <4; oct++){
			if(trace) cout << "oct: " << oct << endl;
			int diry= pow(-1,oct%2);
			int dirz= pow(-1,oct/2);
			m_f0= new Function(v, (diry*v[0]+dirz*v[1]-1));
			ctc0=new CtcFwdBwd(*m_f0,LEQ);
			for (int i=0; i<p; i++){
				m_func[i] = new Function(v,(x->at(i) +v[0]*(y->at(i)- diry*x->at(i))+v[1]*(z->at(i)-dirz*x->at(i))-v[2]-Interval(-epseq,epseq)));
				m_ctc.set_ref(i,(*new CtcFwdBwd(*m_func[i])));
			}

			for (int i=0; i<p; i++){
				/* We must be on the plane defined by v */
				m_ctc1.set_ref(i,m_ctc[i]);
				linfun[i][0]=x->at(i);
				linfun[i][1]=y->at(i)-diry*x->at(i);
				linfun[i][2]=z->at(i)-dirz*x->at(i);
			}

			double _box[3][2];

			if (diry==1) {
				_box[0][0]= 0;
				_box[0][1]= 1;
			}else{
				_box[0][0]= -1;
				_box[0][1]=  0;
			}
			if(dirz==1){
				_box[1][0]= 0;
				_box[1][1]= 1;
			}else{
				_box[1][0]= -1;
				_box[1][1]=  0;
			}
			_box[2][0]= -100;
			_box[2][1]=  100;

			IntervalVector box(3,_box);
			Vector prec(3);
			
			prec[0]=precbc;
			prec[1]=precbc;
			prec[2]=precd;

			if(trace) cout << "precision bisection " <<  prec[0] << " " << prec[1] << " " << prec[2] << endl;
		
			Vector proba(3);
			proba[0]=0.33;
			proba[1]=0.33;
			proba[2]=0.34;

			CellHeapQInter buff;
			BeamSearch str(buff);

			Bsc * bs;
			if (bisect == "rr")
				bs = new RoundRobin(prec,0.5);
			else if (bisect == "rr2")
				bs = new RoundRobinNvar(2,prec,0.5);

			/** 
			 * Optimization zone 
			 **/
			CtcQInterAffPlane ctcq(n,p,m_ctc1,linfun,epseq,Qoct,QINTERPROJ);
						
			CtcCompo ctcqf0(*ctc0,ctcq);
			CtcFixPoint ctcf(ctcqf0,1);

			SolverOptQInter s(ctcf,*bs,str,ctcq);

			s.str.with_oracle = false;
			s.str.with_storage = true;
			s.timeout = 3600;
			s.trace = trace;
			s.gaplimit = gaplimit;
			s.feasible_tries = nbrand;

			s.bestsolpointnumber = bestsolpointnumber;
			s.bestsolpoint=bestsol;
			if(trace) cout << "\tavant resolution q:" << ctcq.q << endl;
			IntervalVector res=s.solve(box);

			if(trace) cout << "\tNumber of branches : " << s.nb_cells << endl;
			nb_cells += s.nb_cells;
			cputime += s.time;
			if (s.bestsolpointnumber){
				Qoct=s.bestsolpointnumber+s.epsobj;
				bestsol=s.bestsolpoint;
				bestsolpointnumber=s.bestsolpointnumber;
			}
			if(trace){
				s.report_possible_inliers();
				s.report_solution();
			}

			end = clock();
			totaltime += ((double)(end)-(double)(start))/CLOCKS_PER_SEC;
			start= clock();
			delete bs;
		}

		for (int i=0; i<p; i++){
			delete m_func[i];
			delete &m_ctc[i];
		}
		
		end = clock();
		cout << "Shape extraction : OK. Time : " << ((double)(end)-(double)(start0))/CLOCKS_PER_SEC << " seconds" << endl;
		cout << "\ttotal time " << totaltime << endl;
		cout << "\tcpu time " << cputime << endl;
		cout << "\ttotal branch number " << nb_cells << endl;

		for (int i=0; i<p; i++)
			delete [] linfun[i];
		delete [] linfun;
	}catch(ibex::SyntaxError& e){
		cout << e << endl;
	}
}
