#include "ibex.h"
#include "ibex_CtcAllPropag.h"

using namespace std;
using namespace ibex;


int main(int argc, char** argv) {
    // Intervals
	Interval x(1,2);
	Interval y(3,4);
	cout << "x+y=" << x+y << endl; // display x+y=[4, 6]

    
    // System (argv)
    if (argc!=2) {
        cerr << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }

    const char* filename = argv[1];
    System sys(filename);

   	double eqeps=0.0;
    ExtendedSystem sys2(sys,eqeps);
    //cout << sys << endl;


    LinearizerXTaylor lr(sys2, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);

    CtcAllPropag ctc(sys2, lr);
    cout << "ctc created" << endl;

    cout << sys2.box << endl;
    ctc.contract(sys2.box);
    cout << sys2.box << endl;

    cout << "Hello, World!" << endl;

}