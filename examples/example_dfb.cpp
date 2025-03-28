#include "ibex.h"
#include "ibex_CtcAllPropag.h"

using namespace std;
using namespace ibex;


int main(int argc, char** argv) {
    // System (argv)
    if (argc<2) {
        cerr << "Usage: " << argv[0] << " example_standar_dfb_1" << endl;
        return 1;
    }

    const char* filename = argv[1];
    System sys(filename);

   	double eqeps=0.0;
    ExtendedSystem sys2(sys,eqeps);
    //cout << sys << endl;


    LinearizerXTaylor lr(sys2, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);
    cout << sys2 << endl;

    CtcHC4 hc4(sys2.ctrs,0.01,true);
    hc4.contract(sys2.box); 

    if(argc>2){ //dfb or ph
        if (strcmp(argv[2], "dfb") == 0) {
            CtcAllPropag ctc(sys2, lr);
            cout << "Using DFB" << endl;
            ctc.contract(sys2.box);
        } else if (strcmp(argv[2], "ph") == 0) {
            CtcPolytopeHull ctc_ph(lr);
            cout << "Using PolyHull" << endl; 
            hc4.contract(sys2.box); 
            ctc_ph.contract(sys2.box); //linearize and contract

            //ctc_ph.optimizer(sys2.box); //only contract
        } else {
            cout << "Unknown option. Use 'dfb' or 'ph'." << endl;
            return 1;   
        }
    }

    cout << "xFinal = " << sys2.box << endl;
}