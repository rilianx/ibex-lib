#include "ibex.h"
#include "ibex_CtcDFBPropag.h"
#include "ibex_CtcDFBManager.h"

using namespace std;
using namespace ibex;

std::streambuf* cout_backup = nullptr;
std::ofstream null_stream;

void disableCout() {
    if (!null_stream.is_open())
        null_stream.open("/dev/null"); // En Windows: "NUL"
    cout_backup = std::cout.rdbuf();  // Guardás el original
    std::cout.rdbuf(null_stream.rdbuf()); // Redirigís a la nada
}

void enableCout() {
    std::cout.rdbuf(cout_backup); // Restaurás el original
}

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

    //disableCout();
    if(argc>2){ //dfb or ph

        if (strcmp(argv[2], "hc4") == 0) {
            cout << "Using HC4" << endl;
            hc4.contract(sys2.box); 
        } else if (strcmp(argv[2], "dfb_hc4") == 0) {
            CtcDFBPropag ctc(sys2, lr, 0.01, true ,true);
            cout << "Using DFB (HC4)" << endl;
            ctc.contract(sys2.box);
        } else if (strcmp(argv[2], "dfb") == 0) {
            CtcDFBPropag ctc(sys2, lr);
            cout << "Using DFB(HC4 only)" << endl;
            ctc.contract(sys2.box);
        } else if (strcmp(argv[2], "dfb_only") == 0) {
            CtcDFBPropag ctc(sys2, lr, 0.01, true ,false, true);
            CtcDFBPropag hc4(sys2, lr, 0.001, true ,true);
            cout << "Using DFB(no HC4)" << endl;
            ctc.contract(sys2.box);
            hc4.contract(sys2.box);

        } else if (strcmp(argv[2], "ph") == 0) {
            CtcPolytopeHull ctc_ph(lr);
            CtcDFBPropag hc4(sys2, lr, 0.001, true ,true);
            cout << "Using PolyHull" << endl; 
            //ctc_ph.contract(sys2.box); //linearize and contract
            ctc_ph.contract(sys2.box); //linearize and contract
            //cout << "xInter = " << sys2.box << endl;
            cout << "xFinal per = " << sys2.box.perimeter() << endl;
            hc4.contract(sys2.box); 
                cout << "xFinal per = " << sys2.box.perimeter() << endl;
            ctc_ph.primal_sol_found.clear();
            ctc_ph.optimizer(sys2.box); //only contract
        } else if (strcmp(argv[2], "acid_dfb") == 0) {

            CtcDFBPropag dfb(sys2, lr, 0.01, false);
            CtcAcid acid_dfb = CtcAcid(sys2, dfb, true);
            CtcDFBManager dfb_manager(dfb, acid_dfb, lr);

            dfb_manager.contract(sys2.box);

        }else if (strcmp(argv[2], "acid_ph") == 0) {
            CtcAcid acid_dfb = CtcAcid(sys2, hc4, true);
            CtcPolytopeHull ctc_ph(lr);
            cout << "Using PolyHull" << endl; 
            acid_dfb.contract(sys2.box); 
            ctc_ph.contract(sys2.box); //linearize and contract
            cout << "xInter = " << sys2.box << endl;
            acid_dfb.contract(sys2.box); 
            ctc_ph.optimizer(sys2.box); //only contract

        }
    }
    //enableCout();

    cout << "xFinal = " << sys2.box << endl;
    cout << "xFinal per = " << sys2.box.perimeter() << endl;
}