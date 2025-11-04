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

    LinearizerXTaylor lr(sys2, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);
    //cout << sys2 << endl;

    disableCout();
    CtcDFBPropag hc4(sys2, lr, 0.001, true ,true);
    hc4.contract(sys2.box);  // ¿POR QUÉ EXISTE UNA CONTRACCIÓN PREVIA?
    
    enableCout();
    //cout <<"x,y"<<endl;
    //cout << "0," << sys2.box.perimeter() << endl;

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
            //cout << "Using DFB" << endl;
            disableCout();
            CtcDFBPropag dfb(sys2, lr, 0.01); //standalone
            dfb.contract(sys2.box);
            enableCout();
            //cout << "count_dfb=" << dfb.count_dfb << " count_hc4=" << dfb.count_hc4 << endl;
        //    for(auto it : dfb.history)
//                cout << it.first<<"," << it.second << endl;

        }  else if (strcmp(argv[2], "dfb2") == 0) {
            //cout << "Using DFB" << endl;
            CtcDFBPropag::b_contraction = true;
            disableCout();
            CtcDFBPropag dfb(sys2, lr, 0.001); //standalone
            dfb.contract(sys2.box);
            CtcDFBPropag dfb2(sys2, lr, 0.001); //standalone
            dfb2.contract(sys2.box);
            enableCout();
            //cout << "count_dfb=" << dfb.count_dfb << " count_hc4=" << dfb.count_hc4 << endl;
  //          for(auto it : dfb.history)
//                cout << it.first<<"," << it.second << endl;
      //      for(auto it : dfb2.history)
    //            cout << it.first + dfb.count_dfb<<"," << it.second << endl;

        } else if (strcmp(argv[2], "dfb_only") == 0) {
            disableCout();
            
            CtcDFBPropag dfb(sys2, lr, 0.01, true ,false, true);
            // CtcDFBPropag hc4(sys2, lr, 0.001, true ,true);
            //cout << "Using DFB(no HC4)" << endl;
            dfb.contract(sys2.box);
            enableCout();

         //   hc4.contract(sys2.box); 
        } else if (strcmp(argv[2], "ph") == 0) {
            CtcPolytopeHull ctc_ph(lr);
            //cout << "Using PolyHull" << endl; 
            ctc_ph.n_soplex_iterations = 0;
            ctc_ph.contract(sys2.box); //linearize and contract

          //  for(auto it : ctc_ph.history)
            //    cout << it.first<<"," << it.second << endl;


        } else if (strcmp(argv[2], "it_ph") == 0) { //iterated ph
            CtcPolytopeHull ctc_ph(lr);
            ctc_ph.n_soplex_iterations = 0;
            //cout << "Using Iterated PolyHull (fixpoint)" << endl; 
            double old_per = sys2.box.perimeter();
            ctc_ph.contract(sys2.box); //linearize and contract

            int iter = 0;
            while(true){ 

                disableCout();
                hc4.contract(sys2.box); 
                enableCout();
                double new_per = sys2.box.perimeter();
                ctc_ph.primal_sol_found.clear();
                ctc_ph.optimizer(sys2.box); //only contract
                if (fabs(new_per-old_per)<0.1) break;
                old_per = new_per;
                iter++;
            }
            
            //for(auto it : ctc_ph.history)
              //  cout << it.first<<"," << it.second << endl;
            

        } else if (strcmp(argv[2], "it_ph2") == 0) { //iterated ph
            CtcPolytopeHull ctc_ph(lr);
            ctc_ph.n_soplex_iterations = 0;
            //cout << "Using Iterated PolyHull (fixpoint)" << endl; 
            double old_per = sys2.box.perimeter();
            ctc_ph.contract(sys2.box); //linearize and contract

            int iter = 0;
            while(true){ 

                disableCout();
                hc4.contract(sys2.box); 
                enableCout();
                double new_per = sys2.box.perimeter();
                ctc_ph.primal_sol_found.clear();
                ctc_ph.contract(sys2.box); 
                if (fabs(new_per-old_per)<0.1) break;
                old_per = new_per;
                iter++;
            }
            
        //    for(auto it : ctc_ph.history)
          //      cout << it.first<<"," << it.second << endl;
            

        } else if (strcmp(argv[2], "acid_dfb") == 0) {
            disableCout();

            CtcDFBPropag dfb(sys2, lr, 0.01, false);
            CtcAcid acid_dfb = CtcAcid(sys2, dfb, true);
            CtcDFBManager dfb_manager(dfb, acid_dfb, lr);

            dfb_manager.contract(sys2.box);
            enableCout();
        }else if (strcmp(argv[2], "acid_ph") == 0) {
            disableCout();

            CtcAcid acid_dfb = CtcAcid(sys2, hc4, true);
            CtcPolytopeHull ctc_ph(lr);
            cout << "Using PolyHull" << endl; 
            acid_dfb.contract(sys2.box); 
            ctc_ph.contract(sys2.box); //linearize and contract
            cout << "xInter = " << sys2.box << endl;
            acid_dfb.contract(sys2.box); 
            ctc_ph.optimizer(sys2.box); //only contract
            enableCout();

        }
    }
    //enableCout();

    cout << "xFinal = " << sys2.box << endl;
    
    //cout << "xFinal per = " << sys2.box.perimeter() << endl;
}