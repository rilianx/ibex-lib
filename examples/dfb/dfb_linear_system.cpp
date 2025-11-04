#include "ibex.h"
#include "ibex_CtcDFBPropag.h"
#include "ibex_CtcDFBManager.h"

using namespace std;
using namespace ibex;

#include <fstream>
#include <iostream>

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

    //System (argv)

 //   std::ofstream file("resultado_oct25_06.txt"); // crea o sobrescribe el archivo
   // std::streambuf* coutbuf = std::cout.rdbuf(); // guarda el buffer original
    //std::cout.rdbuf(file.rdbuf()); // redirige cout al archivo



    if (argc<2) {
        cerr << "Usage: " << argv[0] << " example_standar_dfb_1" << endl;
        return 1;
    }
    
    disableCout();
    const char* filename = argv[1];
    System sys_hc4(filename);
    System sys_dfb_hc4(filename);
    System sys_dfb(filename);
    System sys_dfb2(filename);
    System sys_dfb_only(filename);
    System sys_ph(filename);
    System sys_it_ph(filename);
    System sys_it_ph2(filename);
    System sys_acid_dfb(filename);
    System sys_acid_ph(filename);


   	double eqeps01=0.0;
    ExtendedSystem sys2_hc4(sys_hc4,eqeps01);

    double eqeps02=0.0;
    ExtendedSystem sys2_dfb_hc4(sys_dfb_hc4,eqeps02);

    double eqeps03=0.0;
    ExtendedSystem sys2_dfb(sys_dfb,eqeps03);

    double eqeps04=0.0;
    ExtendedSystem sys2_dfb2(sys_dfb2,eqeps04);

    double eqeps05=0.0;
    ExtendedSystem sys2_dfb_only(sys_dfb_only,eqeps05);

    double eqeps06=0.0;
    ExtendedSystem sys2_ph(sys_ph,eqeps06);
 
    double eqeps07=0.0;
    ExtendedSystem sys2_it_ph(sys_it_ph,eqeps07);

    double eqeps08=0.0;
    ExtendedSystem sys2_it_ph2(sys_it_ph2,eqeps08);

    double eqeps09=0.0;
    ExtendedSystem sys2_acid_dfb(sys_acid_dfb,eqeps09);

    double eqeps10=0.0;
    ExtendedSystem sys2_acid_ph(sys_acid_ph,eqeps10);


    LinearizerXTaylor lr_hc4(sys2_hc4, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);

    LinearizerXTaylor lr_dfb_hc4(sys2_dfb_hc4, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);

    LinearizerXTaylor lr_dfb(sys2_dfb, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);

    LinearizerXTaylor lr_dfb2(sys2_dfb2, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);

    LinearizerXTaylor lr_dfb_only(sys2_dfb_only, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);

    LinearizerXTaylor lr_ph(sys2_ph, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);

    LinearizerXTaylor lr_it_ph(sys2_it_ph, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);

    LinearizerXTaylor lr_it_ph2(sys2_it_ph2, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);

    LinearizerXTaylor lr_acid_dfb(sys2_acid_dfb, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);

    LinearizerXTaylor lr_acid_ph(sys2_acid_ph, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);

    CtcDFBPropag hc4(sys2_hc4, lr_hc4, 0.001, true ,true); // HC4 standalone
    CtcDFBPropag dfb_hc4(sys2_dfb_hc4, lr_dfb_hc4, 0.001, true ,true); // DFB + HC4

    CtcDFBPropag dfb(sys2_dfb, lr_dfb, 0.01); //DFB standalone
    CtcDFBPropag dfb2(sys2_dfb2, lr_dfb2, 0.01); //DFB2 standalone

    CtcDFBPropag dfb_only(sys2_dfb_only, lr_dfb_only, 0.01); //DFB_only standalone

    CtcPolytopeHull ctc_ph(lr_ph);
    ctc_ph.n_soplex_iterations = 0;

    CtcPolytopeHull ctc_it_ph(lr_it_ph);

    CtcDFBPropag hc4_it_ph2(sys2_hc4, lr_hc4, 0.001, true ,true); 
    CtcPolytopeHull ctc_it_ph2(lr_it_ph2);


    CtcDFBPropag dfb_acid(sys2_acid_dfb, lr_acid_dfb, 0.01, false);
    CtcAcid acid_dfb = CtcAcid(sys2_acid_dfb, dfb_acid, true);
    CtcDFBManager dfb_acid_manager(dfb_acid, acid_dfb, lr_acid_dfb);



    CtcDFBPropag hc4_acid_ph(sys2_acid_ph, lr_acid_ph, 0.001, true ,true);

    CtcAcid acid_ph_dfb = CtcAcid(sys2_acid_ph, hc4_acid_ph, true);
    CtcPolytopeHull ctc_acid_ph(lr_acid_ph);

    //cout <<"x,y"<<endl;
    //cout << "0," << sys2.box.perimeter() << endl;



    hc4.contract(sys2_hc4.box); 
    dfb_hc4.contract(sys2_dfb_hc4.box);
    dfb.contract(sys2_dfb.box);

    dfb2.contract(sys2_dfb2.box);
    CtcDFBPropag dfb2_2(sys2_dfb2, lr_dfb2, 0.001); //standalone
    dfb2_2.contract(sys2_dfb2.box);

    dfb_only.contract(sys2_dfb_only.box);


    ctc_ph.contract(sys2_ph.box); //linearize and contract


    ctc_it_ph.n_soplex_iterations = 0;
    //cout << "Using Iterated PolyHull (fixpoint)" << endl; 
    double old_per = sys2_it_ph.box.perimeter();
    ctc_it_ph.contract(sys2_it_ph.box); //linearize and contract
    int iter_it_ph = 0;
    while(true){ 
        hc4.contract(sys2_it_ph.box); 
        double new_per = sys2_it_ph.box.perimeter();
        ctc_ph.primal_sol_found.clear();
        ctc_ph.optimizer(sys2_it_ph.box); //only contract
        if (fabs(new_per-old_per)<0.1) break;
        old_per = new_per;
        iter_it_ph++;
    }


    ctc_it_ph2.n_soplex_iterations = 0;
    //cout << "Using Iterated PolyHull (fixpoint)" << endl; 
    double old_per02 = sys2_it_ph.box.perimeter();
    ctc_it_ph2.contract(sys2_it_ph2.box); //linearize and contract
    int iter_it_ph2 = 0;
    while(true){ 
        hc4_it_ph2.contract(sys2_it_ph2.box);
        double new_per02 = sys2_it_ph2.box.perimeter();
        ctc_ph.primal_sol_found.clear();
        ctc_ph.contract(sys2_it_ph2.box);
        if (fabs(new_per02-old_per02)<0.1) break;
        old_per02 = new_per02;
        iter_it_ph2++;
    }


    dfb_acid_manager.contract(sys2_acid_dfb.box);


    acid_ph_dfb.contract(sys2_acid_ph.box); 
    ctc_acid_ph.contract(sys2_acid_ph.box); 

    acid_ph_dfb.contract(sys2_acid_ph.box); 
    ctc_acid_ph.optimizer(sys2_acid_ph.box); 

   enableCout();
    std::cout << "\n xFinal HC4 = \n" << sys2_hc4.box << std::endl;
    std::cout << "\n xFinal DFB + HC4= \n" << sys2_dfb_hc4.box << std::endl;

    std::cout << "\n xFinal DFB = \n" << sys2_dfb.box << std::endl;
    std::cout << "\n xFinal DFB2 = \n" << sys2_dfb2.box << std::endl;

    std::cout << "\n xFinal DFB_ONLY = \n" << sys2_dfb_only.box << std::endl;

    std::cout << "\nxFinal PH = \n" << sys2_ph.box << std::endl;

    std::cout << "\nxFinal IT_PH = \n" << sys2_it_ph.box << std::endl;
    std::cout << "\nxFinal IT_PH2 = \n" << sys2_it_ph2.box << std::endl;

    std::cout << "\nxFinal ACID_DFB = \n" << sys2_acid_dfb.box << std::endl;
    std::cout << "\nxFinal ACID_PH = \n" << sys2_acid_ph.box << std::endl;
    std::cout << std::endl;


    // 🔽 Cuando termines, restauras cout
//    std::cout.rdbuf(coutbuf);
//    std::cout << "Se guardó resultado.txt correctamente." << std::endl;
}