#include "ibex.h"
#include "ibex_CtcAllPropag.h"

using namespace std;
using namespace ibex;


int main(int argc, char** argv) {
    // Intervals
    Interval x_0 = Interval(-1.565, 2.880);
    Interval x_1 = Interval(-0.478, 4.463);
    Interval x_2 = Interval(-1.038, 6.032);
    Interval x_3 = Interval(-0.048, 3.615);
    Interval x_4 = Interval(-1.076, 2.647);
//	cout << "x_0 + x_1 + x_2 + x_3 + x_4 = " << x_0 + x_1 + x_2 + x_3 + x_4 << endl;

    IntervalVector x = IntervalVector({x_0, x_1, x_2, x_3, x_4});
	cout << "xOG = " << x << endl;

    // System (argv)
    if (argc!=2) {
        cerr << "Usage: " << argv[0] << "standar_system_dfb" << endl;
        return 1;
    }

    const char* filename = argv[1];
    System sys(filename);

   	double eqeps=0.0;
    ExtendedSystem sys2(sys,eqeps);
    //cout << sys << endl;


    LinearizerXTaylor lr(sys2, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM, LinearizerXTaylor::HANSEN);
    cout << sys2 << endl;

    CtcAllPropag ctc(sys2, lr);
    cout << "ctc created" << endl;
//    cout << sys2.box << endl;
    ctc.contract(x);
 //   cout << sys2.box << endl;

    cout << "xFinal = " << sys2.box << endl;
}