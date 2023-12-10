#include "ibex.h"
#include "affine_arithm.h"
#include "ibex_AffineEval.h"



using namespace std;
using namespace ibex;


int main() {
    // Intervals
	Interval x(1,2);
	Interval y(3,4);
	cout << "x+y=" << x+y << endl; // display x+y=[4, 6]

    //Functions
    Function f("x","y", "((sin(x))^2 + 3*y)-1");
    IntervalVector box(2,Interval(0,2));
    cout << f << endl; // display _f_0:(x,y)->((sin(x)^2+(3*y))-1)
    cout << "f("<< box << ")=" << f.eval(box) << endl; // display f(([0, 2] ; [0, 2]))=[-1, 6]

    //Partial Derivatives and gradients
    cout << "gradient of f at " << box.mid() << " is " << f.gradient(box.mid()) << endl;
    //display gradient of f at (1 ; 1) is ([0.9092974268256812, 0.9092974268256822] ; <3, 3>)

    Function g = f.diff(); 
    cout << "g(x,y)=" << g << endl; // display g(x,y)=((2*sin(x)*cos(x)), 3)
    cout << g.eval_vector(IntervalVector(2,Interval(1,2))) << endl;
    // display ([-0.832293673094285, 1.08060461173628] ; <3, 3>)

    //Contractor
    CtcFwdBwd ctc(f, LEQ);

    IntervalVector x1(2);
    x1[0] = Interval(-1,1);
    x1[1] = Interval(0,6);

    ctc.contract(x1);
    cout <<"The resulting box:" << x1 << endl;
    // display The resulting box:([-1, 1] ; [0, 0.3333333333333334])

    IntervalVector x2(1, Interval(1,2));
    Affine af(x2,0);

    Affine af2 = 5.*pow(af,2) - af;
    cout << af2 << endl;
    cout << af2.ev << endl;

    Function f3("x", "5*x^2-x");
    AffineEval af_eval(f, x2.size());
    cout <<"The resulting box:" << af_eval.eval(x2) << endl;
    
}