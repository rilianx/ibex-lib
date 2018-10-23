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
        cout << "maxsat" << endl;
    }catch(ibex::SyntaxError& e){
		cout << e << endl;
	}
}