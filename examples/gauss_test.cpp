#include "ibex.h"

using namespace std;
using namespace ibex;

int main() {
	/*Input matrix define*/
	IntervalMatrix A(2,4);
	IntervalVector x(4);
	IntervalVector b(2);
	A[0][0] = Interval(4,4);
	A[0][1] = Interval(-0.1,-0.1);
	A[0][2] = Interval(0.2,0.2);
	A[0][3] = Interval(-0.6,-0.6);
	A[1][0] = Interval(0.1,0.1);
	A[1][1] = Interval(7,7);
	A[1][2] = Interval(0.3,0.3);
	A[1][3] = Interval(-5,-5);
	x[0] = Interval(0,8);
	x[1] = Interval(0,4);
	x[2] = Interval(-8,8);
	x[3] = Interval(0,1);
	b[0] = Interval(-1,3);
	b[1] = Interval(-16,-7);
	std::vector< std::vector <int> > comb_piv;
	/*get all possible combinations of pivots*/
//	gauss_seidel(A,b,x,1e-7);
//	cout << x << endl;
	vector<Matrix> permutations;
	gauss_jordan_all(A,permutations,1e-8);
	cout << "Matrix A: " << endl;
	cout << A.mid() << endl << endl;
	cout << "Vector x: " << endl;
	cout << x << endl << endl;
	cout << "Matrix b: " << endl;
	cout << b << endl << endl;
	cout << "#P matrices: " << permutations.size() << endl;
	bool change = true;
//	while (!(x.is_empty()) && (change) ){
		change = false;
		IntervalVector aux = x;
		cout << "&&&&&&&&&&&&&&&&&& NEW ITERATION &&&&&&&&&&&&&&&&&" << endl;
		for (int i = 0 ; i < permutations.size() ; i++){
			Matrix P = permutations[i];
			IntervalMatrix PA = P*A;
			IntervalVector Pb = P*b;
			cout << "Using the Matrix P: " << i+1 << endl;
			cout<< "%%%%%%%%"<< endl;
			cout << "PA:" << endl;
			cout << PA.mid() << endl;
			cout << "before:" << endl;
			cout << x << endl;
//			gauss_seidel(PA,Pb,x,1e8);
			bwd_mul(Pb,PA,x,1e8);
			cout << "after:" << endl;
			cout << x << endl;
			cout << endl;
		}
//		if (aux != x) change = true;
//	}
}
