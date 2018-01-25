#include "ibex.h"

using namespace std;
using namespace ibex;

int main() {

	/*EXAMPLE 1*/
//	IntervalMatrix A(2,4);
//	IntervalVector x(4);
//	IntervalVector b(2);
//	A[0][0] = Interval(1,7);
//	A[0][1] = Interval(4,5);
//	A[1][0] = Interval(3,4);
//	A[1][1] = Interval(9,10);
//	x[0] = Interval(0,8);
//	x[1] = Interval(0,4);
//	b[0] = Interval(1,5);
//	b[1] = Interval(-16,-7);
	/*EXAMPLE 2*/
//	IntervalMatrix A(2,4);
//	IntervalVector x(4);
//	IntervalVector b(2);
//	A[0][0] = Interval(5,6);
//	A[0][1] = Interval(-8,-7);
//	A[0][2] = Interval(0,5);
//	A[0][3] = Interval(-0.9,-0.1);
//	A[1][0] = Interval(-0.9,3);
//	A[1][1] = Interval(9,10);
//	A[1][2] = Interval(8,9);
//	A[1][3] = Interval(-4,0);
//	x[0] = Interval(0,8);
//	x[1] = Interval(-1,2);
//	x[2] = Interval(1,5);
//	x[3] = Interval(-4,4);
//	b[0] = Interval(-1,1);
//	b[1] = Interval(-4,-1);
	/*EXAMPLE 3*/
//	IntervalMatrix A(2,6);
//	IntervalVector x(6);
//	IntervalVector b(2);
//	A[0][0] = Interval(1,1);
//	A[0][1] = Interval(2,2);
//	A[0][2] = Interval(3,3);
//	A[0][3] = Interval(4,4);
//	A[1][4] = Interval(-1,-1);
//	A[0][5] = Interval(0,0);
//	A[1][0] = Interval(5,5);
//	A[1][1] = Interval(6,6);
//	A[1][2] = Interval(7,7);
//	A[1][3] = Interval(8,8);
//	A[1][4] = Interval (0,0);
//	A[1][5] = Interval(-1,-1);
//	x[0] = Interval(-15,9);
//	x[1] = Interval(-6,33);
//	x[2] = Interval(0,18);
//	x[3] = Interval(-6,10);
//	x[4] = Interval(-1,3);
//	x[5] = Interval(-16,-7);
//	b[0] = Interval(0,0);
//	b[1] = Interval(0,0);
	/*EXAMPLE 4*/
//	IntervalMatrix A(4,5);
//	IntervalVector x(5);
//	IntervalVector b(4);
//	A[0][0] = Interval(2,2);
//	A[0][1] = Interval(1,1);
//	A[0][2] = Interval(1,1);
//	A[0][3] = Interval(1,1);
//	A[0][4] = Interval(1,1);
//	A[1][0] = Interval(1,1);
//	A[1][1] = Interval(2,2);
//	A[1][2] = Interval(1,1);
//	A[1][3] = Interval(1,1);
//	A[1][4] = Interval(1,1);
//	A[2][0] = Interval(1,1);
//	A[2][1] = Interval(1,1);
//	A[2][2] = Interval(2,2);
//	A[2][3] = Interval(1,1);
//	A[2][4] = Interval(1,1);
//	A[3][0] = Interval(1,1);
//	A[3][1] = Interval(1,1);
//	A[3][2] = Interval(1,1);
//	A[3][3] = Interval(2,2);
//	A[3][4] = Interval(1,1);
//	x[0] = Interval(-100,100);
//	x[1] = Interval(-100,100);
//	x[2] = Interval(-100,100);
//	x[3] = Interval(-100,100);
//	x[4] = Interval(-100,100);
//	b[0] = Interval(6,6);
//	b[1] = Interval(6,6);
//	b[2] = Interval(6,6);
//	b[3] = Interval(6,6);
	/*EXAMPLE 5*/
	IntervalMatrix A(2,3);
	IntervalVector x(3);
	IntervalVector b(2);
	A[0][0] = Interval(4,4);
	A[0][1] = Interval(-0.1,-0.1);
	A[0][2] = Interval(0.2,0.2);
	A[1][0] = Interval(0.1,0.1);
	A[1][1] = Interval(7,7);
	A[1][2] = Interval(4,4);
	x[0] = Interval(-10,10);
	x[1] = Interval(-4,3);
	x[2] = Interval(-6,6);
	b[0] = Interval(-1,3);
	b[1] = Interval(-16,-7);

	std::vector< std::vector <int> > comb_piv;
	/*get all possible combinations of pivots*/
//	bwd_mul(b,A,x,1e-8);
//		gauss_seidel(A,b,x,1e-8);
//	cout << x << endl;
	IntervalVector y(x);
	vector<Matrix> permutations;
	gauss_jordan_all(A,permutations,1e-8);
	cout << "Matrix A: " << endl;
	cout << A.mid() << endl << endl;
	cout << "Vector x: " << endl;
	cout << x << endl << endl;
	cout << "Matrix b: " << endl;
	cout << b << endl << endl;
	cout << "#P matrices: " << permutations.size() << endl;
		IntervalVector aux = x;
		for (int i = 0 ; i < permutations.size()  ; i++){
			Matrix P = permutations[i];
			IntervalMatrix PA = P*A;
			IntervalVector Pb = P*b;
			cout << "Using the Matrix P: " << i+1 << endl;
			cout<< "%%%%%%%%"<< endl;
			cout << "PA:" << endl;
			cout << PA.mid() << endl;
			cout << "before:" << endl;
			cout << x << endl;
//			gauss_seidel(PA,Pb,x,0.001);
			bwd_mul(Pb,PA,x,1e8);
			cout << "after:" << endl;
			cout << x << endl;
			cout << endl;
		}
		cout << "%%&&&&%%&%&%&%&&%&&%" << endl;
		Matrix P(best_gauss_jordan (A,1e-8));
		IntervalMatrix PA = P*A;
		IntervalVector Pb = P*b;
		cout << PA.mid() << endl << endl;
		bwd_mul(Pb,PA,y,1e-8);
//		gauss_seidel(PA,Pb,y,0.001);
		cout << y << endl << endl;
}
