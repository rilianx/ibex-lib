#include "ibex.h"

using namespace std;
using namespace ibex;

bool compare(const std::pair<double, std::pair<int,int> >&i, const std::pair<double, std::pair<int,int> >&j){
  return i.first < j.first;
}

void interval_comp(int eq, int var, IntervalMatrix A, IntervalVector & x, IntervalVector b){

	Interval value(b[eq]);
	for (int i = 0; i < A.nb_cols(); i++){
		if (i != var){
			value = value +(-A[eq][i]*x[i]);
		}
	}
	value = value/A[eq][var];
	x[var] = value&=x[var];

}

void list_strategy1(vector<pair <int,int> >& order, IntervalMatrix A, IntervalVector x){

	order.clear();
	vector< pair< double, pair<int,int> > > order_aux;
	for (int j = 0 ; j < A.nb_cols() ; j++){
		double max = -1;
		int pos = -1;
		for (int i = 0 ; i < A.nb_rows() ; i++){
			double calc = std::abs(A[i][j].mid())*x[j].diam();
			if (calc > max){
				max = calc;
				pos = i;
			}
		}
		order_aux.push_back(make_pair(max,make_pair(pos,j)));
	}
	std::sort(order_aux.begin(), order_aux.end(), compare);
	for (int i = 0 ; i < order_aux.size() ; i++)
		order.push_back(order_aux[i].second);
	for (int i = 0 ; i < order_aux.size() ; i++)
		cout << order_aux[i].second.first << "     "  << order_aux[i].second.second <<  endl;
}
/*
 * pos : equation
 * j : variable
 */
void list_strategy2(vector< pair<int,int> >& order, IntervalMatrix A, IntervalVector x){

	order.clear();
	double impact_values[A.nb_rows()][A.nb_cols()];
	double impact_values_rel[A.nb_rows()][A.nb_cols()];
	vector< pair< double, pair<int,int> > > order_aux;


	for (int j = 0 ; j < A.nb_cols() ; j++)
		for (int i = 0 ; i < A.nb_rows() ; i++)
			impact_values[i][j] = std::abs(A[i][j].mid())*x[j].diam();
	for (int i = 0 ; i < A.nb_rows() ; i++){
		double sum = 0;
		for (int j = 0 ; j < A.nb_cols() ; j++)
			sum = sum+ impact_values[i][j] ;
		for (int j = 0 ; j < A.nb_cols() ; j++)
			impact_values[i][j] = impact_values[i][j]/sum;
	}

	for (int j  = 0; j < A.nb_cols() ; j++){
		double max = -1;
		int pos = -1;
		for (int i = 0 ; i < A.nb_rows() ; i++){
			if (impact_values[i][j] > max){
				pos = i;
				max = impact_values[i][j];
			}
		}
		order_aux.push_back(make_pair(max,make_pair(pos,j)));
	}
	std::sort(order_aux.begin(), order_aux.end(), compare);
	for (int i = 0 ; i < order_aux.size() ; i++)
		order.push_back(order_aux[i].second);
	for (int i = 0 ; i < order_aux.size() ; i++)
		cout << order_aux[i].second.first << "     "  << order_aux[i].second.second <<  endl;
}

void list_strategy3(vector< pair<int,int> >& order, IntervalMatrix A, IntervalVector x){

	order.clear();
	double impact_values[A.nb_rows()][A.nb_cols()];
	double impact_values_rel[A.nb_rows()][A.nb_cols()];
	vector< pair< double, pair<int,int> > > order_aux;


	for (int j = 0 ; j < A.nb_cols() ; j++)
		for (int i = 0 ; i < A.nb_rows() ; i++)
			impact_values[i][j] = std::abs(A[i][j].mid())*x[j].diam();
	for (int i = 0 ; i < A.nb_rows() ; i++){
		double sum = 0;
		for (int j = 0 ; j < A.nb_cols() ; j++)
			sum = sum+ impact_values[i][j] ;
		for (int j = 0 ; j < A.nb_cols() ; j++)
			impact_values[i][j] = impact_values[i][j]/sum;
	}

	for (int j  = 0; j < A.nb_cols() ; j++){
		double max = -1;
		double sum = 0;
		int pos = -1;
		for (int i = 0 ; i < A.nb_rows() ; i++){
			sum = sum + impact_values[i][j];
			if (impact_values[i][j] > max){
				pos = i;
				max = impact_values[i][j];
			}
		}
		order_aux.push_back(make_pair(sum,make_pair(pos,j)));
	}
	std::sort(order_aux.begin(), order_aux.end(), compare);
	for (int i = 0 ; i < order_aux.size() ; i++)
		order.push_back(order_aux[i].second);
	for (int i = 0 ; i < order_aux.size() ; i++)
		cout << order_aux[i].second.first << "     "  << order_aux[i].second.second <<  endl;


}

int main() {
	/*EXAMPLE 1*/
	IntervalMatrix A(2,2);
	IntervalVector x(2);
	IntervalVector b(2);
	A[0][0] = Interval(4,7);
	A[0][1] = Interval(3,5);
	A[1][0] = Interval(0,4);
	A[1][1] = Interval(1,2);
	x[0] = Interval(0,8);
	x[1] = Interval(0,4);
	b[0] = Interval(1,5);
	b[1] = Interval(-16,-7);
	cout << A << endl << endl;
	cout << A.mid() << endl << endl;
	Matrix B(A.mid());
	real_inverse(B,B);
	cout << (B*A) << endl << endl;
	cout << (B*A).mid() << endl << endl;
	/*EXAMPLE 2*/
//	IntervalMatrix A(2,3);
//	IntervalVector x(3);
//	IntervalVector b(2);
//	A[0][0] = Interval(4,4);
//	A[0][1] = Interval(-0.1,-0.1);
//	A[0][2] = Interval(0.2,0.2);
//	A[1][0] = Interval(0.1,0.1);
//	A[1][1] = Interval(7,7);
//	A[1][2] = Interval(4,4);
//	x[0] = Interval(-10,10);
//	x[1] = Interval(-4,3);
//	x[2] = Interval(-6,6);
//	b[0] = Interval(-1,3);
//	b[1] = Interval(-16,-7);
	/*EXAMPLE 2 extended*/
//	IntervalMatrix A(2,5);
//	IntervalVector x(5);
//	IntervalVector b(2);
//	A[0][0] = Interval(4,4);
//	A[0][1] = Interval(-0.1,-0.1);
//	A[0][2] = Interval(0.2,0.2);
//	A[0][3] = Interval(-1,-1);
//	A[0][4] = Interval(0,0);
//	A[1][0] = Interval(0.1,0.1);
//	A[1][1] = Interval(7,7);
//	A[1][2] = Interval(4,4);
//	A[1][3] = Interval(0,0);
//	A[1][4] = Interval(-1,-1);
//	x[0] = Interval(-10,10);
//	x[1] = Interval(-4,3);
//	x[2] = Interval(-6,6);
//	x[3] = Interval(-1,3);
//	x[4] = Interval(-16,-7);
//	b[0] = Interval(0,0);
//	b[1] = Interval(0,0);
//	cout << A.mid() << endl << endl;
//	cout << x << endl<< endl;
//	bool cambio = true;
//	int i = 1;
//	vector< pair<int,int> > order;
//	while(cambio){
//		cout << "iteracion " << i << endl;
//		i++;
//		IntervalVector aux_aux=x;
//		list_strategy2(order, A, x);
//		for (int j = 0 ; j < order.size() ; j++)
//			interval_comp(order[j].first,order[j].second,A,x,b);
//		cout << x << endl;
//		if (aux_aux == x) cambio = false;
//	}
//	cout << x << endl;
}
