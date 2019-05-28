/*
 * ibex_Conditioners.cpp
 *
 *  Created on: 13-10-2017
 *      Author: victor
 */

#include "ibex_Conditioners.h"

using namespace std;

namespace ibex {

	bool compare(const std::pair<double, std::pair<int,int> >&i, const std::pair<double, std::pair<int,int> >&j){
		return i.first > j.first;
	}

	bool compare2(const std::pair<double, Vector >&i, const std::pair<double, Vector >&j){
		return i.first < j.first;
	}


	void combinatorial(IntervalMatrix A, int cols,int rows,std::vector< std::vector <int> > & comb_piv){
		vector<int> pivots;
		/*Initialize the possible pivot combination*/
		for (int i = 0; i < rows ; i++)
			pivots.push_back(i);
		int last = pivots.size()-1;
		bool end = true;
		while (end){
			comb_piv.push_back(pivots);
			if (pivots[last] < cols-1)
				pivots[last]++;
			else{
				int k=0;
				for (int i = last-1 ; i>=0 ; i--){
					if (pivots[i] != cols-2-k){
						end = true;
						pivots[i]++;
						int aux = pivots[i];
						k=1;
						for (int j = i+1 ; j < last+1 ; j++){
							pivots[j] = aux + k;
							k++;
						}
						break;
					}
					else
						end =false;
					k++;
				}

			}
		}
	}

	Matrix gauss_jordan_test (IntervalMatrix& A, vector<pair <double, int > > order_cols, vector<pair <int, int> > & pair_contr,double prec){

		set <int> rows_checked;
		Matrix B(1,1);
		B.resize(A.nb_rows(),A.nb_cols());
		Matrix perm(1,1);
		perm.resize(B.nb_rows(),B.nb_rows());
		B = A.mid();

		for (int i = 0; i<A.nb_rows() ; i++)
			for (int j = 0; j<A.nb_rows() ; j++){
				if (i == j) perm[i][j] = 1;
				else perm[i][j] = 0;
			}
		for (int j = order_cols.size()-A.nb_rows(); j < order_cols.size() ; j++){
			int col_pos = order_cols[j].second;
			int temp_piv = -1;
				for (int i = 0; i < B.nb_rows() ; i++)
					if (( (B[i][col_pos] < -prec) || (B[i][col_pos] > prec)) && (rows_checked.count(i) != 1)){
						rows_checked.insert(i);
						temp_piv = i;
						pair_contr.push_back(make_pair(i,col_pos));
						break;
					}
					if (temp_piv !=-1){
						double coef = B[temp_piv][col_pos];
						Matrix aux_perm(1,1);
						aux_perm.resize(A.nb_rows(),A.nb_rows());
						for (int k = 0; k<A.nb_rows() ; k++)
							for (int l = 0; l<A.nb_rows() ; l++){
								if (k == l) aux_perm[k][l] = 1;
								else aux_perm[k][l] = 0;
							}

						for (int k = 0 ; k < B.nb_rows() ; k++){
							if ((k != temp_piv) &&( (B[k][col_pos] < -prec) || (B[k][col_pos] > prec))) {
								double factor = B[k][col_pos];
								aux_perm[k][temp_piv] = -factor/coef;
								for (int l = 0 ; l < B.nb_cols() ; l++)
									B[k][l]	= B[k][l]-(B[temp_piv][l]*factor/coef);
							}
						}
						/*new:dejar con 1 test*/
						for (int i = 0 ; i < B.nb_cols() ; i++)
							B[temp_piv][i] = B[temp_piv][i]/coef;
						aux_perm[temp_piv][temp_piv] = 1/coef;
						/*end: dejar con 1 test*/
						perm = aux_perm*perm;
					}
				}
		return perm;
	}

	bool max_piv(Matrix B, set<int> & ban_rows, set<int> & ban_cols, pair<int,int> & max_values, double prec){
		bool exist = false;
		for (int i = 0 ; i < B.nb_rows() ; i++){
			double max = 0;
			if (ban_rows.count(i) != 1){
				for (int j = 0 ; j < B.nb_cols() ; j++){
					if (ban_cols.count(j) != 1){
						if ((B[i][j] > prec) || (B[i][j] < -prec)){
							if (std::abs(B[i][j]) > max){
								exist = true;
								max = B[i][j];
								max_values.first = i ;
								max_values.second = j;
							}

						}
					}
				}
			}
		}
		ban_rows.insert(max_values.first);
		ban_cols.insert(max_values.second);
		return exist;
	}

	Matrix best_gauss_jordan (IntervalMatrix& A,  double prec){
		Matrix B(1,1);
		B.resize(A.nb_rows(),A.nb_cols());
		Matrix perm(1,1);
		set<int> ban_rows;
		set<int> ban_cols;
		perm.resize(B.nb_rows(),B.nb_rows());
		pair<int,int> max_values;
		/*Initialize the permutation matrix*/
		for (int i = 0; i<A.nb_rows() ; i++)
			for (int j = 0; j<A.nb_rows() ; j++){
				if (i == j) perm[i][j] = 1;
				else perm[i][j] = 0;
			}
		/*Initialize B*/
		B = A.mid();
		bool exist = true;
		while (exist){
			exist = max_piv(B,ban_rows,ban_cols,max_values, prec);
			if (exist){
				int temp_piv = max_values.first;
				int actual_col = max_values.second;
				double coef = B[temp_piv][actual_col];
				Matrix aux_perm(1,1);
				aux_perm.resize(A.nb_rows(),A.nb_rows());
				for (int m = 0; m<A.nb_rows() ; m++)
				for (int l = 0; l<A.nb_rows() ; l++){
					if (m == l) aux_perm[m][l] = 1;
					else aux_perm[m][l] = 0;
				}
				for (int m = 0 ; m < B.nb_rows() ; m++){
					if (m != temp_piv){
						double factor = B[m][actual_col];
						aux_perm[m][temp_piv] = -factor/coef;
						for (int l = 0 ; l < B.nb_cols() ; l++)
							B[m][l]	= B[m][l]-(B[temp_piv][l]*factor/coef);
					}
				}
				/*new:dejar con 1 test*/
					for (int i = 0 ; i < B.nb_cols() ; i++)
						B[temp_piv][i] = B[temp_piv][i]/coef;
					aux_perm[temp_piv][temp_piv] = 1/coef;
				/*end: dejar con 1 test*/
				perm = aux_perm*perm;
			}
		}
		A=perm*A;
		IntervalMatrix C = A;
		Matrix P_aux = perm;
		int rows_A_aux = A.nb_rows();
		A.resize(ban_rows.size(),A.nb_cols());
		perm.resize(ban_rows.size(),perm.nb_cols());

		int flag = 0;
		for (int i = 0; i < rows_A_aux ; i++){
			if (ban_rows.count(i) == 1){
				A[i-flag] = C[i];
				perm[i-flag] = P_aux[i];
			}
			else flag++;
		}
		return perm;
	}

	void gauss_jordan_all (IntervalMatrix& A, vector<Matrix>& permutations,vector < vector < pair <int, int> > > & pair_contr_all , double prec){
			int temp_piv;
			set <int> rows_checked;
			std::vector< std::vector <int> > comb_piv;
			/*get all possible combinations of pivots*/
			combinatorial(A,A.nb_cols(),A.nb_rows(),comb_piv);
			vector <pair<int,int> > aux_list;

			Matrix B(1,1);
			B.resize(A.nb_rows(),A.nb_cols());

			Matrix perm(1,1);
			perm.resize(B.nb_rows(),B.nb_rows());
			while(comb_piv.size() > 0){
				aux_list.clear();
				/*Initialize the permutation matrix*/
				for (int i = 0; i<A.nb_rows() ; i++)
					for (int j = 0; j<A.nb_rows() ; j++){
						if (i == j) perm[i][j] = 1;
						else perm[i][j] = 0;
				}
				/*Initialize B*/
				B = A.mid();
				rows_checked.clear();
				for (int k = 0 ; k < comb_piv[comb_piv.size()-1].size() ; k++){
					int actual_col = comb_piv[comb_piv.size()-1][k];
					temp_piv = -1;
					for (int i = 0; i < B.nb_rows() ; i++)
						if (( (B[i][actual_col] < -prec) || (B[i][actual_col] > prec)) && (rows_checked.count(i) != 1)){
							rows_checked.insert(i);
							temp_piv = i;
							aux_list.push_back(make_pair(temp_piv,actual_col));
							break;
						}
					if (temp_piv==-1)
						break;
					else{
						double coef = B[temp_piv][actual_col];
						Matrix aux_perm(1,1);
						aux_perm.resize(A.nb_rows(),A.nb_rows());
						for (int m = 0; m<A.nb_rows() ; m++)
							for (int l = 0; l<A.nb_rows() ; l++){
								if (m == l) aux_perm[m][l] = 1;
								else aux_perm[m][l] = 0;
							}
						for (int m = 0 ; m < B.nb_rows() ; m++){
							if (m != temp_piv){
								double factor = B[m][actual_col];
								aux_perm[m][temp_piv] = -factor/coef;
								for (int l = 0 ; l < B.nb_cols() ; l++)
									B[m][l]	= B[m][l]-(B[temp_piv][l]*factor/coef);
							}
						}
						/*new:dejar con 1 test*/
						for (int i = 0 ; i < B.nb_cols() ; i++)
							B[temp_piv][i] = B[temp_piv][i]/coef;
						aux_perm[temp_piv][temp_piv] = 1/coef;
						/*end: dejar con 1 test*/
						perm = aux_perm*perm;
					}
				}

				/*If gauss is perform complete, add the matrix perm to the list permutation*/
				if (temp_piv!=-1){
					pair_contr_all.push_back(aux_list);
					permutations.push_back(perm);
				}
				comb_piv.pop_back();
			}
		}




		Matrix gauss_jordan (IntervalMatrix& A, double prec){

				set <int> rows_checked;
				Matrix B(1,1);
				B.resize(A.nb_rows(),A.nb_cols());
				B = A.mid();
				Matrix perm(1,1);
				perm.resize(B.nb_rows(),B.nb_rows());

				for (int i = 0; i<A.nb_rows() ; i++)
					for (int j = 0; j<A.nb_rows() ; j++){
						if (i == j) perm[i][j] = 1;
						else perm[i][j] = 0;
					}

				/*gauss*/
				int temp_piv;
					for (int j = 0; j < B.nb_cols() ; j++){
					temp_piv = -1;
					for (int i = 0; i < B.nb_rows() ; i++)
						if (( (B[i][j] < -prec) || (B[i][j] > prec)) && (rows_checked.count(i) != 1)){
							rows_checked.insert(i);
							temp_piv = i;
							break;
						}
					if (temp_piv !=-1){
						double coef = B[temp_piv][j];
						Matrix aux_perm(1,1);
						aux_perm.resize(A.nb_rows(),A.nb_rows());
						for (int k = 0; k<A.nb_rows() ; k++)
							for (int l = 0; l<A.nb_rows() ; l++){
								if (k == l) aux_perm[k][l] = 1;
								else aux_perm[k][l] = 0;
							}

						for (int k = 0 ; k < B.nb_rows() ; k++){
							if ((k != temp_piv) &&( (B[k][j] < -prec) || (B[k][j] > prec))) {
	//							if ((k != temp_piv)){
								double factor = B[k][j];
								aux_perm[k][temp_piv] = -factor/coef;
								for (int l = 0 ; l < B.nb_cols() ; l++)
									B[k][l]	= B[k][l]-(B[temp_piv][l]*factor/coef);
							}
						}
						perm = aux_perm*perm;
					}
				}

				//TODO: dejar afuera
				A=perm*A;
				//TODO: modularizar
			    IntervalMatrix C = A;
			    Matrix P_aux = perm;
			    int rows_A_aux = A.nb_rows();
			    A.resize(rows_checked.size(),A.nb_cols());
			    perm.resize(rows_checked.size(),perm.nb_cols());

			    int flag = 0;
			    for (int i = 0; i < rows_A_aux ; i++){
			        if (rows_checked.count(i) == 1){
			            A[i-flag] = C[i];
			            perm[i-flag] = P_aux[i];
			        }
			        else flag++;
			    }
			    return perm;
			}

		void new_pseudoinverse(Array <const ExprNode>& xn, Array <const ExprNode>& bn, Matrix PA_aux, IntervalMatrix & A, double prec){
			int ind_col;
			set <int> li_cols;
			/*look for li columns*/
			for (int i = 0 ; i < PA_aux.nb_cols() ; i++){
				int count = 0;
				for (int j = 0; j < PA_aux.nb_rows() ; j++){
					if ((PA_aux[j][i] < -prec) || (PA_aux[j][i] > prec)){
						ind_col = j;
						count++;
					}
					if (count > 1){
						ind_col = -1;
						break;
					}
				}

				if ((ind_col != -1) && (li_cols.count(i) != 1)){
					li_cols.insert(i);
				}
			}
			int start= A.nb_rows();
			A.resize(A.nb_cols(), A.nb_cols());
			for (int k = 0 ; k < A.nb_rows() ; k++){
				if (li_cols.count(k) != 1){
					for (int j = 0; j < A.nb_cols() ; j++){
						if (j == k) A[start][j] = 1;
						else A[start][j] = 0;

					}
					start++;
					bn.add(xn[k]);
				}
			}
		}

		bool pseudoinverse(Matrix A, IntervalMatrix& P){
			Matrix tmp=A.transpose()*A;
			Matrix tmp2(A.nb_cols(),A.nb_cols());
			try{
				real_inverse(tmp, tmp2);
			}catch(ibex::Exception &e){
				return false;
			}
			P = tmp2*A.transpose();
			return true;
		}

	pair<int,int> find_next_pivot(Matrix & PA, IntervalVector x,set<int> & ban_rows, set<int> & ban_cols){

		vector<pair <double, pair<int,int> > > order_cols; /*impact value,variable,equation*/
		Matrix impact_values(PA.nb_rows(),PA.nb_cols());
		for (int j = 0 ; j < PA.nb_cols() ; j++){
			bool exist_value = false;
			if (ban_cols.count(j) != 1){
				/*aca*/
				for (int i = 0 ; i < PA.nb_rows() ; i++){
					if ((ban_rows.count(i) != 1) && (std::abs(PA[i][j])>1e-8)){
						impact_values[i][j] = std::abs(PA[i][j]);
						exist_value = true;
					}
				}
				if (exist_value)
								order_cols.push_back(make_pair(-1,make_pair(j,-1)));
			}
		}
		if (order_cols.empty()) return make_pair(-1,1);
		for (int j  = 0; j < order_cols.size() ; j++){
				for (int i = 0 ; i < PA.nb_rows() ; i++){
					if ((ban_rows.count(i) != 1)  ){
						int col = order_cols[j].second.first;
						if (impact_values[i][col] > order_cols[j].first){
							order_cols[j].first = impact_values[i][col];
							order_cols[j].second.second = i;
						}
					}
				}
		}

		std::sort(order_cols.begin(), order_cols.end(), compare);


		if (std::abs(order_cols[0].first) < 1e-8) return make_pair(-1,-1);
		else{
			ban_rows.insert(order_cols[0].second.second);
			ban_cols.insert(order_cols[0].second.first);
			return order_cols[0].second;
		}
	}

	void best_gauss_jordan_coll (IntervalMatrix A, IntervalVector x, vector<IntervalMatrix> & perm_list,
						vector <vector <pair <int,int> > > & proj_vars, double prec=1e-8){

			vector <pair<int,int> > aux_list;
			Matrix B(1,1);
			B.resize(A.nb_rows(),A.nb_cols());
			Matrix perm(1,1);
			set<int> ban_rows;
			set<int> ban_cols;
			perm.resize(B.nb_rows(),B.nb_rows());
			pair<int,int> max_values;
			bool available_cols = true;
			while (available_cols){

				aux_list.clear();
				ban_rows.clear();
				/*Initialize B*/
				B = A.mid();
				/*Initialize the permutation matrix*/
				for (int i = 0; i<A.nb_rows() ; i++)
					for (int j = 0; j<A.nb_rows() ; j++){
						if (i == j) perm[i][j] = 1;
						else perm[i][j] = 0;
					}
				while (ban_rows.size() != A.nb_rows()){
					if (ban_cols.size() == A.nb_cols()){
						ban_cols.clear();
						available_cols = false;
					}
					pair<int,int> var_eq = find_next_pivot(B, x, ban_rows, ban_cols);
					if (var_eq.first !=-1){
						aux_list.push_back(make_pair(var_eq.first,var_eq.second));
						double coef = B[var_eq.second][var_eq.first];
						Matrix aux_perm(1,1);
						aux_perm.resize(A.nb_rows(),A.nb_rows());

						for (int k = 0; k<A.nb_rows() ; k++)
							for (int l = 0; l<A.nb_rows() ; l++){
								if (k == l) aux_perm[k][l] = 1;
								else aux_perm[k][l] = 0;
							}
						for (int k = 0 ; k < B.nb_rows() ; k++){
							if ((k != var_eq.second) && (  std::abs(B[k][var_eq.first]) > prec)) {
								double factor = B[k][var_eq.first];
								aux_perm[k][var_eq.second] = -factor/coef;
								for (int l = 0 ; l < B.nb_cols() ; l++){
									if (l == var_eq.first) B[k][l] = 0;
									else B[k][l] = B[k][l]-(B[var_eq.second][l]*factor/coef);

								}
							}
						}
						/*make the pivot position 1*/
						for (int i = 0 ; i < B.nb_cols() ; i++){
							if (i == var_eq.first) B[var_eq.second][i] = 1;
							else B[var_eq.second][i] = B[var_eq.second][i]/coef;
						}
						aux_perm[var_eq.second][var_eq.second] = 1/coef;
						perm = aux_perm*perm;

					}
					else{
						available_cols = false;
						break;
					}
				}
				if (!(aux_list.empty())){
					proj_vars.push_back(aux_list);
					perm_list.push_back(perm);
				}
				if(A.nb_cols()==ban_cols.size()) available_cols = false;
			}
		}
	double impact_heuristic_new(int row, int column,IntervalMatrix PA, IntervalVector x,int heuristic){

			if (heuristic == 1){
				return std::abs(PA[row][column].mag());
			}

			else if (heuristic == 2) {
				return std::abs(PA[row][column].mag()*x[column].diam());
			}

			else if (heuristic == 3) {
				return std::abs(PA[row][column].mag());
			}
			else if (heuristic == 4) {
				return std::abs(PA[row][column].mag()*x[column].diam());
			}


			else return 0;
		}

	pair<int,int> find_pivot_new(IntervalMatrix & PA, IntervalVector x,set<int> & ban_rows, set<int> & ban_cols,int heuristic){
		pair<int,int> position;
		position.first = -1;
		position.second = -1;
		Matrix impact_values(PA.nb_rows(),PA.nb_cols());
		for (int i = 0 ; i < impact_values.nb_rows(); i++)
			for (int j = 0 ; j < impact_values.nb_cols() ; j++)
				impact_values[i][j] = 0;
		for (int j = 0 ; j < PA.nb_cols() ; j++){
			if (ban_cols.count(j) != 1){
				for (int i = 0 ; i < PA.nb_rows() ; i++){
					if ((ban_rows.count(i) != 1) &&  !(PA[i][j].contains(0))){
						impact_values[i][j] = impact_heuristic_new(i,j,PA,x,heuristic);
					}
				}
			}
		}
		/*search for the maximum value*/
		double max_value = 1e-8;
		if ((heuristic == 1) || (heuristic == 2))
			for (int i = 0 ; i < impact_values.nb_rows(); i++)
				for (int j = 0 ; j < impact_values.nb_cols() ; j++){
					if (impact_values[i][j] > max_value){
						max_value = impact_values[i][j];
						position.first = j;
						position.second = i;

					}
				}
		if ((heuristic == 3) || (heuristic == 4)){
			int best_col = -1;
			for (int i = 0 ; i < impact_values.nb_rows();i++){
				double sum = 0;
				for (int j =0 ; j < impact_values.nb_cols();j++)
					sum = sum + impact_values[i][j];
				if (sum != 0){
					for (int j =0 ; j < impact_values.nb_cols();j++)
						impact_values[i][j] = (double)(impact_values[i][j]/sum);
				}
			}
			Matrix impact_relative(1,PA.nb_cols());
			for (int j = 0 ; j < PA.nb_cols() ; j++)
				for (int i = 0 ; i < PA.nb_rows() ; i++)
					impact_relative[0][j] = impact_relative[0][j]+impact_values[i][j];
			for (int j = 0 ; j < PA.nb_cols() ; j++){
				if (impact_relative[0][j] > max_value){
					max_value = impact_relative[0][j];
					best_col = j;
				}
			}

			if (best_col == -1) return position;
			else{
				max_value = 1e-8 ;
				for (int i = 0 ; i < impact_values.nb_rows() ; i++){
					if (impact_values[i][best_col] > max_value){
						max_value = impact_values[i][best_col];
						position.first = best_col;
						position.second = i;

					}
				}
			}
		}




		if(position.first != -1){
			ban_rows.insert(position.second);
			ban_cols.insert(position.first);
		}
		return position;
	}


	void best_gauss_jordan_new (IntervalMatrix A, IntervalVector x, vector<IntervalMatrix> & perm_list,
						vector <vector <pair <int,int> > > & proj_vars, double prec,int heuristic){

			vector <pair<int,int> > aux_list;
			IntervalMatrix B(1,1);
			B.resize(A.nb_rows(),A.nb_cols());
			IntervalMatrix perm(1,1);
			set<int> ban_rows;
			set<int> ban_cols;
			perm.resize(B.nb_rows(),B.nb_rows());
			pair<int,int> max_values;
			bool available_cols = true;
			B = A;
			while (available_cols){

				aux_list.clear();
				ban_rows.clear();
				/*Initialize B*/
				B = A;
				/*Initialize the permutation matrix*/
				for (int i = 0; i<A.nb_rows() ; i++)
					for (int j = 0; j<A.nb_rows() ; j++){
						if (i == j) perm[i][j] = 1;
						else perm[i][j] = 0;
					}
				while (ban_rows.size() != A.nb_rows()){
					if (ban_cols.size() == A.nb_cols()){
						ban_cols.clear();
						available_cols = false;
					}
	//				pair<int,int> var_eq = find_next_pivot(B, x, ban_rows, ban_cols,heuristic);
					pair<int,int> var_eq = find_pivot_new(B, x, ban_rows, ban_cols,heuristic);
					if (var_eq.first !=-1){
						/*dddd*/
						aux_list.push_back(make_pair(var_eq.first,var_eq.second));
						Interval coef = B[var_eq.second][var_eq.first];
						IntervalMatrix aux_perm(1,1);
						aux_perm.resize(A.nb_rows(),A.nb_rows());
						for (int k = 0; k<A.nb_rows() ; k++)
							for (int l = 0; l<A.nb_rows() ; l++){
								if (k == l) aux_perm[k][l] = 1;
								else aux_perm[k][l] = 0;
							}
						for (int k = 0 ; k < B.nb_rows() ; k++){
							if ((k != var_eq.second) && (  (B[k][var_eq.first].mag() > prec))) {
								Interval factor = B[k][var_eq.first];
								aux_perm[k][var_eq.second] = -factor/coef;
								for (int l = 0 ; l < B.nb_cols() ; l++){
									if (l == var_eq.first) B[k][l] = 0;
									else B[k][l] = B[k][l]-(B[var_eq.second][l]*factor/coef);

								}
							}
						}
						/*make the pivot position 1*/
						for (int i = 0 ; i < B.nb_cols() ; i++){
							if (i == var_eq.first) B[var_eq.second][i] = 1;
							else B[var_eq.second][i] = B[var_eq.second][i]/coef;
						}
						aux_perm[var_eq.second][var_eq.second] = 1/coef;
						perm = aux_perm*perm;

					}
					else{
						available_cols = false;
						break;
					}
				}
				if (aux_list.size()>0){
					proj_vars.push_back(aux_list);
					perm_list.push_back(perm);
					 available_cols = false;
//					cout << perm.nb_cols() << " " << perm.nb_rows() << endl;
				}
				if(A.nb_cols()==ban_cols.size()) available_cols = false;
			}
		}


	Matrix best_P (IntervalMatrix& A, IntervalVector x, double prec){
		Matrix perm(A.nb_cols(),A.nb_rows());
		double max =0;

		bool found;
		IntervalVector box2(2*A.nb_cols()+A.nb_rows());
		vector <pair<double,Vector > > P_rows;
		//Vector row(box2.size());
		double suma[A.nb_cols()];
		for (int i = 0 ; i < A.nb_cols() ; i++){
			max=0;
			for (int j = 0 ; j < A.nb_rows(); j++)
				if (A[j][i].mag()>max){
				suma[i] = max;
				max = A[j][i].mag();
				}
		}
		int nb_eq=0;
		for (int i = 0 ; i < A.nb_cols() ; i++){
			LPSolver lp_solver(box2.size()); // p, s y w
			lp_solver.clean_ctrs();
			//lp_solver.set_bounds(IntervalVector(row.size()));

			//initializing domains and objective function
			for (int j = 0 ; j < box2.size() ; j++){
				if (j < A.nb_rows()){
					//initialize the variables p_i
					box2[j] = Interval(-1e7,1e7);
					lp_solver.set_obj_var(j,0);
				}
				else if (j >= A.nb_rows() && j < A.nb_rows() + A.nb_cols()){
					//initialize variables s_i
					box2[j]=Interval(-1e7, 1e7);
					lp_solver.set_obj_var(j,0);
				}
				else{
					//initialize auxiliary variables w_i
					box2[j]=Interval(-1e7, 1e7);
//					lp_solver.set_obj_var(j,1);
					lp_solver.set_obj_var(j,(x[j- A.nb_rows() -A.nb_cols()].diam()));
//					lp_solver.set_obj_var(j,(suma[j- A.nb_rows() -A.nb_cols()])*(x[j- A.nb_rows() -A.nb_cols()].diam()));
				}
			}
			lp_solver.set_bounds(box2);

			Vector row(box2.size());
			//s_i = 1
			row[A.nb_rows()+i] = 1;
			lp_solver.add_constraint(row,GEQ,1);
			lp_solver.add_constraint(row,LEQ,1);

			// w_i - s_i > ; w_i + s_i > 0
			for (int j = 0 ; j < A.nb_cols() ; j++){
				Vector row(box2.size());
				row[A.nb_rows() + A.nb_cols()+j] = 1;
				row[A.nb_rows()+j]=-1;
				lp_solver.add_constraint(row,GEQ,0);

				Vector row2(box2.size());
				row2[A.nb_rows() + A.nb_cols()+j] = 1;
				row2[A.nb_rows()+j]=1;
				lp_solver.add_constraint(row2,GEQ,0);
			}
			//s_j - sum_{j=1} A[i][j] p_i = 0
			for (int j = 0 ; j < A.nb_cols() ; j++){
				Vector row(box2.size());
				for (int ii = 0 ; ii < A.nb_rows() ; ii++)
					row[ii] = -A[ii][j].mid();
				row[A.nb_rows()+j]=1;
				lp_solver.add_constraint(row,GEQ,0);
				lp_solver.add_constraint(row,LEQ,0);
			}

			LPSolver::Status_Sol stat = lp_solver.solve();
			Vector v(box2.size());
			lp_solver.get_primal_sol(v);
			v.resize(A.nb_rows());

			found = false;
			for (int j = 0 ; j < nb_eq ; j++)
				if (v==perm[j]) {found = true; break;}


			if (found == false){
				perm.put(nb_eq,0,v,true);
				nb_eq++;
			}

		}
		perm.resize(nb_eq,A.nb_rows());
		A = perm*A;
		return perm;
	}


} /* namespace ibex */
