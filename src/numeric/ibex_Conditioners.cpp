/*
 * ibex_Conditioners.cpp
 *
 *  Created on: 13-10-2017
 *      Author: victor
 */

#include "ibex_Conditioners.h"

using namespace std;

namespace ibex {

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

	void gauss_jordan_all (IntervalMatrix& A, vector<Matrix>& permutations, double prec){
		int temp_piv;
		set <int> rows_checked;
		std::vector< std::vector <int> > comb_piv;
		/*get all possible combinations of pivots*/
		combinatorial(A,A.nb_cols(),A.nb_rows(),comb_piv);

//		cout << comb_piv.size() << endl;
		Matrix B(1,1);
		B.resize(A.nb_rows(),A.nb_cols());

		Matrix perm(1,1);
		perm.resize(B.nb_rows(),B.nb_rows());
		/*perform the gauss elimination for each comb_piv element*/
//		comb_piv.clear();
//		std::vector<int> aux;
//		aux.push_back(3); aux.push_back(1);
//		comb_piv.push_back(aux);
		while(comb_piv.size() > 0){
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
						break;
					}
				if (temp_piv==-1) break;
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
//					for (int i = 0 ; i < B.nb_cols() ; i++)
//						B[temp_piv][i] = B[temp_piv][i]/coef;
//					aux_perm[temp_piv][temp_piv] = 1/coef;
					/*end: dejar con 1 test*/
					perm = aux_perm*perm;
				}
			}

			/*If gauss is perform complete, add the matrix perm to the list permutation*/
			if (temp_piv!=-1) permutations.push_back(perm);
			comb_piv.pop_back();
		}
	}

	bool max_piv(Matrix B, set<int> & ban_rows, set<int> & ban_cols, pair<int,int> & max_values,double prec){
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

	Matrix best_gauss_jordan (IntervalMatrix& A, double prec){
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
			exist = max_piv(B,ban_rows,ban_cols,max_values,prec);
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
		return perm;
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
} /* namespace ibex */
