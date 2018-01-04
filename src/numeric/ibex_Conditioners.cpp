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

	void gauss_jordan_all (IntervalMatrix& A, double prec){

		cout << A << endl;
		int temp_piv;
		set <int> rows_checked;
		vector<Matrix> permutations;
		std::vector< std::vector <int> > comb_piv;
		/*get all possible combinations of pivots*/
		combinatorial(A,A.nb_cols(),A.nb_rows(),comb_piv);
		cout << comb_piv.size() << endl;
		Matrix B(1,1);
		B.resize(A.nb_rows(),A.nb_cols());

		Matrix perm(1,1);
		perm.resize(B.nb_rows(),B.nb_rows());
		/*perform the gauss elimination for each comb_piv element*/
		while(comb_piv.size() != 0){
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
					perm = aux_perm*perm;
				}
			}

			/*If gauss is perform complete, add the matrix perm to the list permutation*/
			if (temp_piv!=-1) permutations.push_back(perm);
			comb_piv.pop_back();
		}
		cout << permutations.size() << endl;
		for (int i = 0 ; i < permutations.size() ; i++)
			cout << permutations[i]*A.mid() << endl;
		exit(0);
//		return permutations;
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
			gauss_jordan_all (A,prec);
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
						if (k != temp_piv){
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
