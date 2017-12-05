/*
 * ibex_Conditioners.cpp
 *
 *  Created on: 13-10-2017
 *      Author: victor
 */

#include "ibex_Conditioners.h"

using namespace std;

namespace ibex {


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

	void gauss_jordan_test (IntervalMatrix& A, double prec){

		Matrix perm(1,1);
		Matrix B(1,1);
		int row_pos,col_pos;
		B.resize(A.nb_rows(),A.nb_cols());
		B = A.mid();
		perm.resize(B.nb_cols(),B.nb_cols());
		for (int i = 0 ; i < perm.nb_rows() ; i++)
			for (int j = 0; j<perm.nb_rows() ; j++){
				if (i == j) perm[i][j] = 1;
				else perm[i][j] = 0;
			}
		for (int k = 0 ; k < B.nb_rows() ; k++){
			double max = -1;
			for (int i = k; i < B.nb_rows() ; i++)
				for (int j = k ; j < B.nb_cols() ; j++)
					if (std::abs(B[i][j]) > max){
						max = std::abs(B[i][j]);
						row_pos = i ;
						col_pos = j ;
					}

		/*pivot exchange*/


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
