/*
 * ibex_CtcSystems.cpp
 *
 *  Created on: 13-10-2017
 *      Author: victor
 */

#include "ibex_CtcSystems.h"

using namespace std;

namespace ibex {

/*
 *
 *
 */
	void LinearSystem::contract(IntervalVector& box){

		if(ctc_type==SIMPLEX){
			ctc.contract(box);
			if(box.is_empty()) return;
		}

		if(ctc_type==PSEUDOINVERSE || ctc_type==GAUSS_JORDAN || ctc_type == GAUSS_PSEUDOINV ){
			IntervalMatrix tmp(PA);
			if(!bwd_mul(Pb, tmp, box, 1e8)){
				box.set_empty();
				return;
			}
		}
	}


	virtual int LinearSystem::linearization(const IntervalVector& x, LinearSolver& lp_solver){
		int num=0;
		for (int i=0; i<A.nb_rows(); i++) {


			Interval error = (IntervalVector(A[i].diam())* IntervalVector(x.mag())).ub();
			Vector row=A[i].mid();
			try {
			    Interval eval=(row*x); // left side image
				if(eval.lb() < (b[i] - error).lb()){
				   lp_solver.add_constraint(row,GEQ, (b[i] - error).lb());
				   num++;
				}
				if(eval.ub() > (b[i] + error).ub()){
					lp_solver.add_constraint(row,LEQ, (b[i] + error).ub());
					num++;

				}
			} catch (LPException&) { }
		}
		return num;
	}

	EmbeddedLinearSystem::EmbeddedLinearSystem(const IntervalMatrix& A, const IntervalMatrix& P, const IntervalMatrix& PA,
				Array<const ExprNode>& xn, Array<const ExprNode>& bn, map<const ExprNode*, int>& node2i,
				ExprDomain& d, bool is_mult=false, bool extended=false) :
					LinearSystem(A, P, PA), xn(xn), bn(bn), node2i(node2i), d(d),
						is_mult(is_mult), x(xn.size()), extended(extended)  {

		b.resize(A.nb_rows());
		if(extended)
			b=Vector::zeros(A.nb_rows());
	}

	EmbeddedLinearSystem::EmbeddedLinearSystem(const IntervalMatrix& AA, Array<const ExprNode>& xn,
				Array<const ExprNode>& bn, map<const ExprNode*, int>& node2i,
				ExprDomain& d, bool is_mult=false, int ctc_type = SIMPLEX, bool extended=false) :

		LinearSystem(AA.nb_cols()+((extended)? bn.size():0), ctc_type), xn(xn), bn(bn), node2i(node2i), d(d),
		is_mult(is_mult), x(xn.size()), extended(extended) {

		A=AA;
		b.resize(bn.size());
		if(extended){
			Matrix I = Matrix::diag(-Vector::ones(bn.size()));
			A.resize(AA.nb_rows(), AA.nb_cols()+bn.size());
			A.put(0,AA.nb_cols(),I);
			this->xn.add(bn);
			x.resize(this->xn.size());
			b=Vector::zeros(bn.size());
		}
	   if(ctc_type == PSEUDOINVERSE){
		   pseudoinverse(A.mid(), P);
		   PA=P*A;
	   }
	   else if(ctc_type == GAUSS_PSEUDOINV ){
		   //A is the matrix A after performing gauss elimination
		   //A <-- G*A : A is a diagonal matrix
		   //A = G*b
		   bool exist;
		   exist = pseudoinverse(A.mid(), P);
		   if (exist)
				PA=P*A;
		   else{
			   PA=A;
			   P = gauss_jordan(PA);
		   }
		   // after this: PA*x = P*b
	   }else if(ctc_type == GAUSS_JORDAN ){
		   PA=A;
		   P = gauss_jordan(PA);
	   }
	   cout.precision(3);
	}

	virtual void EmbeddedLinearSystem::contract(IntervalVector& box){
		IntervalVector logx(xn.size()), absx(xn.size());
		for(int i=0;i<xn.size(); i++){
		     x[i] = d[ node2i[ &xn[i]] ].i();
		     if(is_mult){
		    	 absx[i] = abs(x[i]);
		    	 logx[i] = log(absx[i]);
		     }
		}
		if(!extended)
			for (int i=0; i<b.size(); i++) {
				if(is_mult)
					b[i] = log(abs(d[ node2i[ &bn[i]] ].i()));
				else
					b[i] = d[ node2i[ &bn[i]] ].i();
			}
		if(ctc_type == PSEUDOINVERSE || ctc_type == GAUSS_JORDAN || ctc_type == GAUSS_PSEUDOINV ) Pb = P*b;
		LinearSystem::contract((is_mult)? logx:x);
		if(((is_mult)? logx:x).is_empty()){ box.set_empty(); return;}
		for(int i=0;i<xn.size(); i++){
			if(is_mult){
				bwd_log(logx[i], absx[i]);
				bwd_abs(absx[i], x[i]);
			}
			d[ node2i[ &xn[i]] ].i() &= x[i];
		}
		d.read_arg_domains(box);
	}

	virtual int EmbeddedLinearSystem::linearization(const IntervalVector& x, LinearSolver& lp_solver){
		return LinearSystem::linearization(x,lp_solver);
	}



	void EmbeddedLinearSystemBuilder::set_P(Matrix* P2){
		if(P) delete P;
		if(P2)
			P= new IntervalMatrix(*P2);
		else P=NULL;
	}

	void EmbeddedLinearSystemBuilder::set_PA(IntervalMatrix* PA2){
		if(PA) delete PA;
		if(PA2)
			PA= new IntervalMatrix(*PA2);
		else PA=NULL;
	}


	void create_subsystems(list<EmbeddedLinearSystemBuilder *> &ls_list, IntervalMatrix& A, Array<const ExprNode> &x, Array<const ExprNode> &b,
			vector<pair <set <int>,set <int> > >& subsets){

		set<int>::iterator it;
		set<int>::iterator it2;
		int j,k;
		for (int i = 0 ; i < subsets.size() ; i++){

			IntervalMatrix temp_A(1,1);
			Array<const ExprNode> temp_x; // size:nb_cols
			Array<const ExprNode> temp_b; // size:nb_rows
			j= 0; k=0;
			temp_A.resize(subsets[i].second.size(),subsets[i].first.size());
			if(subsets[i].second.size() > 1){
				for (it=subsets[i].second.begin(); it!=subsets[i].second.end(); it++){
					temp_b.add(b[*it]);
					for (it2=subsets[i].first.begin(); it2!=subsets[i].first.end(); it2++){
						temp_A[j][k] = A[*it][*it2];
						k++;
					}
					k = 0;
					j++;
				}
				for (it2=subsets[i].first.begin(); it2!=subsets[i].first.end(); it2++)
					temp_x.add(x[*it2]);
					EmbeddedLinearSystemBuilder*  system_A = new EmbeddedLinearSystemBuilder(temp_A,temp_x,temp_b);
					ls_list.push_back(system_A);

			}
		}
	}

	void find_subsystems(list<EmbeddedLinearSystemBuilder *> &ls_list, IntervalMatrix& A, Array<const ExprNode> &x,
			Array<const ExprNode> &b, int& nb_rows, int& nb_cols){
		vector<pair <set <int>,set <int> > > subsets;   //each subset contains the rows that are related

		vector< pair < set<int>, set<int> > > rows_sub;

		std::pair<std::set<int>::iterator,bool> ret;
		set<int>::iterator it;
		set<int>::iterator it2;

		for (int i=0; i< nb_rows ; i++){
			set<int> temp1;
			set<int> temp2;
			temp2.insert(i);
			for (int j=0; j< nb_cols ; j++)
				if (A[i][j] != 0)  temp1.insert(j);
			rows_sub.push_back( make_pair (temp1,temp2) );
			temp1.clear();
			temp2.clear();
		}

		while (!rows_sub.empty()){
			subsets.push_back(rows_sub[rows_sub.size()-1]);
			rows_sub.pop_back();
			int TAM = subsets.size()-1;
			int i=0;
			while(i < rows_sub.size()){
				for (it=rows_sub[i].first.begin(); it!=rows_sub[i].first.end(); ++it){
					if (subsets[TAM].first.count(*it)){
						/*insert the row in the subset*/
						subsets[TAM].second.insert(rows_sub[i].second.begin(),rows_sub[i].second.end());
						/*insert all columns in the subset*/
						subsets[TAM].first.insert(rows_sub[i].first.begin(),rows_sub[i].first.end());
						rows_sub.erase(rows_sub.begin()+i);
						i = -1;
						break;
					}
				}
				i++;
			}
		}

		create_subsystems(ls_list, A,x,b,subsets);
	}

	template<typename T>
	void add_row(IntervalMatrix& A, Array<const ExprNode> &b, const T* bi, map<const ExprNode*, int>& xmap,
			int& nb_rows, int& nb_cols){

		nb_rows++;
		A.resize(nb_rows,nb_cols);
		for(int i=0;i<nb_cols;i++) A[nb_rows-1][i]=0.0;
		b.add(*bi);

		map< const ExprNode*, Interval, lt > children;

		if(CAST(ExprAdd, bi))
			ExprNode2Dag().get_add_children(*bi,children);
		else
			ExprNode2Dag().get_mul_children(*bi,children);

		for(map< const ExprNode*, Interval, lt >::iterator it =children.begin(); it!=children.end();it++){
			if(xmap.find(it->first)==xmap.end()) {
				nb_cols++;
				xmap[it->first] = nb_cols-1;
				A.resize(nb_rows,nb_cols);
				for(int i=0;i<nb_rows-1;i++) A[i][nb_cols-1]=0.0;
			}
			A[nb_rows-1][ xmap[it->first] ]=it->second;
		}
	}

	template<typename T>
	Array<Ctc> getEmbeddedLinearSystems(CtcDag& dag_ctc, bool is_mult, int ctc_type, bool extended, int nb_nodes=-1){

		Function& f = dag_ctc.get_f();
		ExprDomain& d  = dag_ctc.getDomains();

		map<const ExprNode*, int> node2i;

		list<EmbeddedLinearSystemBuilder *> ls_list;
		Array<Ctc> linear_systems;

		IntervalMatrix A(1,1);
		Array<const ExprNode> x; // size:nb_cols
		Array<const ExprNode> b; // size:nb_rows
		int nb_rows=0;
		int nb_cols=0;

		map<const ExprNode*, int> xmap;
		set<const T*> bmap;

	    if(nb_nodes==-1) nb_nodes=f.nb_nodes();

		for(int i=0;i<f.nb_nodes(); i++){
			const ExprNode* n = &f.node(i);
			node2i[n]=i;


			if(CAST(T, n) && i<=nb_nodes){

				for(int i=0; i<n->fathers.size(); i++){
					if(CAST(T, &n->fathers[i])) continue;

					const T* bi;
					if(! CAST(ExprConstant, &CAST(T,n)->left) )
						bi=CAST(T,n);
					else if (! CAST(ExprConstant,& CAST(T,n)->right) &&
							CAST(T, & CAST(T,n)->right) )
						bi= CAST(T, & CAST(T,n)->right) ;
					else break;

	                if(bmap.find(bi)==bmap.end()   ){
	                   bmap.insert(bi);
					   add_row(A, b, bi, xmap, nb_rows , nb_cols);
	                }

					break;
				}
			}
		}

		x.resize(xmap.size());

		for(map<const ExprNode*, int>::iterator it =xmap.begin(); it!=xmap.end();it++){
			x.set_ref(it->second,*node2i.find(it->first)->first);
		}

		find_subsystems(ls_list, A, x, b, nb_rows , nb_cols);

		list<EmbeddedLinearSystemBuilder *>::iterator it;
		for(it=ls_list.begin(); it!=ls_list.end(); it++){
			(*it)->set_node2i(&node2i);
			(*it)->set_ctc_type(ctc_type);
			(*it)->set_domain(&d);
			(*it)->set_is_mult(is_mult);
			(*it)->set_extended(extended);


	        if(ctc_type==LinearSystem::MULT_GAUSS_JORDAN){

	    		if(extended){
	    			Matrix I = Matrix::diag(-Vector::ones((*it)->get_bn().size()));
	    			int nb_rows=(*it)->getA().nb_rows();
	    			int nb_cols=(*it)->getA().nb_cols();
	    			(*it)->getA().resize(nb_rows, nb_cols+(*it)->get_bn().size());
	    			(*it)->getA().put(0,nb_cols,I);
	    			(*it)->get_xn().add((*it)->get_bn());
	    		}

	    		list<Matrix> Ps;
	    		list<IntervalMatrix> PAs;

	    		gauss_jordan_collection((*it)->getA(),Ps, PAs);

	    		list<Matrix>::iterator itP=Ps.begin();
	    		list<IntervalMatrix>::iterator itPA=PAs.begin();

	    		for(;itP!=Ps.end();itP++, itPA++){
	    			(*it)->set_P(&(*itP));
	    			(*it)->set_PA(&(*itPA));
	    			linear_systems.add(*(*it)->create());
	    		}
				(*it)->set_P(NULL);
				(*it)->set_PA(NULL);

	        }else{

	        	try{
	        		linear_systems.add(*(*it)->create());
	        	}catch(SingularMatrixException &e) {cout << "SingularMatrixException" << endl;}
	        }

		}
		return linear_systems;
	}

}

