/*
 * ibex_CtcSystems.h
 *
 *  Created on: 13-10-2017
 *      Author: victor
 */

#ifndef SRC_CONTRACTOR_IBEX_CTCSYSTEMS_H_
#define SRC_CONTRACTOR_IBEX_CTCSYSTEMS_H_
#include <list>
#include <map>
#include <set>

#include "ibex_Ctc.h"
#include "ibex_Linearizer.h"
#include "ibex_CtcPolytopeHull.h"
#include "ibex_CtcDAG.h"
#include "ibex_Conditioners.h"
#include "ibex_ExprNode2Dag.h"


#define CAST(m,n) dynamic_cast<const m*>(n)

using namespace std;

namespace ibex {

/**
 * \brief A contractor using a linear system A.x=b
 *
 * Linear system of the form A.x=b.
 *
 * If P exists, the contractor use it for contracting x by using the system PA.x=Pb
 * and a simple projection based method
 * If P does not exist the contractor perform at most 2n simplex for contracting x
 */

class LinearSystem : public Linearizer, public Ctc {
public:


	enum {SIMPLEX, PSEUDOINVERSE, GAUSS_JORDAN, GAUSS_PSEUDOINV, MULT_GAUSS_JORDAN, NEW_PSEUDOINV};

	CtcPolytopeHull ctc;
	int ctc_type;


	/**
	 * \brief Create the linear inequalities Ax<=b.
	 */
	LinearSystem(int nb_var, int ctc_type = SIMPLEX) :
		Linearizer(nb_var),
		ctc(*this, LPSolver::default_max_iter,
		LPSolver::default_max_time_out, LPSolver::default_eps, Interval (1e-14, 1e10)),
		ctc_type(ctc_type), Ctc(nb_var), A(1,1), b(1), P(1,1), PA(1,1), Pb(1) {

	}

	/**
	 * \brief Create the linear inequalities Ax<=b.
	 */
	LinearSystem(const IntervalMatrix& A, const IntervalMatrix& P, const IntervalMatrix& PA,  int ctc_type = GAUSS_JORDAN) :
		Linearizer(A.nb_cols()), ctc(*this, LPSolver::default_max_iter,
		LPSolver::default_max_time_out, LPSolver::default_eps, Interval (1e-14, 1e10)),
		ctc_type(ctc_type), Ctc(A.nb_cols()), A(A), b(1), P(P), PA(PA), Pb(1) {

	}

	/**
	 * Performs the contraction of the current box using a linear system
	 */
	void contract(IntervalVector& box);

	/**
	 * \brief Add the inequalities in the solver (does not perform a linearization)
	 *
	 * This method is required by the PolytopeHull
	 */
	int linearization(const IntervalVector& x, LPSolver& lp_solver);




protected:

	IntervalMatrix A;
	IntervalVector b;

	IntervalMatrix P;
	IntervalMatrix PA;
	IntervalVector Pb;

};

/** \ingroup contractor
 *
 * \brief EmbeddedLinearSystem (contractor).
 *
 *
 */
class EmbeddedLinearSystem : public LinearSystem {
public:

	/**
	 * All but SIMPLEX
	 */
	EmbeddedLinearSystem(const IntervalMatrix& A, const IntervalMatrix& P, const IntervalMatrix& PA,
					Array<const ExprNode>& xn, Array<const ExprNode>& bn, map<const ExprNode*, int>& node2i,
					ExprDomain& d, bool is_mult, bool extended);

	EmbeddedLinearSystem(const IntervalMatrix& AA, Array<const ExprNode>& xn,
					Array<const ExprNode>& bn, map<const ExprNode*, int>& node2i,
					ExprDomain& d, bool is_mult, int ctc_type, bool extended);

	//TODO: implement this
	//static Array<Ctc> get_ctcs(CtcDag& dag_ctc, bool is_mult, int ctc_type, bool extended, int nb_nodes=-1);

	/**
	 * \brief Maps intervals from the nodes, then contracts them by using
	 * PA.x=b and then maps the reductions back into the DAG
	 * \param box
	 */
	void contract(IntervalVector& box);

	/**
	 * TODO: add comments
	 */
	virtual int linearize(const IntervalVector& x, LPSolver& lp_solver){
		return LinearSystem::linearization(x,lp_solver);
	}

	/**
	 * TODO: why public?
	 */
	bool is_mult;
	bool extended;


private:

	/**
	 * TODO: add comments
	 */
	Array<const ExprNode> xn; // size:nb_cols
	Array<const ExprNode> bn; // size:nb_rows
	IntervalVector x;

	map<const ExprNode*, int> node2i;
	ExprDomain& d;

};

/**
 * TODO: add comments
 */

class EmbeddedLinearSystemBuilder {

public:

	EmbeddedLinearSystemBuilder(const IntervalMatrix& A, Array<const ExprNode>& xn,
			Array<const ExprNode>& bn) : A(A), xn(xn), bn(bn), PA(NULL), P(NULL), node2i(NULL), d(NULL), is_mult(false), ctc_type(-1), extended(false) { }

	void set_P(Matrix* P2);
	void set_PA(IntervalMatrix* PA2);

	void set_node2i(map<const ExprNode*, int>* n) { node2i=n; }
	void set_domain(ExprDomain* dd) { d=dd; }
	void set_is_mult(bool m) { is_mult=m; }
	void set_ctc_type(int t) { ctc_type=t; }
	void set_extended(bool e) { extended=e; }


	EmbeddedLinearSystem* create();

	inline Array<const ExprNode>& get_bn() { return bn;}
	inline Array<const ExprNode>& get_xn() { return xn;}


    inline IntervalMatrix& getA() { return A; }


private:
	IntervalMatrix A;
	IntervalMatrix* P;
	IntervalMatrix* PA;

	Array<const ExprNode> xn; // size:nb_cols
	Array<const ExprNode> bn; // size:nb_rows

	map<const ExprNode*, int>* node2i;
	ExprDomain* d;

	bool is_mult;
	int ctc_type;
	bool extended;
};

	/**
	 * \brief This function creates the subsystems once is found by the function find_subsystems.
	 */
	void create_subsystems(list<EmbeddedLinearSystemBuilder *> &ls_list, IntervalMatrix& A, Array<const ExprNode> &x, Array<const ExprNode> &b,
		vector<pair <set <int>,set <int> > >& subsets);

	/**
	 * \brief This function finds all the independent subsystems from the matrix A, in order to separate them. Note that
	 * no information is lost by doing this process.
	 */
	void find_subsystems(list<EmbeddedLinearSystemBuilder *> &ls_list, IntervalMatrix& A, Array<const ExprNode> &x,
		Array<const ExprNode> &b, int& nb_rows, int& nb_cols);


	/**
	 * \brief This function adds a new row to the matrix A once a new subexpression is found.
	 */
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
//		cout<< ls_list.size() << endl;
		list<EmbeddedLinearSystemBuilder *>::iterator it;
		for(it=ls_list.begin(); it!=ls_list.end(); it++){
			(*it)->set_node2i(&node2i);
			(*it)->set_ctc_type(ctc_type);
			(*it)->set_domain(&d);
			(*it)->set_is_mult(is_mult);
			(*it)->set_extended(extended);



	        try{
	        	linear_systems.add(*(*it)->create());
	        }catch(SingularMatrixException &e) {cout << "SingularMatrixException" << endl;}


		}
		return linear_systems;
	}



}


#endif /* SRC_CONTRACTOR_IBEX_CTCSYSTEMS_H_ */
