//============================================================================
//                                  I B E X
// File        : Composition of contractors
// Author      : Ignacio Araya
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Mar 23, 2018
// Last Update :
//============================================================================

#include "ibex_CtcAdaptive.h"

using namespace std;

namespace ibex {



CtcAdaptive::CtcAdaptive(const Array<Ctc>& list, int m, int L, bool bf, bool fp ) :
		Ctc(list), list(list),iter(0), nb_ctr(m), L(L), bf(bf), fp(fp) {
	assert(check_nb_var_ctc_list(list));

	gcalls = new int[list.size()];
	calls = new int[list.size()];
	effective_calls = new int[list.size()];
	nb_input_ctr = new int[list.size()];
	nb_act_ctr = new int[list.size()];

	for(int i=0; i<list.size(); i++){
		gcalls[i]=0;
		calls[i]=0;
		effective_calls[i]=0;
		nb_input_ctr[i]=0;
		nb_act_ctr[i]=0;
	}

}

AdaptCell* CtcAdaptive::closest_informed_ancestor(Cell& cell, int ctc, int c){

	AdaptCell* anc=&cell.get<AdaptCell>();

	while(anc->p){
		anc=anc->p;
		if(anc->T.find(make_pair(ctc,c)) != anc->T.end())	 return anc;
	}


	return anc; // <-- root

}



void CtcAdaptive::contract(Cell& c) {
  //cout << "CONTRACT" << endl;
	if(!c.data.used(typeid(AdaptCell).name())) c.add<AdaptCell>(); //root node

	//ctc.contract(c.box);

	BitSet flags(BitSet::empty(Ctc::NB_OUTPUT_FLAGS));
	BitSet impact(BitSet::all(c.box.size())); // always set to "all" for the moment (to be improved later)

    int k=0;
    int nb_succ_ctc = list.size();
		bool success[nb_succ_ctc];
		for(int i=0;i<list.size();i++)
		  success[i]=true;

   do{

	 for(int j=0; j<list.size();j++){

		 if(k>0 && nb_succ_ctc<2) break;

		//array of informed nodes for each constraint i
		AdaptCell* a[nb_ctr];

		if(list[j].input_ctr){
			list[j].input_ctr->clear();

			for(int i=0; i<nb_ctr; i++){
				a[i] = closest_informed_ancestor(c, j, i);
				if(a[i]->T.find(make_pair(j,i)) == a[i]->T.end()) a[i]->T[make_pair(j,i)]=1;


				if( j==0 ||(k==0 && iter % a[i]->T[make_pair(j,i)] == 0) ||
					(	(a[i]->T[make_pair(j,i)]==1 && a[i]->F[make_pair(j,i)]==0)) || L==10000 )
					list[j].input_ctr->add(i);

			}

			gcalls[j]++;
			if(list[j].input_ctr->size()==0) {
				if(k==0){
					success[j]=false;
					nb_succ_ctc--;
				}
				continue;
			}
		}



		if(list[j].active_ctr) list[j].active_ctr->fill(0,nb_ctr-1);
		list[j].contract(c.box,impact,flags);


		//active constraints
		BitSet ca(nb_ctr);
		if(list[j].active_ctr) ca=*list[j].active_ctr;
		if(!list[j].active_ctr || (ca.size()==0 && c.box.is_empty()) ) ca.fill(0,nb_ctr-1);

		calls[j]++;
		if(list[j].active_ctr && ca.size()>0) {
			effective_calls[j]++;
		}else if(success[j]){
			success[j]=false;
			nb_succ_ctc--;
		}

		if(list[j].input_ctr) nb_input_ctr[j]+=list[j].input_ctr->size();
		if(list[j].active_ctr) nb_act_ctr[j]+=list[j].active_ctr->size();



    if (c.box.is_empty()) cout << j << ": empty_box" << endl;
		if(list[j].input_ctr)
			cout << j << ":" << list[j].input_ctr->size() << "-->"<< ca.size() << endl;
		else
			cout << j << ":" << ca.size() << endl;
    //cout << c.box << endl;

		if(list[j].input_ctr){
			for(int i=0; i<nb_ctr; i++){

				if((*list[j].input_ctr)[i] && ca[i]){

					a[i]->T[make_pair(j,i)] = 1;
					a[i]->F[make_pair(j,i)] = 0;

					AdaptCell* p = c.get<AdaptCell>().p;
					if(p && bf){
						p->T[make_pair(j,i)] = 1;
						p->F[make_pair(j,i)] = 0;
					}


				}else if((*list[j].input_ctr)[i]){
					a[i]->F[make_pair(j,i)]++;

					if(a[i]->F[make_pair(j,i)]==L){
						a[i]->T[make_pair(j,i)] *= 2;
						a[i]->F[make_pair(j,i)] = 0;
						//cout  << a[i] <<  j << "," << i << ":" << a[i]->T[make_pair(j,i)] << endl;
					}

				}

			}
		}

		if (c.box.is_empty()){
			iter++;
			return;
		}

	}
		k++;

	}while(fp && nb_succ_ctc>=2);

	iter++;

}

} // end namespace ibex
