/*
 * ibex_CtcCompoHiper.cpp
 *
 *  Created on: 06-11-2017
 *      Author: victor
 */

#include "ibex_CtcCompoHiper.h"

namespace ibex {


void CtcCompoHiper::init_impacts() {
	for(int i = 0 ; i < list.size() ; i++){
			ctc_calls.push_back(0);
			ctc_contractions.push_back(0);
		}
		all_ctc_contract = 0;
		llamadas = 0;
}





CtcCompoHiper::CtcCompoHiper(const Array<Ctc>& list) :
		Ctc(list), list(list){
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompoHiper::CtcCompoHiper(Ctc& c1, Ctc& c2) :
		Ctc(c1.nb_var),list(Array<Ctc>(c1,c2)){
	assert(check_nb_var_ctc_list(list));
	init_impacts();
}

CtcCompoHiper::CtcCompoHiper(Ctc& c1, Ctc& c2, Ctc& c3) :
		Ctc(c1.nb_var),list(Array<Ctc>(c1,c2,c3)) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}


CtcCompoHiper::CtcCompoHiper(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4):
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4)) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompoHiper::CtcCompoHiper(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5)) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompoHiper::CtcCompoHiper(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6)) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompoHiper::CtcCompoHiper(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7)){
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}



CtcCompoHiper::~CtcCompoHiper() {
//	delete[] impacts;
}


void CtcCompoHiper::contract(IntervalVector& box) {

	bool contraction[list.size()];
	bool flag = true;
	for (int i=0; i<list.size(); i++)
		contraction[i] = false;
	for (int i=0; i<list.size(); i++) {
		ctc_calls[i]++;
		IntervalVector box_aux = box;
		list[i].contract(box);
		/*primera forma*/
		if(box.is_strict_subset(box_aux)){
			ctc_contractions[i]++;
			contraction[i] = true;
		}
		/*segunda forma*/
//		for (int j = 0; j < box.size(); j++){
//			if (box_aux[j].diam()-box[j].diam() >= 0.01) {
//				ctc_contractions[i]++;
//				contraction[i] = true;
//				break;
//			}
//		}
		/*tercera forma*/
//		for (int j = 0; j < box.size(); j++){
//			if (box[j].diam()/box_aux[j].diam() <= 0.99) {
//				ctc_contractions[i]++;
//				contraction[i] = true;
//				break;
//			}
//		}

		if (box.is_empty()) {
			for (int i=0; i<list.size(); i++)
				if (contraction[i] == false){
					flag = false;
					break;
				}
			if (flag){
				all_ctc_contract++;
			}
			return;
		}
	}
	for (int i=0; i<list.size(); i++)
		if (contraction[i] == false){
			flag = false;
			break;
		}
	if (flag)
		all_ctc_contract++;
}

} // end namespace ibex
