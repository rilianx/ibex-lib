/*
 * ibex_TreeCellOpt.cpp
 *
 *  Created on: 26 feb. 2018
 *      Author: iaraya
 */


#include "ibex_BufferPSO.h"

using namespace std;

namespace ibex {
	System* Part::orig_sys=orig_sys;



  /*
	* Removes every node with no children in the last_node branch
	*/
	void BufferPSO::trim() {

    //particles are moved into the subboxes (diving)
    for(int i=0;i<popsize; i++){
        if(tree->contains(last_node->box, particles[i]->position)){
            double d1 = (last_node->get<CellPSO>().left)?
                    tree->manhattan_distance(particles[i]->position, last_node->get<CellPSO>().left->box) : POS_INFINITY;
            double d2 = (last_node->get<CellPSO>().right)?
                    tree->manhattan_distance(particles[i]->position, last_node->get<CellPSO>().right->box) : POS_INFINITY;
            if(d1 < d2 && d1 > 0.0) particles[i]->reinitialize(last_node->get<CellPSO>().left->box.random());
            else if(d2 != POS_INFINITY && d2 > 0.0) particles[i]->reinitialize(last_node->get<CellPSO>().right->box.random());
            else if(d1>0.0 && d2 > 0.0) particles[i] = NULL;
        }
    }

		tree->trim(last_node);
		last_node=NULL;
	}

	Cell* BufferPSO::top(double current_ub) const{

		for(int i=0;i<popsize; i++){
		  if(!particles[i]){
			//TODO: seleccionar nodos prometedores de zonas altas del árbol
			Cell* c = tree->minlb(tree->root);
			particles[i] = new Part(c->box.random());
		  }
		//  cout << i << ":" << -particles[i]->fitness(current_ub) << endl;
		}

		/*
			Seleccionar 2 partículas
			h = p2 + rand(0.0,2.0)*(p1-p2)
			Si h es box-feasible y h>p2, p2=h
			Retornar mejor de las 3 partículas
			*/

		int i1 = rand()%popsize;
		int i2 = rand()%popsize;
		while(i2==i1) i2=rand()%popsize;
	  //

		Part* p1=particles[i1];
		Part* p2=particles[i2];
		if(p1->fitness(current_ub) < p1->fitness(current_ub)){
		  p1=p2;
		  p2=particles[i1];
		}
		double r = (2.0*rand())/RAND_MAX;
		Vector h = p2->position + r*(p1->position-p2->position);

		if(tree->search(h)){
		  Part p(h);
		  if (p.fitness(current_ub) > p2->fitness(current_ub)){
			 *p2=p;
		//	 cout << i1 << "," << i2 << endl;
		//	 cout << "improved:" << -p2->fitness(current_ub) << endl;
		   }
		}

		Part* pbest=particles[0];
		for(int i=1;i<popsize; i++)
			if(particles[i]->fitness(current_ub) > pbest->fitness(current_ub))
				pbest=particles[i];

		Cell* aux=tree->search(pbest->position);
		last_node = aux;
		return aux;

	}
}
