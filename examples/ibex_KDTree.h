#ifndef __KDTree__
#define __KDTree__

#include "ibex_IntervalVector.h"
#include "ibex_Vector.h"
#include <iostream>
#include <queue>
#include <stack>
#include <list>
#include <fstream>
#include<array>
#include <set>
#include "ibex_Function.h"

#include <math.h>


using namespace std;
using namespace ibex;
namespace kdtree {



/////////////////////////////////////////////////////////////////////////////////


    class Node  : public Vector {
        public:

            //Rango de cada dimension, con k el numero de dimensiones
            IntervalVector range;

            //Rango de evaluacion de la funcion en el nodo (se actualiza al buscar un punto)
            //TODO: evaluar si conviene que sea un mapa f->f.eval(range)
            Interval eval_range;

            int height;
            int dimension;

            //Nodos hijos y el que lo contiene, null si es raíz
            Node *left, *right, *padre;

            //Constructor
            Node(Vector&);
            
            //Destructor
            ~Node();

            Interval eval_r(const Function& function){
                eval_range = function.eval(range);
                return eval_range;
            }

            //Actualiza los rangos a partir del ultimo nodo insertado
            static void update_ranges(Node*);

            //Busca un punto en el árbol.
            Node* search(Vector&);

            //Agrega un hijo al nodo
            void add_child(Node* n);

            //Retorna true si tiene al menos un hijo, false en caso contrario.
            bool children();
            
            //Retorna la dimensión en la que el punto no está contenido en los rangos del nodo
            //  Si no está contenido en ninguna de ellas, retorna -1.
            int max_disjoint_dimension(IntervalVector&);
            

            void clear(){
                if(left) {
                    left->clear();
                    delete left;
                }

                if(right){
                    right->clear();
                    delete right;
                }
            }


    };


   
/////////////////////////////////////////////////////////////////////////////////////////////////////


    class IKDTree{           
        public:        

            Node* raiz;
            IKDTree();
            ~IKDTree(){if(raiz) raiz->clear();}
            
            //basic operators
            void insert(Vector&);
            bool remove(Vector& punto);   
            Vector* search(Vector& punto);   

            //generic operator
            set<pair <double, const Vector*>>  get_kmin(const Function* function, const list<Function>& constraints, int kmin=1, const IntervalVector& bounds=IntervalVector(0));  
            
            static bool satisfiesConstraints(const IntervalVector& box, const list<Function>& constraints);
            void updateMinSet(set<pair<double, const Vector*>>& min_set, int kmin, double point_eval, const Vector* punto);

            int height(){
                if(raiz == nullptr) return 0;
                return raiz->height;
            }

            int size(){
                return npoints;
            }

        private:
            int npoints;

    };
   
}

#endif //KDTreeRange2
