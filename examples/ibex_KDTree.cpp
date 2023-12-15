
#include "ibex_KDTree.h"
#include <bits/stdc++.h>
#include <queue>
#include <stack>

//
#include "ibex_Function.h"
#include "set"
//

using namespace ibex ;
using namespace std;
namespace kdtree {


Node::Node(Vector& punto_) : Vector(punto_), range(punto_){
        /* Se inicializa todos las 'caracteristicas' desde cero*/
        height = 0;
        dimension = -1;
        left = nullptr;
        right = nullptr;
        padre = nullptr;
    }
    
    Node::~Node(){

    }


set<pair<double, const Vector*>> IKDTree::get_kmin(const Function* function, const list<Function>& constraints, int k, const IntervalVector& bounds) {
    set<pair<double, const Vector*>> min_set;  // Conjunto mínimo a devolver

    int n_evals = 0;
    int n_visited = 0;
    cout << "bounds: " << bounds << endl;

    if (raiz == nullptr) {
        cout << "No existe el árbol" << endl; 
        return min_set;
    }

    std::stack<Node*> stack;
    stack.push(raiz);

    while (!stack.empty()) {
        Node* nCurrent = stack.top();  // Nodo actual en la parte superior de la pila
        stack.pop();

        // Si el conjunto mínimo ya tiene al menos 'k' elementos 
        // verifica que algún punto del subárbol podría ser mejor que el peor punto del conjunto mínimo
        // y que el subárbol intersecta con el rango de búsqueda.
        // y que algún punto dentro del subárbol podría satisfacer las restricciones
        if ( (min_set.size() == k && nCurrent->eval_range.lb() > (--min_set.end())->first) ||
                !nCurrent->range.intersects(bounds) ||
                !satisfiesConstraints(nCurrent->range, constraints)) continue;
        
        Vector* punto = nCurrent; 

        // Evalúa la función en el punto del nodo.
        double point_eval = 0.0;
        if (function != NULL) point_eval=function->eval(*punto).ub(); 
        
        n_visited ++;
        n_evals ++;
        
        // Verifica si el punto del nodo satisface las restricciones.
        if (satisfiesConstraints(*punto, constraints) && bounds.contains(*punto)) 
            // Actualiza el conjunto mínimo si el punto es uno de los k mejores.
            updateMinSet(min_set, k, point_eval, punto);

        if (function == NULL && min_set.size() == k) 
            break; //ya se juntaron k puntos que satisfacen las restricciones

        // Agrega los hijos del nodo actual a la pila para explorar más nodos.
        if (nCurrent->left != nullptr && nCurrent->right != nullptr) {
            if (function) n_evals += 2;

            if (!function || nCurrent->left->eval_r(*function).lb() < nCurrent->right->eval_r(*function).lb()) {
                stack.push(nCurrent->left);  
                stack.push(nCurrent->right); // Agrega el hijo derecho al final si es más prometedor.
            } else {
                stack.push(nCurrent->right);  
                stack.push(nCurrent->left); // Agrega el hijo izquierdo al final si es más prometedor.
            }
        } else {
            if (nCurrent->left != nullptr) 
                stack.push(nCurrent->left);  // Agrega el hijo izquierdo si existe.
            else if (nCurrent->right != nullptr) 
                stack.push(nCurrent->right);  // Agrega el hijo derecho si existe.
            
        }
    }

    cout << "Número de evaluaciones: " << n_evals << endl;
    cout << "Número de nodos visitados: " << n_visited << endl;
    return min_set;  // Devuelve el conjunto mínimo resultante.
}

bool IKDTree::satisfiesConstraints(const IntervalVector& box, const list<Function>& constraints) {
    for (const auto& constraint : constraints)
        if (constraint.eval(box).lb() > 0) 
            return false;
    return true;
}

void IKDTree::updateMinSet(set<pair<double, const Vector*>>& min_set, int k, double point_eval, const Vector* punto) {

    if (min_set.size() == k) {
        auto lastSet = --min_set.end();
        if (point_eval < lastSet->first) {
            min_set.erase(lastSet);
            min_set.insert(make_pair(point_eval, punto));
        }
    }else {
        min_set.insert(make_pair(point_eval, punto));
    }
}

    void IKDTree::insert(Vector& punto){
        /* Se crea un nodo con el punto que se desea insertar*/
        Node* aux = new Node(punto);

        /* Si no existe la raiz, se inserta el punto como raiz*/
        if(raiz == nullptr){
            raiz = aux;
            npoints++;
            return;
        }

        // Se crea un nodo para ir recorriendo el arbol
        Node *nCurrent = raiz;
        while(true){

            /* Si el nodo actual no tiene hijos*/
            if(nCurrent->left == nullptr && nCurrent->right == nullptr){
                /* Se inserta el nuevo punto como hijo derecho */
                nCurrent->add_child(aux); 
                break;
            }

            /* Si solo tiene un hijo,
                se verifica si tiene al menos un rango disjunto con ese hijo*/
            else if(nCurrent->left == nullptr || nCurrent->right == nullptr){
                Node* child = (nCurrent->left)? nCurrent->left:nCurrent->right;

                /* Si el hijo contiene al punto, se actualiza el nodo actual al hijo*/
                if (child->range.contains(punto)){
                    nCurrent = child;
                    continue;
                }

                /* Se inserta el hijo,
                    se actualiza la dimension disjunta y ordena los hijos*/
                nCurrent->add_child(aux); 

                break; 
            }

            /* Si el nodo actual tiene ambos hijos*/
            else{
                
                int dim = nCurrent->dimension;
                /* Se busca la caja mas cercana al punto,
                    en relacion a la dimension disjunta, y se actualiza el nodo actual*/
                if(punto[dim] - nCurrent->left->range[dim].ub() <= nCurrent->right->range[dim].lb() - punto[dim] )
                    nCurrent = nCurrent->left;     
                else
                    nCurrent = nCurrent->right;
                
            }   
        }

        /* Se asigna el padre al nodo insertado,
            y se actualiza los rangos de los ancestros*/
        aux->padre = nCurrent;
        Node::update_ranges(aux);
        npoints++;

    }
   

    void Node::update_ranges(Node* nodo){
        Node* nodoP = nodo->padre;
        bool ranges_ready=false;

        /* Se recorren todos los ancestros del nodo*/
        while(nodoP != nullptr){         

            int h = nodoP->height;

            int hleft = -1;
            int hright= -1;

            if(nodoP->left) hleft = nodoP->left->height;
            if(nodoP->right) hright = nodoP->right->height;

            nodoP->height = std::max(hleft,hright) + 1;

            /* Si los rangos no se han actualizado*/
            if(!ranges_ready){  
                //se actualizan los rangos.
                if (nodo->range.is_subset(nodoP->range)) ranges_ready=true;
                else nodoP->range |= nodo->range; 
            }

            if(ranges_ready && nodoP->height == h) break;

            nodoP = nodoP->padre;
        }
    }


    Node* Node::search(Vector& punto){
        Node* actual = this;
        int visitados = 1;
        
        while(actual != nullptr){
            // Si el punto actual es el punto buscado se retorna
            if( *actual == punto){
                return actual;     
            }

            /* Si el nodo actual tiene ambos hijos*/
            if(actual->left && actual->right) 
                if(actual->left->range[actual->dimension].intersects(actual->right->range[actual->dimension])){
                    cout << actual->left->range[actual->dimension] << " -- " << actual->right->range[actual->dimension] << endl;
                    exit(0);
                }

            /* Si tiene un hijo derecho 
                y el punto buscado esta contenido en la caja de ese hijo*/
            if(actual->right && (actual->dimension==-1 ||
                         (actual->right && punto[actual->dimension] >= actual->right->range[actual->dimension].lb()))){

                actual = actual->right;
                visitados++;
                continue;

            } else {
                actual = actual->left;
                visitados++;
                continue;
            }

        }

        return nullptr;
    }


    void Node::add_child(Node* nodo){
        nodo->padre = this;

        /* Si el nodo actual no tiene hijos,
            se inserta en el hijo derecho*/
        if(left == nullptr && right == nullptr)
            right = nodo;

        /* Si posee un hijo*/
        else{
            Node* current_child;
            /* Si el hijo izquierdo no existe lo inserta ahi,
                y se guarda el hijo actual para sacar su max dimension disjunta*/
            if(left == nullptr){
                left = nodo;
                current_child = right;
            }

            /* Si el hijo derecho no existe, lo inserta ahi,
                y se guarda el hijo actual para sacar su max dimension disjunta*/
            else {
                right = nodo;
                current_child = left;
            }

            // Disjoint dimension
            /* Se actualiza la dimension disjunta del nodo*/
            int dim = current_child->max_disjoint_dimension(nodo->range);
            dimension = dim;

            /* Se ordenan los hijos acorde a la dimension disjunta del nodo*/
            if(right->range[dim].ub() <= left->range[dim].lb()){
                Node* aux = right;
                right = left;
                left = aux;
            }
        }
    }

    bool IKDTree::remove(Vector& punto){
        // Se busca el nodo que contiene el punto buscado
        Node* nodoPunto = raiz->search(punto);

        //int altura = nodoPunto->height;
        //cout << "heiGHTT:" << altura << endl;
        
        /* Si el punto que se desea eliminar no existe*/
        if(nodoPunto == nullptr) {
            cout << "not found:" << punto << endl;
            //exit(0); -> Solo para testear la eliminación; descomentar despues
            return false;
        }

        npoints--;

        Node* nodoPoint = nodoPunto;
        
        while(true){
            /* Si el nodo tiene dos hijos*/
            if(nodoPoint->left != nullptr && nodoPoint->right != nullptr){
                /* Si la altura del nodo izquierdo es mayor que la del derecho,
                    se continua con ese nodo*/
                if(nodoPoint->left->height > nodoPoint->right->height)
                    nodoPoint = nodoPoint->left;
                /* Si la altura del nodo derecho es mayor que la del izquierdo,
                    se continua con ese nodo*/
                else
                    nodoPoint = nodoPoint->right;

            /* Si tiene uno o cero hijos*/
            }else{
                nodoPunto->Vector::operator=(*nodoPoint); 
                
                /* Se elimina el nodo que contiene el punto,
                    del nodo padre*/
                Node* nodoPadre = nodoPoint->padre;
                if(nodoPadre){
                    if(nodoPadre->left == nodoPoint) nodoPadre->left=nullptr;
                    else nodoPadre->right=nullptr;
                    nodoPadre->dimension=-1;
                }
                
                /* Se verifica si el nodo a eliminar tiene
                    un hijo izquierdo, derecho o si no tiene hijo*/
                Node* nodoHijo = (nodoPoint->left)? nodoPoint->left: (nodoPoint->right)? nodoPoint->right: nullptr;
                
                /* Si el nodo si tiene un hijo,
                    se inserta el hijo al nodo padre*/
                if(nodoHijo != nullptr){
                    if(nodoPadre) nodoPadre->add_child(nodoHijo);
                    else raiz=nodoHijo;

                    /* Se actualizan los ancestros desde el nodo hijo*/
                    Node::update_ranges(nodoHijo);
                
                /* Si el nodo que se desea eliminar no tiene hijo*/
                }else {
                    /* Si existe el padre se actualizan los rangos*/
                    if(nodoPadre){
                        nodoHijo = (nodoPadre->left)? nodoPadre->left : (nodoPadre->right)? nodoPadre->right : nullptr;
                        if(nodoHijo) Node::update_ranges(nodoHijo);
                        else{
                            nodoPadre->range=*nodoPadre;
                            Node::update_ranges(nodoPadre);
                        }
                    }else raiz=NULL;
                }

                delete nodoPoint;
                return true;
            }
        }

        return false;
    }


    bool Node::children(){
        if(this->left != nullptr || this->right != nullptr)
            return true;
        return false;
    }


    //Retorna la dimensión en la que el punto es disjunto con el rango del nodo
    int Node::max_disjoint_dimension(IntervalVector& range2){

        double max_dist = 0;
        int dim = -1;

        for(int i=0 ; i<range.size() ; i++){
            /*Si el punto no está contenido en el rango, 
                se retorna la dimensión*/
            double dist = std::max(range[i].lb()-range2[i].ub(), range2[i].lb()-range[i].ub());
            if (dist > max_dist){
                max_dist=dist;
                dim = i;
            }
        }

        if(dim==-1) {
            cout << "no disjoint dimension" << endl;
            exit(0);
        }

        return dim;
    }  
    
    IKDTree::IKDTree(){
            raiz = nullptr;
    }


    Vector* IKDTree::search(Vector& punto){
        return raiz->search(punto);
    }

}