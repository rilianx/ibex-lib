#include "ibex.h"
#include "../src/bin/parse_args.h"
#include "time.h"
#include "ibex_Function.h"
#include "ibex_KDTree.h"

#include <cstdlib>
#include <cctype>

#include "iostream"
#include <string>
#include <math.h>

#include <ctime>
#include <chrono>

#include <sstream>

using namespace std;
using namespace ibex ;
using namespace kdtree ;

std::vector<Vector> points_in_hypersphere(int n, double r, int cantidadPts) {
    std::vector<Vector> puntos;
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    for (int i = 0; i < cantidadPts; ++i) {
        Vector punto(n);
        double norma = 0;

        // Generar puntos en un espacio de dimensión superior
        for (int j = 0; j < n; ++j) {
            punto[j] = (2.0 * std::rand() / RAND_MAX - 1.0) * r;
            norma += punto[j] * punto[j];
        }

        // Normalizar para proyectar sobre la hiperesfera
        norma = std::sqrt(norma);
        for (int j = 0; j < n; ++j) {
            punto[j] *= r / norma;
        }

        puntos.push_back(punto);
    }

    return puntos;
}

std::vector<Vector> boxes_in_hypersphere(int n, double r, double boxsize, int cantidadPts) {
    std::vector<Vector> boxes_as_points;
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    for (int i = 0; i < cantidadPts; ++i) {
        double norma = 0;
        Vector punto(n);

        // Generar puntos en un espacio de dimensión superior
        for (int j = 0; j < n; ++j) {
            punto[j] = (2.0 * std::rand() / RAND_MAX - 1.0) * r;
            norma += punto[j] * punto[j];
        }

        // Normalizar para proyectar sobre la hiperesfera
        norma = std::sqrt(norma);
        for (int j = 0; j < n; ++j) {
            punto[j] *= r / norma;
        }

        IntervalVector box(punto);
        box.inflate(boxsize);
        Vector box_v(n*2);

        for (int k=0;k<2*n;k+=2){
            box_v[k] = box[k/2].lb();
            box_v[k+1] = box[k/2].ub();
        }
        
        boxes_as_points.push_back(box_v);
    }

    return boxes_as_points;
}

IntervalVector input_bounds(int n){
    IntervalVector box(n);
    string input;
    std::getline(std::cin, input);
    for (int i=0;i<n;i++){
        
        double a=NEG_INFINITY,b=POS_INFINITY;
        cout << "ingrese los límtes para la dimensión " << i << " (ej: 0 1): ";
        
        std::getline(std::cin, input);
        if (input != "") {
            std::istringstream iss(input);
            iss >> a >> b;
        }
        box[i] = Interval(a,b);
    }
    return box;
}

Function* input_function(int n) {
    // Pedir al usuario que ingrese la función
    std::string input;
    std::string function_str;
    std::cout << "Ingrese la función a minimzar o presione enter en para satisfacción de restricciones: ";
    std::getline(std::cin, input);
    
    if (input != "") function_str = input;
    else return NULL;
    
    // Definir la función a partir de la entrada del usuario
    Function* f = new Function(string("x[" + std::to_string(n) + "]").c_str(), function_str.c_str());
    return f;
}

list<Function> input_constraints(int n) {
    list<Function> ctrs;

    while (true){
        // Pedir al usuario que ingrese una restricción g(x) <= 0
        std::string input;
        std::string function_str;
        std::cout << "Ingrese una nueva restricción g(x)<=0 o presione enter para terminar (e.g., 'x[0] + x[1]'): ";
        std::getline(std::cin, input);
        
        if (input != "") function_str = input;
        else break;

        // Definir la función a partir de la entrada del usuario
        Function f(string("x[" + std::to_string(n) + "]").c_str(), function_str.c_str());
        ctrs.push_back(f);
    }
    return ctrs;
}

list<Function> collision_contraints(IntervalVector& y){
    list<Function> ctrs;
    int dim = y.size();   

    for (int i=0;i<dim;i++){
        ctrs.push_back(Function(("x["+std::to_string(dim*2)+"]").c_str(), ("x["+std::to_string(2*i)+"]-" + std::to_string(y[i].ub())).c_str()));
        ctrs.push_back(Function(("x["+std::to_string(dim*2)+"]").c_str(), (std::to_string(y[i].lb()) + "-x["+std::to_string(2*i+1)+"]").c_str()));
    }
    return ctrs;
}

Function* knn_function(int n){
    string input;
    Vector y(n, 0.0);
    cout << "Ingrese el punto de referencia: (default: 0 0 0): ";
    std::getline(std::cin, input);
    std::getline(std::cin, input);
    if (input != "") {
        std::istringstream iss(input);
        for (int i=0;i<n;i++) iss >> y[i];
    }            

    string dist_funct="";
    for (int i=0;i<n;i++){
        dist_funct += "(x["+std::to_string(i)+"]-"+std::to_string(y[i])+")^2";
        if (i<n-1) dist_funct += "+";
    }
    return new Function(string("x[" + std::to_string(n) + "]").c_str(), dist_funct.c_str());
}

int main(int argc, char** argv){

    args::ArgumentParser parser("********* KD-Tree *********.", "Tests for the interval kd-tree structure");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Flag _generic_search(parser, "generic", "Búsqueda Genérica", {"generic"});
    args::Flag _knn(parser, "int", "K vecinos más cercano", {"knn"});
    args::Flag _collision_detection(parser, "int", "Detección de colisiones", {"collision"});
    args::ValueFlag<int> _seed(parser, "int", "Seed for random generator", {"seed"});

	try	{ parser.ParseCLI(argc, argv); 	}
	catch (args::Help&) { std::cout << parser; return 0; }
	catch (args::ParseError& e)	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}catch (args::ValidationError& e){
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}

    if(_seed){
        srand(_seed.Get());
    }

    IKDTree kdtree;

    int n;
    cout << "Ingrese número de dimensiones: "; 
    cin >> n;

    if(_generic_search || _knn || _collision_detection){
        if(_collision_detection) n *= 2;
        

        Function* f = NULL;
        if(_generic_search || _collision_detection){
            f = input_function(n);
        }else if(_knn){
            f = knn_function(n);
        }

        IntervalVector bounds = input_bounds(n);
        IntervalVector y(n/2, 0.0); //para colisiones

        list<Function> ctrs;
        if(_generic_search || _knn)
            ctrs = input_constraints(n);
        else if(_collision_detection){
            cout << "Ingrese caja de referencia: (ejemplo: 0 1 0 1 0 1): ";
            for (int i=0;i<n/2;i++){
                double lb,ub;
                cin >> lb >> ub;
                y[i] = Interval(lb,ub);
            }

            ctrs = collision_contraints(y); 
        }

        if (f != NULL)
            cout << "Objective function (minimization): " << *f << endl;  

        for (auto it = ctrs.begin(); it != ctrs.end(); ++it){
            cout << "Constraint: " << *it << " <= 0" << endl;
        }

        //input the number of points
        cout << "Ingrese la cantidad de elementos a generar (hypershpere): ";
        int cantidadPts; cin >> cantidadPts;

        vector<Vector> points;
        if (_collision_detection){
            points = boxes_in_hypersphere(n/2,1.0,0.05,cantidadPts);
        }else 
            points = points_in_hypersphere(n,1.0,cantidadPts);

        if(f!=NULL){
            double min_eval = POS_INFINITY;
            double min_i = -1;

            for (int i=0;i<cantidadPts;i++){
                if (IKDTree::satisfiesConstraints(points[i], ctrs) && bounds.contains(points[i])){
                    double eval = f->eval(points[i]).mid();
                    if (eval < min_eval){
                        min_eval = eval;
                        min_i = i;
                    }
                }
            }

            if(min_i != -1){
                cout << "Punto mínimo: " << points[min_i] << endl;
                cout << "Valor mínimo (f): " << min_eval << endl;
            }else{
                cout << "No se encontró un punto que satisfaga las restricciones" << endl;
            }
        }

        
        for (int i=0;i<cantidadPts;i++)
            kdtree.insert(points[i]);
        
        cout << "\nUsando Interval KD-Tree" << endl;
        cout << "Ingrese la máxima cantidad de puntos a buscar (default: inf): ";
        int k=100000; string input;
        std::getline(std::cin, input);
        std::getline(std::cin, input);
        if (input != "") {
            std::istringstream iss(input);
            iss >> k;
        }

        int altura = kdtree.height();
        cout << "Altura: " << altura/log2(cantidadPts) << "*log(n)" << endl;

        set<pair <double, const Vector*>> points_ = kdtree.get_kmin(f, ctrs, k, bounds);
        if (points_.size() == 0){
            cout << "No se encontró un punto que satisfaga las restricciones" << endl;
            return 0;
        }else{
            for (auto it = points_.begin(); it != points_.end(); ++it){
                if (_collision_detection){
                    //transform to box
                    IntervalVector box(n/2);
                    for (int i=0;i<n/2;i++){
                        box[i] = Interval((*it->second)[2*i], (*it->second)[2*i+1]);
                    }
                    cout << "Box: " << box << " .. " << box.intersects(y) ;
                }else 
                    cout << "Punto: " << *it->second << " .. ";

                if (f != NULL)
                    cout << "Eval (f): " << it->first << endl;
                else 
                    cout << endl;
            }
        }
    }

    

    
}
