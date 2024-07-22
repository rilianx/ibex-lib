/* ============================================================================
 * I B E X - Function basic evaluation
 * ============================================================================
 * Copyright   : IMT Atlantique (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Ignacio Araya
 * Created     : 
 * Last update : 
 * ---------------------------------------------------------------------------- */

#ifndef __IBEX_AFFINE_ARITH_H__
#define __IBEX_AFFINE_ARITH_H__

#include <iostream>
#include "ibex.h"
#include "ibex_ExprDomain.h"
#include "ibex_Function.h"
#include "ibex_Interval.h"
#include <vector>
#include <typeinfo>

using namespace std;

namespace ibex {
    class Affine {
        public:

        IntervalVector a;
        Interval err;
        Interval ev; //auxiliar
        bool is_constant; //err=ev 
        const IntervalVector* box;


        Affine(int n) : a(n,Interval(0,0)), err(0), ev(0), box(NULL), is_constant(true) {
        }

        Affine(const IntervalVector& box, int i=-1) : a(box.size(),Interval(0,0)), err(0), ev(0), box(&box) {
            if (i!=-1){
                a[i] = Interval(1.0);
                ev = box[i];
                is_constant = false;
            }else 
                is_constant = true;
        }


        //operator=
        Affine& operator=(const Affine& af) {
            a = af.a;
            err = af.err;
            ev = af.ev;
            box = af.box;
            is_constant = af.is_constant;
            return *this;
        }

     //operator<<
        friend std::ostream& operator<<(std::ostream& os, const Affine& af) {
            os << "Affine(";
            for (int i=0; i<af.a.size(); i++) {
                if (af.a[i] == 0) continue;
                else os << af.a[i] << "*x_" << i+1;
                os << " + ";
            }
            os << af.err << ")";
            return os;
        }

    //add overload based on python code
        friend Affine operator+(const Affine& a, const Affine& b) {
            if (a.box != b.box) {
                throw std::invalid_argument("Boxes must be the same");
            }
            Affine c(a);
            c.a += b.a;
            c.err += b.err;
            c.ev += b.err;
            c.is_constant = a.is_constant && b.is_constant;
            return c;
        }

        friend Affine operator+(const Affine& a, const Interval& b) {
            Affine c(a);
            c.err += b; c.ev += b;
            return c;
        }

        friend Affine operator+(const Interval& b, const Affine& a) {
            return a + b;
        }

        friend Affine operator+(const Affine& a, const double& b) {
            return a + Interval(b);
        }

        friend Affine operator+(const double& b, const Affine& a) {
            return a + Interval(b);
        }

        friend Affine operator*(const Affine& a, const Affine& b) {
            if (a.box != b.box) {
                throw std::invalid_argument("Boxes must be the same");
            }

            if (a.is_constant) 
                return b*a.ev;
            else if (b.is_constant) 
                return a*b.ev;
            else
                //not implemented
                throw std::invalid_argument("Not implemented yet multiplication of affine functions");
        }

        friend Affine operator*(const Affine& a, const Interval& b) {
            Affine c(a);
            c.a *= b; c.err *= b; c.ev *= b;
            return c;
        }

        friend Affine operator*(const Interval& b, const Affine& a) {
            return a*b;
        }

        friend Affine operator*(const Affine& a, const double& b) {
            return a * Interval(b);
        }

        friend Affine operator*(const double& b, const Affine& a) {
            return a * Interval(b);
        }

        friend Affine operator-(const Affine& a) {
            Affine c(a);
            c.a = -c.a; c.err = -c.err; c.ev = -c.ev;
            return c;
        }
        
        friend Affine operator-(const Affine& a, const Affine& b) {
            return a + (-b);
        }

        friend Affine cheby_convex(const Affine& af, const std::function<Interval(Interval)> f, 
                            const std::function<Interval(Interval)> df_proj) {
            Interval x = af.ev;
            Interval m = (f(x.ub())-f(x.lb()))/x.diam();
            Interval xp = df_proj(m);
            Interval A = f(x.lb())-m*x.lb();
            Interval B = f(xp)-m*xp;
            return m*af + (A+B)/2 + Interval((-abs(A-B)/2).lb(),(abs(A-B)/2).ub());
            //return m*af + Interval(A,B)
        }

        friend Affine operator/(const Interval& a, const Affine& b) {
            if(b.ev.lb()>0 || b.ev.ub()<0)
                return a * pow(b,-1);
            else
                throw std::invalid_argument("Division by zero");
        }

        friend Affine operator/(const double& a, const Affine& b) {
            return Interval(a) * pow(b,-1);
        }

        //sqr
        friend Affine sqr(const Affine& a) {
            return pow(a,2);
        }
        
        friend Affine pow(const Affine& a, const int& b) {
            if (b % 2 == 0 && b > 0) {
                std::function<Interval(Interval)> f = [b](Interval x){ return pow(x, b); };
                std::function<Interval(Interval)> df_proj = [b](Interval m){ return m/b; };
                return cheby_convex(a, f, df_proj);
            } else if (b==-1) {
                std::function<Interval(Interval)> f = [](Interval x){ return 1/x; };
                std::function<Interval(Interval)> df_proj = [](Interval m){ return sqrt(-1/m); };
                return cheby_convex(a, f, df_proj);
            } else if (b > 0 && b % 2 != 0) {
                if (a.ev.lb() >= 0 || a.ev.ub() <= 0 ){
                    std::function<Interval(Interval)> f = [b](Interval x){ return pow(x, b); };
                    std::function<Interval(Interval)> df_proj = [b](Interval m){ return m/b; };
                    return cheby_convex(a, f, df_proj);
                } else {
                    return Affine(*a.box, -1) + pow(a.ev,b);
                }
            } else {
                // Exponente negativo distinto a -1.
                throw std::invalid_argument("exponent is not implemented");
            }
        }

        friend Affine exp(const Affine& a) {
            return 1;
            // std::function<Interval(Interval)> f = [](Interval x){ return ibex::exp(Interval::E.lb(), x); };
            // std::function<Interval(Interval)> df_proj = [](Interval m){ return exp(e, m); };
            // return cheby_convex(a, f, df_proj);
        }

        friend Affine log(const Affine& a) {
            if (a.ev.lb() < 0 || a.ev.ub() < 0) {
                throw std::invalid_argument("the function is not defined in the interval");
            } else {
                throw std::invalid_argument("the function is not convex on the interval");
            }
        }

        friend Affine cos(const Affine& a) {
            // Es convexa cuando su segunda derivada (-cos(x)) es no negativa.
            // [pi/2, 3pi/2], [5pi/2, 7pi/2], [9pi/2, 11pi/2]...
            if (a.ev.diam() <= Interval::PI.lb()) {
                if (std::cos(a.ev.ub()) <= 0 && std::cos(a.ev.lb()) <= 0) {
                    std::function<Interval(Interval)> f = [](Interval x){ return ibex::cos(x); };
                    std::function<Interval(Interval)> df_proj = [](Interval m){ return ibex::sin(m); };
                    return cheby_convex(a, f, df_proj);
                }
            } else {
                throw std::invalid_argument("the function is not convex on the interval");  
            }
        }

        friend Affine sin(const Affine& a) {
            // Es convexa cuando su segunda derivada (-sin(x)) es no negativa.
            // [pi, 2pi], [3pi, 4pi], [5pi, 6pi]...
            if (a.ev.diam() <= Interval::PI.lb()) {
                if (std::sin(a.ev.ub()) <= 0 && std::sin(a.ev.lb()) <= 0) {
                    std::function<Interval(Interval)> f = [](Interval x){ return ibex::sin(x); };
                    std::function<Interval(Interval)> df_proj = [](Interval m){ return ibex::cos(m); };
                    return cheby_convex(a, f, df_proj);
            } else {
                throw std::invalid_argument("the function is not convex on the interval");
            }}
        }

        friend Affine tan(const Affine& a) {
            if (a.ev.diam() < Interval::PI.lb() / 2) {
                if (std::tan(a.ev.ub()) >= 0 && std::tan(a.ev.lb()) >= 0) {
                    std::function<Interval(Interval)> f = [](Interval x){ return ibex::tan(x); };
                    // std::function<Interval(Interval)> df_proj = [](Interval m){ return ibex::acos(m); };  // 1/m = cos^2(x) -> despejar x
                    // return cheby_convex(a, f, df_proj);
                    return 1;
                }}

            return 1;
        }

        friend Affine cosh(const Affine& a) {
            std::function<Interval(Interval)> f = [](Interval x){ return ibex::cosh(x); };
            std::function<Interval(Interval)> df_proj = [](Interval m){ return ibex::sinh(m); };
            return cheby_convex(a, f, df_proj);
        }

        friend Affine sinh(const Affine& a) {
            if (a.ev.lb() >= 0) {
                std::function<Interval(Interval)> f = [](Interval x){ return ibex::sinh(x); };
                std::function<Interval(Interval)> df_proj = [](Interval m){ return ibex::cosh(m); };
                return cheby_convex(a, f, df_proj);
            }
        }

        friend Affine tanh(const Affine& a) {
            if (a.ev.ub() <= 0) {
                std::function<Interval(Interval)> f = [](Interval x){ return ibex::tanh(x); };
                // std::function<Interval(Interval)> df_proj = [](Interval m){ return ibex::sinh(m); };  // 1/m = cosh^2(x) -> despejar x
                // return cheby_convex(a, f, df_proj);
            } else {
                throw std::invalid_argument("the function is not convex on the interval");
            }
        }

        friend Affine acos(const Affine& a) {
            if (a.ev.ub() <= 0 && a.ev.lb() >= -1) {
                std::function<Interval(Interval)> f = [](Interval x){ return ibex::acos(x); };
                std::function<Interval(Interval)> df_proj = [](Interval m){ return -(1/sqrt(1-pow(m, 2))); };
                return cheby_convex(a, f, df_proj);
            } else {
                throw std::invalid_argument("the function is not convex on the interval");
            }
        }

        friend Affine asin(const Affine& a) {
            if (a.ev.lb() >= 0 && a.ev.ub() <= 1) {
                std::function<Interval(Interval)> f = [](Interval x){ return ibex::asin(x); };
                std::function<Interval(Interval)> df_proj = [](Interval m){ return 1/sqrt(1-sqr(m)); };
                return cheby_convex(a, f, df_proj);
            } else {
                throw std::invalid_argument("the function is not convex on the interval");
            }
        }

        friend Affine atan(const Affine& a) {
            if (a.ev.ub() <= 0) {
                std::function<Interval(Interval)> f = [](Interval x){ return ibex::atan(x); };
                std::function<Interval(Interval)> df_proj = [](Interval m){ return 1/(sqr(m)+1); };  
                return cheby_convex(a, f, df_proj);
            } else {
                throw std::invalid_argument("the function is not convex on the interval");
            }
        }

        friend Affine acosh(const Affine& a) {
            // Es siempre concava?
            throw std::invalid_argument("the function is not convex on the interval");
        }

        friend Affine asinh(const Affine& a) {
            if (a.ev.ub() <= 0) {
                std::function<Interval(Interval)> f = [](Interval x){ return ibex::asinh(x); };
                std::function<Interval(Interval)> df_proj = [](Interval m){ return ibex::acosh(m); };  
                return cheby_convex(a, f, df_proj);
            } else {
                throw std::invalid_argument("the function is not convex on the interval");
            }
        }

        friend Affine atanh(const Affine& a) {
            if (a.ev.lb() >= 0) {
                std::function<Interval(Interval)> f = [](Interval x){ return ibex::atanh(x); };
                std::function<Interval(Interval)> df_proj = [](Interval m){ return 1/(1-sqr(m)); };  
                return cheby_convex(a, f, df_proj);
            } else {
                throw std::invalid_argument("the function is not convex on the interval");
            }
        }

        friend Affine floor(const Affine& a) {
            // TODO
            // Como determinar si es convexa o concava?
            return 1;
        }

        friend Affine ceil(const Affine& a) {
            // TODO
            // Como determinar si es convexa o concava?
            return 1;
        }
    };
}

#endif // __IBEX_AFFINE_ARITH_H__