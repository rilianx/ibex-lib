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
        Interval ev;
        const IntervalVector* box;


        Affine(int n) : a(n,Interval(0,0)), err(0), ev(0), box(NULL) {
        }

        Affine(const IntervalVector& box, int i=-1) : a(box.size(),Interval(0,0)), err(0), ev(0), box(&box) {
            if (i!=-1){
                a[i] = Interval(1.0);
                ev = box[i];
            }
        }


        //operator=
        Affine& operator=(const Affine& af) {
            a = af.a;
            err = af.err;
            ev = af.ev;
            box = af.box;
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
            c.ev += c.err;
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

        //aff*aff
        /*
                new_a0 = self.a[0] * other.a[0]
        new_a = [self.a[0] * bi + ai * other.a[0] for ai, bi in zip(self.a[1:], other.a[1:])]
        new_error = self.error * other.to_interval() + other.error * self.to_interval()
        y = self.x - self.x.mid()
        other_y = other.x - other.x.mid()


        new_error += sum(self.a[i] * other.a[j] * y[i-1] * other_y[j-1] for i in range(1, len(self.a)) for j in range(1, len(other.a)))
        return AF([new_a0] + new_a, self.x, new_error)*/
        friend Affine operator*(const Affine& a, const Affine& b) {
            if (a.box != b.box) {
                throw std::invalid_argument("Boxes must be the same");
            }
            Affine c(a);
            c.a[0] = a.err.mid() * b.err.mid();

            for (int i=0; i<a.a.size(); i++) 
                c.a[i] = a.a[0] * b.a[i] + a.a[i] * b.a[0];
            
            c.err = a.err * b.ev + b.err * a.ev;
            for (int i=1; i<a.a.size(); i++) {
                for (int j=0; j<a.a.size(); j++) {
                    c.err += a.a[i] * b.a[j] * (*a.box - a.box->mid()) * (*b.box - b.box->mid());
                }
            }
            c.ev = c.a*(*c.box) + c.err;
            return c;
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
            Affine c(*a.box, -1);
            c.a = -a.a; c.err = -a.err; c.ev = -a.ev;
            return c;
        }
        
        friend Affine operator-(const Affine& a, const Affine& b) {
            return a + (-b);
        }

        Affine cheby_convex(const Affine& af, const std::function<Interval(Interval)> f, 
                            const std::function<Interval(Interval)> df_proj) const {
            Interval x = af.ev;
            Interval m = (f(x.ub())-f(x.lb()))/x.diam();
            Interval xp = df_proj(m);
            Interval A = f(x.lb())-m*x.lb();
            Interval B = f(xp)-m*xp;

            return m*af + (A+B)/2 + Interval((-abs(A-B)/2).lb(),(abs(A-B)/2).ub());
        }

        friend Affine operator/(const Interval& a, const Affine& b) {
            return a * pow(b,-1);
        }

        friend Affine operator/(const double& a, const Affine& b) {
            return Interval(a) * pow(b,-1);
        }

        //sqr
        friend Affine sqr(const Affine& a) {
            return pow(a,2);
        }

        friend Affine pow(const Affine& a, const int& b) {
            if (b==2) {
                std::function<Interval(Interval)> f = [](Interval x){ return sqr(x); };
                std::function<Interval(Interval)> df_proj = [](Interval m){ return m/2; };

                return a.cheby_convex(a, f, df_proj);
            } else if (b==-1) {
                std::function<Interval(Interval)> f = [](Interval x){ return 1/x; };
                std::function<Interval(Interval)> df_proj = [](Interval m){ return sqrt(-1/m); };

                return a.cheby_convex(a, f, df_proj);
            } else {
                throw std::invalid_argument("exponent is not implemented");
            }
        }

    };
}

#endif // __IBEX_AFFINE_ARITH_H__