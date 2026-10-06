//============================================================================
//                                  I B E X
// File        : ibex_LSmearDFB.h
// Author      : (DFB) LSmear con los duales del contractor DFB
//============================================================================

#ifndef __IBEX_LSMEAR_DFB_H__
#define __IBEX_LSMEAR_DFB_H__

#include "ibex_SmearFunction.h"
#include "ibex_ExtendedSystem.h"
#include "ibex_CtcDFBHull.h"

namespace ibex {

/**
 * \brief LSmear con las sensibilidades del LP del objetivo de DFB.
 *
 * `LSmear` resuelve en cada nodo un LP con SoPlex (jacobiana en el punto
 * medio) para obtener los multiplicadores duales del objetivo, y puntua cada
 * variable por |sensibilidad| * diametro. DFB ya resolvio ese LP en el
 * contractor —sobre la relajacion afin, mas apretada— y su certificado
 * `gamma` da la misma sensibilidad: `|gamma_j / gamma_goal|`. Este bisector
 * usa eso y no llama a SoPlex. Si no hay sensibilidades validas para la caja
 * (no se resolvio la cota del objetivo, o la caja no esta contenida en la de
 * la linealizacion), cae a `SmearSumRelative`, igual que `LSmear`.
 *
 * `DFBH_BISEC_CMP=1`: ademas corre un `LSmear` de verdad sobre el mismo nodo
 * y cuenta cuantas veces eligen la misma variable (solo diagnostico).
 */
class LSmearDFB : public SmearSumRelative {
public:
    LSmearDFB(ExtendedSystem& sys, const Vector& prec, LargestFirst& lf, Bsc* comparar = NULL);
    /** Sin objetivo (solver): las sensibilidades son la suma sobre los `2n`
     *  LPs de la pasada (ver `CtcDFBHull`). */
    LSmearDFB(System& sys, double prec);
    virtual BisectionPoint choose_var(const Cell& cell);
    virtual int var_to_bisect(IntervalMatrix& J, const IntervalVector& box) const;
    static long n_llamadas, n_sin_sens, n_fuera, n_dfb, n_iguales;
private:
    System& esys;
    int gv;   //!< variable objetivo, o -1
    Bsc* comparar;
};

} /* namespace ibex */

#endif
