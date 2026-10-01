//============================================================================
//                                  I B E X                                   
// File        : ibex_Optimizer05Config.cpp
// Author      : Bertrand Neveu, Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Dec 11, 2014
// Last Update : Oct 15, 2019
//============================================================================

#include "ibex_Optimizer05Config.h"
#include "ibex_CtcDFBPropag.h"
#include "ibex_CtcDFBManager.h"
#include "ibex_CtcHC4.h"
#include "ibex_CtcAcid.h"
#include "ibex_Ctc3BCid.h"
#include "ibex_CtcPolytopeHull.h"
#include "ibex_CtcDFBHull.h"
#include "ibex_CtcCompo.h"
#include "ibex_CtcNewton.h"
#include "ibex_CtcFixPoint.h"
#include "ibex_CtcKuhnTucker.h"
#include "ibex_Exception.h"
#include "ibex_OptimLargestFirst.h"
#include "ibex_RoundRobin.h"
#include "ibex_SmearFunction.h"
#include "ibex_LSmear.h"
#include "ibex_Random.h"
#include "ibex_NormalizedSystem.h"
#include "ibex_LinearizerXTaylor.h"
#include "ibex_LinearizerCompo.h"
#include "ibex_LinearizerAffine2.h"
#include "ibex_CellHeap.h"
#include "ibex_CellDoubleHeap.h"
#include "ibex_CellBeamSearch.h"
#include "ibex_LoupFinderDefault.h"
#include "ibex_SyntaxError.h"

#include <sstream>
#include <vector>

using namespace std;

namespace ibex {

Optimizer05Config::Optimizer05Config(int argc, char** argv) {

	try {
		if (argc<8) {
			ibex_error("usage: optimizer05 filename filtering linear_relaxation bisection strategy [beamsize] prec goal_prec timelimit randomseed");
		}

		load_sys(argv[1]);

	} catch(SyntaxError& e) {
		ibex_error(e.msg.c_str());
	}

//	if (simpl_level)
//		sys->set_simplification_level(simpl_level.Get());

	filtering = argv[2];
	linearrelaxation = argv[3];
	bisection = argv[4];
	strategy = argv[5];

	int nbinput=5;

	beamsize = -1;
	if (strategy=="bs" || strategy== "beamsearch") { beamsize=atoi(argv[6]); nbinput++; }

	double prec= atof(argv[nbinput+1]);
	double goalprec= atof (argv[nbinput+2]);
	double timelimit= atof(argv[nbinput+3]);
	int randomseed = atoi(argv[nbinput+4]);

	RNG::srand(randomseed);

	// the extended system
	double eqeps=1e-8;
	ext_sys = &rec(new ExtendedSystem(*sys,eqeps));
	norm_sys = &rec(new NormalizedSystem(*sys,eqeps));

	/*
	cout << "file " << argv[1] << endl;
	cout << " filtering " << filtering;
        cout << " linearrelaxation " << linearrelaxation;
	cout << " bisection " << bisection ;
	cout << " strategy " << strategy ;
	cout << " randomseed " << randomseed << endl;
	 */

	if (getenv("IBEX_ARGS")) fprintf(stderr, "[args] eps_x=%.17g eps=%.17g\n", prec, goalprec);
	set_eps_x(Vector(sys->nb_var,prec));
	set_rel_eps_f(goalprec);
	set_abs_eps_f(goalprec);

	// the optimizer : the same precision goalprec is used as relative and absolute precision
	//	Optimizer o(sys->nb_var,*ctcxn,*bs,loupfinder,*buffer,ext_sys.goal_var(),1.e-7,goalprec,goalprec);
	// the trace
	set_trace(getenv("IBEX_TRACE") ? atoi(getenv("IBEX_TRACE")) : 0);   /* IBEX_TRACE=1: loup/uplo */

	// the allowed time for search
	set_timeout(timelimit);
}

unsigned int Optimizer05Config::nb_var() {
	return sys->nb_var;
}

void Optimizer05Config::load_sys(const char* filename) {
	std::size_t found = string(filename).find(".nl");

	if (found!=std::string::npos) {
		cerr << "\n\033[31mAMPL files can only be read with optimizer05 (ibex-opt-extra package).\n\n";
		exit(0);
	} else
		sys = &rec(new System(filename));
}

Linearizer* Optimizer05Config::get_linear_relax() {
	Linearizer* lr;

	/* DOS EJES SEPARADOS.
	 *
	 * `--lr` elegia a la vez QUE linearizador y QUE contractor lineal, asi que
	 * `--lr=art` significaba «afin + PolytopeHull» y no habia forma de pedir
	 * «afin + DFB», que es justamente la celda que interesa. `DFBH_LR` elige la
	 * relajacion con independencia del contractor; si no esta, se conserva el
	 * comportamiento anterior leyendo el modo de `--lr`. */
	/* AFIN POR OMISION.
	 *
	 * Medido sobre 139 instancias, 4 semillas, 3892 corridas: en las instancias
	 * donde produccion pasa de 1 s, cambiar SOLO el linearizador da celdas
	 * x0.322 y cpu x0.294 con PolytopeHull, y x0.322 / x0.229 con DFB. Los dos
	 * ejes son independientes y se multiplican. `compo` poda todavia mas
	 * (celdas x0.281) pero no lo paga (cpu x0.342), asi que no es el default.
	 *
	 * `DFBH_LR=xn` restaura xtaylor, para reproducir mediciones anteriores. */
	const char* env_lr = getenv("DFBH_LR");
	std::string modo = env_lr ? std::string(env_lr)
	                 : ((linearrelaxation=="art" || linearrelaxation=="compo")
	                    ? linearrelaxation : std::string("art"));

	if (modo=="art")
		lr = &rec(new LinearizerAffine2(*ext_sys));
	else if (modo=="compo")
		/* Compo escribe las restricciones de AMBOS en el MISMO LPSolver y suma
		 * las filas: el poliedro es la interseccion, no un reemplazo. */
		lr = &rec(new LinearizerCompo(rec(new LinearizerAffine2(*ext_sys)),
		                              rec(new LinearizerXTaylor(*ext_sys))));
	else //if (modo=="xn")
		lr = &rec(new LinearizerXTaylor(*ext_sys));
	/*	else {
			stringstream ss;
			ss << "[optimizer05] " << linearrelaxation << " is not an implemented relaxation mode ";
			ibex_error(ss.str().c_str());
		}
*/
	return lr;
}

Ctc& Optimizer05Config::get_ctc() {

	// the first contractor called
	CtcHC4* hc4 = &rec(new CtcHC4(ext_sys->ctrs,0.01,true));
	// hc4 inside acid and 3bcid : incremental propagation beginning with the shaved variable
	CtcHC4* hc44cid = &rec(new CtcHC4(ext_sys->ctrs,0.1,true));
	// hc4 inside xnewton loop
	CtcHC4* hc44xn = &rec(new CtcHC4(ext_sys->ctrs,0.01,false));

	/* DFBH_ENACID=1: el test barato de DFB se agrega al contractor que usan
	 * ACID y 3BCID para rebanar. Se construye ANTES de los dos para que compartan
	 * el mismo objeto `CtcDFBHull` y vean sus `gamma`.
	 *
	 * 3BCID es la prueba de una hipotesis: si el test pierde dentro de ACID
	 * porque perturba una heuristica ADAPTATIVA —ACID decide cuantas variables
	 * rebanar segun cuanto le rinde—, entonces en 3BCID, que no adapta nada,
	 * la poda extra deberia aparecer como ganancia neta. */
	Ctc* sub_shave = hc44cid;
	CtcDFBHull* hull_acid = NULL;
	if (getenv("DFBH_ENACID")) {
		hull_acid = &rec(new CtcDFBHull(*get_linear_relax(), ext_sys->goal_var()));
		sub_shave = &rec(new CtcCompo(*hc44cid, *&rec(new CtcDFBGamma(*hull_acid))));
	}

	// The 3BCID contractor on all variables (component of the contractor when filtering == "3bcidhc4")
	Ctc3BCid* c3bcidhc4 = &rec(new Ctc3BCid(*sub_shave));
	// hc4 followed by 3bcidhc4 : the actual contractor used when filtering == "3bcidhc4"
	CtcCompo* hc43bcidhc4 = &rec(new CtcCompo(*hc4, *c3bcidhc4));

	// The ACID contractor (component of the contractor  when filtering == "acidhc4")
	/* DFBH_ENACID=1: el test barato de DFB entra DENTRO de ACID, junto al HC4
	 * que ACID usa para rebanar. Cuesta 2n productos punto por trozo, sin
	 * pivotear ni re-linealizar, y usa los `gamma` del nodo anterior —validos
	 * solo si la caja de ahora esta contenida en aquella, cosa que el test
	 * verifica—. Medido sobre los trozos: mata el 11.9 % en `alkylbis` contra el
	 * 5.5 % de HC4, y son complementarios (§28). */
	CtcAcid* acidhc4 = &rec(new CtcAcid(*ext_sys, *sub_shave, true));
	// hc4 followed by acidhc4 : the actual contractor used when filtering == "acidhc4"
	CtcCompo* hc4acidhc4 = &rec(new CtcCompo(*hc4, *acidhc4));

    Linearizer* lr = get_linear_relax();
	CtcDFBPropag* dfb = &rec(new CtcDFBPropag(*ext_sys, *lr, 0.1, true));
	CtcDFBPropag* dfb_hc4 = &rec(new CtcDFBPropag(*ext_sys, *lr, 0.1, true, true)); //hc4

    /* DFB puro, sin los HC4 internos: es el analogo fiel de CtcPolytopeHull
     * para la sustitucion directa (--lr=dfb). PolyHull es un contractor lineal
     * a secas y el HC4 del bucle lo aporta el envoltorio, asi que meter aca un
     * CtcDFBPropag con HC4 adentro contaria HC4 dos veces. */
    CtcDFBPropag* dfb_solo = &rec(new CtcDFBPropag(*ext_sys, *lr, 0.01, true, false, true));

    CtcDFBPropag* propag_dfb = &rec(new CtcDFBPropag(*ext_sys, *lr, 0.01, false));
	CtcAcid* acid_dfb = &rec(new CtcAcid(*ext_sys, *propag_dfb, true));
	CtcDFBManager* dfb_manager = &rec(new CtcDFBManager(*propag_dfb, *acid_dfb, *lr));

	Ctc* ctc;
	if (filtering == "hc4")
		ctc = hc4;
	else if (filtering == "dfb")
	    ctc = dfb;
	else if (filtering == "dfb_hc4")
	    ctc = dfb_hc4;
    else if (filtering == "acid_dfb")
        ctc = dfb_manager;
	else if (filtering =="acidhc4")
		ctc = hc4acidhc4;
	else if (filtering =="3bcidhc4")
		ctc = hc43bcidhc4;
	else {
		stringstream ss;
		ss << "[optimizer05] " << filtering << " is not an implemented  contraction mode ";
		ibex_error(ss.str().c_str());
	}

	
	/* DFBH_RATIO=r: umbral del punto fijo (omision 0.2). Como el nodo con DFB
	 * sale ~15 % mas barato (§38) y poda ~4.9 % menos, lo natural es gastar el
	 * ahorro en iterar mas el punto fijo: mas poda, pagada con el ahorro. */
	const double ratio_fp = getenv("DFBH_RATIO") ? atof(getenv("DFBH_RATIO"))
	                                             : relax_ratio;
	Ctc* cxn;
	CtcPolytopeHull* cxn_poly;
	CtcCompo* cxn_compo;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn") {
		cxn_poly = &rec(new CtcPolytopeHull(*lr));
		cxn_compo = &rec(new CtcCompo(*cxn_poly, *hc44xn));
	        if (sys->nb_ctr==0)
		  cxn = cxn_poly;
		else if (getenv("DFBH_NOFIX"))
		  cxn = cxn_compo;          /* sin iteracion, para comparar */
		else
		  cxn = &rec(new CtcFixPoint (*cxn_compo, ratio_fp));
	}
	/* --lr=dfb: la sustitucion directa que el objetivo pide. DFB ocupa el lugar
	 * exacto de CtcPolytopeHull, con el mismo envoltorio
	 * CtcFixPoint(CtcCompo(., hc44xn)) y el mismo relax_ratio, de modo que
	 * `--filtering=acidhc4 --lr=dfb` es la configuracion de produccion con lo
	 * unico cambiado siendo el contractor lineal. Las configuraciones
	 * --filtering=dfb y --filtering=acid_dfb cambian ademas la composicion, asi
	 * que no aislan la sustitucion. */
	/* `--lr=dfbhull`: la estrategia del §14 de MEDICIONES_DUELO.md en el hueco
	 * EXACTO de CtcPolytopeHull, con el mismo envoltorio y una sola pasada por
	 * llamada. `--filtering=acidhc4 --lr=dfbhull` es produccion con lo unico
	 * cambiado que se quiere medir: quien resuelve el LP. */
	else if (linearrelaxation=="dfbhull") {
		/* Si el test barato esta dentro de ACID, el contractor lineal tiene que
		 * ser el MISMO objeto, para que ACID vea los `gamma` que este calcula. */
		CtcDFBHull* cxn_dfbh = hull_acid ? hull_acid
		    : &rec(new CtcDFBHull(*lr, ext_sys->goal_var(),
		                          LPSolver::default_tolerance, hc4));
		cxn_compo = &rec(new CtcCompo(*cxn_dfbh, *hc44xn));
		/* SIN PUNTO FIJO. Sirve para saber si la propagacion intercalada dentro
		 * de la pasada es redundante con la iteracion del CtcFixPoint: si lo
		 * es, quitando la iteracion deberia pasar a rendir. `DFBH_NOFIX=1`. */
		if (sys->nb_ctr==0 || getenv("DFBH_NOFIX"))
		  cxn = (sys->nb_ctr==0) ? (Ctc*)cxn_dfbh : (Ctc*)cxn_compo;
		else
		  cxn = &rec(new CtcFixPoint (*cxn_compo, ratio_fp));
	}
	else if (linearrelaxation=="dfb") {
		cxn_compo = &rec(new CtcCompo(*dfb_solo, *hc44xn));
		if (sys->nb_ctr==0)
		  cxn = dfb_solo;
		else
		  cxn = &rec(new CtcFixPoint (*cxn_compo, ratio_fp));
	}
	//  the actual contractor  ctc + linear relaxation
	Ctc* ctcxn;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn"
	    || linearrelaxation=="dfb" || linearrelaxation=="dfbhull")
		/* DFBH_LRFIRST=1: el contractor lineal ANTES de ACID, para que ACID vea
		 * los `gamma` de SU nodo. Con el orden normal el test barato se rechaza por
		 * contencion el 85-98 % de las veces —los `gamma` disponibles vienen de otro
		 * nodo del recorrido— y queda en puro costo. OJO: esto cambia la
		 * composicion, asi que ya no es la sustitucion pura del §19. */
		if (getenv("DFBH_LRFIRST"))
			ctcxn = &rec(new CtcCompo  (*cxn, *ctc));
		else
			ctcxn = &rec(new CtcCompo  (*ctc, *cxn));
	else
		ctcxn = ctc;
	Ctc* ctckkt;

	if (sys->nb_ctr == 0){
	  ctckkt = &rec(new CtcKuhnTucker(*norm_sys, true));
	  ctcxn = &rec(new CtcCompo  (*ctcxn, *ctckkt));
	  }

	return *ctcxn;
}

Bsc& Optimizer05Config::get_bsc() {
	Bsc* bs;

	const Vector& eps_x=get_eps_x();
	Vector prec(ext_sys->nb_var);

	if (eps_x.size()==1) // not really initialized
		ext_sys->write_ext_vec(Vector(ext_sys->nb_var,eps_x[0]), prec);
	else
		ext_sys->write_ext_vec(eps_x, prec);

	// TODO: should the following value be set to "abs_eps_f" instead?
	// This question is probably related to the discussion #400
	prec[ext_sys->goal_var()] = OptimizerConfig::default_eps_x;


	if (bisection=="roundrobin")
		bs = &rec(new RoundRobin (prec,0.5));
	else if (bisection== "largestfirst")
                bs = &rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.5));
	else if (bisection== "largestfirstnoobj")
                bs = &rec(new OptimLargestFirst(ext_sys->goal_var(),false,prec,0.5));
	else if (bisection=="smearsum")
		bs = &rec(new SmearSum(*ext_sys,prec,
				       rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.5))));
	else if (bisection=="smearmax")
		bs = &rec(new SmearMax(*ext_sys,prec,
				       rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.5))));
	else if (bisection=="smearsumrel")
                bs = &rec(new SmearSumRelative(*ext_sys,prec,
					       rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.5))));
	else if (bisection=="smearmaxrel")
		bs = &rec(new SmearMaxRelative(*ext_sys,prec,
					       rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.5))));
	else if  (bisection=="lsmear")
                bs = &rec (new LSmear(*ext_sys,prec,
				      rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.5)),
				      LSMEAR));
	else if (bisection=="lsmearmg")
	        bs = &rec (new LSmear(*ext_sys,prec,
				      rec(new OptimLargestFirst(ext_sys->goal_var(),true,prec,0.5))));
	else {
		stringstream ss;
		ss << "[optimizer05] " << bisection << " is not an implemented  bisection mode ";
		ibex_error(ss.str().c_str());
	}

	return *bs;
}


LoupFinder& Optimizer05Config::get_loup_finder() {
	return rec(new LoupFinderDefault(*norm_sys, true));
	//LoupFinderDefault loupfinder (norm_sys,false);
}

CellBufferOptim& Optimizer05Config::get_cell_buffer() {
	CellBufferOptim* buffer;

	CellHeap* futurebuffer = &rec(new CellHeap(*ext_sys));
	CellHeap* currentbuffer = &rec(new CellHeap(*ext_sys));

	//cout << "strategy=[" << strategy << "]\n";
	if (strategy=="bfs")
		buffer = &rec(new CellHeap(*ext_sys));
	else if (strategy=="dh")
		buffer = &rec(new CellDoubleHeap(*ext_sys));
	else if (strategy=="bs")
		buffer = &rec(new CellBeamSearch(*currentbuffer, *futurebuffer, *ext_sys, beamsize));

	return *buffer;
}

int Optimizer05Config::goal_var() {
	return ext_sys->goal_var();
}

} // namespace ibex
