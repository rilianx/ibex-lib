// Aplica UN contractor (CtcSC o CtcPolytopeHull) sobre una caja leida en
// hexadecimal (%a) e informa si la vacia y, para CtcSC, cuantos LPs dieron
// multiplicadores no finitos (NaN).
//   ./sc_una sistema.mbx caja.txt art|xt sc|ph
#include "ibex.h"
#include "ibex_CtcSC.h"
#include "ibex_LinearizerAffine2.h"
#include <cstdio>
#include <string>
using namespace ibex;
int main(int argc, char** argv) {
	if (argc < 5) { fprintf(stderr, "uso: sc_una sistema.mbx caja.txt art|xt sc|ph\n"); return 1; }
	System sys(argv[1]); int n = sys.nb_var; IntervalVector b(n);
	FILE* f = fopen(argv[2], "r");
	for (int i = 0; i < n; i++) { double lo, hi; if (fscanf(f, "%la %la", &lo, &hi) != 2) return 1; b[i] = Interval(lo, hi); }
	fclose(f);
	std::string lin = argv[3], ctc = argv[4];
	Linearizer* l = lin == "art" ? (Linearizer*) new LinearizerAffine2(sys) : (Linearizer*) new LinearizerXTaylor(sys);
	CtcSC* sc = ctc == "sc" ? new CtcSC(*l, CtcSC::DUAL, 1e-9) : NULL;
	Ctc* c = sc ? (Ctc*) sc : (Ctc*) new CtcPolytopeHull(*l);
	IntervalVector b0 = b;
	try { c->contract(b); } catch (...) { b.set_empty(); }
	int ch = 0;
	if (!b.is_empty()) for (int i = 0; i < n; i++) if (!(b[i] == b0[i])) ch++;
	printf("%s %s: %s (%d variables contraidas)", argv[3], argv[4], b.is_empty() ? "VACIA" : "no vacia", ch);
	if (sc) printf(", LPs=%ld, certificados con NaN=%ld, filas descartadas=%ld", sc->nb_lps(), sc->nb_nan_certificates(), sc->nb_dropped_rows());
	printf("\n");
	return 0;
}
