#include "ibex.h"
#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <set>

#include "SAT_Clauses.h"
#include "SAT_Solver.h"

using namespace std;
using namespace ibex;
//main(int argc, char *argv[]) {
int main(int argc, char *argv[]) {
    try{
        cout << "maxsat" << endl;

        //SAT_Clauses *cl = new SAT_Clauses();
        SAT_Solver *s = new SAT_Solver();

        char saved_input_file[WORD_LENGTH];
        int i,  var; 
        long begintime, endtime, mess;
        struct tms *a_tms;
        FILE *fp_time;

        if (argc<2) {
            printf("Using format: maxsatz input_instance [upper_bound]\n");
            return FALSE;
        }
        //for (i=0; i<WORD_LENGTH; i++) saved_input_file[i]=argv[1][i];
        //for (i=0; i<WORD_LENGTH; i++) cout << saved_input_file[i] << endl;

        a_tms = ( struct tms *) malloc( sizeof (struct tms));
        mess=times(a_tms); begintime = a_tms->tms_utime;

        switch (s.build_simple_sat_instance(argv[1])) {
            case FALSE:
                printf("Input file error\n");
                return FALSE;
            case TRUE:
                if (argc>2)
                    UB=atoi(argv[2]);
                else
                    UB=NB_CLAUSE;
                s->init();
                s->dpl();
                break;
            case NONE:
                printf("An empty resolvant is found!\n"); break;
        }
        mess=times(a_tms); endtime = a_tms->tms_utime;

        printf("s OPTIMUM FOUND\nc Optimal Solution (minimum number of unsatisfied clauses) = %d\n", UB);
        printf("v");
        for (i = 0; i < NB_VAR; i++) {
            if (s->var_best_value[i] == FALSE)
                printf(" -%i", i + 1);
            else
                printf(" %i", i + 1);
        }
        printf(" 0\n");
        printf("NB_MONO= %ld, NB_UNIT= %ld, NB_BRANCHE= %ld, NB_BACK= %ld \n", NB_MONO, NB_UNIT, NB_BRANCHE, NB_BACK);
                    
        printf ("Program terminated in %5.3f seconds.\n", ((double)(endtime-begintime)/CLOCKS_PER_SEC));

        fp_time = fopen("timetable", "a");
        fprintf(fp_time, "maxsatz14bis+fl %s %5.3f %ld %ld %d %d %d %d\n", 
            saved_input_file, ((double)(endtime-begintime)/CLOCKS_PER_SEC), 
            NB_BRANCHE, NB_BACK,  
            UB, NB_VAR, INIT_NB_CLAUSE, NB_CLAUSE-INIT_NB_CLAUSE);
        printf("maxsatz14bis+fl %s %5.3f %ld %ld %d %d %d %d\n", 
            saved_input_file, ((double)(endtime-begintime)/CLOCKS_PER_SEC), 
            NB_BRANCHE, NB_BACK,
            UB, NB_VAR, INIT_NB_CLAUSE, NB_CLAUSE-INIT_NB_CLAUSE);
        fclose(fp_time);
        return TRUE;
    }catch(ibex::SyntaxError& e){
		cout << e << endl;
	}
}