#include "ibex.h"
#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <set>

// #define WORD_LENGTH 100
// #define TRUE 1
// #define FALSE 0
// #define NONE -1

// #define WEIGHT 4
// #define WEIGHT1 25
// #define WEIGHT2 5
// #define WEIGHT3 1
// #define T 10

// #define NO_CONFLICT -3
// #define NO_REASON -3
static const int WORD_LENGTH = 100;
static const int TRUE = 1;
static const int FALSE = 0;
static const int NONE = -1;

static const int WEIGHT = 4;
static const int WEIGHT1 = 25;
static const int WEIGHT2 = 5;
static const int WEIGHT3 = 1;
static const int T = 10;

static const int NO_CONFLICT = -3;
static const int NO_REASON = -3;

typedef signed char my_type;
typedef unsigned char my_unsigned_type;

#include "ibex_SAT_Solver.h"
#include "ibex_SAT_Clauses.h"

using namespace std;
using namespace ibex;

int main(int argc, char *argv[]) {
    try{
        cout << "maxsat" << endl;

        SAT_Clauses *cl = new SAT_Clauses();
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

        switch (cl->build_simple_sat_instance(argv[1])){
            case FALSE:
                s->INIT_NB_CLAUSE = s->NB_CLAUSE;
                printf("Input file error\n");
                return FALSE;
            case TRUE:
                s->INIT_NB_CLAUSE = s->NB_CLAUSE;
                if (argc>2)
                    s->UB=atoi(argv[2]);
                else
                    s->UB=s->NB_CLAUSE;
                s->init();
                s->dpl();
                break;
            case NONE:
                s->INIT_NB_CLAUSE = s->NB_CLAUSE;
                printf("An empty resolvant is found!\n"); break;
        }
        mess=times(a_tms); endtime = a_tms->tms_utime;

        printf("s OPTIMUM FOUND\nc Optimal Solution (minimum number of unsatisfied clauses) = %d\n", s->UB);
        printf("v");
        for (i = 0; i < s->NB_VAR; i++) {
            if (s->var_best_value[i] == FALSE)
                printf(" -%i", i + 1);
            else
                printf(" %i", i + 1);
        }
        printf(" 0\n");
        printf("NB_MONO= %ld, NB_UNIT= %ld, NB_BRANCHE= %ld, NB_BACK= %ld \n", s->NB_MONO, s->NB_UNIT, s->NB_BRANCHE, s->NB_BACK);

        printf ("Program terminated in %5.3f seconds.\n", ((double)(endtime-begintime)/CLOCKS_PER_SEC));

        fp_time = fopen("timetable", "a");
        fprintf(fp_time, "maxsatz14bis+fl %s %5.3f %ld %ld %d %d %d %d\n",
            saved_input_file, ((double)(endtime-begintime)/CLOCKS_PER_SEC),
            s->NB_BRANCHE, s->NB_BACK,
            s->UB, s->NB_VAR, s->INIT_NB_CLAUSE, s->NB_CLAUSE-s->INIT_NB_CLAUSE);
        printf("maxsatz14bis+fl %s %5.3f %ld %ld %d %d %d %d\n",
            saved_input_file, ((double)(endtime-begintime)/CLOCKS_PER_SEC),
            s->NB_BRANCHE, s->NB_BACK,
            s->UB, s->NB_VAR, s->INIT_NB_CLAUSE, s->NB_CLAUSE-s->INIT_NB_CLAUSE);
        fclose(fp_time);
        return TRUE;
    }catch(ibex::SyntaxError& e){
		cout << e << endl;
	}
}