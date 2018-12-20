#include "ibex.h"
#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <set>

#include "ibex_SAT_Solver.h"
#include "ibex_SAT_Clauses.h"

#define tab_variable_size  10000
#define tab_clause_size 40000
#define tab_unitclause_size \
 ((tab_clause_size/4<2000) ? 2000 : tab_clause_size/4)
#define my_tab_variable_size \
 ((tab_variable_size/2<1000) ? 1000 : tab_variable_size/2)
#define my_tab_clause_size \
 ((tab_clause_size/2<2000) ? 2000 : tab_clause_size/2)
#define my_tab_unitclause_size \
 ((tab_unitclause_size/2<1000) ? 1000 : tab_unitclause_size/2)
#define tab_literal_size 2*tab_variable_size
#define double_tab_clause_size 2*tab_clause_size
#define positive(literal) literal<NB_VAR
#define negative(literal) literal>=NB_VAR
#define get_var_from_lit(literal) \
  ((literal<NB_VAR) ? literal : literal-NB_VAR)
#define complement(lit1, lit2) \
 ((lit1<lit2) ? lit2-lit1 == NB_VAR : lit1-lit2 == NB_VAR)

#define inverse_signe(signe) \
 (signe == POSITIVE) ? NEGATIVE : POSITIVE
#define unsat(val) (val==0)?"UNS":"SAT"
#define pop(stack) stack[--stack ## _fill_pointer]
#define push(item, stack) stack[stack ## _fill_pointer++] = item
#define satisfiable() CLAUSE_STACK_fill_pointer == NB_CLAUSE

static const int NO_CONFLICT = -3;
static const int NO_REASON = -3;
static const int NEGATIVE = 0;
static const int POSITIVE = 1;
static const int PASSIVE = 0;
static const int ACTIVE = 1;

static const int WORD_LENGTH = 100;
static const int TRUE = 1;
static const int FALSE = 0;
static const int NONE = -1;

static const int WEIGHT = 4;
static const int WEIGHT1 = 25;
static const int WEIGHT2 = 5;
static const int WEIGHT3 = 1;
static const int T = 10;

typedef signed char my_type;
typedef unsigned char my_unsigned_type;

/* the tables of variables and clauses are statically allocated. Modify the 
   parameters tab_variable_size and tab_clause_size before compilation if 
   necessary */

int *neg_in[tab_variable_size];
int *pos_in[tab_variable_size];
int neg_nb[tab_variable_size];
int pos_nb[tab_variable_size];
my_type var_current_value[tab_variable_size];
my_type var_rest_value[tab_variable_size];
my_type var_state[tab_variable_size];

int saved_clause_stack[tab_variable_size];
int saved_reducedclause_stack[tab_variable_size];
int saved_unitclause_stack[tab_variable_size];
int saved_nb_empty[tab_variable_size];
my_unsigned_type nb_neg_clause_of_length1[tab_variable_size];
my_unsigned_type nb_pos_clause_of_length1[tab_variable_size];
my_unsigned_type nb_neg_clause_of_length2[tab_variable_size];
my_unsigned_type nb_neg_clause_of_length3[tab_variable_size];
my_unsigned_type nb_pos_clause_of_length2[tab_variable_size];
my_unsigned_type nb_pos_clause_of_length3[tab_variable_size];

float reduce_if_negative[tab_variable_size];
float reduce_if_positive[tab_variable_size];

int *sat[tab_clause_size];
int *var_sign[tab_clause_size];
my_type clause_state[tab_clause_size];
my_type clause_length[tab_clause_size];

int VARIABLE_STACK_fill_pointer = 0;
int CLAUSE_STACK_fill_pointer = 0;
int UNITCLAUSE_STACK_fill_pointer = 0;
int REDUCEDCLAUSE_STACK_fill_pointer = 0;

int VARIABLE_STACK[tab_variable_size];
int CLAUSE_STACK[tab_clause_size];
int UNITCLAUSE_STACK[tab_unitclause_size];
int REDUCEDCLAUSE_STACK[tab_clause_size];

int PREVIOUS_REDUCEDCLAUSE_STACK_fill_pointer = 0;

int NB_VAR;
int NB_CLAUSE;
int INIT_NB_CLAUSE;
int REAL_NB_CLAUSE;

long NB_UNIT=1, NB_MONO=0, NB_BRANCHE=0, NB_BACK = 0;
int NB_EMPTY=0, UB;

int reason[tab_variable_size];
int REASON_STACK[tab_variable_size];
int REASON_STACK_fill_pointer=0;

int MY_UNITCLAUSE_STACK[tab_variable_size];
int MY_UNITCLAUSE_STACK_fill_pointer=0;
int CANDIDATE_LITERALS[2*tab_variable_size];
int CANDIDATE_LITERALS_fill_pointer=0;
int NEW_CLAUSES[tab_clause_size][7];
int NEW_CLAUSES_fill_pointer=0;
int lit_to_fix[tab_clause_size];
int *SAVED_CLAUSE_POSITIONS[tab_clause_size];
int SAVED_CLAUSES[tab_clause_size];
int SAVED_CLAUSES_fill_pointer=0;
int lit_involved_in_clause[2*tab_variable_size];
int INVOLVED_LIT_STACK[2*tab_variable_size];
int INVOLVED_LIT_STACK_fill_pointer=0;
int fixing_clause[2*tab_variable_size];
int saved_nb_clause[tab_variable_size];
int saved_saved_clauses[tab_variable_size];
int saved_new_clauses[tab_variable_size];

using namespace std;
using namespace ibex;

main(int argc, char *argv[]) {
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

  SAT_Clauses *sat_clauses = new SAT_Clauses();

  switch (sat_clauses->build_simple_sat_instance(argv[1])) {
    case FALSE: printf("Input file error\n"); return FALSE;
    case TRUE:
      if (argc>2)
        UB=atoi(argv[2]);
      else
        UB=NB_CLAUSE;
      init();
      dpl();
      break;
    case NONE: printf("An empty resolvant is found!\n"); break;
  }
  mess=times(a_tms); endtime = a_tms->tms_utime;

  printf("s OPTIMUM FOUND\nc Optimal Solution (minimum number of unsatisfied clauses) = %d\n", UB);
  printf("v");
  for (i = 0; i < NB_VAR; i++) {
    if (var_best_value[i] == FALSE)
      printf(" -%i", i + 1);
    else
      printf(" %i", i + 1);
  }
  printf(" 0\n");
  printf("NB_MONO= %ld, NB_UNIT= %ld, NB_BRANCHE= %ld, NB_BACK= %ld \n", 
	 NB_MONO, NB_UNIT, NB_BRANCHE, NB_BACK);
	        
  printf ("Program terminated in %5.3f seconds.\n",
	  ((double)(endtime-begintime)/CLOCKS_PER_SEC));

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
}