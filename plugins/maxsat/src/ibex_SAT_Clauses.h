//============================================================================
//                                  I B E X
// File        : ibex_SAT_Clauses.h
// Author      : Cristian Duartes
// Copyright   : See README
// Created     : Oct 30, 2018
//============================================================================

#ifndef __IBEX_SAT_CLAUSES_H__
#define __IBEX_SAT_CLAUSES_H__

#define pop_clause(stack) stack[--stack ## _fill_pointer]
#define push_clause(item, stack) stack[stack ## _fill_pointer++] = item
#define positive(literal) (literal < NB_VAR)
#define negative(literal) (literal >= NB_VAR)

typedef signed char my_type;
typedef unsigned char my_unsigned_type;

using namespace std;
namespace ibex {
    class SAT_Clauses{
        public:
            int smaller_than(int lit1, int lit2);
            my_type redundant(int *new_clause, int *old_clause, int NB_VAR);
            void remove_passive_clauses(int NB_CLAUSE, int *clause_state, int *sat, int *var_sign, int *clause_length);
            void remove_passive_vars_in_clause(int clause, int **var_sign, int *var_state);
            int clean_structure(int NB_CLAUSE, int NB_VAR, int **neg_nb, int **pos_nb, int **var_sign, int **pos_in, int **neg_in, int *clause_state, int *sat, int *clause_length, int *var_state);
            void lire_clauses(FILE *fp_in, int NB_CLAUSE, int NB_VAR, int **sat, int *clause_length, int *clause_state);
            void build_structure(int NB_VAR, int NB_CLAUSE, int **neg_nb, int **pos_nb, int *clause_length, int **sat, int **var_sign, int **neg_in, int **pos_in, int *var_state);
            void eliminate_redundance(int NB_CLAUSE, int *clause_state, int *clause_length, int UNITCLAUSE_STACK);
            my_type build_simple_sat_instance(char *input_file, int NB_CLAUSE, int NB_VAR, int **sat, int **neg_nb, int **pos_nb,
                        int *clause_length, int *clause_state, int **var_sign, int **neg_in, int **pos_in, int *var_state, int UNITCLAUSE_STACK);

            static bool complement(int lit1, int lit2, int NB_VAR){if(lit1<lit2) return lit2-lit1 == NB_VAR; else return lit1-lit2 == NB_VAR;}

            static const int NEGATIVE = 0;
            static const int POSITIVE = 1;
            static const int PASSIVE = 0;
            static const int ACTIVE = 1;

            static const int OLD_CLAUSE_REDUNDANT = -77;
            static const int NEW_CLAUSE_REDUNDANT = -7;

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
    };
}

#endif // __IBEX_SAT_CLAUSES_H__