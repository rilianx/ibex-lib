//============================================================================
//                                  I B E X
// File        : ibex_SAT_Solver.h
// Author      : Cristian Duartes
// Copyright   : See README
// Created     : Oct 30, 2018
//============================================================================

#ifndef __IBEX_SAT_SOLVER_H__
#define __IBEX_SAT_SOLVER_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <limits.h>

#include "SAT_Clauses.h"
using namespace std;
namespace ibex {
    class SAT_Solver{
        public:
            void remove_clauses(int var);
            int reduce_clauses(int var);
            int my_reduce_clauses(int var);
            int my_reduce_clauses_for_fl(int var);
            void print_values(int nb_var);
            int backtracking();
            int verify_solution();
            void reset_context(int saved_clause_stack_fill_pointer, int saved_reducedclause_stack_fill_pointer, int saved_unitclause_stack_fill_pointer, int saved_variable_stack_fill_pointer);
            int replace_clause(int newclause, int clause_to_replace, int *clauses);
            void create_binaryclause(int var1, int sign1, int var2, int sign2, int clause1, int clause2);
            int verify_binary_clauses(int *varssigns, int var1, int sign1, int var2, int sign2);
            int create_clause_from_conflict_clauses(int clause1, int clause2, int clause3);
            int search_linear_reason1(int var);
            int search_linear_reason2(int var);
            int check_reason(int *varssigns, int clause, int clause1, int clause2);
            int create_complementary_binclause(int clause, int clause1, int clause2);
            int get_satisfied_literal(int clause);
            void create_ternary_clauses(int var1, int sign1, int var2, int sign2, int var3, int sign3, int clause1, int clause2, int clause3);
            int non_linear_conflict(int empty_clause, int var1, int sign1, int var2, int sign2);
            int linear_conflict(int clause);
            void remove_linear_reasons();
            int there_is_unit_clause( int var_to_check );
            int assign_and_unitclause_process( int var, int value, int starting_point );
            int store_reason_clauses( int clause, int starting );
            void remove_reason_clauses();
            int failed_literal( int conflict );
            int lookahead();
            int satisfy_unitclause(int unitclause);
            int my_unitclause_process(int starting_point);
            int get_complement(int lit);
            void create_unitclause(int lit, int subsumedclause, int *clauses);
            int verify_resolvent(int lit, int clause1, int clause2);
            int searching_two_clauses_to_fix_neglit(int clause, int lit);
            int simple_get_neg_clause_nb(int var);
            int simple_get_pos_clause_nb(int var);
            int get_neg_clause_nb(int var);
            int searching_two_clauses_to_fix_poslit(int clause, int lit);
            int get_pos_clause_nb(int var);
            int satisfy_literal(int lit);
            int assign_value(int var, int current_value, int rest_value);
            int unitclause_process();
            int choose_and_instantiate_variable();
            int dpl();
            void init();

            /** 
             * Custom type
             * */
            typedef signed char my_type;
            typedef unsigned char my_unsigned_type;
            static const int tab_variable_size = 10000;
            static const int tab_clause_size = 40000;
            int NB_EMPTY=0, UB;
            int REAL_NB_CLAUSE;
            int UNITCLAUSE_STACK_fill_pointer = 0;
            int VARIABLE_STACK_fill_pointer = 0;
            int CLAUSE_STACK_fill_pointer = 0;
            int REDUCEDCLAUSE_STACK_fill_pointer = 0;
            my_type var_current_value[tab_variable_size];
            my_type clause_state[tab_clause_size];
            my_type clause_length[tab_clause_size];
            // static const int WORD_LENGTH = 100 ;
            // static const int TRUE = 1;
            // static const int FALSE = 0;
            // static const int NONE = -1;

            // static const int WEIGHT = 4;
            // static const int WEIGHT1 = 25;
            // static const int WEIGHT2 = 5;
            // static const int WEIGHT3 = 1;
            // static const int T = 10;

            // static const int NEGATIVE = 0;
            // static const int POSITIVE = 1;
            // static const int PASSIVE = 0;
            // static const int ACTIVE = 1;

            // static const int NO_CONFLICT = -3;
            // static const int NO_REASON = -3;

            int *neg_in[tab_variable_size];
            int *pos_in[tab_variable_size];
            int neg_nb[tab_variable_size];
            int pos_nb[tab_variable_size];

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

            int VARIABLE_STACK[tab_variable_size];
            int CLAUSE_STACK[tab_clause_size];
            int UNITCLAUSE_STACK[tab_unitclause_size];
            int REDUCEDCLAUSE_STACK[tab_clause_size];

            int PREVIOUS_REDUCEDCLAUSE_STACK_fill_pointer = 0;

            int NB_VAR;
            int NB_CLAUSE;
            int INIT_NB_CLAUSE;
            long NB_UNIT=1, NB_MONO=0, NB_BRANCHE=0, NB_BACK = 0;
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

            my_type var_best_value[tab_variable_size]; // Best assignment of variables
                        
            
            
    };
}

#endif // __IBEX_SAT_SOLVER_H__