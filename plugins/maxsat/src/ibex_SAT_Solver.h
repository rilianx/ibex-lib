//============================================================================
//                                  I B E X
// File        : ibex_SAT_Solver.h
// Author      : Cristian Duartes
// Copyright   : See README
// Created     : Oct 30, 2018
//============================================================================

#include <stdio.h>
// #include <stdlib.h>
#include <time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <limits.h>

#include "ibex_SAT_Clauses.h"
#ifndef __IBEX_SAT_SOLVER_H__
#define __IBEX_SAT_SOLVER_H__

using namespace std;
namespace ibex {
    class SAT_Solver{
        public:
            int tab_unitclause_size;
            int my_tab_variable_size;
            int my_tab_clause_size;
            int my_tab_unitclause_size;
            int tab_literal_size;
            double double_tab_clause_size;
            int *UNITCLAUSE_STACK;
            SAT_Solver(){
                if(tab_clause_size/4<2000)
                    tab_unitclause_size = 2000;
                else
                    tab_unitclause_size = tab_clause_size/4;

                if(tab_variable_size/2<1000)
                    my_tab_variable_size = 1000;
                else
                    my_tab_variable_size = tab_variable_size/2;

                if(tab_clause_size/2<2000)
                    my_tab_clause_size = 2000;
                else
                    my_tab_clause_size = tab_clause_size/2;

                if(tab_unitclause_size/2<1000)
                    my_tab_unitclause_size = 1000;
                else   
                    my_tab_unitclause_size = tab_unitclause_size/2;

                tab_literal_size = 2*tab_variable_size;
                double_tab_clause_size = 2*tab_clause_size;
                UNITCLAUSE_STACK = new int[tab_unitclause_size];
            }
            inline bool positive(int literal) { return literal<NB_VAR;}
            inline bool negative(int literal) { return literal>=NB_VAR;}

            inline int get_var_from_lit(int literal) {if(literal<NB_VAR) return literal; else return literal-NB_VAR;}

            // static bool complement(int lit1, int lit2){if(lit1<lit2) return lit2-lit1 == NB_VAR; else return lit1-lit2 == NB_VAR;}

            inline int inverse_sign(int signe){ if(signe == SAT_Clauses::POSITIVE) return SAT_Clauses::NEGATIVE; else return SAT_Clauses::POSITIVE;}

            inline int unsat(int val){if (val==0) return 0; else return 1;} //inline string "UNS" = 0; "SAT" = 1

            inline bool satisfiable() { return CLAUSE_STACK_fill_pointer == NB_CLAUSE; }

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
            int assign_and_unitclause_process( int var, int value, int starting_point);
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
            int init();

            //typedef signed char my_type;
            //typedef unsigned char my_unsigned_type;
            int NB_VAR;

            static const int tab_variable_size = 10000;
            static const int tab_clause_size = 40000;
            
            int NB_EMPTY=0;
            int UB;
            int REAL_NB_CLAUSE;
            int UNITCLAUSE_STACK_fill_pointer = 0;
            int VARIABLE_STACK_fill_pointer = 0;
            int CLAUSE_STACK_fill_pointer = 0;
            int REDUCEDCLAUSE_STACK_fill_pointer = 0;
            my_type var_current_value[tab_variable_size];
            my_type clause_state[tab_clause_size];
            my_type clause_length[tab_clause_size];

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
            int REDUCEDCLAUSE_STACK[tab_clause_size];

            int PREVIOUS_REDUCEDCLAUSE_STACK_fill_pointer = 0;

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

            int CLAUSES_TO_REMOVE[tab_clause_size];
            int CLAUSES_TO_REMOVE_fill_pointer=0;
            
            
    };
}

#endif // __IBEX_SAT_SOLVER_H__