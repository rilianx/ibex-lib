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
                // NUMBER BY REFERENCES
                int *NB_CLAUSE, *NB_EMPTY, *UB, *REAL_NB_CLAUSE;
                long *NB_BACK;
                // fill_pointers
                int *SAVED_CLAUSES_fill_pointer, *NEW_CLAUSES_fill_pointer, *CLAUSE_STACK_fill_pointer, *UNITCLAUSE_STACK_fill_pointer,
                        *REDUCEDCLAUSE_STACK_fill_pointer, *REDUCEDCLAUSE_STACK_fill_pointer, *MY_UNITCLAUSE_STACK_fill_pointer, *UNITCLAUSE_STACK_fill_pointer,
                        *VARIABLE_STACK_fill_pointer, *CLAUSES_TO_REMOVE_fill_pointer, *LINEAR_REASON_STACK1_fill_pointer, *LINEAR_REASON_STACK2_fill_pointer;

                // MAIN STACKS
                int *SAVED_CLAUSES, *CLAUSE_STACK, *REDUCEDCLAUSE_STACK, *MY_UNITCLAUSE_STACK, *UNITCLAUSE_STACK, *VARIABLE_STACK, *CLAUSES_TO_REMOVE, *LINEAR_REASON_STACK1,
                        *LINEAR_REASON_STACK2;
                // OTHER STACKS
                int *saved_reducedclause_stack, *saved_new_clauses, *saved_unitclause_stack, *saved_clause_stack, *saved_nb_empty, *saved_nb_clause, *saved_saved_clauses,
                        *clause_involved, *reason;;
                // DUAL
                int **var_sign, **pos_in, **neg_in, **NEW_CLAUSES, **SAVED_CLAUSE_POSITIONS;

                // MY TYPE
                my_type *clause_state, *clause_length, *var_current_value, *var_state, *var_rest_value;
                // CellMOP() : depth(0), id(0), a(0.0), w_lb(POS_INFINITY), ub_distance(POS_INFINITY) {}
            SAT_Solver(): NB_CLAUSE(NB_CLAUSE), NB_EMPTY(NB_EMPTY), REAL_NB_CLAUSE(REAL_NB_CLAUSE), NB_BACK(NB_BACK) {}

            void remove_clauses(int var, int **pos_in, int **neg_in, my_type *var_current_value, my_type *clause_state, int *CLAUSE_STACK, int *CLAUSE_STACK_fill_pointer);
            int reduce_clauses(int var, int *NB_EMPTY, int *UB, int **neg_in, int **pos_in, int *REDUCEDCLAUSE_STACK, int *UNITCLAUSE_STACK,
                    int *UNITCLAUSE_STACK_fill_pointer, int *REDUCEDCLAUSE_STACK_fill_pointer, my_type *var_current_value, my_type *clause_state, my_type *clause_length);
            int my_reduce_clauses(int var, int **neg_in, int **pos_in, int *REDUCEDCLAUSE_STACK, int *MY_UNITCLAUSE_STACK, int *REDUCEDCLAUSE_STACK_fill_pointer,
                    int *MY_UNITCLAUSE_STACK_fill_pointer, my_type *var_current_value ,my_type *clause_state, my_type *clause_length);
            int my_reduce_clauses_for_fl(int var, int **neg_in, int **pos_in, int *REDUCEDCLAUSE_STACK, int *UNITCLAUSE_STACK, int *REDUCEDCLAUSE_STACK_fill_pointer,
                    int *UNITCLAUSE_STACK_fill_pointer, my_type *var_current_value, my_type *clause_state, my_type *clause_length);
            void print_values(int nb_var, my_type *var_current_value);
            int backtracking(int *saved_reducedclause_stack, int *saved_new_clauses, int *NB_EMPTY, int *NB_CLAUSE, int *CLAUSE_STACK, int *CLAUSE_STACK_fill_pointer,
                    int *REDUCEDCLAUSE_STACK_fill_pointer, int *REDUCEDCLAUSE_STACK, int *UNITCLAUSE_STACK, int *UNITCLAUSE_STACK_fill_pointer, int *saved_unitclause_stack,
                    int *VARIABLE_STACK, int *NEW_CLAUSES_fill_pointer, int *saved_clause_stack, int *saved_nb_empty, int *saved_nb_clause, int *saved_saved_clauses,
                    int **SAVED_CLAUSE_POSITIONS, int *SAVED_CLAUSES, int *SAVED_CLAUSES_fill_pointer, int *VARIABLE_STACK_fill_pointer, long *NB_BACK, int **pos_in,
                    int **neg_in, int *UB, my_type *var_state, my_type *var_rest_value, my_type *var_current_value, my_type *clause_state, my_type *clause_length);
            int verify_solution(int *REAL_NB_CLAUSE, int **var_sign, my_type *var_current_value);
            void reset_context(int saved_clause_stack_fill_pointer, int saved_reducedclause_stack_fill_pointer, int saved_unitclause_stack_fill_pointer,
                    int saved_variable_stack_fill_pointer, int *CLAUSE_STACK_fill_pointer, int *REDUCEDCLAUSE_STACK_fill_pointer, int *REDUCEDCLAUSE_STACK,
                    int *VARIABLE_STACK_fill_pointer, int *VARIABLE_STACK, int *CLAUSE_STACK, int *reason, int *UNITCLAUSE_STACK_fill_pointer,
                    my_type *clause_length, my_type *clause_state, my_type *var_state);
            int replace_clause(int newclause, int clause_to_replace, int *clauses, int **SAVED_CLAUSE_POSITIONS, int *SAVED_CLAUSES_fill_pointer, int *SAVED_CLAUSES);
            void create_binaryclause(int var1, int sign1, int var2, int sign2, int clause1, int clause2, int **NEW_CLAUSES, int *NB_CLAUSE, int **var_sign,
                    int **SAVED_CLAUSE_POSITIONS, int *SAVED_CLAUSES_fill_pointer, int *SAVED_CLAUSES,
                    my_type *clause_state, my_type *clause_length, int *NEW_CLAUSES_fill_pointer, int **pos_in, int **neg_in);
            int verify_binary_clauses(int *varssigns, int var1, int sign1, int var2, int sign2);
            int create_clause_from_conflict_clauses(int clause1, int clause2, int clause3, int **var_sign, int **NEW_CLAUSES, int *NB_CLAUSE,
                    int **SAVED_CLAUSE_POSITIONS, int *SAVED_CLAUSES_fill_pointer, int *SAVED_CLAUSES, int *NEW_CLAUSES_fill_pointer, int **pos_in, int **neg_in,
                    int *CLAUSES_TO_REMOVE, int *CLAUSES_TO_REMOVE_fill_pointer, my_type *clause_state, my_type *clause_length, my_type *var_state);
            int search_linear_reason1(int var, int *reason, int **var_sign, int *LINEAR_REASON_STACK1, int *clause_involved);
            int search_linear_reason2(int var, int *reason, int **var_sign, int *LINEAR_REASON_STACK1, int *LINEAR_REASON_STACK2, int *LINEAR_REASON_STACK1_fill_pointer,
                    int *LINEAR_REASON_STACK2_fill_pointer, int *clause_involved);
            int check_reason(int *varssigns, int clause, int clause1, int clause2, int *reason, int **var_sign, my_type *var_current_value);
            int create_complementary_binclause(int clause, int clause1, int clause2, int **var_sign, int *reason, my_type *var_current_value, int **NEW_CLAUSES,
                    int *NB_CLAUSE, int **SAVED_CLAUSE_POSITIONS, int *SAVED_CLAUSES_fill_pointer, int *SAVED_CLAUSES, int *NEW_CLAUSES_fill_pointer, int **pos_in,
                    int **neg_in, my_type *clause_state, my_type *clause_length);
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

    };
}

#endif // __IBEX_SAT_SOLVER_H__