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

            SAT_Solver(){
                // if(tab_clause_size/4<2000)
                //     tab_unitclause_size = 2000;
                // else
                //     tab_unitclause_size = tab_clause_size/4;

                // if(tab_variable_size/2<1000)
                //     my_tab_variable_size = 1000;
                // else
                //     my_tab_variable_size = tab_variable_size/2;

                // if(tab_clause_size/2<2000)
                //     my_tab_clause_size = 2000;
                // else
                //     my_tab_clause_size = tab_clause_size/2;

                // if(tab_unitclause_size/2<1000)
                //     my_tab_unitclause_size = 1000;
                // else   
                //     my_tab_unitclause_size = tab_unitclause_size/2;

                // tab_literal_size = 2*tab_variable_size;
                // double_tab_clause_size = 2*tab_clause_size;
                // UNITCLAUSE_STACK = new int[tab_unitclause_size];
            }

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

    };
}

#endif // __IBEX_SAT_SOLVER_H__