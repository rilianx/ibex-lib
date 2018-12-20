//============================================================================
//                                  I B E X
// File        : ibex_SAT_Clauses.h
// Author      : Cristian Duartes
// Copyright   : See README
// Created     : Oct 30, 2018
//============================================================================

#ifndef __IBEX_SAT_CLAUSES_H__
#define __IBEX_SAT_CLAUSES_H__

typedef signed char my_type;
typedef unsigned char my_unsigned_type;

using namespace std;
namespace ibex {
    class SAT_Clauses{
        public:
            SAT_Clauses(){
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
            int smaller_than(int lit1, int lit2);
            my_type redundant(int *new_clause, int *old_clause);
            void remove_passive_clauses();
            void remove_passive_vars_in_clause(int clause);
            int clean_structure();
            void lire_clauses(FILE *fp_in);
            void build_structure();
            void eliminate_redundance();
            my_type build_simple_sat_instance(char *input_file, int *NB_VAR, int *NB_CLAUSE, int *INIT_NB_CLAUSE);

    };
}

#endif // __IBEX_SAT_CLAUSES_H__