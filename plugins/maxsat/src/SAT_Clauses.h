//============================================================================
//                                  I B E X
// File        : ibex_SAT_Clauses.h
// Author      : Cristian Duartes
// Copyright   : See README
// Created     : Oct 30, 2018
//============================================================================

#ifndef __IBEX_SAT_CLAUSES_H__
#define __IBEX_SAT_CLAUSES_H__

#define OLD_CLAUSE_REDUNDANT -77
#define NEW_CLAUSE_REDUNDANT -7

using namespace std;
namespace ibex {
    class SAT_Clauses{
        public:
            /** 
             * Custom type
             * */
            typedef signed char my_type;
            typedef unsigned char my_unsigned_type;
            /* */
            int smaller_than(int lit1, int lit2);
            my_type redundant(int *new_clause, int *old_clause);
            void remove_passive_clauses();
            void remove_passive_vars_in_clause(int clause);
            int clean_structure();
            void lire_clauses(FILE *fp_in);
            void build_structure();
            void eliminate_redundance();
            my_type build_simple_sat_instance(char *input_file);
    };
}
#endif // __IBEX_SAT_CLAUSES_H__