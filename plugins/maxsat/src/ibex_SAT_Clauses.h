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
            int smaller_than(int lit1, int lit2);
            my_type redundant(int *new_clause, int *old_clause);
            void remove_passive_clauses();
            void remove_passive_vars_in_clause(int clause);
            int clean_structure();
            void lire_clauses(FILE *fp_in);
            void build_structure();
            void eliminate_redundance();
            my_type build_simple_sat_instance(char *input_file);

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