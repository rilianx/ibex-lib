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
            int smaller_than(int lit1, int lit2, int *NB_VAR);
            my_type redundant(int *new_clause, int *old_clause, int *NB_VAR);
            void remove_passive_clauses(int *NB_CLAUSE, my_type *clause_state, int **sat, int *var_sign, my_type *clause_length);
            void remove_passive_vars_in_clause(int clause, int **var_sign, my_type *var_state);
            int clean_structure(int *NB_CLAUSE, int *NB_VAR, my_type *clause_state, int **sat, int **var_sign, my_type *clause_length, my_type *var_state, int **pos_in, int *pos_nb, int **neg_in, int *neg_nb);
            void lire_clauses(FILE *fp_in, int *NB_CLAUSE, int *NB_VAR, int **sat, my_type *clause_length, my_type *clause_state);
            void build_structure(int *NB_VAR, int *NB_CLAUSE, int *neg_nb, int **neg_in, int *pos_nb, int **pos_in, my_type *clause_length, int **sat, int **var_sign, my_type *var_state);
            void eliminate_redundance(int *NB_CLAUSE, my_type *clause_state, my_type *clause_length, int *UNITCLAUSE_STACK, int *UNITCLAUSE_STACK_fill_pointer);

            my_type build_simple_sat_instance(char *input_file, int *NB_VAR, int *NB_CLAUSE, int *INIT_NB_CLAUSE, my_type *clause_state, int **sat, int **var_sign, my_type *clause_length, my_type *var_state, int **pos_in, int *pos_nb, int **neg_in, int *neg_nb, int *UNITCLAUSE_STACK, int *UNITCLAUSE_STACK_fill_pointer);

            void push_clause(int item, int *UNITCLAUSE_STACK, int *UNITCLAUSE_STACK_fill_pointer){
                (*UNITCLAUSE_STACK_fill_pointer)++;
                UNITCLAUSE_STACK[(*UNITCLAUSE_STACK_fill_pointer)] = item;
            }

            // #define complement(lit1, lit2) \ ((lit1<lit2) ? lit2-lit1 == NB_VAR : lit1-lit2 == NB_VAR)
            int complement(int lit1, int lit2, int NB_VAR){ // TODO: revisar
                if (lit2-lit1 == NB_VAR)
                    return 1;
                else if (lit1-lit2 == NB_VAR)
                    return 1;
                else
                    return 0;
            }

            // #define negative(literal) literal>=NB_VAR

            inline int negative_clause(int literal, int NB_VAR) { return literal >= NB_VAR; }

            // #define positive(literal) literal<NB_VAR

            inline int positive_clause(int literal, int NB_VAR) { return literal < NB_VAR; }

            // #define get_var_from_lit(literal) \ ((literal<NB_VAR) ? literal : literal-NB_VAR)

            int get_var_from_lit_clause(int literal, int NB_VAR){
                if (literal < NB_VAR)
                    return literal;
                else
                    return literal - NB_VAR;
            }

    };
}

#endif // __IBEX_SAT_CLAUSES_H__