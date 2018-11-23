//============================================================================
//                                  I B E X
// File        : ibex_SAT_Solver.cpp
// Author      : Cristian Duartes
// Copyright   : See README
// Created     : Oct 30, 2018
//============================================================================

#include "ibex_SAT_Solver.h"
using namespace std;
namespace ibex {

    void SAT_Solver::remove_clauses(int var) {
        register int clause;
        register int *clauses;
        if (var_current_value[var] == SAT_Clauses::POSITIVE)
            clauses = pos_in[var];
        else
            clauses = neg_in[var];
        for(clause=*clauses; clause!=SAT_Clauses::NONE; clause=*(++clauses)) {
            if (clause_state[clause] == SAT_Clauses::ACTIVE) {
                clause_state[clause] = SAT_Clauses::PASSIVE;
                push_clause(clause, CLAUSE_STACK);
            }
        }
    }

    int SAT_Solver::reduce_clauses(int var) {
        register int clause;
        register int *clauses;
        if (var_current_value[var] == SAT_Clauses::POSITIVE)
            clauses = neg_in[var];
        else
            clauses = pos_in[var];
        for(clause=*clauses; clause!=SAT_Clauses::NONE; clause=*(++clauses)) {
            if (clause_state[clause] == SAT_Clauses::ACTIVE) {
                clause_length[clause]--;
                push_clause(clause, REDUCEDCLAUSE_STACK);
                switch (clause_length[clause]) {
                    case 0:
                        NB_EMPTY++;
                        if (UB<=NB_EMPTY) return SAT_Clauses::NONE;
                        break;
                    case 1: 
                        push_clause(clause, UNITCLAUSE_STACK);
                        break;
                }
            }
        }
        return SAT_Clauses::TRUE;
    }

    int SAT_Solver::my_reduce_clauses(int var) {
        register int clause;
        register int *clauses;
        if (var_current_value[var] == SAT_Clauses::POSITIVE)
            clauses = neg_in[var];
        else
            clauses = pos_in[var];
        for(clause=*clauses; clause!=SAT_Clauses::NONE; clause=*(++clauses)) {
            if (clause_state[clause] == SAT_Clauses::ACTIVE) {
                clause_length[clause]--;
                push_clause(clause, REDUCEDCLAUSE_STACK);
                    switch (clause_length[clause]) {
                        case 0:
                            return clause;
                        case 1: 
                            push_clause(clause, MY_UNITCLAUSE_STACK);
                    break;
                }
            }
        }
        return SAT_Clauses::NO_CONFLICT;
    }

    int SAT_Solver::my_reduce_clauses_for_fl(int var) {
        register int clause;
        register int *clauses;
        if (var_current_value[var] == SAT_Clauses::POSITIVE)
            clauses = neg_in[var];
        else
            clauses = pos_in[var];
        for(clause=*clauses; clause!=SAT_Clauses::NONE; clause=*(++clauses)) {
            if (clause_state[clause] == SAT_Clauses::ACTIVE) {
                clause_length[clause]--;
                push_clause(clause, REDUCEDCLAUSE_STACK);
                switch (clause_length[clause]) {
                    case 0:
                        return clause;
                    case 1:
                        push_clause(clause, UNITCLAUSE_STACK);
                        break;
                }
            }
        }
        return SAT_Clauses::NO_CONFLICT;
    }

    void SAT_Solver::print_values(int nb_var) {
        FILE* fp_out;
        int i;
        fp_out = fopen("satx.sol", "w");
        for (i=0; i<nb_var; i++) {
            if (var_current_value[i] == 1) 
                fprintf(fp_out, "%d ", i+1);
            else
                fprintf(fp_out, "%d ", 0-i-1);
        }
        fprintf(fp_out, "\n");
        fclose(fp_out);			
    } 

    int SAT_Solver::backtracking() {
        int var, index,clause, *position, saved;
            
        NB_BACK++;

        do {
            var = pop_clause(VARIABLE_STACK);
            if (var_rest_value[var] == SAT_Clauses::NONE){
                var_state[var] = SAT_Clauses::ACTIVE;
            }else{
                for (index = saved_clause_stack[var]; 
                        index < CLAUSE_STACK_fill_pointer;
                        index++)
                    clause_state[CLAUSE_STACK[index]] = SAT_Clauses::ACTIVE;
                CLAUSE_STACK_fill_pointer = saved_clause_stack[var];

                for (index = saved_reducedclause_stack[var];
                        index < REDUCEDCLAUSE_STACK_fill_pointer;
                        index++) {	
                    clause = REDUCEDCLAUSE_STACK[index];
                    clause_length[REDUCEDCLAUSE_STACK[index]]++;
                }
                REDUCEDCLAUSE_STACK_fill_pointer = saved_reducedclause_stack[var];
                UNITCLAUSE_STACK_fill_pointer=saved_unitclause_stack[var];
                NB_EMPTY=saved_nb_empty[var];
                NB_CLAUSE=saved_nb_clause[var];
                NEW_CLAUSES_fill_pointer=saved_new_clauses[var];
                
                saved=saved_saved_clauses[var];
                for (index = SAVED_CLAUSES_fill_pointer-1 ;
                        index >= saved;
                        index--)
                    *SAVED_CLAUSE_POSITIONS[index]=SAVED_CLAUSES[index];
                SAVED_CLAUSES_fill_pointer=saved;

                if (NB_EMPTY<UB) {
                    var_current_value[var] = var_rest_value[var];
                    var_rest_value[var] = SAT_Clauses::NONE;
                    push_clause(var, VARIABLE_STACK);
                    if (reduce_clauses(var)==SAT_Clauses::NONE)
                        return SAT_Clauses::NONE;
                    remove_clauses(var);
                    return SAT_Clauses::TRUE;
                }else
                    var_state[var] = SAT_Clauses::ACTIVE;
            }
        } while (VARIABLE_STACK_fill_pointer > 0);
        return SAT_Clauses::FALSE;
    }

    int SAT_Solver::verify_solution() {
        int i, nb=0, var, *vars_signs, clause_truth,cpt;

        for (i=0; i<REAL_NB_CLAUSE; i++) {
            clause_truth = SAT_Clauses::FALSE;
            vars_signs = var_sign[i];
            for(var=*vars_signs; var!=SAT_Clauses::NONE; var=*(vars_signs+=2))
                if (*(vars_signs+1) == var_current_value[var] ) {
                    clause_truth = SAT_Clauses::TRUE;
                    break;
            }
            if (clause_truth == SAT_Clauses::FALSE) nb++;
        }
        return nb;
    }

    void SAT_Solver::reset_context(int saved_clause_stack_fill_pointer, 
            int saved_reducedclause_stack_fill_pointer,
            int saved_unitclause_stack_fill_pointer,
            int saved_variable_stack_fill_pointer) {
        int index, var, clause;
        for (index = saved_clause_stack_fill_pointer; 
                index < CLAUSE_STACK_fill_pointer;
                index++)
            clause_state[CLAUSE_STACK[index]] = SAT_Clauses::ACTIVE;
        CLAUSE_STACK_fill_pointer = saved_clause_stack_fill_pointer;

        for (index = saved_reducedclause_stack_fill_pointer;
                index < REDUCEDCLAUSE_STACK_fill_pointer;
                index++) {	
            clause = REDUCEDCLAUSE_STACK[index];
            clause_length[REDUCEDCLAUSE_STACK[index]]++;
        }
        REDUCEDCLAUSE_STACK_fill_pointer = saved_reducedclause_stack_fill_pointer;

        for(index=saved_variable_stack_fill_pointer;
                index<VARIABLE_STACK_fill_pointer;
                index++) {
            var=VARIABLE_STACK[index];
            reason[var]=SAT_Clauses::NO_REASON;
            var_state[var]=SAT_Clauses::ACTIVE;
        }
        VARIABLE_STACK_fill_pointer=saved_variable_stack_fill_pointer;

        UNITCLAUSE_STACK_fill_pointer=saved_unitclause_stack_fill_pointer;
    }

    int SAT_Solver::replace_clause(int newclause, int clause_to_replace, int *clauses) {
        int clause, flag=SAT_Clauses::FALSE;
        
        for(clause=*clauses; clause!=SAT_Clauses::NONE; clause=*(++clauses)) {
            if (clause==clause_to_replace) {
                *clauses=newclause;
                SAVED_CLAUSE_POSITIONS[SAVED_CLAUSES_fill_pointer]=clauses;
                push_clause(clause_to_replace, SAVED_CLAUSES);
                flag=SAT_Clauses::TRUE;
                break;
            }
        }
        if (flag==SAT_Clauses::FALSE)
            printf("problem...");
        return flag;
    }

    void SAT_Solver::create_binaryclause(int var1, int sign1, int var2, int sign2, 
            int clause1, int clause2) {
        int clause, *vars_signs, flag=SAT_Clauses::FALSE, *clauses1, *clauses2;
        if (sign1==SAT_Clauses::POSITIVE)
            clauses1=pos_in[var1];
        else
            clauses1=neg_in[var1];
        if (sign2==SAT_Clauses::POSITIVE)
            clauses2=pos_in[var2];
        else
            clauses2=neg_in[var2];
        vars_signs=NEW_CLAUSES[NEW_CLAUSES_fill_pointer++];
        if (var1<var2) {
            vars_signs[0]=var1; vars_signs[1]=sign1;
            vars_signs[2]=var2; vars_signs[3]=sign2;
        }else{
            vars_signs[0]=var2; vars_signs[1]=sign2;
            vars_signs[2]=var1; vars_signs[3]=sign1;
        }
        vars_signs[4]=SAT_Clauses::NONE;
        var_sign[NB_CLAUSE]=vars_signs;
        clause_state[NB_CLAUSE]=SAT_Clauses::ACTIVE;
        clause_length[NB_CLAUSE]=2;
        // if (NB_CLAUSE==305)
        // printf("aaa...");
        replace_clause(NB_CLAUSE, clause1, clauses1);
        replace_clause(NB_CLAUSE, clause2, clauses2);
        NB_CLAUSE++;
    }

    int SAT_Solver::verify_binary_clauses(int *varssigns, int var1, int sign1, int var2, int sign2) {
        int nb=0;

        if (var1==*varssigns) {
            if ((*(varssigns+1)!=1-sign1) || (var2!=*(varssigns+2)) || (*(varssigns+3)!=1-sign2)) {
                printf("problem..");
                return SAT_Clauses::FALSE;
            }
        }else{
            if ((var2 != *varssigns) || (*(varssigns+1)!=1-sign2) || (var1!=*(varssigns+2)) || (*(varssigns+3)!=1-sign1)) {
                printf("problem..");
                return SAT_Clauses::FALSE;
            }
        }
        return SAT_Clauses::TRUE;
    }

    int SAT_Solver::create_clause_from_conflict_clauses(int clause1, int clause2, int clause3) {
        int var3, sign3, var2, sign2,*clauses2, *clauses3, *vars_signs, 
            varssigns[4], i=0, var;

        if ((clause_state[clause1]==SAT_Clauses::ACTIVE) && (clause_length[clause1]==2) &&
                (clause_state[clause2]==SAT_Clauses::ACTIVE) && (clause_length[clause2]==1) &&
                (clause_state[clause3]==SAT_Clauses::ACTIVE) && (clause_length[clause3]==1)) {
            vars_signs = var_sign[clause1];
            for(var=*vars_signs; var!=SAT_Clauses::NONE; var=*(vars_signs+=2)) {
                if (var_state[var]==SAT_Clauses::ACTIVE) {
                    varssigns[i++]=var;
                    varssigns[i++]=*(vars_signs+1);
                }
            }
            /*
            vars_signs = var_sign[clause2];
            for(var=*vars_signs; var!=SAT_Clauses::NONE; var=*(vars_signs+=2)) {
            if (var_state[var]==SAT_Clauses::ACTIVE) {
            var2=var; sign2=*(vars_signs+1);
            }
            }
            vars_signs = var_sign[clause3];
            for(var=*vars_signs; var!=SAT_Clauses::NONE; var=*(vars_signs+=2)) {
            if (var_state[var]==SAT_Clauses::ACTIVE) {
            var3=var; sign3=*(vars_signs+1);
            }
            }
            verify_binary_clauses(varssigns, var2, sign2, var3, sign3);
            */
            var2=varssigns[0]; sign2=1-varssigns[1];
            var3=varssigns[2]; sign3=1-varssigns[3];
            create_binaryclause(var2, sign2, var3, sign3, clause2, clause3);
            push_clause(clause1, CLAUSES_TO_REMOVE);
            push_clause(clause2, CLAUSES_TO_REMOVE);
            push_clause(clause3, CLAUSES_TO_REMOVE);
            return SAT_Clauses::TRUE;
        }
        else {
            return SAT_Clauses::FALSE;
        }
    }

    int LINEAR_REASON_STACK1[SAT_Solver::tab_clause_size];
    int LINEAR_REASON_STACK1_fill_pointer=0;
    int LINEAR_REASON_STACK2[SAT_Solver::tab_clause_size];
    int LINEAR_REASON_STACK2_fill_pointer=0;
    int clause_involved[SAT_Solver::tab_clause_size];

    int SAT_Solver::search_linear_reason1(int var) {
        int *vars_signs, clause, fixed_var, index_var, new_fixed_var;

        for(fixed_var=var; fixed_var!=SAT_Clauses::NONE; fixed_var=new_fixed_var) {
            clause=reason[fixed_var];
            vars_signs = var_sign[clause]; new_fixed_var=SAT_Clauses::NONE;
            push_clause(clause, LINEAR_REASON_STACK1);
            clause_involved[clause]=SAT_Clauses::TRUE;
            for(index_var=*vars_signs; index_var!=SAT_Clauses::NONE; index_var=*(vars_signs+=2)) {
                if ((index_var!=fixed_var) && (reason[index_var]!=SAT_Clauses::NO_REASON)) {
                    if (new_fixed_var==SAT_Clauses::NONE)
                        new_fixed_var=index_var;
                    else
                        return SAT_Clauses::FALSE;
                }
            }
        }
        return SAT_Clauses::TRUE;
    }

    #define SIMPLE_NON_LINEAR_CASE 2

    int SAT_Solver::search_linear_reason2(int var) {
        int *vars_signs, clause, fixed_var, index_var, new_fixed_var;

        for(fixed_var=var; fixed_var!=SAT_Clauses::NONE; fixed_var=new_fixed_var) {
            clause=reason[fixed_var];
            if (clause_involved[clause]==SAT_Clauses::TRUE) {
                if ( LINEAR_REASON_STACK2_fill_pointer == 2 && LINEAR_REASON_STACK1_fill_pointer > 2 && LINEAR_REASON_STACK1[ 2 ] == clause ) 
                    return SIMPLE_NON_LINEAR_CASE;
                else
                    return SAT_Clauses::FALSE;
                }else 
                push_clause(clause, LINEAR_REASON_STACK2);
            vars_signs = var_sign[clause]; new_fixed_var=SAT_Clauses::NONE;
            for(index_var=*vars_signs; index_var!=SAT_Clauses::NONE; index_var=*(vars_signs+=2)) {
                if ((index_var!=fixed_var) && (reason[index_var]!=SAT_Clauses::NO_REASON)) {
                    if (new_fixed_var==SAT_Clauses::NONE)
                        new_fixed_var=index_var;
                    else
                        return SAT_Clauses::FALSE;
                }
            }
        }
        return SAT_Clauses::TRUE;
    }

    // clause1 is l1->l2, clause is l2->l3, clause3 is ((not l3) or (not l4))
    // i.e., the reason of l2 is clause1, the reason of l3 is clause
    int SAT_Solver::check_reason(int *varssigns, int clause, int clause1, int clause2) {
        int var, *vars_signs, var1, var2, flag;

        if ((reason[varssigns[0]]!=clause1) || (reason[varssigns[2]]!=clause)) 
            return SAT_Clauses::FALSE;
        vars_signs = var_sign[clause2]; flag=SAT_Clauses::FALSE;
        for(var=*vars_signs; var!=SAT_Clauses::NONE; var=*(vars_signs+=2)) {
            if ((varssigns[2]==var) && (reason[var]!=SAT_Clauses::NO_REASON) && (*(vars_signs+1) != var_current_value[var])) {
                flag=SAT_Clauses::TRUE;
            }
        }
        return flag;
    }

    int SAT_Solver::create_complementary_binclause(int clause, int clause1, int clause2) {
        int var, *vars_signs, i=0, varssigns[4], sign, j=0;
        vars_signs = var_sign[clause];
        for(var=*vars_signs; var!=SAT_Clauses::NONE; var=*(vars_signs+=2)) {
            if (reason[var]!=SAT_Clauses::NO_REASON) {
                varssigns[i++]=var; varssigns[i++]=*(vars_signs+1); 
            }
        }
        if (reason[varssigns[2]]==clause1) {
            var=varssigns[2]; sign=varssigns[3];
            varssigns[2]=varssigns[0]; varssigns[3]=varssigns[1];
            varssigns[0]=var; varssigns[1]=sign;
        }
        if ((i!=4) || (check_reason(varssigns, clause, clause1, clause2)==SAT_Clauses::FALSE))
            printf("problem...");
        create_binaryclause(varssigns[0], 1-varssigns[1],
            varssigns[2], 1-varssigns[3], clause1, clause2);
        return SAT_Clauses::TRUE;
    }

    int SAT_Solver::get_satisfied_literal(int clause) {
        int var, *vars_signs;
        vars_signs = var_sign[clause];
        for(var=*vars_signs; var!=SAT_Clauses::NONE; var=*(vars_signs+=2)) {
            if (*(vars_signs+1) == var_current_value[var])
                return var;
        }
        printf("erreur");
        return SAT_Clauses::NONE;
    }

    void SAT_Solver::create_ternary_clauses(int var1, int sign1, int var2, int sign2, int var3, int sign3, int clause1, int clause2, int clause3) {
        int clause, *vars_signs, flag=SAT_Clauses::FALSE, *clauses1, *clauses2, *clauses3;
        if (sign1==SAT_Clauses::POSITIVE) clauses1=pos_in[var1]; else clauses1=neg_in[var1];
        if (sign2==SAT_Clauses::POSITIVE) clauses2=pos_in[var2]; else clauses2=neg_in[var2];
        if (sign3==SAT_Clauses::POSITIVE) clauses3=pos_in[var3]; else clauses3=neg_in[var3];
        vars_signs=NEW_CLAUSES[NEW_CLAUSES_fill_pointer++];
        vars_signs[0]=var1; vars_signs[1]=sign1;
        vars_signs[2]=var2; vars_signs[3]=sign2;
        vars_signs[4]=var3; vars_signs[5]=sign3;
        vars_signs[6]=SAT_Clauses::NONE;
        var_sign[NB_CLAUSE]=vars_signs;
        clause_state[NB_CLAUSE]=SAT_Clauses::ACTIVE;
        clause_length[NB_CLAUSE]=3;
        // if (NB_CLAUSE==305)
        // printf("aaa...");
        replace_clause(NB_CLAUSE, clause1, clauses1);
        replace_clause(NB_CLAUSE, clause2, clauses2);
        replace_clause(NB_CLAUSE, clause3, clauses3);
        NB_CLAUSE++;
    }  

    int SAT_Solver::non_linear_conflict(int empty_clause, int var1, int sign1, int var2, int sign2) {
        int var, sign, j;
        // driving unit clause is LINEAR_REASON_STACK1[2] (propagate
        // it resulting the empty_clause by simple non-linear derivation
        // var1, sign1, var2, and sign2 are the two literals of empty_clause
        var=get_satisfied_literal(LINEAR_REASON_STACK1[2]);
        sign=var_current_value[var];
        for(j=2; j<LINEAR_REASON_STACK1_fill_pointer-1; j++) {
            create_complementary_binclause(LINEAR_REASON_STACK1[j],
                LINEAR_REASON_STACK1[j+1],
                LINEAR_REASON_STACK1[j-1]);
            push_clause(LINEAR_REASON_STACK1[j], CLAUSES_TO_REMOVE);
        }
        push_clause(LINEAR_REASON_STACK1[j], CLAUSES_TO_REMOVE);
        create_ternary_clauses(var, sign, var1, sign1, var2, sign2,
            LINEAR_REASON_STACK1[2],
            empty_clause, empty_clause);
        create_ternary_clauses(var, 1-sign, var1, 1-sign1, var2, 1-sign2,
            LINEAR_REASON_STACK2[1],
            LINEAR_REASON_STACK1[1],
            LINEAR_REASON_STACK2[1]);
        push_clause(empty_clause, CLAUSES_TO_REMOVE);
        push_clause( LINEAR_REASON_STACK1[1], CLAUSES_TO_REMOVE);
        push_clause( LINEAR_REASON_STACK2[1], CLAUSES_TO_REMOVE);
        return SAT_Clauses::TRUE;
    }

        
    int SAT_Solver::linear_conflict(int clause) {
        int var, *vars_signs, i=0, varssigns[6], j=0, res;
        vars_signs = var_sign[clause];
        for(var=*vars_signs; var!=SAT_Clauses::NONE; var=*(vars_signs+=2)) {
            if (reason[var]!=SAT_Clauses::NO_REASON) {
                varssigns[i++]=var; varssigns[i++]=*(vars_signs+1); 
                if (i>4)
                    return SAT_Clauses::FALSE;
            }
        }
        if (i>4) return SAT_Clauses::FALSE;
        if (i==0){
            printf("bizzar...\n");
        }else {
            for(j=0; j<LINEAR_REASON_STACK1_fill_pointer; j++) 
                clause_involved[LINEAR_REASON_STACK1[j]]=SAT_Clauses::NONE;
            LINEAR_REASON_STACK1_fill_pointer=1; LINEAR_REASON_STACK2_fill_pointer=1;
            LINEAR_REASON_STACK1[0]=clause; LINEAR_REASON_STACK2[0]=clause;
            if (search_linear_reason1(varssigns[0])==SAT_Clauses::FALSE){
                return SAT_Clauses::FALSE;
            }else {
                if (i==4) {
                    res=search_linear_reason2(varssigns[2]);
                    if (res==SAT_Clauses::FALSE){
                        return SAT_Clauses::FALSE;
                    }else if (res==SIMPLE_NON_LINEAR_CASE) {
                        // printf("zskjehrz  \n");
                        return non_linear_conflict(clause, varssigns[0], varssigns[1], 
                            varssigns[2], varssigns[3]);
                    }
                    create_binaryclause(varssigns[0], 1-varssigns[1], 
                        varssigns[2], 1-varssigns[3], 
                        LINEAR_REASON_STACK1[1], LINEAR_REASON_STACK2[1]);
                    for(j=1; j<LINEAR_REASON_STACK2_fill_pointer-1; j++) {
                        create_complementary_binclause(LINEAR_REASON_STACK2[j],
                            LINEAR_REASON_STACK2[j+1],
                            LINEAR_REASON_STACK2[j-1]);
                        push_clause(LINEAR_REASON_STACK2[j], CLAUSES_TO_REMOVE);
                    }
                    push_clause(LINEAR_REASON_STACK2[j], CLAUSES_TO_REMOVE);
                }
                push_clause(clause, CLAUSES_TO_REMOVE);
                for(j=1; j<LINEAR_REASON_STACK1_fill_pointer-1; j++) {
                    create_complementary_binclause(LINEAR_REASON_STACK1[j],
                        LINEAR_REASON_STACK1[j+1],
                        LINEAR_REASON_STACK1[j-1]);
                    push_clause(LINEAR_REASON_STACK1[j], CLAUSES_TO_REMOVE);
                }
                push_clause(LINEAR_REASON_STACK1[j], CLAUSES_TO_REMOVE);
            }
            return SAT_Clauses::TRUE;
        }
    }

    void SAT_Solver::remove_linear_reasons() {
        int i, clause;
        for(i=0; i<LINEAR_REASON_STACK1_fill_pointer; i++) {
            clause=LINEAR_REASON_STACK1[i];
            clause_state[clause]=SAT_Clauses::PASSIVE;
            push_clause(clause, CLAUSE_STACK);
        }
        for(i=1; i<LINEAR_REASON_STACK2_fill_pointer; i++) {
            clause=LINEAR_REASON_STACK2[i];
            clause_state[clause]=SAT_Clauses::PASSIVE;
            push_clause(clause, CLAUSE_STACK);
        }
    }      

    int SAT_Solver::there_is_unit_clause( int var_to_check ) {
        int unitclause_position, unitclause, var, *vars_signs;

        for( unitclause_position = 0;
            unitclause_position < UNITCLAUSE_STACK_fill_pointer;
            unitclause_position++) {
            unitclause = UNITCLAUSE_STACK[ unitclause_position ];
            if ((clause_state[unitclause] == SAT_Clauses::ACTIVE)  && (clause_length[unitclause]>0)){
                vars_signs = var_sign[unitclause];
                for(var=*vars_signs; var!=SAT_Clauses::NONE; var=*(vars_signs+=2)) {
                    if ( var == var_to_check && var_state[var] == SAT_Clauses::ACTIVE ) {
                        return SAT_Clauses::TRUE;
                    }
                }
            }
        }
        return SAT_Clauses::FALSE;
    }

    int SAT_Solver::assign_and_unitclause_process( int var, int value, int starting_point ) {
        int clause;
        var_current_value[var] = value;
        var_rest_value[var] = SAT_Clauses::NONE;
        var_state[var] = SAT_Clauses::PASSIVE;
        push_clause(var, VARIABLE_STACK);
        if ((clause=my_reduce_clauses_for_fl(var))==SAT_Clauses::NO_CONFLICT) {
            remove_clauses(var);
            return my_unitclause_process( starting_point );
        }else {
            return clause;
        }
    }

    int SAT_Solver::store_reason_clauses( int clause, int starting ) {
        int *vars_signs, var, i;
        push_clause(clause, REASON_STACK);
        for(i=starting; i<REASON_STACK_fill_pointer; i++) {
            clause=REASON_STACK[i];
            vars_signs = var_sign[clause];
            for(var=*vars_signs; var!=SAT_Clauses::NONE; var=*(vars_signs+=2)) {
                if (reason[var]!=SAT_Clauses::NO_REASON) {
                    push_clause(reason[var], REASON_STACK);
                    reason[var]=SAT_Clauses::NO_REASON;
                }
            }
        }
        return i;
    }

    void SAT_Solver::remove_reason_clauses() {
        int i, clause;
        for(i=0; i<REASON_STACK_fill_pointer; i++) {
            clause=REASON_STACK[i];
            clause_state[clause]=SAT_Clauses::PASSIVE;
            push_clause(clause, CLAUSE_STACK);
        }
        REASON_STACK_fill_pointer=0;
    }

    int SAT_Solver::failed_literal( int conflict ) {
        int clause, var, la = 0;
        int saved_clause_stack_fill_pointer, saved_reducedclause_stack_fill_pointer,
            saved_unitclause_stack_fill_pointer, saved_variable_stack_fill_pointer,
            my_saved_clause_stack_fill_pointer, saved_reason_stack_fill_pointer;

        saved_clause_stack_fill_pointer= CLAUSE_STACK_fill_pointer;
        saved_reducedclause_stack_fill_pointer = REDUCEDCLAUSE_STACK_fill_pointer;
        saved_unitclause_stack_fill_pointer = UNITCLAUSE_STACK_fill_pointer;
        saved_variable_stack_fill_pointer=VARIABLE_STACK_fill_pointer;
        my_saved_clause_stack_fill_pointer= CLAUSE_STACK_fill_pointer;

        for( var=0; var < NB_VAR && la+conflict+NB_EMPTY<UB; var++ ) {
            if ( var_state[ var ] == SAT_Clauses::ACTIVE &&
                !there_is_unit_clause( var )) {
            simple_get_pos_clause_nb(var); simple_get_neg_clause_nb(var);
            if (nb_neg_clause_of_length2[ var ] > 1 &&  nb_pos_clause_of_length2[ var ] > 1 ) {
            //do {
                if ((clause=assign_and_unitclause_process(var, SAT_Clauses::FALSE, saved_unitclause_stack_fill_pointer))!=SAT_Clauses::NO_CONFLICT) {
            //  printf("One conflict found\n");
                saved_reason_stack_fill_pointer = store_reason_clauses( clause, 0 );
                reset_context(my_saved_clause_stack_fill_pointer,
                                saved_reducedclause_stack_fill_pointer,
                                saved_unitclause_stack_fill_pointer,
                                saved_variable_stack_fill_pointer);
                //remove_reason_clauses();
                //saved_reason_stack_fill_pointer = 0;
                my_saved_clause_stack_fill_pointer=CLAUSE_STACK_fill_pointer;
                if ((clause=assign_and_unitclause_process(var, SAT_Clauses::TRUE, saved_unitclause_stack_fill_pointer))>=0) {
                    la++;
                    store_reason_clauses( clause, saved_reason_stack_fill_pointer );
                    reset_context(my_saved_clause_stack_fill_pointer,
                                saved_reducedclause_stack_fill_pointer,
                                saved_unitclause_stack_fill_pointer,
                                saved_variable_stack_fill_pointer);
                    remove_reason_clauses();
                    my_saved_clause_stack_fill_pointer=CLAUSE_STACK_fill_pointer;
                } else {
                    REASON_STACK_fill_pointer = 0;
                    reset_context(my_saved_clause_stack_fill_pointer,
                                saved_reducedclause_stack_fill_pointer,
                                saved_unitclause_stack_fill_pointer,
                                saved_variable_stack_fill_pointer);
                }
                } else {
                reset_context(my_saved_clause_stack_fill_pointer,
                                saved_reducedclause_stack_fill_pointer,
                                saved_unitclause_stack_fill_pointer,
                                saved_variable_stack_fill_pointer);
                }
                //!!There could be more conflicts than just one
                //} while( clause != SAT_Clauses::NO_CONFLICT );
            }
            }
        }
        //if ( satisfiable() ) printf("Empty set\n");
        reset_context(saved_clause_stack_fill_pointer,
                        saved_reducedclause_stack_fill_pointer,
                        saved_unitclause_stack_fill_pointer,
                        saved_variable_stack_fill_pointer);
        return la;
    }

    int SAT_Solver::lookahead() {
        int saved_clause_stack_fill_pointer, saved_reducedclause_stack_fill_pointer,
            saved_unitclause_stack_fill_pointer, saved_variable_stack_fill_pointer,
            my_saved_clause_stack_fill_pointer,
            clause, conflict=0, var, *vars_signs, i, unitclause;

        // if (NB_BACK==160)
        //  printf("sqhvdzhj");

        CLAUSES_TO_REMOVE_fill_pointer=0;
        saved_clause_stack_fill_pointer= CLAUSE_STACK_fill_pointer;
        saved_reducedclause_stack_fill_pointer = REDUCEDCLAUSE_STACK_fill_pointer;
        saved_unitclause_stack_fill_pointer = UNITCLAUSE_STACK_fill_pointer;
        saved_variable_stack_fill_pointer=VARIABLE_STACK_fill_pointer;
        my_saved_clause_stack_fill_pointer= CLAUSE_STACK_fill_pointer;
        while ((clause=my_unitclause_process(0))!=SAT_Clauses::NO_CONFLICT) {
            conflict++;
            if (conflict+NB_EMPTY>=UB) break;
            
            if (linear_conflict(clause)==SAT_Clauses::TRUE) {
                conflict--; NB_EMPTY++;
                reset_context(my_saved_clause_stack_fill_pointer, 
                        saved_reducedclause_stack_fill_pointer,
                        saved_unitclause_stack_fill_pointer,
                        saved_variable_stack_fill_pointer);
                remove_linear_reasons();
                my_saved_clause_stack_fill_pointer=CLAUSE_STACK_fill_pointer;
            }else{
                push_clause(clause, REASON_STACK);
                for(i=0; i<REASON_STACK_fill_pointer; i++) {
                    clause=REASON_STACK[i]; vars_signs = var_sign[clause];
                    for(var=*vars_signs; var!=SAT_Clauses::NONE; var=*(vars_signs+=2)) {
                        if (reason[var]!=SAT_Clauses::NO_REASON) {
                            push_clause(reason[var], REASON_STACK);
                            reason[var]=SAT_Clauses::NO_REASON;
                        }
                    }
                }
                reset_context(my_saved_clause_stack_fill_pointer, 
                        saved_reducedclause_stack_fill_pointer,
                        saved_unitclause_stack_fill_pointer,
                        saved_variable_stack_fill_pointer);
                for(i=0; i<REASON_STACK_fill_pointer; i++) {
                    clause=REASON_STACK[i];
                    clause_state[clause]=SAT_Clauses::PASSIVE; push_clause(clause, CLAUSE_STACK);
                }
                REASON_STACK_fill_pointer=0;
                my_saved_clause_stack_fill_pointer=CLAUSE_STACK_fill_pointer;
            }
        }
        if ( conflict+NB_EMPTY < UB ) {
            reset_context(my_saved_clause_stack_fill_pointer, 
                saved_reducedclause_stack_fill_pointer,
                saved_unitclause_stack_fill_pointer,
                saved_variable_stack_fill_pointer); 
            conflict += failed_literal( conflict );
        }
            
        reset_context(saved_clause_stack_fill_pointer, 
                saved_reducedclause_stack_fill_pointer,
                saved_unitclause_stack_fill_pointer,
                saved_variable_stack_fill_pointer);
        if (conflict+NB_EMPTY>=UB) 
            return SAT_Clauses::NONE;
        for (i=0; i<CLAUSES_TO_REMOVE_fill_pointer; i++) {
            clause=CLAUSES_TO_REMOVE[i];
            push_clause(clause, CLAUSE_STACK); clause_state[clause]=SAT_Clauses::PASSIVE;
        }
        CLAUSES_TO_REMOVE_fill_pointer=0;
        return conflict;
    }

    int SAT_Solver::satisfy_unitclause(int unitclause) {
        int *vars_signs, var, clause;

        vars_signs = var_sign[unitclause];
        for(var=*vars_signs; var!=SAT_Clauses::NONE; var=*(vars_signs+=2)) {
            if (var_state[var] == SAT_Clauses::ACTIVE ){
                var_current_value[var] = *(vars_signs+1);
                var_rest_value[var] = SAT_Clauses::NONE;
                reason[var]=unitclause;
                var_state[var] = SAT_Clauses::PASSIVE;
                push_clause(var, VARIABLE_STACK);
                if ((clause=my_reduce_clauses(var))==SAT_Clauses::NO_CONFLICT) {
                    remove_clauses(var);
                    return SAT_Clauses::NO_CONFLICT;
                }else 
                    return clause;
            }
        }
        return SAT_Clauses::NO_CONFLICT;
    }
    
    int SAT_Solver::my_unitclause_process(int starting_point) {
        int unitclause, var, *vars_signs, unitclause_position,clause,
            my_unitclause_position, my_unitclause;

        for (unitclause_position = starting_point;unitclause_position < UNITCLAUSE_STACK_fill_pointer;unitclause_position++) {
            unitclause = UNITCLAUSE_STACK[unitclause_position];
            if ((clause_state[unitclause] == SAT_Clauses::ACTIVE)  && (clause_length[unitclause]>0)) {
                MY_UNITCLAUSE_STACK_fill_pointer=0;
                if ((clause=satisfy_unitclause(unitclause)) != SAT_Clauses::NO_CONFLICT){
                    return clause;
                }else{
                    for (my_unitclause_position = 0; 
                        my_unitclause_position < MY_UNITCLAUSE_STACK_fill_pointer;
                        my_unitclause_position++) {
                    my_unitclause = MY_UNITCLAUSE_STACK[my_unitclause_position];
                    if ((clause_state[my_unitclause] == SAT_Clauses::ACTIVE)  
                        && (clause_length[my_unitclause]>0)) {
                        if ((clause=satisfy_unitclause(my_unitclause)) != SAT_Clauses::NO_CONFLICT)
                        return clause;
                    }     
                    }
                }
            }
        }
        return SAT_Clauses::NO_CONFLICT;
    }

    int SAT_Solver::get_complement(int lit) {
        if (positive(lit))
            return lit+NB_VAR;
        else
            return lit-NB_VAR;
    }

    void SAT_Solver::create_unitclause(int lit, int subsumedclause, int *clauses) {
        int clause, *vars_signs, flag=SAT_Clauses::FALSE;

        vars_signs=NEW_CLAUSES[NEW_CLAUSES_fill_pointer++];
        if (lit<NB_VAR) {
            vars_signs[0]=lit;
            vars_signs[1]=SAT_Clauses::POSITIVE;
        }else{
            vars_signs[0]=lit-NB_VAR;
            vars_signs[1]=SAT_Clauses::NEGATIVE;
        }
        vars_signs[2]=SAT_Clauses::NONE;
        var_sign[NB_CLAUSE]=vars_signs;
        clause_state[NB_CLAUSE]=SAT_Clauses::ACTIVE;
        clause_length[NB_CLAUSE]=1;
        push_clause(NB_CLAUSE, UNITCLAUSE_STACK);

        for(clause=*clauses; clause!=SAT_Clauses::NONE; clause=*(++clauses)) {
            if (clause==subsumedclause) {
                *clauses=NB_CLAUSE;
                SAVED_CLAUSE_POSITIONS[SAVED_CLAUSES_fill_pointer]=clauses;
                push_clause(subsumedclause, SAVED_CLAUSES);
                flag=SAT_Clauses::TRUE;
                break;
            }
        }
        if (flag==SAT_Clauses::FALSE)
            printf("erreur ");
        NB_CLAUSE++;
    }

    int SAT_Solver::verify_resolvent(int lit, int clause1, int clause2) {
        int *vars_signs1, *vars_signs2, lit1, lit2, temp, flag=SAT_Clauses::FALSE, var, nb=0;

        if ((clause_state[clause1]!=SAT_Clauses::ACTIVE) || (clause_state[clause2]!=SAT_Clauses::ACTIVE))
            printf("erreur ");
        if ((clause_length[clause1]!=2) || (clause_length[clause2]!=2))
            printf("erreur ");
        vars_signs1=var_sign[clause1];
        vars_signs2=var_sign[clause2];
        for(var=*vars_signs1; var!=SAT_Clauses::NONE; var=*(vars_signs1+=2)) {
            if (var_state[var] == SAT_Clauses::ACTIVE ) {
                nb++;
                if (*(vars_signs1+1)==SAT_Clauses::POSITIVE) 
                    temp=var;
                else
                    temp=var+NB_VAR;
                if (temp==lit)
                    flag=SAT_Clauses::TRUE;
                else
                    lit1=temp;
            }
        }
        if ((nb!=2) || (flag==SAT_Clauses::FALSE))
            printf("erreur ");
        nb=0; flag=SAT_Clauses::FALSE;
        for(var=*vars_signs2; var!=SAT_Clauses::NONE; var=*(vars_signs2+=2)) {
            if (var_state[var] == SAT_Clauses::ACTIVE ) {
                nb++;
                if (*(vars_signs2+1)==SAT_Clauses::POSITIVE) 
                    temp=var;
                else
                    temp=var+NB_VAR;
                if (temp==lit)
                    flag=SAT_Clauses::TRUE;
                else
                    lit2=temp;
            }
        }
        if ((nb!=2) || (flag==SAT_Clauses::FALSE))
            printf("erreur ");
        if (!complement(lit1, lit2))
            printf("erreur ");
    }

    int SAT_Solver::searching_two_clauses_to_fix_neglit(int clause, int lit) {
        int lit1, clause1, var1, opp_lit1;
        if (lit_to_fix[clause]==SAT_Clauses::NONE) {
            lit_to_fix[clause]=lit;
        }else{
            lit1=lit_to_fix[clause];
            var1=get_var_from_lit(lit1);
            //  if (var_state[var1]!=SAT_Clauses::ACTIVE)
            //    printf("erreur2  ");
            opp_lit1=get_complement(lit1);
            clause1=fixing_clause[opp_lit1];
            if ((clause1!= SAT_Clauses::NONE) && (clause_state[clause1]==SAT_Clauses::ACTIVE)) {
                fixing_clause[opp_lit1]=SAT_Clauses::NONE;
                lit_involved_in_clause[opp_lit1]=SAT_Clauses::NONE;
                // verify_resolvent(lit, clause1, clause);
                push_clause(clause1, CLAUSE_STACK);
                clause_state[clause1]=SAT_Clauses::PASSIVE;
                push_clause(clause, CLAUSE_STACK);
                clause_state[clause]=SAT_Clauses::PASSIVE;
                create_unitclause(lit, clause1, neg_in[lit-NB_VAR]);
                var1=get_var_from_lit(lit1);
                nb_neg_clause_of_length2[var1]--;
                nb_pos_clause_of_length2[var1]--;
                return SAT_Clauses::TRUE;
            }else {
                fixing_clause[lit1]=clause;
                push_clause(lit1, CANDIDATE_LITERALS);
                lit_involved_in_clause[lit1]=clause;
                push_clause(lit1, INVOLVED_LIT_STACK);
            }
        }
        return SAT_Clauses::FALSE;
    }

    int SAT_Solver::simple_get_neg_clause_nb(int var) {
        my_type neg_clause1_nb=0,neg_clause3_nb = 0, neg_clause2_nb = 0;
        int *clauses, clause, i;
        clauses = neg_in[var]; MY_UNITCLAUSE_STACK_fill_pointer=0;

        for(clause=*clauses; clause!=SAT_Clauses::NONE; clause=*(++clauses))
            if ((clause_state[clause] == SAT_Clauses::ACTIVE) && (clause_length[clause]==2))
                neg_clause2_nb++;
            nb_neg_clause_of_length2[var] = neg_clause2_nb;
            return neg_clause2_nb;
    }

    int SAT_Solver::simple_get_pos_clause_nb(int var) {
        my_type pos_clause1_nb=0,pos_clause3_nb = 0, pos_clause2_nb = 0;
        int *clauses, clause, i;
        clauses = pos_in[var]; MY_UNITCLAUSE_STACK_fill_pointer=0;

        for(clause=*clauses; clause!=SAT_Clauses::NONE; clause=*(++clauses))
            if ((clause_state[clause] == SAT_Clauses::ACTIVE) && (clause_length[clause]==2))
                pos_clause2_nb++;
            nb_pos_clause_of_length2[var] = pos_clause2_nb;
            return pos_clause2_nb;
    }

    int SAT_Solver::get_neg_clause_nb(int var) {
        my_type neg_clause1_nb=0,neg_clause3_nb = 0, neg_clause2_nb = 0;
        int *clauses, clause, i;
        clauses = neg_in[var]; MY_UNITCLAUSE_STACK_fill_pointer=0;

        for(clause=*clauses; clause!=SAT_Clauses::NONE; clause=*(++clauses)) {
            if ((clause_state[clause] == SAT_Clauses::ACTIVE) && (clause_length[clause]>0)) {
                switch(clause_length[clause]) {
                    case 1:
                        neg_clause1_nb++; 
                        push_clause(clause, MY_UNITCLAUSE_STACK); break;
                    case 2:
                        neg_clause2_nb++; 
                        if (searching_two_clauses_to_fix_neglit(clause, var+NB_VAR)==SAT_Clauses::TRUE) {
                            neg_clause2_nb-=2; neg_clause1_nb++; 
                        }
                        break;
                    default:
                        neg_clause3_nb++;
                        break;
                }
            }
        }
        for(i=0; i<CANDIDATE_LITERALS_fill_pointer; i++) 
            fixing_clause[CANDIDATE_LITERALS[i]]=SAT_Clauses::NONE;
        CANDIDATE_LITERALS_fill_pointer=0;
        nb_neg_clause_of_length1[var] = neg_clause1_nb;
        nb_neg_clause_of_length2[var] = neg_clause2_nb;
        nb_neg_clause_of_length3[var] = neg_clause3_nb;
        return neg_clause1_nb+neg_clause2_nb + neg_clause3_nb;
    }

    #define OTHER_LIT_FIXED 1
    #define THIS_LIT_FIXED 2

    int SAT_Solver::searching_two_clauses_to_fix_poslit(int clause, int lit) {
        int lit1, clause1, var1, opp_lit1;
        if (lit_to_fix[clause]==SAT_Clauses::NONE) {
            lit_to_fix[clause]=lit;
        }else{
            lit1=lit_to_fix[clause];
            var1=get_var_from_lit(lit1);
            //   if (var_state[var1]!=SAT_Clauses::ACTIVE)
            //   printf("erreur2  ");
            clause1=lit_involved_in_clause[lit1];
            if ((clause1!=SAT_Clauses::NONE) && (clause_state[clause1]==SAT_Clauses::ACTIVE)) {
                //  verify_resolvent(lit1, clause1, clause);
                push_clause(clause1, CLAUSE_STACK);
                clause_state[clause1]=SAT_Clauses::PASSIVE;
                push_clause(clause, CLAUSE_STACK);
                clause_state[clause]=SAT_Clauses::PASSIVE;
                if (lit1<NB_VAR) {
                    create_unitclause(lit1, clause1, pos_in[lit1]);
                    nb_pos_clause_of_length2[lit1]-=2;
                    nb_pos_clause_of_length1[lit1]++;
                }else{
                    create_unitclause(lit1, clause1, neg_in[lit1-NB_VAR]);
                    nb_neg_clause_of_length2[lit1-NB_VAR]-=2;
                    nb_neg_clause_of_length1[lit1-NB_VAR]++;
                }
                return OTHER_LIT_FIXED;
            }else{
                opp_lit1=get_complement(lit1);
                clause1=fixing_clause[opp_lit1];
                if ((clause1!= SAT_Clauses::NONE) && (clause_state[clause1]==SAT_Clauses::ACTIVE)) {
                fixing_clause[opp_lit1]=SAT_Clauses::NONE;
                //	verify_resolvent(lit, clause1, clause);
                push_clause(clause1, CLAUSE_STACK);
                clause_state[clause1]=SAT_Clauses::PASSIVE;
                push_clause(clause, CLAUSE_STACK);
                clause_state[clause]=SAT_Clauses::PASSIVE;
                create_unitclause(lit, clause1, pos_in[lit]);
                var1=get_var_from_lit(lit1);
                nb_neg_clause_of_length2[var1]--;
                nb_pos_clause_of_length2[var1]--;
                return THIS_LIT_FIXED;
                }
                else {
                fixing_clause[lit1]=clause;
                push_clause(lit1, CANDIDATE_LITERALS);
                }
            }
        }
        return SAT_Clauses::FALSE;
    }

    int SAT_Solver::get_pos_clause_nb(int var) {
        my_type pos_clause1_nb=0, pos_clause3_nb = 0, pos_clause2_nb = 0;
        int *clauses, clause, clause1, i;
        clauses = pos_in[var];
        for(clause=*clauses; clause!=SAT_Clauses::NONE; clause=*(++clauses)) {
            if ((clause_state[clause] == SAT_Clauses::ACTIVE) && (clause_length[clause]>0)) {
                switch(clause_length[clause]) {
                    case 1:
                        if (MY_UNITCLAUSE_STACK_fill_pointer>0) {
                            clause1=pop_clause(MY_UNITCLAUSE_STACK);
                            clause_state[clause]=SAT_Clauses::PASSIVE;
                            push_clause(clause, CLAUSE_STACK);
                            clause_state[clause1]=SAT_Clauses::PASSIVE;
                            push_clause(clause1, CLAUSE_STACK);
                            nb_neg_clause_of_length1[var]--;
                            NB_EMPTY++;
                        }else pos_clause1_nb++; 
                        break;
                    case 2:
                        pos_clause2_nb++; 
                        switch(searching_two_clauses_to_fix_poslit(clause, var)) {
                            case OTHER_LIT_FIXED: nb_neg_clause_of_length2[var]--;
                            pos_clause2_nb--;
                            break;
                            case THIS_LIT_FIXED: pos_clause2_nb-=2;
                            pos_clause1_nb++;
                            break;
                        }
                        break;
                    default:
                        pos_clause3_nb++;
                        break;
                }
            }
        }
        for(i=0; i<CANDIDATE_LITERALS_fill_pointer; i++) 
            fixing_clause[CANDIDATE_LITERALS[i]]=SAT_Clauses::NONE;
        CANDIDATE_LITERALS_fill_pointer=0;
        for(i=0; i<INVOLVED_LIT_STACK_fill_pointer; i++) 
            lit_involved_in_clause[INVOLVED_LIT_STACK[i]]=SAT_Clauses::NONE;
        INVOLVED_LIT_STACK_fill_pointer=0;
        nb_pos_clause_of_length1[var] = pos_clause1_nb;
        nb_pos_clause_of_length2[var] = pos_clause2_nb;
        nb_pos_clause_of_length3[var] = pos_clause3_nb;
        return pos_clause1_nb+pos_clause2_nb + pos_clause3_nb;
    }

    int SAT_Solver::satisfy_literal(int lit) {
        int var;
        if (positive(lit)) {
            if (var_state[lit]==SAT_Clauses::ACTIVE) {
                var_current_value[lit] = SAT_Clauses::TRUE;
                if (reduce_clauses(lit)==SAT_Clauses::FALSE)
                    return SAT_Clauses::NONE;
                var_rest_value[lit]=SAT_Clauses::NONE;
                var_state[lit] = SAT_Clauses::PASSIVE;
                push_clause(lit, VARIABLE_STACK);
                remove_clauses(lit);
            }else if (var_current_value[lit]==SAT_Clauses::FALSE) return SAT_Clauses::NONE;
        }else {
            var = get_var_from_lit(lit);
            if (var_state[var]==SAT_Clauses::ACTIVE) {
                var_current_value[var] = SAT_Clauses::FALSE;
                if (reduce_clauses(var)==SAT_Clauses::FALSE) return SAT_Clauses::NONE;
                var_rest_value[var]=SAT_Clauses::NONE;
                var_state[var] = SAT_Clauses::PASSIVE;
                push_clause(var, VARIABLE_STACK);
                remove_clauses(var);
            }else if (var_current_value[var]==SAT_Clauses::TRUE) return SAT_Clauses::NONE;
        }
        return SAT_Clauses::TRUE;
    }

    int SAT_Solver::assign_value(int var, int current_value, int rest_value) {
        if (var_state[var]==SAT_Clauses::PASSIVE)
            printf("erreur1...\n");
        var_state[var] = SAT_Clauses::PASSIVE;
        push_clause(var, VARIABLE_STACK);
        var_current_value[var] = current_value;
        var_rest_value[var] = rest_value;
        if (reduce_clauses(var)==SAT_Clauses::NONE) 
            return SAT_Clauses::NONE;
        remove_clauses(var);
        return SAT_Clauses::TRUE;
    }

    int SAT_Solver::unitclause_process() {
        int unitclause, var, *vars_signs, unitclause_position,clause;
        
        for (unitclause_position = 0; unitclause_position < UNITCLAUSE_STACK_fill_pointer; unitclause_position++) {
            unitclause = UNITCLAUSE_STACK[unitclause_position];
            if ((clause_state[unitclause] == SAT_Clauses::ACTIVE)  && (clause_length[unitclause]>0)) {
                vars_signs = var_sign[unitclause];
                for(var=*vars_signs; var!=SAT_Clauses::NONE; var=*(vars_signs+=2)) {
                    if (var_state[var] == SAT_Clauses::ACTIVE ){
                        var_current_value[var] = *(vars_signs+1);
                        var_rest_value[var] = SAT_Clauses::NONE;
                        var_state[var] = SAT_Clauses::PASSIVE;
                        push_clause(var, VARIABLE_STACK);
                        if ((clause=reduce_clauses(var)) !=SAT_Clauses::NONE) {
                            remove_clauses(var);
                            break;
                        }
                        else {
                            return SAT_Clauses::NONE;
                        }
                    }
                }
            }     
        }
        return SAT_Clauses::TRUE;
    }

    int SAT_Solver::choose_and_instantiate_variable() {
        int var, nb=0, chosen_var=SAT_Clauses::NONE,cont=0, cont1; 
        int  i;
        float posi, nega;
        int a,b,c,clause;
        float poid, max_poid = -1.0; 
        my_type pos2, neg2, flag=0;
        NB_BRANCHE++;

        // if (NB_BRANCHE==173)
        //   printf("zerza ");

        if (lookahead()==SAT_Clauses::NONE)
            return SAT_Clauses::NONE;

        if (UB-NB_EMPTY==1)
            if (unitclause_process() ==SAT_Clauses::NONE)
                return SAT_Clauses::NONE;

        for (clause=0; clause<NB_CLAUSE; clause++) 
            lit_to_fix[clause]=SAT_Clauses::NONE;

        for (var = 0; var < NB_VAR; var++) {
            if (var_state[var] == SAT_Clauses::ACTIVE) {
                reduce_if_negative[var]=0;
                reduce_if_positive[var]=0;
                if (get_neg_clause_nb(var) == 0) {
                    NB_MONO++;
                    var_current_value[var] = SAT_Clauses::TRUE;
                    var_rest_value[var] = SAT_Clauses::NONE;
                    var_state[var] = SAT_Clauses::PASSIVE;
                    push_clause(var, VARIABLE_STACK);
                    remove_clauses(var);
                }else if (get_pos_clause_nb(var) == 0) {
                    NB_MONO++;
                    var_current_value[var] = SAT_Clauses::FALSE;
                    var_rest_value[var] = SAT_Clauses::NONE;
                    var_state[var] = SAT_Clauses::PASSIVE;
                    push_clause(var, VARIABLE_STACK);
                    remove_clauses(var);
                }else if (nb_neg_clause_of_length1[var]+NB_EMPTY>=UB) {
                    flag++;
                    if (assign_value(var, SAT_Clauses::FALSE, SAT_Clauses::NONE)==SAT_Clauses::NONE)
                        return SAT_Clauses::NONE;
                }else if (nb_pos_clause_of_length1[var]+NB_EMPTY>=UB) {
                    flag++;
                    if (assign_value(var, SAT_Clauses::TRUE, SAT_Clauses::NONE)==SAT_Clauses::NONE)
                        return SAT_Clauses::NONE;
                }else if (nb_neg_clause_of_length1[var] >=
                        nb_pos_clause_of_length1[var]
                        + nb_pos_clause_of_length2[var]
                        + nb_pos_clause_of_length3[var]){
                    flag++;
                    if (assign_value(var, SAT_Clauses::FALSE, SAT_Clauses::NONE)==SAT_Clauses::NONE)
                        return SAT_Clauses::NONE;
                }else if (nb_pos_clause_of_length1[var] >= nb_neg_clause_of_length1[var]
                        + nb_neg_clause_of_length2[var]
                        + nb_neg_clause_of_length3[var]){
                    flag++;
                    if (assign_value(var, SAT_Clauses::TRUE, SAT_Clauses::NONE)==SAT_Clauses::NONE)
                        return SAT_Clauses::NONE;
                }else{
                    if (nb_neg_clause_of_length1[var]>nb_pos_clause_of_length1[var]) {
                    cont+=nb_pos_clause_of_length1[var];
                    }else{
                    cont+=nb_neg_clause_of_length1[var];
                    }
                }
            }
        }

        if (cont+NB_EMPTY>=UB)
            return SAT_Clauses::NONE;
        for (var = 0; var < NB_VAR; var++) {
            if (var_state[var] == SAT_Clauses::ACTIVE) {
                /*   
                if (nb_neg_clause_of_length1[var]>nb_pos_clause_of_length1[var])
                cont1=cont-nb_pos_clause_of_length1[var];
                else cont1=cont-nb_neg_clause_of_length1[var];
                if (nb_neg_clause_of_length1[var]+cont1+NB_EMPTY>=UB) {
                if (assign_value(var, SAT_Clauses::FALSE, SAT_Clauses::NONE)==SAT_Clauses::NONE)
                return SAT_Clauses::NONE;
                }
                else if (nb_pos_clause_of_length1[var]+cont1+NB_EMPTY>=UB) {
                if (assign_value(var, SAT_Clauses::TRUE, SAT_Clauses::NONE)==SAT_Clauses::NONE)
                return SAT_Clauses::NONE;
                }
                else {
                */
                reduce_if_positive[var]=nb_neg_clause_of_length1[var]*2+
                nb_neg_clause_of_length2[var]*4+ 
                nb_neg_clause_of_length3[var];
                reduce_if_negative[var]=nb_pos_clause_of_length1[var]*2+
                nb_pos_clause_of_length2[var]*4+ 
                nb_pos_clause_of_length3[var];
                poid=reduce_if_positive[var]*reduce_if_negative[var]*64+
                reduce_if_positive[var]+reduce_if_negative[var];
                if (poid>max_poid) {
                chosen_var=var;
                max_poid=poid;
                }
                //	 }
            }
        }
        if (chosen_var == SAT_Clauses::NONE) return SAT_Clauses::FALSE;
        //      printf("%d \n",NB_BACK);
        //   printf("Chosen_va %d\n",chosen_var);
        saved_clause_stack[chosen_var] = CLAUSE_STACK_fill_pointer;
        saved_reducedclause_stack[chosen_var] = REDUCEDCLAUSE_STACK_fill_pointer;
        saved_unitclause_stack[chosen_var] = UNITCLAUSE_STACK_fill_pointer;
        saved_nb_empty[chosen_var]=NB_EMPTY;
        // return assign_value(chosen_var, SAT_Clauses::TRUE, SAT_Clauses::FALSE);
        saved_nb_clause[chosen_var]=NB_CLAUSE;
        saved_saved_clauses[chosen_var]=SAVED_CLAUSES_fill_pointer;
        saved_new_clauses[chosen_var]=NEW_CLAUSES_fill_pointer;
        
        if (reduce_if_positive[chosen_var]<reduce_if_negative[chosen_var])
            return assign_value(chosen_var, SAT_Clauses::TRUE, SAT_Clauses::FALSE);
        else
            return assign_value(chosen_var, SAT_Clauses::FALSE, SAT_Clauses::TRUE);
    }

    int SAT_Solver::dpl() {
        int var, nb;
        do {
            if (VARIABLE_STACK_fill_pointer==NB_VAR) {
                UB=NB_EMPTY; 
                nb=verify_solution();
                if (nb!=NB_EMPTY)
                    printf("problem nb...");
                printf("o %d\n", UB);
                for (var = 0; var < NB_VAR; var++)
                    var_best_v alue[var] = var_current_value[var];
                while (backtracking()==SAT_Clauses::NONE);
                    if (VARIABLE_STACK_fill_pointer==0)
                        break;
            }
            if (UB-NB_EMPTY==1)
                if (unitclause_process() ==SAT_Clauses::NONE)
                    while (backtracking()==SAT_Clauses::NONE);
            if (choose_and_instantiate_variable()==SAT_Clauses::NONE)
                while (backtracking()==SAT_Clauses::NONE);
                // else if (lookahead()==SAT_Clauses::NONE) 
                //  while (backtracking()==SAT_Clauses::NONE);
        } while (VARIABLE_STACK_fill_pointer > 0);
    }

    int SAT_Solver::init() {
        int var, clause;
        NB_EMPTY=0; REAL_NB_CLAUSE=NB_CLAUSE;
        UNITCLAUSE_STACK_fill_pointer=0;
        VARIABLE_STACK_fill_pointer=0;
        CLAUSE_STACK_fill_pointer = 0;
        REDUCEDCLAUSE_STACK_fill_pointer = 0;
        for (var=0; var<NB_VAR; var++) {
            reason[var]=SAT_Clauses::NO_REASON;
            fixing_clause[var]=SAT_Clauses::NONE;
            fixing_clause[var+NB_VAR]=SAT_Clauses::NONE;
            lit_involved_in_clause[var]=SAT_Clauses::NONE;
            lit_involved_in_clause[var+NB_VAR]=SAT_Clauses::NONE;
        }
        for (clause=0; clause<NB_CLAUSE; clause++) {
            lit_to_fix[clause]=SAT_Clauses::NONE;
            clause_involved[clause]=SAT_Clauses::NONE;
        }
        return 0;
    }
}