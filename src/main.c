/* 
 * File:   main.c
 * Author: oro
 *
 * Created on 30 / desembre / 2010, 11:35
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "utils.h"
#include "sequence_utils.h"
#include "distance_matrix_utils.h"
#include "guide_tree_utils.h"
#include "parameters_utils.h"

int my_rank;
int np;

//ARREGLAR ELS OPENS DE FITGERS, TRACTAR ERRORS
//ARREGLAR PARAMETRES
//PARAMETRE PER REPES O NO - fet
//PARAMETRE PER IDENTIFICAR ARBRES DE FORA - Fet
//PARAMETRE PER PODER AFEGIR UNA LLISTA DARBRES - fet

char* other_pg_identify_tree(Parameters *P);

/*MTA main function*/

int main(int argc, char** argv) {

    Parameters *P;
    Sequence *S;
    Distance_matrix *DM;
    
    double tini=0.0, tfi=0.0;
    int mpi_mode = 0;
    int other_pg = 0;

    other_pg = search_mode(argc, argv, "-other_pg");
    
    /*Other functionalities as guide tree identification*/
  
    if(other_pg){
        P = read_parameters(argc, argv);
        if(strm(P->other_pg_mode, "-tree_id")){
            char* tree_id;
            tree_id = other_pg_identify_tree(P);
            fprintf(stdout, "\nTree id: %s\n", tree_id);
            vfree(tree_id);
            free_parameters(P);
            exit(EXIT_SUCCESS);
        }
        else {
            fprintf(stderr, "ERROR: BAD OTHER PROGRAM MODE\n");
            exit(EXIT_FAILURE);
        }
        free_parameters(P);
    }
    
    mpi_mode = search_mode(argc, argv, "-mpi");
    tini = get_time();   

#ifdef MPI_FLAG /*DISTRIBUTED VERSION*/
    
    if(mpi_mode == 1){   

        int nseq, max_len, min_len;
        int position=0, tam=0;
        char *buffer;
       
        MPI_Init(&argc, &argv);

        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &np);
        
        P = read_parameters(argc, argv);
        
        tam = 3 * sizeof(int);
        buffer = (char *) malloc(tam);

        if(my_rank == 0){      
            
            fprintf(stdout, "\n-------------------MULTIPLE TREES METHOD---------------------------\n\n");           
            fprintf(stdout, "Input Files:\n");
            fprintf(stdout, "\tInput: (Sequences) %s.%s\n", (P->F)->name, (P->F)->suffix);
            fprintf(stdout, "\tInput: (Method) %s\n", P->align_method);
            fprintf(stdout, "\tInput: (Number of trees) %d\n", P->ntree);
       
            S=read_fasta_sequences(P->seqfile);
            
            if(S==NULL){
                fprintf(stderr, "Error - No Sequences\n");
                return(EXIT_FAILURE);
            }
            
            nseq = S->nseq;
            max_len = S->max_len;
            min_len = S->min_len;

            MPI_Pack(&nseq, 1, MPI_INT, buffer, tam, &position, MPI_COMM_WORLD);
            MPI_Pack(&max_len, 1, MPI_INT, buffer, tam, &position, MPI_COMM_WORLD);
            MPI_Pack(&min_len, 1, MPI_INT, buffer, tam, &position, MPI_COMM_WORLD);

            MPI_Bcast(buffer, position, MPI_PACKED, 0, MPI_COMM_WORLD);

            DM = make_distance_matrix(S, P->dm_method);
        }
        else {

            MPI_Bcast(buffer, tam, MPI_PACKED, 0, MPI_COMM_WORLD);

            MPI_Unpack(buffer, tam, &position, &nseq, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, tam, &position, &max_len, 1, MPI_INT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, tam, &position, &min_len, 1, MPI_INT, MPI_COMM_WORLD);      

            S = declare_sequence(nseq, max_len, min_len);
            DM = declare_distance_matrix(nseq);
        }
        
        mta_program_mpi(P, S, DM, nseq);

        if(my_rank == 0){
            printf("\nOUTPUT RESULTS:\n");
            printf("\tFile Type: Guide Tree / Format: Newick / File Name: %s.dnd\n", (P->F)->name);
            printf("\tFile Type: Alignment / Format: Fasta / File Name:  %s.fasta_aln\n", (P->F)->name);
            printf("\tFile Type: Scores / Format: txt / File Name: %s.scores\n", (P->F)->name);
            tfi = get_time();
            printf("\nMultiples Trees Method Time: %f seconds\n", tfi-tini);
        }
        
        MPI_FINALIZE();
        
    }
    
#elif SEQ_FLAG

    if(mpi_mode == 0) { /*Sequential version*/
        
        double max_score;
        
        P = read_parameters(argc, argv);

        fprintf(stdout, "\n-------------------MULTIPLE TREES METHOD---------------------------\n\n");           
        fprintf(stdout, "Input Files:\n");
        fprintf(stdout, "\tInput: (Sequences) %s.%s\n", (P->F)->name, (P->F)->suffix);
        fprintf(stdout, "\tInput: (Method) %s\n", P->align_method);
        fprintf(stdout, "\tInput: (Number of trees) %d\n", P->ntree);
        
        S=read_fasta_sequences(P->seqfile);

        if(S==NULL){
            fprintf(stderr, "Error - No Sequences\n");
            return(EXIT_FAILURE);
        }
 
        DM = make_distance_matrix(S, P->dm_method);

        if(P->rep){
            max_score = mta_program_repeated_trees(P, S, DM);
        }
        else{
            max_score = mta_program_no_repeated_trees(P, S, DM);
        }
        
        printf("\nOUTPUT RESULTS:\n");
        printf("\tFile Type: Guide Tree / Format: Newick / File Name: %s.dnd\n", (P->F)->name);
        printf("\tFile Type: Alignment / Format: Fasta / File Name:  %s.fasta_aln\n", (P->F)->name);
        printf("\tFile Type: Scores / Format: txt / File Name: %s.scores\n", (P->F)->name);
        
        tfi = get_time();
        printf("\nMultiples Trees Method Time: %f seconds\n", tfi-tini);
        
       
    } 
#endif
    
    free_parameters(P);
    free_distance_matrix(DM);
    free_sequence(S);  
        
    return (EXIT_SUCCESS); 
}

/*Identify tree id*/

 char* other_pg_identify_tree(Parameters *P){

    char *tree_id;
    NT_node **T;
    Sequence *S;
    int tot_node=0, tree_id_lenght=0;

    fprintf(stdout, "\nTREE IDENTIFICATION METHOD:\n");
    S=read_fasta_sequences(P->seqfile);
    tree_id_lenght = (5 * S->nseq);

    T=read_tree(P->tree_file, &tot_node, S->nseq, S->seq_names, S->seq_id, NULL);
    T[3][0]->level = 0;
    tree_id = vcalloc(tree_id_lenght, sizeof(char));
    tree_id = identify_tree(T[3][0], tree_id, S->nseq);

    free_sequence(S);  
    return tree_id;

 }