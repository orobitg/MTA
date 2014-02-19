
#ifndef ALIGNMENT_UTILS_H
#define	ALIGNMENT_UTILS_H

#include "utils.h"
#include "distance_matrix_utils.h"
#include "parameters_utils.h"
#include "guide_tree_utils.h"
#include "evaluation_utils.h"


//SERIAL METHOD

int mta_program_repeated_trees(Parameters *P, Sequence *S, Distance_matrix *DM);
int mta_program_no_repeated_trees(Parameters *P, Sequence *S, Distance_matrix *DM);

//MPI METHOD
#ifdef MPI_FLAG
void mta_program_mpi(Parameters *P, Sequence *S, Distance_matrix *DM, int nseq);
void master_generate_mta(Parameters *P, Sequence *S, Distance_matrix *DM);
void worker_generate_mta(Parameters *P, Sequence *S, Distance_matrix *DM);
#endif

void align_tree(Parameters *P, char *treename, char *libname, char *alnname, int ntree);
char** read_tree_list(char* tree_list, int ntree);

#endif	/* ALIGNMENT_UTILS_H */