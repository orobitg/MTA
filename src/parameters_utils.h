/* 
 * File:   parameters_utils.h
 * Author: oro
 *
 * Created on 5 / gener / 2011, 12:45
 */

#ifndef PARAMETERS_UTILS_H
#define	PARAMETERS_UTILS_H

#include "utils.h"
#include "guide_tree_utils.h"

/*Parameters list*/

struct Parameters{

    char *seqfile;
    int ntree;
    int seq_activated;
    int rep;
    int nattemps;
    int tree_id_activated;
    int random_percentage;
    char *dm_method;
    char *align_method;
    char *align_method_bin;
    char *score_method;
    char *outdir;
    char *tree_id_file;
    char *tree_file;
    int mpi_mode;
    int other_pg;
    char *other_pg_mode;
    
    //SP & NORMD parameters
    char *mat; //matrix
    float gop; //gap opening penalization 
    float gep; //gap extended penalization
    float egap; //end gap penalization
    
    //STRIKE & IRMSD parameters
    char *str_file;
    char *irmsd_file;
        
    Fname *F;
    Paths *PB;
};
typedef struct Parameters Parameters;

Parameters *declare_parameters();
void free_parameters(Parameters *P);
Parameters *default_values(Parameters *P, char *outdir);
Parameters *read_parameters(int argc, char** argv);
Parameters *score_parameters(Parameters *P, int argc, char** argv);
int search_mode(int argc, char** argv, char *mode);
void printhelp();


#endif	/* PARAMETERS_UTILS_H */

