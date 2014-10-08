
#ifndef _GUIDE_TREE_H
#define	_GUIDE_TREE_H

#include <stdarg.h>
#include <signal.h>

#include "utils.h"
#include "distance_matrix_utils.h"

#define NODE 0
#define LEAF 1
#define LEFT 1
#define RIGHT 2
#define INTERVAL 0.50

/*
 * Guide Tree
 */

typedef struct tnode *NT_node;

/**
* Node of a tree
*/
typedef struct tnode{
    int visited;
    char *name;
    char *file;
    int seq_id;

    ///The parent node
    NT_node parent;
    ///Left child node
    NT_node left;
    ///Right child node
    NT_node right;
    NT_node bot;
    /// is leaf?
    int isseq;
    int seq;
    int maxnseq;
    int nseq;
    int nnodes;
    ///contains a list of the sequences
    int *lseq;
    ///contains a coded version of the node: 10010101
    int *lseq2;
    ///contains distances to the root, in nodes
    int *idist;
    ///contains real distances *1000
    int *ldist;
    float dist;
    float bootstrap;
    float dp;
    int order;
    int aligned;
    ///Number of leave below the considered node
    int leaf;
    ///Number of nodes below the considered node
    int node;
    int group;
    float score;
    int align;
    char *seqal;
    int index;
    int fork;
    int level;
    
    
    int iscluster;
    int isorcluster;
    int ant;
    int cl_id;
    //To recalculate a distance matrix from a tree
    //int node_id;
    //int nnode;
   // int sub_root;
    //int sub_tree;
}Treenode;

FILE * create_tree(NT_node ptree, NT_node parent,int *nseq,int  *ntotal,int  *nnodes,NT_node **lu, FILE *fp);
NT_node declare_tree_node(int nseq);
void set_info(NT_node p, NT_node parent, int pleaf, char *pname, float pdist, float bootstrap);
NT_node insert_tree_node(NT_node pptr);
void create_tree_node(NT_node pptr, NT_node parent);
//NT_node make_root_tree ( Alignment *A,Constraint_list *CL,int gop, int gep,Sequence *S,  char *tree_file,int maximise);
NT_node reroot(NT_node ptree, int nseq, int ntotal, int nnodes, NT_node **lu);
float calc_root_mean(NT_node root, float *maxdist, int nseq, NT_node **lu);
float calc_mean(NT_node nptr, float *maxdist, int nseq,NT_node **lu);
NT_node insert_root(NT_node p, float diff);
int tree_file2nseq (char *fname);

NT_node **make_tree(Distance_matrix *DM, Sequence *S, char *tree_file, char *mode, int tree, int nrand, char *method);
NT_node **make_nj_tree(int **distances, char **out_seq, char **out_seq_name, int *out_seq_id, int out_nseq, char *tree_file, char *mode, int tree, int nrand, char *method);
NT_node **int_dist2nj_tree(int **distances, char **out_seq_name, int *out_seq_id, int out_nseq,  char *tree_file, char *mode, int tree, int nrand, char *method);
NT_node **float_dist2nj_tree(float **distances, char **out_seq_name, int *out_seq_id, int out_nseq,  char *tree_file, char *mode, int tree, int nrand, char *method);
NT_node **dist2nj_tree (double **distances, char **out_seq_name, int *out_seq_id, int out_nseq,  char *tree_file, char *mode, int tree, int nrand, char *method);
void guide_tree(char *fname, double **saga_tmat, char **saga_seq_name, int saga_nseq, char *mode, int tree, int nrand, char *method);
void nj_tree(char **tree_description, int nseq, char *mode, int tree, int nrand);
void print_phylip_tree(char **tree_description, FILE *tree, int bootstrap);
void print_phylip_tree_clo(char **tree_description, FILE *tree, int bootstrap);
void two_way_split(char **tree_description, FILE *tree, int start_row, int flag, int bootstrap);
NT_node** read_tree(char *treefile, int *tot_node,int nseq, char **seq_names, int *out_seq_id, char *method);
void slow_nj_tree_mta(char **tree_description, int tree, int nrand);
void fast_nj_tree_mta(char **tree_description, int tree, int nrand);
void slow_nj_tree_bgt(char **tree_description);
void fast_nj_tree_bgt(char **tree_description);
char* identify_tree(NT_node root, char *tree_id, int nseqs);
int search_tree_id(char** tree_id_list, char* tree_id, int ntreeid);
int count_entries_file(char *fname, int nseq);
char** file2tree_id_list(char *fname, int ntrees, int nseqs, int entries);

char* root_unrooted_tree(char *treename, int ntree, char *retree_bin);


#endif	/* _GUIDE_TREE_HPP */
