/* 
 * File:   kmcoffee.h
 * Author: oro
 *
 * Created on 12 / abril / 2012, 19:16
 */

#ifndef KMCOFFEE_H
#define	KMCOFFEE_H

#include "utils.h"
#include "guide_tree_utils.h"
#include "sequence_utils.h"
#include "distance_matrix_utils.h"
#include "parameters_utils.h"

typedef struct clnode *CL_Node;

typedef struct clnode {
    
    char *alnfname;
    char *seqsfname;
    int cl_id;
    int nseq;
    int isprofile;
    int iscluster;
    int isroot;
    
    CL_Node parent;
    CL_Node left;
    CL_Node right;
    
}ClusterNode;



CL_Node declare_cluster_node(int n);
CL_Node free_cluster_node(CL_Node cluster);

int kmcoffee(Parameters *P);
void align_clusters(NT_node root, Sequence *S, Parameters *P, CL_Node **C, int cutlevel, int *other, int *nclusters);
void identify_clusters(NT_node root, Sequence *S, Parameters *P, int cutlevel, int *other, int *nclusters);

void identify_clusters2(NT_node **T, Sequence *S, Parameters *P, int *nclusters);
void correct_cluster(NT_node node, NT_node *new_node, int *mod, int *where);
void loss_cl(NT_node root, NT_node *list_node, int *i);

void align_profiles_main(NT_node root, Sequence *S, Parameters *P, int nclusters);
void align_profiles_root(NT_node root, Sequence *S, Parameters *P, CL_Node **C, int *cl_id);
void align_profiles_tree(NT_node root, Sequence *S, Parameters *P, CL_Node **C, int *cl_id);
CL_Node build_cluster_file(NT_node cl_root, Sequence *S, Parameters *P, int cl_id);
CL_Node build_profile_file(NT_node cl_root, Sequence *S, Parameters *P, CL_Node C, CL_Node child, int cl_id, int type);
int build_cluster_file_runtree(NT_node cl_root, Sequence *S, FILE *fp);
//void run_tree(NT_node root);

#endif	/* KMCOFFEE_H */

