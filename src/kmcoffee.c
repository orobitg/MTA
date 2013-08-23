
#include "kmcoffee.h"

/* 
 * Implementar el sistema de tall
 * Implementar una estructura de cluster (nseq, nom fitxers, ordre tipo arbre?) - Sa de veure
*/

CL_Node declare_cluster_node(int n){
    
    CL_Node cluster;
    
    cluster= (CL_Node) vcalloc(1, sizeof(ClusterNode));
    
    cluster->seqsfname = (char *) vcalloc(ALLPATH, sizeof(char));
    cluster->alnfname = (char *) vcalloc(ALLPATH, sizeof(char));
    cluster->cl_id = n;
    cluster->nseq = 0;
    cluster->isprofile = 0;
    cluster->iscluster = 0;
    cluster->isroot = 0;
    cluster->parent = NULL;
    cluster->left = NULL;
    cluster->right = NULL;
    
    return cluster;
    
}

CL_Node free_cluster_node(CL_Node cluster){
    
    vfree(cluster->seqsfname);
    vfree(cluster->alnfname);
    cluster->parent = NULL;
    cluster->left = NULL;
    cluster->right = NULL;
    
    vfree(cluster);
        
    return NULL;
}

int kmcoffee(Parameters *P){
    
    Sequence *S;
    Distance_matrix *DM;
    NT_node **T;
    int nclusters=0;
    //int cutlevel=0, other=0;
    char *treename;
    char *tree_id;

    S=read_fasta_sequences(P->seqfile);
    if(P->k > S->nseq || P->k == 0){
        P->k = S->nseq / 2;
        fprintf(stderr, "\tERROR: The number of clusters (k) cannot be bigger than th number of sequences or 0. New k to %d\n", P->k);
    }
    if(P->n == 0){
        P->n = S->nseq;
    }

    if(S==NULL){
        fprintf(stderr, "Error - No Sequences\n");
        return(EXIT_FAILURE);
    }
           
    treename = (char *) vcalloc(ALLPATH, sizeof(char));
    tree_id = (char *) vcalloc(S->nseq * 10, sizeof(char));
    sprintf(treename, "%s%s.dnd", P->outdir, (P->F)->name);
    
    DM = make_distance_matrix(S, "ktup");

    T = make_tree(DM, S, treename, "nj_bgt", 0, 0, NULL);
    T[3][0]->level = 0;
    tree_id = identify_tree(T[3][0], tree_id, S->nseq);
    
    if(P->align_method == NULL){
        sprintf(P->align_method, "proba_pair");
    }
    sprintf(P->align_method, "proba_pair");
    //printf("%f\n", log2(P->k));
    fprintf(stdout, "\n---> Identifying the balanced guide tree clusters\n");
    //cutlevel = log2(P->k);
    //other = (P->k) - (cutlevel * cutlevel); 
    //printf("P->k: %d - cutlevel: %d - %d\n", P->k+1, cutlevel, other);
    identify_clusters2(T, S, P, &nclusters);
    fprintf(stdout, "\tFinal number of clusters: %d\n", nclusters);
    fprintf(stdout, "---> DONE\n");
    fprintf(stdout, "\n---> Generating and aligning the balanced guide tree clusters\n");     
    align_profiles_main(T[3][0], S, P, nclusters);

    fprintf(stdout, "---> DONE\n");
    
    vfree(treename);
    vfree(tree_id);

    free_distance_matrix(DM);
    free_sequence(S);  
}


void align_clusters(NT_node root, Sequence *S, Parameters *P, CL_Node **C, int cutlevel, int *other, int *nclusters){
    
    NT_node LL, RL;
    
    if(root->level == cutlevel){
        
        if((*other) > 0 && root->nseq > 2){ 
            //printf("\t\tNEW CLUSTER - %d -- %d\n", root->level, *other);
            C[0][*nclusters] = build_cluster_file(root->left, S, P, *nclusters);
            (*nclusters)++;
            //printf("\t\tNEW CLUSTER - %d -- %d\n", root->level, *other);
            C[0][*nclusters] = build_cluster_file(root->right, S, P, *nclusters); 
            (*nclusters)++;
            (*other)--;  
            
            
/*
            root->visited = 2;
            (*nclusters)++;
*/
            //C[nclusters]
            //alinear profiles
            
        }
        else {
            //printf("\t\tNEW CLUSTER - %d -- %d\n", root->level, *other);
            C[0][*nclusters] = build_cluster_file(root, S, P,  *nclusters);
            (*nclusters)++;
/*
            if((root->parent)->visited == 0){
                (root->parent)->visited = 1;
            }
            else {
                (*nclusters)++;
                (root->parent)->visited += 1;
                //alinear profiles
            }
*/
        }
    }
    
    LL = root->left;
    RL = root->right;
    
    if(!LL->isseq && !RL->isseq){
        //printf("Left Node - %d\n", LL->level);
        align_clusters(LL, S, P, C, cutlevel, other, nclusters);
        //printf("Right Node - %d\n", RL->level);
        align_clusters(RL, S, P, C, cutlevel, other, nclusters);
    }
    else if(!LL->isseq && RL->isseq){
        //printf("Left Node - %d\n", LL->level);
        align_clusters(LL, S, P, C, cutlevel, other, nclusters);
        //printf("\tRight Leaf %d - %s - %d\n", RL->seq_id, RL->name, RL->level);
        if(RL->level == cutlevel || root->level < cutlevel-1){
         
        //    printf("\t\tNEW CLUSTER - %d\n", RL->level);
            
            C[0][*nclusters] = build_cluster_file(RL, S, P, *nclusters);
            (*nclusters)++;
/*
            if(root->visited == 0){
                root->visited = 1;
            }
            else {
                (*nclusters)++;
                root->visited += 1;
                //alinear profiles
            }
*/
        }   
    }
    else if(LL->isseq && !RL->isseq){
        //printf("\tLeft Leaf %d - %s - %d\n", LL->seq_id, LL->name, LL->level);
        if(LL->level == cutlevel || root->level < cutlevel-1){
        //        printf("\t\tNEW CLUSTER - %d\n", LL->level);
           
           C[0][*nclusters] = build_cluster_file(LL, S, P,*nclusters);
            (*nclusters)++;
/*
            if(root->visited == 0){
                root->visited = 1;
            }
            else {
                (*nclusters)++;
                root->visited += 1;
                //alinear profiles
            }
*/
                
                
        }   
        //printf("Right Node - %d\n", RL->level);
        align_clusters(RL, S, P, C, cutlevel, other, nclusters);
    }
    else if(LL->isseq && RL->isseq){   
        //printf("\tLeft Leaf %d - %s - %d\n", LL->seq_id, LL->name, LL->level);
        if(LL->level == cutlevel){
        //        printf("\t\tNEW CLUSTER - %d\n", LL->level);
                
                C[0][*nclusters] = build_cluster_file(LL, S, P, *nclusters);
                (*nclusters)++;
        }  
        //printf("\tRight Leaf %d - %s - %d\n", RL->seq_id, RL->name, RL->level);
        if(RL->level == cutlevel){
        //        printf("\t\tNEW CLUSTER - %d\n", RL->level);
                C[0][*nclusters] = build_cluster_file(RL, S, P, *nclusters);
                (*nclusters)++;
        }   
        
        if((LL->parent)->level < cutlevel-1){
            //printf("\n\nCluster per dalt\n");
        //    printf("\t\tNEW CLUSTER - %d -- %d\n", root->level, *other);
            C[0][*nclusters] = build_cluster_file(LL->parent, S, P, *nclusters);
            (*nclusters)++;
            
        }
/*
        root->visited = 2;
        (*nclusters)++;
*/
        //alinear profiles
    }  
    
    
    
    
    
    
}

void identify_clusters(NT_node root, Sequence *S, Parameters *P, int cutlevel, int *other, int *nclusters){
    
    NT_node LL, RL;
    LL = root->left;
    RL = root->right;
       
    if(!LL->isseq && !RL->isseq){
        //printf("Left Node - %d\n", LL->level);
        identify_clusters(LL, S, P, cutlevel, other, nclusters);
        //printf("Right Node - %d\n", RL->level);
        identify_clusters(RL, S, P, cutlevel, other, nclusters);
    }
    else if(!LL->isseq && RL->isseq){
        //printf("Left Node - %d\n", LL->level);
        identify_clusters(LL, S, P, cutlevel, other, nclusters);
        //printf("\tRight Leaf %d - %s - %d\n", RL->seq_id, RL->name, RL->level);
        if(RL->level == cutlevel || root->level < cutlevel-1){
         
            //printf("\t\tNEW CLUSTER - %d\n", RL->level);
            (*nclusters)++;
            //build_cluster_file(RL, S, P, C, *nclusters);
            RL->aligned = 1;
            RL->iscluster = 1;

        }   
    }
    else if(LL->isseq && !RL->isseq){
        //printf("\tLeft Leaf %d - %s - %d\n", LL->seq_id, LL->name, LL->level);
        if(LL->level == cutlevel || root->level < cutlevel-1){
                //printf("\t\tNEW CLUSTER - %d\n", LL->level);
            (*nclusters)++;
            //build_cluster_file(LL, S, P, C, *nclusters);  
            LL->aligned = 1;
            LL->iscluster = 1;
                
        }   
        //printf("Right Node - %d\n", RL->level);
        identify_clusters(RL, S, P, cutlevel, other, nclusters);
    }
    else if(LL->isseq && RL->isseq){   
        //printf("\tLeft Leaf %d - %s - %d\n", LL->seq_id, LL->name, LL->level);
        if(LL->level == cutlevel){
                //printf("\t\tNEW CLUSTER - %d\n", LL->level);
                (*nclusters)++;
                //build_cluster_file(LL, S, P, C, *nclusters);
                LL->aligned = 1;
                LL->iscluster = 1;
        }  
        //printf("\tRight Leaf %d - %s - %d\n", RL->seq_id, RL->name, RL->level);
        if(RL->level == cutlevel){
                //printf("\t\tNEW CLUSTER - %d\n", RL->level);
                (*nclusters)++;
                //build_cluster_file(RL, S, P, C, *nclusters);
                RL->aligned = 1;
                RL->iscluster = 1;
        }   
        
        if((LL->parent)->level < cutlevel-1){
            //printf("\n\nCluster per dalt\n");
            (*nclusters)++;
            //build_cluster_file(LL->parent, S, P, C, *nclusters);
            LL->parent->aligned = 1;
            LL->parent->iscluster = 1;
            
        }
    }
    
    if(root->level == cutlevel){
        //printf("\t\tNEW CLUSTER - %d -- %d\n", root->level, *other);
        if((*other) > 0 && root->nseq > 2){
            (*nclusters)++;
            //build_cluster_file(root->left, S, P, C, *nclusters);
            root->left->aligned = 1;
            root->left->iscluster = 1;
            (*nclusters)++;
            //build_cluster_file(root->right, S, P *nclusters); 
            root->right->aligned = 1;
            root->right->iscluster = 1;
            (*other)--;  

            
        }
        else {
/*
            if(root->nseq > P->n){
                NT_node RL2, LL2, tmp;
                RL2 = root->right;
                while(RL2->nseq > P->n){
                    tmp = RL2;
                }
                
            }
*/
            (*nclusters)++;
            //build_cluster_file(root, S, P, C, *nclusters);
            root->aligned = 1;
            root->iscluster = 1;
        }
    }
    
}
void print_cluster(NT_node root);

void identify_clusters2(NT_node **T, Sequence *S, Parameters *P, int *nclusters){
    
    NT_node *list_clusters, *tmp_clusters;
    NT_node root, new_node, node;
    int i=0, j=0, k=0;
    int tmp_nclusters = 0;
    int found=0, new=0, mod=0, where=0;
    int seqincl=0, dif=0;
    
    root = T[3][0];
    
    list_clusters = (NT_node *) vcalloc(S->nseq+1, sizeof(NT_node));
    tmp_clusters = (NT_node *) vcalloc(S->nseq+1, sizeof(NT_node));
    
    //printf("Nnodes: %d\n", root->nnodes);
    for(i=0, (*nclusters)=0; i<root->nnodes; i++){      
        if(T[0][i]->isseq){
            //printf("\t- %d - %s\n", T[0][i]->isseq, T[0][i]->name);
            list_clusters[*nclusters] = T[0][i];
            tmp_clusters[*nclusters] = NULL;
            (*nclusters)++;
        }  
    }
   
    tmp_nclusters = *nclusters;
    while((*nclusters) > P->k){
        
        for(i=0; i<(*nclusters); i++){
            found = 0;
            new=0;
            if(tmp_nclusters > P->k){
                new=1;
                if(list_clusters[i]->parent != NULL){
                    node = list_clusters[i]->parent;
                }
                else {
                    node = list_clusters[i];
                }
            }
            else {
                node = list_clusters[i];
            }
           
            for(j=0; j<k; j++){
                if(tmp_clusters[j] == node){
                    found=1;
                    //printf("\t\tRepeated\n");
                    break;
                }
            }
            if(found == 0){
                tmp_clusters[k] = node;
                k++;
                if(new){
                       tmp_nclusters--;
                   //     printf("\t\tNew Node: %d\n", node->nseq);
                }
                
                //printf("\tNO FOUND - TMP-CLUSTERS: %d - %d\n", tmp_nclusters, k);
            }       
        } 
        //tmp_nclusters = k;        
        *nclusters = k;
        tmp_nclusters = k;
        k=0;
        for(i=0; i<S->nseq; i++){
            if(list_clusters[i]!=NULL){
                list_clusters[i]->isorcluster=0;
            }
            if(i < (*nclusters)){
                list_clusters[i] = tmp_clusters[i];  
                list_clusters[i]->isorcluster=1;
                tmp_clusters[i] = NULL;
            }
            else {
                list_clusters[i] = NULL;
                tmp_clusters[i] = NULL;
            }
        }
        k=0;
        NT_node tmp;
        for(i=0; i<(*nclusters)-1; i++){
            for(j=i+1; j<(*nclusters); j++){
                if(list_clusters[j]->nseq < list_clusters[i]->nseq){
                    tmp = list_clusters[j];
                    list_clusters[j] = list_clusters[i];
                    list_clusters[i] = tmp;
                }
            }
        }
/*
        printf("\n");
        for(i=0; i<(*nclusters); i++){
            printf("\tNode: %d - %d - %d - %d\n", list_clusters[i]->index, list_clusters[i]->isseq, list_clusters[i]->node, list_clusters[i]->nseq);
        }
        printf("\n");
*/  
    }
    
   // print_cluster(root);
/*
    if(root->isorcluster){
        printf("root is cluster: %d\n", root->nseq);
    }
*/
    for(i=0; i<(*nclusters); i++){  
            new_node = list_clusters[i];
            new_node->visited = 1;
            correct_cluster(list_clusters[i], &new_node, &mod, &where);
            if(mod!=0){
                list_clusters[i]->isorcluster = 0;
                list_clusters[i]->visited = 0;
            }        
            new_node->iscluster=1;
            new_node->isorcluster=1;
            new_node->aligned=1;
            new_node->visited=0;
            list_clusters[i] = new_node;
            seqincl += list_clusters[i]->nseq;
            mod = 0;     
            where = 0;
    }
/*
    for(i=0; i<(*nclusters); i++){
        printf("\tNode: %d - %d - %d - %d\n", list_clusters[i]->index, list_clusters[i]->isseq, list_clusters[i]->node, list_clusters[i]->nseq);
    }
*/
  //  print_cluster(root);
    if(seqincl < S->nseq){
        dif = S->nseq - seqincl;
        tmp_nclusters = *nclusters;
        loss_cl(root, list_clusters, &tmp_nclusters);
        *nclusters = tmp_nclusters;
    }
/*
    for(i=0; i<(*nclusters); i++){
        printf("\tNode: %d - %d - %d - %d\n", list_clusters[i]->index, list_clusters[i]->isseq, list_clusters[i]->node, list_clusters[i]->nseq);
    }
 */
    //print_cluster(root);
    vfree(list_clusters);
    vfree(tmp_clusters);
}

void correct_cluster(NT_node node, NT_node *new_node, int *mod, int *where){
    
    NT_node RL, LL;
    
    LL = node->left;
    RL = node->right;
   
   //printf("Node seqs: %d LL leaf: %d RL leaf: %d - LL cl: %d RL cl: %d - where: %d\n", (*new_node)->nseq, LL->isseq, RL->isseq, LL->isorcluster, RL->isorcluster, *where);
    if(!LL->isseq && !RL->isseq){
        if(*mod == 0){
            if(!LL->isorcluster && !RL->isorcluster){
                if(*where == 0 || node->visited == 1){
                    *where = 1;
                }
                correct_cluster(LL, new_node, mod, where);
                if(*where == 0 || node->visited == 1){
                    *where = 2;
                }
                correct_cluster(RL, new_node, mod, where);
            }
        }
           
        if(*mod == 0){
            if(LL->isorcluster){
/*
                (or_node->right)->iscluster = 1;
                (or_node->right)->aligned = 1;
*/
                (*new_node)->visited = 0;
              if(*where == 1){
                    *new_node = (*new_node)->right;
                }
                else if(*where == 2){
                    *new_node = (*new_node)->left;
                }
                else if(*where == 0){

                   *new_node = (*new_node)->right;   
                   *where = 2;
                }
                  
                if(!((*new_node)->isseq)){
                    *where = 0;
                    (*new_node)->visited=1;
                    correct_cluster(*new_node, new_node, mod, where);     
                }
                *mod = 2;
            }
/*
            else {
                if(*where == 0){
                    *where = 1;
                }
                correct_cluster(LL, new_node, mod, where);
            }
*/          
        }
        
        if(*mod == 0){
            if(RL->isorcluster){
                (*new_node)->visited = 0;
                if(*where == 1){
                    *new_node = (*new_node)->right;
                }
                else if(*where == 2){
                    *new_node = (*new_node)->left;

                }
 
                else if(*where == 0){
                   *new_node = (*new_node)->left; 
                   *where = 1;
                }
                  
                if(!((*new_node)->isseq)){
                    *where = 0;
                    (*new_node)->visited=1;
                    correct_cluster(*new_node, new_node, mod, where); 
                }
                *mod = 1;
            }
/*
            else{
                if(*where == 0){
                    *where = 2;
                }
                correct_cluster(RL, new_node, mod, where);
            }
*/
        }
    }
    else if(LL->isseq && !RL->isseq){
        if(*mod == 0){
            if(!LL->isorcluster && !RL->isorcluster){
                if(*where == 0 || node->visited == 1){
                    *where = 2;
                }
                correct_cluster(RL, new_node, mod, where);
            }

            if(LL->isorcluster){
                (*new_node)->visited = 0;
                if(*where == 1){
                    *new_node = (*new_node)->right;
                }
                else if(*where == 2){
                    *new_node = (*new_node)->left;
                }
                else if(*where == 0){
                   *new_node = (*new_node)->right; 
                   *where = 2;
                }
                  
                if(!((*new_node)->isseq)){
                    *where = 0;
                    (*new_node)->visited=1;
                    correct_cluster(*new_node, new_node, mod, where);   
                }
                *mod = 2;

            }
            if(RL->isorcluster){
                (*new_node)->visited = 0;
                if(*where == 1){
                    *new_node = (*new_node)->right;
                }
                else if(*where == 2){
                    *new_node = (*new_node)->left;
                }
                else if(*where == 0){
                   *new_node = (*new_node)->left;
                   *where = 1;
                }
                  
                if(!((*new_node)->isseq)){
                    *where = 0;
                    (*new_node)->visited=1;
                    correct_cluster(*new_node, new_node, mod, where); 
                }
                *mod = 1;
            }
        }     
    }
    else if(!LL->isseq && RL->isseq){
        if(*mod == 0){
            if(!LL->isorcluster && !RL->isorcluster){
                if(*where == 0 || node->visited == 1){
                    *where = 1;
                }
                correct_cluster(LL, new_node, mod, where);
            }

            if(LL->isorcluster){
                (*new_node)->visited = 0;
                if(*where == 1){
                    *new_node = (*new_node)->right;
                }
                else if(*where == 2){
                    *new_node = (*new_node)->left;
                }
                else if(*where == 0){
                   *new_node = (*new_node)->right;
                   *where = 2;
                }
                  
                if(!((*new_node)->isseq)){
                    *where = 0;
                    (*new_node)->visited=1;
                    correct_cluster(*new_node, new_node, mod, where);               
                }
                *mod = 2;
            }
            if(RL->isorcluster){
                (*new_node)->visited = 0;
               if(*where == 1){
                    *new_node = (*new_node)->right;
                }
                else if(*where == 2){
                    *new_node = (*new_node)->left;
                }
                else if(*where == 0){
                   *new_node = (*new_node)->left;
                   *where = 1;
                }
                  
                if(!((*new_node)->isseq)){
                    *where = 0;
                    (*new_node)->visited=1;
                    correct_cluster(*new_node, new_node, mod, where);               
                }
                *mod = 1;
            }
        }
    }
    else if(LL->isseq && RL->isseq){
        if(*mod == 0){
            if(LL->isorcluster){
                (*new_node)->visited = 0;
                if(*where == 1){
                    *new_node = (*new_node)->right;
                }
                else if(*where == 2){
                    *new_node = (*new_node)->left;
                }
                else if(*where == 0){
                   *new_node = (*new_node)->right;
                   *where = 2;
                }
                  
                if(!((*new_node)->isseq)){
                    *where = 0;
                    (*new_node)->visited=1;
                    correct_cluster(*new_node, new_node, mod, where);
                }               
                *mod = 2;
            }
            else if(RL->isorcluster){
                (*new_node)->visited = 0;
                if(*where == 1){
                    *new_node = (*new_node)->right;
                }
                else if(*where == 2){
                    *new_node = (*new_node)->left;
                }
                else if(*where == 0){
                   *new_node = (*new_node)->left;
                   *where = 1;
                }
                  
                if(!((*new_node)->isseq)){
                    *where = 0;
                    (*new_node)->visited=1;
                    correct_cluster(*new_node, new_node, mod, where);   
                }
                *mod = 1;
            }
        }
    }
    
}

void loss_cl(NT_node root, NT_node *list_node, int *tmp_nclusters){
    
    NT_node RL, LL;
    
    LL = root->left;
    RL = root->right;
    //printf("i: %d - nseqs: %d\n", *tmp_nclusters, root->nseq);
    if(!LL->isseq && !RL->isseq){        
        if(!LL->iscluster){
            loss_cl(LL, list_node, tmp_nclusters);
        }      
        if(!RL->iscluster){
            loss_cl(RL, list_node, tmp_nclusters);
        }
        
    }
    else if(LL->isseq && !RL->isseq){
        if(!LL->iscluster){
            LL->iscluster = 1;
            LL->aligned = 1;
            LL->isorcluster = 1;
            list_node[*tmp_nclusters] = LL;
            (*tmp_nclusters)++;
            //printf("lost_cl: %d\n", *tmp_nclusters);
        }
        if(!RL->iscluster){
            loss_cl(RL, list_node, tmp_nclusters);
        }
        
    }
    else if(!LL->isseq && RL->isseq){
        if(!LL->iscluster){
            loss_cl(LL, list_node, tmp_nclusters);
        }
        if(!RL->iscluster){
            RL->iscluster = 1;
            RL->aligned = 1;
            RL->isorcluster = 1;
            list_node[*tmp_nclusters] = RL;
            (*tmp_nclusters)++;    
            //printf("lost_cl: %d\n", *tmp_nclusters);
        }
        
    }
    else if(LL->isseq && RL->isseq){
        if(!LL->iscluster && !RL->iscluster) {
            root->iscluster = 1;
            root->aligned = 1;
            root->isorcluster = 1;
            list_node[*tmp_nclusters] = root;
            (*tmp_nclusters)++;    
            //printf("lost_cl: %d\n", *tmp_nclusters);
        }
        else if(!LL->iscluster){
            LL->iscluster = 1;
            LL->aligned = 1;
            LL->isorcluster = 1;
            list_node[*tmp_nclusters] = LL;
            (*tmp_nclusters)++;    
            //printf("lost_cl: %d\n", *tmp_nclusters);
        }
        else if(!RL->iscluster){
            RL->iscluster = 1;
            RL->aligned = 1;
            RL->isorcluster = 1;
            list_node[*tmp_nclusters] = RL;
           (*tmp_nclusters)++; 
            //printf("lost_cl: %d\n", *tmp_nclusters);
        }
        
    }
    
}

//tmps al temporal en futures versions

void align_profiles_main(NT_node root, Sequence *S, Parameters *P, int nclusters){
    
    CL_Node **C;
    FILE *prf_F, *st;
    char *profilename, *command, *st_name, *tmp_aln_name;
    int cl_id=0, i=0, nprofiles=0;
    
    tmp_aln_name = (char *) vcalloc(ALLPATH, sizeof(char));
    
    C = (CL_Node **) vcalloc(1, sizeof(CL_Node *));
    C[0] = (CL_Node *) vcalloc(((P->k * 2) +1), sizeof(CL_Node));  
    
    if(strm(P->cl_mode, "rootcl")){
        fprintf(stdout, "\tAlign profiles in rootmode\n");
        profilename = (char *) vcalloc(ALLPATH, sizeof(char));
        command = (char *) vcalloc(ALLPATH, sizeof(char));
        
        align_profiles_root(root, S, P, C, &cl_id);
        
        fprintf(stdout, "\tFinal number of clusters: %d\n", cl_id);
        sprintf(profilename, "%s/tmp_%s_prf%d.fa", TMP, P->F->name, P->k);
        prf_F = openfile(profilename, "w");
        for(i=0; i<cl_id; i++){		
            printf("\t\t\t%s\n", C[0][i]->alnfname);
            fprintf(prf_F, "%s\n", C[0][i]->alnfname);
        }
        fclose(prf_F);
        fprintf(stdout, "---> DONE\n");
        fprintf(stdout, "\n---> Aligning the profiles\n");
        sprintf(command, "%s -output fasta_aln -quiet -outfile %s%s.fasta_aln -method %s -profile FILE::%s >/dev/null 2>/dev/null", (P->PB)->tcoffee_bin, P->outdir, P->F->name, P->align_method, profilename);
        if (system(command))
                printf("%s\n",command);
         
        nprofiles=1;

        remove(profilename);
        
        sprintf(tmp_aln_name, "rm ./tmp_aln_%s_*.dnd", P->F->name);
        i = system(tmp_aln_name);
        
        vfree(command);
        vfree(profilename);
    
    }
    else if(strm(P->cl_mode, "treecl")){
        fprintf(stdout, "\tAlign profiles in treemode\n");
        align_profiles_tree(root, S, P, C, &cl_id);
        nprofiles = cl_id - nclusters;        
    }
    else {
        fprintf(stderr, "ERROR: Wrong cluster alignment mode\n");       
        exit(-1);
    }
    
    if(P->stats){
        st_name = (char *) vcalloc(ALLPATH, sizeof(char));
        sprintf(st_name, "%s%s_statistics.txt", P->outdir, P->F->name);
        st = openfile(st_name, "w");
        fprintf(st, "%s;%d;%d\n", P->F->name, nclusters, nprofiles);
        for(i=0; i<cl_id; i++){
            if(C[0][i]->iscluster){
                fprintf(st, "%d;%d\n", i+1, C[0][i]->nseq);
            }
        }
        fclose(st);
        vfree(st_name);
    }
    
    
    sprintf(tmp_aln_name, "rm ./tmp_%s_*.dnd", P->F->name);
    i = system(tmp_aln_name);



    vfree(tmp_aln_name);
    for(i=0; i<cl_id; i++){
        if(C[0][i]->isroot ==0){
            remove(C[0][i]->alnfname);
        }
        //if(C[0][i]->nseq > 1){
        remove(C[0][i]->seqsfname);                 
    }

  
    for(i=1; i<nclusters; i++){
        C[0][i] = free_cluster_node(C[0][i]);
    }
    vfree(C);
    
}

void align_profiles_root(NT_node root, Sequence *S, Parameters *P, CL_Node **C, int *cl_id){
    NT_node LL, RL;
    LL = root->left;
    RL = root->right;
    
    if(!LL->isseq && !RL->isseq){
        //printf("Left node %d\n", LL->aligned);
        align_profiles_root(LL, S, P, C, cl_id);
        //printf("Right node %d\n", RL->aligned);
        align_profiles_root(RL, S, P, C, cl_id);
    }
    else if(!LL->isseq && RL->isseq){
        //printf("Left node %d\n", LL->aligned);
        align_profiles_root(LL, S, P, C, cl_id);
        //printf("Right Leaf %s - %d\n", RL->name, RL->aligned);
        if(RL->aligned == 1){
            if(RL->iscluster){
                C[0][*cl_id] = build_cluster_file(RL, S, P, *cl_id);
                RL->cl_id = *cl_id;
                (*cl_id)++;
            }
        }
    }
    else if(LL->isseq && !RL->isseq){
        //printf("Left Leaf %s - %d\n", LL->name, LL->aligned);
        if(LL->aligned == 1){
            //printf("\tCluster Leaf\n");
            if(LL->iscluster){
                C[0][*cl_id] = build_cluster_file(LL, S, P, *cl_id);
                LL->cl_id = *cl_id;
                (*cl_id)++;
            }
        }
        align_profiles_root(RL, S, P, C, cl_id);
    }
    else if(LL->isseq && RL->isseq){
       // printf("Left Leaf %s - %d\n", LL->name, LL->aligned);
        if(LL->aligned == 1){
            //printf("\tCluster Leaf\n");
            if(LL->iscluster){
                C[0][*cl_id] = build_cluster_file(LL, S, P, *cl_id);
                LL->cl_id = *cl_id;
                (*cl_id)++;
            }
        }
        if(RL->aligned == 1){
            //printf("\tCluster Leaf\n");         
            if(RL->iscluster){
                C[0][*cl_id] = build_cluster_file(RL, S, P, *cl_id);
                RL->cl_id = *cl_id;
                (*cl_id)++;
            }
            
        }
    }
    
    if(root->aligned){      
        if(root->iscluster){          
            C[0][*cl_id] = build_cluster_file(root, S, P,*cl_id);
            root->cl_id = *cl_id;
            (*cl_id)++;
        }
    }

}

void align_profiles_tree(NT_node root, Sequence *S, Parameters *P, CL_Node **C, int *cl_id){
    
    NT_node LL, RL;
    LL = root->left;
    RL = root->right;
    int i, child;
    
    if(!LL->isseq && !RL->isseq){
        //printf("Left node %d\n", LL->aligned);
        align_profiles_tree(LL, S, P, C, cl_id);
        //printf("Right node %d\n", RL->aligned);
        align_profiles_tree(RL, S, P, C, cl_id);
    }
    else if(!LL->isseq && RL->isseq){
        //printf("Left node %d\n", LL->aligned);
        align_profiles_tree(LL, S, P, C, cl_id);
        //printf("Right Leaf %s - %d\n", RL->name, RL->aligned);
        if(RL->aligned == 1){
            if(RL->iscluster){
                C[0][*cl_id] = build_cluster_file(RL, S, P, *cl_id);
                RL->cl_id = *cl_id;
                (*cl_id)++;
            }
            //printf("\tCluster Leaf\n");
            RL->parent->visited += 1;
            if((RL->parent)->visited == 1){
                (RL->parent)->cl_id=*cl_id;
                //printf("\t\tCHILD ID: %d\n", RL->cl_id);
                C[0][*cl_id] = build_profile_file(RL->parent, S, P, C[0][*cl_id], C[0][RL->cl_id], *cl_id, (RL->parent)->visited);
                (RL->parent)->ant=*cl_id;
                C[0][RL->cl_id]->parent = C[0][*cl_id];
                (*cl_id)++;
            }
            if((RL->parent)->visited == 2){
                //printf("\t\tCHILD ID: %d\n", RL->cl_id);
                i = (RL->parent)->ant;
                C[0][i] = build_profile_file(RL->parent, S, P, C[0][i], C[0][RL->cl_id], i, (RL->parent)->visited); 
                C[0][RL->cl_id]->parent = C[0][*cl_id];
                RL->parent->aligned = 1;
            }
        }
    }
    else if(LL->isseq && !RL->isseq){
        //printf("Left Leaf %s - %d\n", LL->name, LL->aligned);
        if(LL->aligned == 1){
            //printf("\tCluster Leaf\n");
            if(LL->iscluster){
                C[0][*cl_id] = build_cluster_file(LL, S, P, *cl_id);
                LL->cl_id = *cl_id;
                (*cl_id)++;
            }
            LL->parent->visited += 1;
            if((LL->parent)->visited == 1){
                (LL->parent)->cl_id=*cl_id;
                //printf("\t\tCHILD ID: %d\n", LL->cl_id);
                C[0][*cl_id] = build_profile_file(LL->parent, S, P, C[0][*cl_id], C[0][LL->cl_id], *cl_id, (LL->parent)->visited);
                (LL->parent)->ant=*cl_id;
                C[0][LL->cl_id]->parent = C[0][*cl_id];
                (*cl_id)++;
            }
            if((LL->parent)->visited == 2){ 
               // printf("\t\tCHILD ID: %d\n", LL->cl_id);
                i = (LL->parent)->ant;
                C[0][i] = build_profile_file(LL->parent, S, P, C[0][i], C[0][LL->cl_id], i, (LL->parent)->visited); 
                C[0][LL->cl_id]->parent = C[0][*cl_id];
                LL->parent->aligned = 1;
            }
        }
        //printf("Right node %d\n", RL->aligned);
        align_profiles_tree(RL, S, P, C, cl_id); 
    }
    else if(LL->isseq && RL->isseq){
       // printf("Left Leaf %s - %d\n", LL->name, LL->aligned);
        if(LL->aligned == 1){
            //printf("\tCluster Leaf\n");
            if(LL->iscluster){
                C[0][*cl_id] = build_cluster_file(LL, S, P, *cl_id);
                LL->cl_id = *cl_id;
                (*cl_id)++;
            }
            LL->parent->visited += 1;
            if((LL->parent)->visited == 1){
                (LL->parent)->cl_id=*cl_id;
                //printf("\t\tCHILD ID: %d\n", LL->cl_id);
                C[0][*cl_id] = build_profile_file(LL->parent, S, P, C[0][*cl_id], C[0][LL->cl_id], *cl_id, (LL->parent)->visited);            
                (LL->parent)->ant=*cl_id;
                 C[0][LL->cl_id]->parent = C[0][*cl_id];
                (*cl_id)++;
            }
            if((LL->parent)->visited == 2){
                //printf("\t\tCHILD ID: %d\n", LL->cl_id);
                i = (LL->parent)->ant;
                C[0][i] = build_profile_file(LL->parent, S, P, C[0][i], C[0][LL->cl_id], i, (LL->parent)->visited);    
                C[0][LL->cl_id]->parent = C[0][*cl_id];
                LL->parent->aligned = 1;
            }
        }
        //printf("Right Leaf %s - %d\n", RL->name, RL->aligned);
        if(RL->aligned == 1){
            //printf("\tCluster Leaf\n");         
            if(RL->iscluster){
                C[0][*cl_id] = build_cluster_file(RL, S, P, *cl_id);
                RL->cl_id = *cl_id;
                (*cl_id)++;
            }
            RL->parent->visited += 1;
            if((RL->parent)->visited == 1){
                (RL->parent)->cl_id=*cl_id;
                //printf("\t\tCHILD ID: %d\n", RL->cl_id);
                C[0][*cl_id] = build_profile_file(RL->parent, S, P, C[0][*cl_id], C[0][RL->cl_id], *cl_id, (RL->parent)->visited);
                (RL->parent)->ant=*cl_id;
                C[0][RL->cl_id]->parent = C[0][*cl_id];
                (*cl_id)++;
            }
            if((RL->parent)->visited == 2){   
                //printf("\t\tCHILD ID: %d\n", RL->cl_id);
                i = (RL->parent)->ant;
                C[0][i] = build_profile_file(RL->parent, S, P, C[0][i], C[0][RL->cl_id], i, (RL->parent)->visited);
                C[0][RL->cl_id]->parent = C[0][*cl_id];
                RL->parent->aligned = 1;
            }
        }
    }
    if(root->aligned){
       
        if(root->iscluster){          
            C[0][*cl_id] = build_cluster_file(root, S, P, *cl_id);
            root->cl_id = *cl_id;
            (*cl_id)++;
            
        }
        //printf("\tCluster Node\n");
        if(root->parent!=NULL){
            root->parent->visited += 1;
            if((root->parent)->visited == 1){
                root->parent->cl_id = *cl_id;
                //printf("\t\tCHILD ID: %d\n", root->cl_id);             
                C[0][*cl_id] = build_profile_file(root->parent, S, P, C[0][*cl_id], C[0][root->cl_id], *cl_id, (root->parent)->visited);
                (root->parent)->ant = *cl_id;
                C[0][root->cl_id]->parent = C[0][*cl_id];
                (*cl_id)++;
            }
            else if((root->parent)->visited == 2){
                //printf("\t\tCHILD ID: %d\n", root->cl_id);
                i = (root->parent)->ant;
                C[0][i] = build_profile_file(root->parent, S, P, C[0][i], C[0][root->cl_id], i, (root->parent)->visited);
                C[0][root->cl_id]->parent = C[0][*cl_id];
                root->parent->aligned = 1;
            }
        }
    }
}
    
CL_Node build_cluster_file(NT_node cl_root, Sequence *S, Parameters *P, int cl_id){
    
    FILE *fp;
    char command[ALLPATH];
    
    CL_Node C;
    
    C = (CL_Node) declare_cluster_node(cl_id);
    
    sprintf(C->seqsfname, "%s/tmp_%s_%d.fa", TMP, P->F->name, C->cl_id);
    sprintf(C->alnfname, "%s/tmp_aln_%s_%d.fa", TMP, P->F->name, C->cl_id);
    //printf("\t\t%d - %d - %s - %s\n", C->cl_id, C->nseq, C->seqsfname, C->alnfname);
    C->iscluster = 1;
    C->isprofile = 0;
    C->isroot = 0;
    if(!(cl_root->isseq)){   
        C->nseq = cl_root->nseq;
        fp = fopen(C->seqsfname, "w");
        build_cluster_file_runtree(cl_root, S, fp);  
        fclose(fp);
        sprintf(command, "%s -in %s -output fasta_aln -outfile %s -method %s -quiet >/dev/null 2>/dev/null", (P->PB)->tcoffee_bin, C->seqsfname, C->alnfname, P->align_method);
        
        if (system(command))
                printf("%s\n",command);
    }
    else {
        C->nseq = 1;
        fp = fopen(C->alnfname, "w");
        build_cluster_file_runtree(cl_root, S, fp);  
        fclose(fp);     
    }
    //printf("\t\tNseq: %d\n", C->nseq); 
    return C;
    //cl_root->aligned = 1; 
}

CL_Node build_profile_file(NT_node cl_root, Sequence *S, Parameters *P, CL_Node C, CL_Node child, int cl_id, int type){
    
    FILE *prf_F;
    char command[ALLPATH];
    
    //printf("\t\tBuild Profile - %d\n", cl_id);
 
    if(type == 1){
        C = (CL_Node) declare_cluster_node(cl_id);
        C->iscluster = 0;
        C->isprofile = 1;
        C->nseq = cl_root->nseq;
        if(cl_root->parent != NULL){
            C->isroot = 0;
            sprintf(C->seqsfname, "%s/tmp_prf_%s_%d.fa", TMP, P->F->name, C->cl_id);
            sprintf(C->alnfname, "%s/tmp_prfaln_%s_%d.fa", TMP, P->F->name, C->cl_id);
        }
        else {
            C->isroot = 1;
            sprintf(C->seqsfname, "%s/tmp_prf_%s_%d.fa", TMP, P->F->name, C->cl_id);
            sprintf(C->alnfname, "%s%s.fasta_aln", P->outdir, P->F->name);
        }
        
        C->left = child;      
        //printf("\t\tLeft child: %d\n", child->cl_id);
    }
    else if(type == 2){
        C->right = child;
        //printf("\t\tLeft Child: %d -- Rigth child: %d\n", C->left->cl_id, C->right->cl_id);
        //printf("\t\tAlign\n");


        prf_F = openfile(C->seqsfname, "w");   	
        fprintf(prf_F, "%s\n", (C->left)->alnfname);
        fprintf(prf_F, "%s\n", (C->right)->alnfname);  
        fclose(prf_F); 
        //printf("P->align: %s\n", P->align_method);

        sprintf(command, "%s -output fasta_aln -quiet -outfile %s -method %s -profile FILE::%s >/dev/null 2>/dev/null", (P->PB)->tcoffee_bin, C->alnfname, P->align_method, C->seqsfname);
        //printf("\t%s%s -output fasta_aln -quiet -outfile %s -method %s -profile FILE::%s\n", MT_HOME, TCOFFEE_BIN, C->alnfname, P->align_method, C->seqsfname);
        //sprintf(command, "cat %s", C->seqsfname);      
        if (system(command))
            printf("%s\n",command); 

    }
   //printf("\t\t%d - %d - %s - %s\n", C->cl_id, C->nseq, C->seqsfname, C->alnfname);
    return C;
}

int build_cluster_file_runtree(NT_node cl_root, Sequence *S, FILE *fp){
    
    NT_node LL, RL;
    int nseq;
    
    if(cl_root->isseq){
        fprintf(fp, ">%s\n%s\n", S->seq_names[cl_root->seq_id], S->seq[cl_root->seq_id]);
        return;
    }
    
    LL = cl_root->left;
    RL = cl_root->right;
    
    if(!LL->isseq && !RL->isseq){
        build_cluster_file_runtree(LL, S, fp);
        build_cluster_file_runtree(RL, S, fp);
    }
    else if(!LL->isseq && RL->isseq){
        build_cluster_file_runtree(LL, S, fp);
        fprintf(fp, ">%s\n%s\n", S->seq_names[RL->seq_id], S->seq[RL->seq_id]);
    }
    else if(LL->isseq && !RL->isseq){
        fprintf(fp, ">%s\n%s\n", S->seq_names[LL->seq_id], S->seq[LL->seq_id]);
        build_cluster_file_runtree(RL, S, fp);
    }
    else if(LL->isseq && RL->isseq){
        fprintf(fp, ">%s\n%s\n", S->seq_names[LL->seq_id], S->seq[LL->seq_id]);
        fprintf(fp, ">%s\n%s\n", S->seq_names[RL->seq_id], S->seq[RL->seq_id]);
    }   
    
}

void print_cluster(NT_node root){
    
    NT_node RL, LL;
    
    LL = root->left;
    RL = root->right;

    if(!LL->isseq && !RL->isseq){
        print_cluster(LL);
        if(LL->isorcluster){
            printf("LL is cluster: %d\n", LL->nseq);
        }    
        
        print_cluster(RL);
        if(RL->isorcluster){
            printf("RL is cluster: %d\n", RL->nseq);
        } 
        
    }
    else if(LL->isseq && !RL->isseq){
        if(LL->isorcluster){
            printf("LL is cluster: %d - %s\n", LL->nseq, LL->name);
        }
        print_cluster(RL);
          if(RL->isorcluster){
            printf("RL is cluster: %d\n", RL->nseq);
        } 
   
        
    }
    else if(!LL->isseq && RL->isseq){
        print_cluster(LL);
        if(LL->isorcluster){
            printf("LL is cluster: %d\n", LL->nseq);
        }    
        if(RL->isorcluster){
            printf("RL is cluster: %d - %s\n", RL->nseq, RL->name);
        }
        
    }
    else if(LL->isseq && RL->isseq){
        if(LL->isorcluster){
            printf("LL is cluster: %d - %s\n", LL->nseq, LL->name);
        }
        if(RL->isorcluster){
            printf("RL is cluster: %d - %s\n", RL->nseq, RL->name);
        }
        
    }
/*
    if(root->isorcluster){
        printf("Node is cluster: %d\n", root->nseq);
    }    
*/
    
}
