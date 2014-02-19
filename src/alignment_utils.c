
#include <stdio.h>
#include <unistd.h>
#include "alignment_utils.h"
#include "utils.h"

extern int my_rank;
extern int np;


/*MTA sequential method with trees repeated*/

#ifdef SEQ_FLAG

int mta_program_repeated_trees(Parameters *P, Sequence *S, Distance_matrix *DM){

    FILE *fp;
    FILE *fptmp;
    int i=0, j=0, k=0, ret, max_ntree=-1;
    double *sc_list;
    double max_score=-99999.99999;
    char *treename, *alnname, *libname=NULL, *bestalignment, *besttree, *scoresfile;
    char *tmpname;
    char **list;
    NT_node **T;
    struct timeval tim;

    fprintf(stdout, "\n---> Generating the Multiples Trees\n");
    
    /*Declare and initialize variables*/
    treename = (char *) vcalloc(FILENAMELEN, sizeof(char));
    alnname = (char *) vcalloc(FILENAMELEN, sizeof(char));
    scoresfile = (char *) vcalloc(FILENAMELEN, sizeof(char));
    bestalignment = (char *) vcalloc(FILENAMELEN, sizeof(char));
    besttree = (char *) vcalloc(FILENAMELEN, sizeof(char));
    sc_list = (double *) vcalloc(2, sizeof(double));
    
    if(P->only_tree == 0){
        if(P->only_aln == 0){
            sprintf(scoresfile, "%s%s.scores", P->outdir, (P->F)->name);
            fp=openfile(scoresfile, "w");
        }

        if(strm(P->align_method, "tcoffee")){
            libname = (char *) vcalloc(FILENAMELEN, sizeof(char));
            sprintf(libname, "%s%s.lib", P->outdir, (P->F)->name);      
        }

        if(strm(P->align_method, "clustalo") || strm(P->align_method, "mafft")){
            tmpname = (char *) vcalloc(FILENAMELEN, sizeof(char));
            sprintf(tmpname, "./outtree");
            fptmp = openfile(tmpname, "w");
            fclose(fptmp);
            vfree(tmpname);
        }
    }
    //gettimeofday(&tim, NULL);
    srand(1985);
    if(P->treelist == 1){
        list = read_tree_list(P->tree_list, P->ntree);
    }
    
    for(i=0; i<P->ntree; i++){ 
        sprintf(treename, "%s%s_%d.dnd", P->outdir, (P->F)->name, i);
        sprintf(alnname, "%s%s_%d.fasta_aln", P->outdir, (P->F)->name, i);
        if(P->treelist == 1){
            sprintf(treename, "%s", list[i]);
        }
        else{
            T = make_tree(DM, S, treename, "nj", i, P->random_percentage, P->align_method); 
        }
        if(P->only_tree == 0){
            align_tree(P, treename, libname, alnname, i);
            if(P->only_aln == 0){
                sc_list[1] = calculate_score(P, alnname, i, sc_list);

                fprintf(fp, "%d;%s_%d;%lf\n", i, (P->F)->name, i, sc_list[1]);
                fprintf(stdout, "\t- %d - Tree: %s_%d.dnd - Aln: %s_%d.fasta_aln  - Score: %lf\n", i, (P->F)->name, i, (P->F)->name, i, sc_list[1]);

                if(sc_list[1] > max_score){
                    max_score = sc_list[1];
                    max_ntree=i;
                    if(i!=0){
                        if(P->treelist == 0){
                                remove_file(besttree);
                                remove_file(bestalignment);
                        }
                    }
                    sprintf(besttree, "%s", treename);
                    sprintf(bestalignment, "%s", alnname);
                }
                else { 
                    if(P->treelist == 0){
                        remove_file(treename);
                        remove_file(alnname);
                    }

                }
            }
        
        }
    }
    
    if(P->only_tree == 0){
        if(strm(P->align_method, "tcoffee")){
            remove_file(libname);
            vfree(libname);
        }
        if(strm(P->align_method, "clustalo") || strm(P->align_method, "mafft")){
            remove(tmpname);
        }
        if(P->only_aln == 0){
            fprintf(stdout, "---> Done\n");
            fprintf(stdout, "\n\t- BEST TREE: %s_%d - Score: %lf\n", (P->F)->name, max_ntree, max_score);
            if(P->treelist == 0){
                sprintf(treename, "%s%s.dnd", P->outdir, (P->F)->name);
                rename(besttree, treename);
                sprintf(alnname, "%s%s.fasta_aln", P->outdir, (P->F)->name);
                rename(bestalignment, alnname);
            }

            fprintf(fp, "Best;%s_%d;%lf\n", (P->F)->name, max_ntree, max_score);
            fclose(fp);
        }
    }
    
    
    vfree(scoresfile);
    vfree(alnname);
    vfree(treename);
    vfree(bestalignment);
    vfree(besttree);
    vfree(sc_list);
    vfree(list);

    return max_score;
}

/*MTA sequential method with no trees repeated*/

int mta_program_no_repeated_trees(Parameters *P, Sequence *S, Distance_matrix *DM){

    FILE *fp;
    FILE *fptmp;
    int i=0, j=0, k=0, ntreeid=0, ret=0, entries=0, max_ntree=-1;   
    int tree_id_lenght;
    double max_score=-99999.99999;
    double *sc_list;
    char *tree_id;
    char *treename, *alnname, *libname=NULL, *scoresfile, *bestalignment, *besttree;
    char *tmpname;
    char **tree_id_list;
    char **list;
    NT_node **T;

    struct timeval tim;
    
    fprintf(stdout, "\n---> Generating the Multiples Trees\n");

    treename = (char *) vcalloc(FILENAMELEN, sizeof(char));
    alnname = (char *) vcalloc(FILENAMELEN, sizeof(char));
    scoresfile = (char *) vcalloc(FILENAMELEN, sizeof(char));
    bestalignment = (char *) vcalloc(FILENAMELEN, sizeof(char));
    besttree = (char *) vcalloc(FILENAMELEN, sizeof(char));
    sc_list = (double *) vcalloc(2, sizeof(double));
    
    
   
    if(P->only_tree == 0){
        if(P->only_aln == 0){
            sprintf(scoresfile, "%s%s.scores", P->outdir, (P->F)->name);
            fp=openfile(scoresfile, "w");
        }

        if(strm(P->align_method, "tcoffee")){
            libname = (char *) vcalloc(FILENAMELEN, sizeof(char));
            sprintf(libname, "%s%s.lib", P->outdir, (P->F)->name);
        }

        if(strm(P->align_method, "clustalo") || strm(P->align_method, "mafft")){
            tmpname = (char *) vcalloc(FILENAMELEN, sizeof(char));
            sprintf(tmpname, "./outtree");
            fptmp = openfile(tmpname, "w");
            fclose(fptmp);
            vfree(tmpname);
        }
    }
    
    if(P->treelist == 1){
        list = read_tree_list(P->tree_list, P->ntree);
    }
    else {
        tree_id_lenght = (5*S->nseq);
        if(!P->tree_id_activated){
            tree_id_list = declare_char(P->ntree, tree_id_lenght);
            ntreeid =0;
        }
        else {
            entries = count_entries_file(P->tree_id_file, S->nseq);
            tree_id_list = file2tree_id_list(P->tree_id_file, P->ntree, S->nseq, entries);
            ntreeid = entries;
        }
    }
     
    //gettimeofday(&tim, NULL);
    srand(1985);
  
    for(i=0, j=0; i<P->ntree && j<P->nattemps; i++, j++){
        
        tree_id = vcalloc(tree_id_lenght, sizeof(char));
        sprintf(treename, "%s%s_%d.dnd", P->outdir, (P->F)->name, i);
        sprintf(alnname, "%s%s_%d.fasta_aln", P->outdir, (P->F)->name, i);
        
        if(P->treelist == 1){
            sprintf(treename, "%s", list[i]);
        }
        else{
            T = make_tree(DM, S, treename, "nj", i, P->random_percentage, P->align_method); 
        }
       
        T[3][0]->level = 0;
        tree_id = identify_tree(T[3][0], tree_id, S->nseq);
        if(search_tree_id(tree_id_list, tree_id, ntreeid)){
            i--;
        }
        else {         
            j=0;
            sprintf(tree_id_list[ntreeid], "%s", tree_id);         
            ntreeid++;
            if(P->only_tree == 0){
                align_tree(P, treename, libname, alnname, i);
                if(P->only_aln == 0){
                    sc_list[1] = calculate_score(P, alnname, i, sc_list);

                    fprintf(fp, "%d;%s_%d;%lf\n", i, (P->F)->name, i, sc_list[1]);
                    fprintf(stdout, "\t- %d - Tree: %s_%d.dnd - Aln: %s_%d.fasta_aln - Score: %lf\n", i, (P->F)->name, i, (P->F)->name, i, sc_list[1]);
                    if(sc_list[1] > max_score){
                        max_score = sc_list[1];
                        max_ntree=i;
                        if(i!=0){
                            if(P->treelist == 0){
                                remove_file(besttree);
                                remove_file(bestalignment);
                            }
                        }
                        sprintf(besttree, "%s", treename);
                        sprintf(bestalignment, "%s", alnname);
                    }
                    else {
                        if(P->treelist == 0){
                            remove_file(treename);
                            remove_file(alnname);
                        }
                    }
                }
            }
        }
        
        vfree(tree_id);    
    }
    
    if(P->only_tree == 0){
        if(strm(P->align_method, "tcoffee")){
            remove_file(libname);
            vfree(libname);
        }
        if(strm(P->align_method, "clustalo") || strm(P->align_method, "mafft")){
            remove(tmpname);
        }
        if(P->only_aln == 0){
            fprintf(stdout, "---> Done\n");
            fprintf(stdout, "\n\t- BEST TREE: %s_%d - Score: %lf\n", (P->F)->name, max_ntree, max_score);
            if(P->treelist == 0){
                sprintf(treename, "%s%s.dnd", P->outdir, (P->F)->name);
                rename(besttree, treename);
                sprintf(alnname, "%s%s.fasta_aln", P->outdir, (P->F)->name);
                rename(bestalignment, alnname);
            }

            fprintf(fp, "Best;%s_%d;%lf\n", (P->F)->name, max_ntree, max_score);

            fclose(fp);
        }
    }
    
    free_char(tree_id_list, -1);
    vfree(scoresfile);
    vfree(alnname);
    vfree(treename);
    vfree(bestalignment);
    vfree(besttree);
    vfree(sc_list);  
    vfree(list);
    
    return max_score;
}


/* MPI METHOD*/
#elif MPI_FLAG

/*MTA Distributed version*/

 void mta_program_mpi(Parameters *P, Sequence *S, Distance_matrix *DM, int nseq){
     
    int i=0;   
    char *bufferDM, *bufferSeq;
    int tamDM=0, tamSeq=0, positionDM=0, positionSeq=0;
    
    /*Pack and Broadcast Distance Matrix*/
    tamDM = 3 * ((nseq * nseq) * sizeof(int));
    bufferDM = (char *) malloc(tamDM);
    
    tamSeq = (nseq * (S->max_len + 1) * sizeof(char)) + (nseq * (MAXNAMES + 1) * sizeof(char)) + (nseq * sizeof(int)) + (nseq * sizeof(int)) + (30 * sizeof(char));
    bufferSeq = (char *) malloc(tamSeq);
      
    if(my_rank == 0){ //master
        //srand(1985);
        for(i=0; i<nseq; i++){
            MPI_Pack((DM->similarity_matrix[i]), nseq, MPI_INT, bufferDM, tamDM, &positionDM, MPI_COMM_WORLD);
        }
        for(i=0; i<nseq; i++){
            MPI_Pack((DM->score_similarity_matrix[i]), nseq, MPI_INT, bufferDM, tamDM, &positionDM, MPI_COMM_WORLD);
        }
        for(i=0; i<nseq; i++){
            MPI_Pack((DM->distance_matrix[i]), nseq, MPI_INT, bufferDM, tamDM, &positionDM, MPI_COMM_WORLD);
        }   
        
        MPI_Bcast(bufferDM, positionDM, MPI_PACKED, 0, MPI_COMM_WORLD);
   
        /*Pack and broadcast the Sequences*/        
        MPI_Pack((S->seq_len), nseq, MPI_INT, bufferSeq, tamSeq, &positionSeq, MPI_COMM_WORLD);
        MPI_Pack((S->seq_id), nseq, MPI_INT, bufferSeq, tamSeq, &positionSeq, MPI_COMM_WORLD);
        MPI_Pack((S->type), 30, MPI_CHAR, bufferSeq, tamSeq, &positionSeq, MPI_COMM_WORLD);

        for(i=0; i<nseq; i++){
            MPI_Pack((S->seq[i]), (S->max_len + 1), MPI_CHAR, bufferSeq, tamSeq, &positionSeq, MPI_COMM_WORLD);
        }
        for(i=0; i<nseq; i++){
            MPI_Pack((S->seq_names[i]), (MAXNAMES + 1), MPI_CHAR, bufferSeq, tamSeq, &positionSeq, MPI_COMM_WORLD);
        } 

        MPI_Bcast(bufferSeq, positionSeq, MPI_PACKED, 0, MPI_COMM_WORLD);
        master_generate_mta(P, S, DM);
    }
    else {

        struct timeval tim;
        
        //gettimeofday(&tim, NULL);
        srand(1985+my_rank);
        
        
        /*Receive and unpack Distance Matrix*/        
        MPI_Bcast(bufferDM, tamDM, MPI_PACKED, 0, MPI_COMM_WORLD);
        
        for(i=0; i<nseq; i++){
            MPI_Unpack(bufferDM, tamDM, &positionDM, (DM->similarity_matrix[i]), nseq, MPI_INT, MPI_COMM_WORLD);
        } 
        for(i=0; i<nseq; i++){
            MPI_Unpack(bufferDM, tamDM, &positionDM, (DM->score_similarity_matrix[i]), nseq, MPI_INT, MPI_COMM_WORLD);
        } 
        for(i=0; i<nseq; i++){
            MPI_Unpack(bufferDM, tamDM, &positionDM, (DM->distance_matrix[i]), nseq, MPI_INT, MPI_COMM_WORLD);
        }

        /*Receive and unpack Sequences*/
        MPI_Bcast(bufferSeq, tamSeq, MPI_PACKED, 0, MPI_COMM_WORLD);
        MPI_Unpack(bufferSeq, tamSeq, &positionSeq, (S->seq_len), nseq, MPI_INT, MPI_COMM_WORLD);    
        MPI_Unpack(bufferSeq, tamSeq, &positionSeq, (S->seq_id), nseq, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack(bufferSeq, tamSeq, &positionSeq, (S->type), 30, MPI_CHAR, MPI_COMM_WORLD);

        for(i=0; i<nseq; i++){
            MPI_Unpack(bufferSeq, tamSeq, &positionSeq, (S->seq[i]), (S->max_len+1), MPI_CHAR, MPI_COMM_WORLD);
        }
        for(i=0; i<nseq; i++){
            MPI_Unpack(bufferSeq, tamSeq, &positionSeq, (S->seq_names[i]), (MAXNAMES+1), MPI_CHAR, MPI_COMM_WORLD);
        }             
        worker_generate_mta(P, S, DM);
    }

    free(bufferDM);
    free(bufferSeq);
 }
 

 
/*MTA distributed version master*/
 
void master_generate_mta(Parameters *P, Sequence *S, Distance_matrix *DM){
    
    FILE *fp;
    FILE*fptmp;
    int i=0, ret, max_ntree=-1, tam=0, position=0, tree=0;
    int work=0, from_where=0, m=0;
    int *dest;
    double max_score=-999999.99999;  
    double *sc_list;
    char *treename, *alnname, *libname, *treatf, *scoresfile;
    char *bestalignment, *besttree, *buffer;
    char *tmpname;
    NT_node **T;
   
    MPI_Status status;
    MPI_Request req;
    
    treename = (char *) vcalloc(FILENAMELEN, sizeof(char));
    alnname = (char *) vcalloc(FILENAMELEN, sizeof(char));
    besttree = (char *) vcalloc(FILENAMELEN, sizeof(char));
    bestalignment = (char *) vcalloc(FILENAMELEN, sizeof(char));
    scoresfile = (char *) vcalloc(FILENAMELEN, sizeof(char));
    treatf = (char *) vcalloc(ALLPATH, sizeof(char));

    sc_list = (double *) vcalloc(2, sizeof(double));
    tam = (1 * sizeof(int)) + (1* sizeof(double));
      
    fprintf(stdout, "\n---> Generating the Multiples Trees\n");
    
    dest = (int *)malloc((np)*sizeof(int));
    assert(dest);
    dest[0]=(np-1);
    for (i=1;i<np;i++){
        dest[i]=i;
    }
    if(P->only_tree == 0){
        if((strm(P->align_method, "tcoffee"))){
            printf("\tBuiding lib....\n");
            libname = (char *) vcalloc(FILENAMELEN, sizeof(char));
            sprintf(libname, "%s/%s.lib", P->outdir, (P->F)->name);
            sprintf(treatf, "%s %s -lib_only -out_lib %s > /dev/null 2>&1", P->align_method_bin, (P->F)->full, libname);
            ret = system(treatf);
            printf("\tDone....\n");
        }
        if(strm(P->align_method, "clustalo") || strm(P->align_method, "mafft")){
            tmpname = (char *) vcalloc(FILENAMELEN, sizeof(char));
            sprintf(tmpname, "./outtree");
            fptmp = openfile(tmpname, "w");
            fclose(fptmp);
            vfree(tmpname);
        }
        if(P->only_aln == 0){
            sprintf(scoresfile, "%s%s.scores", P->outdir, (P->F)->name);
            fp = fopen(scoresfile, "w");
        }
    }
    
   
    
    for(i=0; i<P->ntree; i++){
        buffer = malloc(tam);
        position=0;
        if(dest[0]){          
            MPI_Send(&i, 1, MPI_INT, dest[1], TAG_MSTREE, MPI_COMM_WORLD);  
            work++;
            dest[0]--;
            memmove((dest+1),(dest+2),dest[0]*sizeof(int));      
        }
        else {
            MPI_Recv(sc_list, NSCORE_TOTAL, MPI_DOUBLE, MPI_ANY_SOURCE, TAG_WSSCORE, MPI_COMM_WORLD, &status);

            tree = (int) sc_list[0];  
            if(P->only_tree == 0 && P->only_aln == 0){
                sprintf(treename, "%s%s_%d.dnd", P->outdir, (P->F)->name, tree);
                sprintf(alnname, "%s%s_%d.fasta_aln", P->outdir, (P->F)->name, tree);
                fprintf(fp, "%d;%s_%d;%lf\n", tree, (P->F)->name, tree, sc_list[1]);
                fprintf(stdout, "\t- %d - Tree: %s_%d.dnd - Aln: %s_%d.fasta_aln  - Score: %lf\n", tree, (P->F)->name, tree, (P->F)->name, (int) tree, sc_list[1]);
                if(sc_list[1] > max_score){
                    max_score = sc_list[1];
                    max_ntree=tree;
                    if(P->treelist == 0){
                        remove_file(besttree);
                        remove_file(bestalignment);
                    }

                    sprintf(besttree, "%s", treename);
                    sprintf(bestalignment, "%s", alnname);
                }
                else {
                    if(P->treelist == 0){
                        remove_file(treename);
                        remove_file(alnname);
                    }
                }
            }

            work--;
            dest[dest[0]+1]=status.MPI_SOURCE;
            dest[0] += 1;

            MPI_Send(&i, 1, MPI_INT, dest[1], TAG_MSTREE, MPI_COMM_WORLD); 
            work++;
            dest[0]--;
            memmove((dest+1),(dest+2),dest[0]*sizeof(int));
        }
        free(buffer);    
    }
    
    for(i=work; work>0; work--){
        buffer = malloc(tam);
        position=0;
        
        MPI_Recv(sc_list, NSCORE_TOTAL, MPI_DOUBLE, MPI_ANY_SOURCE, TAG_WSSCORE, MPI_COMM_WORLD, &status);
        
        if(P->only_tree == 0 && P->only_aln == 0){
            sprintf(treename, "%s%s_%d.dnd", P->outdir, (P->F)->name, tree);
            sprintf(alnname, "%s%s_%d.fasta_aln", P->outdir, (P->F)->name, tree);
            fprintf(fp, "%d;%s_%d;%lf\n", tree, (P->F)->name, tree, sc_list[1]);
            fprintf(stdout, "\t- %d - Tree: %s_%d.dnd - Aln: %s_%d.fasta_aln  - Score: %lf\n", tree, (P->F)->name, tree, (P->F)->name, tree, sc_list[1]);
            if(sc_list[1] > max_score){
                printf("%f - %f\n", sc_list[1], max_score);
                max_score = sc_list[1];
                max_ntree=tree;
                if(P->treelist == 0){
                    remove_file(besttree);
                    remove_file(bestalignment);
                }
                #printf("removed_best: %s\n", alnname);

                sprintf(besttree, "%s", treename);
                sprintf(bestalignment, "%s", alnname);
            }
            else {
                if(P->treelist == 0){
                    remove_file(treename);
                    remove_file(alnname);
                }
            }
        }
 
        dest[dest[0]+1]=from_where;
        dest[0] += 1;
        free(buffer);
    }
    
    if(work == 0){
        int fin = FIN_MSG;
        for(i=0; i<=dest[0]; i++){
            MPI_Isend(&fin, 1, MPI_INT, i, TAG_MSTREE, MPI_COMM_WORLD, &req);
        }  
    }
    if(P->only_tree == 0){
        if(strm(P->align_method, "tcoffee")){
            remove_file(libname);
            vfree(libname);
        }  
         if(strm(P->align_method, "clustalo") || strm(P->align_method, "mafft")){
            remove(tmpname);
        }
        if(P->only_aln == 0){
            if(P->treelist == 0){
                sprintf(treename, "%s%s.dnd", P->outdir, (P->F)->name);
                rename(besttree, treename);
                sprintf(alnname, "%s%s.fasta_aln", P->outdir, (P->F)->name);
                rename(bestalignment, alnname);
            }

            fprintf(stdout, "\n\t- BEST TREE: %s_%d - Score: %lf\n", (P->F)->name, max_ntree, max_score);
            fprintf(fp, "Best;%s_%d;%lf\n", (P->F)->name, max_ntree, max_score);

            fprintf(stdout, "---> Done\n");

            fclose(fp);
        }
    }
    
    vfree(treename);
    vfree(alnname);
    vfree(besttree);
    vfree(bestalignment);
    vfree(scoresfile);
    vfree(treatf);
    vfree(sc_list);
    free(dest);    
}

/*MTA distributed version master*/

void worker_generate_mta(Parameters *P, Sequence *S, Distance_matrix *DM){
    
    int tam=0, i=0, tree=0;
    double *sc_list;
    char *treename, *alnname, *libname=NULL;
    char *buffer;
    NT_node **T;
    char **list;
    
    MPI_Status status;
    MPI_Request req;
      
    tam = (1 * sizeof(int)) + (1 * sizeof(double));
    sc_list = (double *) vcalloc(NSCORE_TOTAL, sizeof(double));
    
    for(i=0; i<NSCORE_TOTAL; i++){
        sc_list[i]=0.0;
    }
    if(P->only_tree == 0){
        if(strm(P->align_method, "tcoffee")){
            libname = (char *) vcalloc(FILENAMELEN, sizeof(char));
            sprintf(libname, "%s/%s.lib", P->outdir, (P->F)->name);
        }
    }
    
    if(P->treelist == 1){
        list = read_tree_list(P->tree_list, P->ntree);
    }
    
    while(tree != FIN_MSG){
        buffer = (char *) malloc(tam);
        //position=0;
        MPI_Recv(&tree, 1, MPI_INT, 0, TAG_MSTREE, MPI_COMM_WORLD, &status);
        
        if(tree != FIN_MSG){
            treename = (char *) vcalloc(FILENAMELEN, sizeof(char));
            alnname = (char *) vcalloc(FILENAMELEN, sizeof(char));
            
            if(P->treelist == 1){
                sprintf(treename, "%s, list[i]);
            }
            else {
                sprintf(treename, "%s/%s_%d.dnd", P->outdir, (P->F)->name, tree);
            }
            sprintf(alnname, "%s/%s_%d.fasta_aln", P->outdir, (P->F)->name, tree);

            T = make_tree(DM, S, treename, "nj", tree, P->random_percentage, P->align_method);
            if(P->only_tree == 0){
                align_tree(P, treename, libname, alnname, tree);
                if(P->only_aln == 0){
                        sc_list[1] = calculate_score(P, alnname, tree, sc_list);
                }
            }

            MPI_Isend(sc_list, NSCORE_TOTAL, MPI_DOUBLE, 0, TAG_WSSCORE, MPI_COMM_WORLD, &req);   
              
            vfree(alnname); 
            vfree(treename);
        }
        free(buffer);      
    }
    if(P->only_tree == 0){
        if(strm(P->align_method, "tcoffee")){
            vfree(libname);
        }
    }
    
    vfree(sc_list);
    vfree(list);
}

#endif
 
/*Align a Tree with the chosen MSA method*/
void align_tree(Parameters *P, char *treename, char *libname, char *alnname, int ntree){

    char *align_app;
    int ret;
    align_app = (char *) vcalloc(ALLPATH, sizeof(char));

    if(strm(P->align_method, "tcoffee")){
        if((P->PB)->tcoffee_bin != NULL){
            if(P->mpi_mode == 0 && libname != NULL){
                if(ntree == 0){              
                    sprintf(align_app, "%s %s -usetree %s -out_lib %s -output=fasta -n_core=1 -outfile=%s > /dev/null 2>&1", (P->PB)->tcoffee_bin, (P->F)->full, treename, libname,alnname);
                }
                else {
                    sprintf(align_app, "%s %s -usetree %s -lib %s -output=fasta -n_core=1 -outfile=%s > /dev/null 2>&1", (P->PB)->tcoffee_bin, (P->F)->full, treename, libname, alnname);
                }
            }
            else {           
                sprintf(align_app, "%s %s -usetree %s -lib %s -output=fasta -n_core=1 -outfile=%s > /dev/null 2>&1", (P->PB)->tcoffee_bin, (P->F)->full, treename, libname, alnname);
            }
            ret=system(align_app);
        }
        else {
            fprintf(stderr, "Error - T-Coffee not found. Correct T-Coffee path\n");
            exit(EXIT_FAILURE);
        }
    }
    else if(strm(P->align_method, "clustalw")){
        if((P->PB)->clustalw_bin != NULL){
            sprintf(align_app, "%s -infile=%s -usetree=%s -output=fasta -outfile=%s > /dev/null 2>&1", (P->PB)->clustalw_bin, (P->F)->full, treename, alnname);
            ret=system(align_app);
        }
        else {
            fprintf(stderr, "Error - ClustalW not found. Correct ClustalW path\n");
            exit(EXIT_FAILURE);
        }
    }  
    else if(strm(P->align_method, "clustalo")){
        if((P->PB)->clustalo_bin != NULL){
            root_unrooted_tree(treename, ntree, (P->PB)->retree_bin);
            sprintf(align_app, "%s -i %s --guidetree-in=%s --outfmt=fa -o %s", (P->PB)->clustalo_bin, (P->F)->full, treename, alnname);
            ret=system(align_app);
        }
        else {
            fprintf(stderr, "Error - ClustalO not found. Correct ClustalO path\n");
            exit(EXIT_FAILURE);
        }
    }
    else if(strm(P->align_method, "mafft")){
        if((P->PB)->mafft_bin != NULL && (P->PB)->nwtomafft != NULL){
            char *treename2;
            treename2 = (char *) vcalloc(ALLPATH, sizeof(char));
            
            root_unrooted_tree(treename, ntree, (P->PB)->retree_bin);

            sprintf(treename2, "%s%s_%d.mafft", P->outdir, (P->F)->name, ntree);
            sprintf(align_app, "%s %s > %s", (P->PB)->nwtomafft, treename, treename2);
            ret=system(align_app);

            sprintf(align_app, "%s --auto --treein %s %s > %s 2> /dev/null", (P->PB)->mafft_bin, treename2, (P->F)->full, alnname);
            ret=system(align_app);

            remove_file(treename2);
            vfree(treename2);
        }
        else {
            fprintf(stderr, "Error - Mafft not found. Correct Mafft path\n");
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "Error - Bad MSA Method\n");
        exit(EXIT_FAILURE);
    }
    
    vfree(align_app);
}

char** read_tree_list(char* tree_list, int ntree){
    FILE *fp;
    char **list;
    int i=0;
    
    list = (char **) vcalloc(ntree, sizeof(char *));
    for(i=0; i<ntree; i++){
        list[i] = (char *) vcalloc(ALLPATH, sizeof(char));
    }
    
    fp = openfile(tree_list, "r");
    
    for(i=0; (i<ntree) && (fscanf(fp, "%s", list[i]) != EOF); i++);
       
       
    if (i != ntree){
        fprintf(stderr, "ERROR - Incorrect number of trees (ntree) %d != (list) %d\n", ntree, i);
    }
    
    fclose(fp);
    
    return list;
    
}