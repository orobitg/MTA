
#include <limits.h>

#include "parameters_utils.h"

/*Parameters declaration*/

Parameters *declare_parameters(){

    Parameters *P;

    P = vcalloc(1, sizeof(Parameters));
    P->dm_method = (char *) vcalloc(NAME_MAX, sizeof(char));
    P->align_method = (char *) vcalloc(NAME_MAX, sizeof(char));
    P->align_method_bin = (char *) vcalloc(NAME_MAX, sizeof(char));
    P->score_method = (char *) vcalloc(NAME_MAX, sizeof(char));
    P->outdir = (char *) vcalloc(FILENAMELEN, sizeof(char));
    P->workdir = (char *) vcalloc(FILENAMELEN, sizeof(char));
    P->seqfile = (char *) vcalloc(FILENAMELEN, sizeof(char));
    P->tree_id_file = (char *) vcalloc(FILENAMELEN, sizeof(char));
    P->tree_file = (char *) vcalloc(FILENAMELEN, sizeof(char));
    P->tree_list = (char *) vcalloc(FILENAMELEN, sizeof(char));
    P->mat = (char *) vcalloc(STRING, sizeof(char));
    P->str_file = (char *) vcalloc(ALLPATH, sizeof(char));
    P->irmsd_file = (char *) vcalloc(ALLPATH, sizeof(char));
    P->other_pg_mode = (char *) vcalloc(FILENAMELEN, sizeof(char));
    
    return P;
}

void free_parameters(Parameters *P){

    vfree(P->dm_method);
    vfree(P->align_method);
    vfree(P->align_method_bin);
    vfree(P->score_method);
    vfree(P->outdir);
    vfree(P->workdir);
    vfree(P->seqfile);
    vfree(P->tree_id_file);
    vfree(P->tree_file);
    vfree(P->tree_list);
    vfree(P->mat);
    vfree(P->str_file);
    vfree(P->irmsd_file);
    vfree(P->other_pg_mode);
    free_fname(P->F);
    free_paths(P->PB);
    
    vfree(P);
}


Parameters *default_values(Parameters *P, char *outdir){

    sprintf(P->align_method, "t_coffee");
    sprintf(P->align_method_bin, "%s", (P->PB)->tcoffee_bin);
    sprintf(P->outdir, "./results/");
    sprintf(P->workdir, "%swork/", P->outdir);
    P->mpi_mode = 0;
    
    if(!P->other_pg){
        P->other_pg = 0;
        P->ntree = 10;
        P->seq_activated = 0;
        P->rep = 1;
        P->nattemps = 500;
        P->tree_id_activated = 0;
        P->random_percentage = 100;
        P->treelist = 0;
        P->gop = -11;
        P->gep = -1;
        P->only_tree = 0;
        P->only_aln = 0;
        
        sprintf(P->mat, "blosum62mt");
        sprintf(P->dm_method, "ktup");      
        sprintf(P->score_method, "sp");    
    }
    
    return P;
}

/*Function to read the input parameters*/

Parameters *read_parameters(int argc, char** argv){

    Parameters *P;
    int i;
    int rep_activated=0;
    int mpi_activated=0;
    
    if(argc < 2 || argc > 26){
        fprintf(stderr, "ERROR - Number of parameters wrong\n\n");
        printhelp();
        exit(EXIT_FAILURE);
    }

    P = declare_parameters();
    P->PB = parse_paths(argv[0]); 

    if(argc == 2 && strm(argv[1], "-help")){
        printhelp();
        exit(EXIT_SUCCESS);
    }
    else if(argc != 2 && strm(argv[1], "-help")){
            fprintf(stderr, "ERROR -help: Number of parameters wrong\n");
            printhelp();
            exit(EXIT_FAILURE);
    }

    if(strm(argv[1], "-other_pg")){
        P->other_pg = 1;
        if(strm(argv[2], "-tree_id")){
            sprintf(P->other_pg_mode, "%s", argv[2]);
            for(i=3; i<argc; i++){
                if(i<argc && strm(argv[i], "-seq")){
                    i++;
                    if(argv[i][0] != '-'){
                        sprintf(P->seqfile, "%s", argv[i]);
                        P->F = parse_fname(P->seqfile);
                        P = default_values(P, (P->F)->path);
                    }
                    else {
                        fprintf(stderr, "ERROR - Parameter -seq: No Sequences\n");
                        exit(EXIT_FAILURE);
                    }
                }
                else if(strm(argv[i], "-tree")){
                    i++;
                    if(i<argc && argv[i][0] != '-'){
                        sprintf(P->tree_file, "%s", argv[i]);
                    }
                    else {
                        fprintf(stderr, "ERROR - Parameter -tree: No Tree\n");
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                     fprintf(stderr, "ERROR - Parameter -tree_id: Bad parameter\n");
                        exit(EXIT_FAILURE);
                }
            }
        }   
    }

    if(!P->other_pg){
        for(i=1; i<argc; i++){
            if(strm(argv[i], "-seq")){
                i++;
                if(i<argc && argv[i][0] != '-'){
                    sprintf(P->seqfile, "%s", argv[i]);
                    P->F = parse_fname(P->seqfile);
                    P = default_values(P, (P->F)->path);
                    P->seq_activated = 1;
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -seq: No Sequences\n");
                    exit(EXIT_FAILURE);
                }
            }
            else if(strm(argv[i], "-rep")){
                i++;
                if(i<argc && argv[i][0] != '-'){
                    rep_activated = 1;
                    P->rep = atoi(argv[i]);
                    if((P->rep != 0) && (P->rep != 1)){
                        fprintf(stderr, "ERROR - Parameter -rep: Wrong parameter value (0 or 1)\n");
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -rep: No value (0 or 1)\n");
                    exit(EXIT_FAILURE);
                }
            }
            else if(strm(argv[i], "-mpi")){
                mpi_activated = 1;
                P->mpi_mode = 1;              
            }
        }

        if(P->seq_activated == 0){
            fprintf(stderr, "ERROR - No Sequences\n\n");
            printhelp();
            exit(EXIT_FAILURE);
        }

        for(i=1; i<argc; i++){          
            if(strm(argv[i], "-seq")){
                i++;
            }
            else if(strm(argv[i], "-ntree")){
                i++;
                if(i<argc && argv[i][0] != '-'){
                    P->ntree = atoi(argv[i]);
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -ntree: No Value (Number of trees)\n");
                    exit(EXIT_FAILURE);
                }
            }
            else if(strm(argv[i], "-dm")){
                i++;
                if(i<argc && argv[i][0] != '-'){
                    if(strm(argv[i], "ktup")){ /*I put this parameter if we want to use another dm methods. At the moment, ktup DM*/
                        sprintf(P->dm_method, "%s", argv[i]);
                    }
                    else {
                        fprintf(stderr, "ERROR - Parameter -dm: Wrong distance matrix values (ktup)\n");
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -dm: No Value (ktup)\n");
                    exit(EXIT_FAILURE);

                }
            }
            else if(strm(argv[i], "-msa")){
                i++;
                if(i<argc && argv[i][0] != '-'){  //Aixo es pot unir
                    if(strm(argv[i], "t_coffee") || strm(argv[i], "clustalw") || strm(argv[i], "clustalo") || strm(argv[i], "mafft")){
                        sprintf(P->align_method, "%s", argv[i]);  
                    }
                    else {
                        fprintf(stderr, "ERROR - Parameter -msa: Wrong MSA method (t_coffee, clustalw, clustalo or mafft)\n");
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -msa: No Value (t_coffee, clustalw, clustalo)\n");
                    exit(EXIT_FAILURE);
                }
            }
            else if(strm(argv[i], "-score")){
                i++;
                if(i<argc && argv[i][0] != '-'){
                    if(strm(argv[i], "sp") || strm(argv[i], "normd") || strm(argv[i], "tcs") || strm(argv[i], "triplet") || strm(argv[i], "coffee") || strm(argv[i], "strike") || strm(argv[i], "irmsd")){
                        sprintf(P->score_method, "%s", argv[i]);
                        P = score_parameters(P, argc, argv);
                    }
                    else {
                        fprintf(stderr, "ERROR - Parameter -score: Wrong score method (sp, normd)\n");
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -score: No Value (sp, normd)\n");
                    exit(EXIT_FAILURE);
                }
            }
            else if(strm(argv[i], "-only_tree")){
                P->only_tree = 1;
            }
             else if(strm(argv[i], "-only_aln")){
                P->only_aln = 1;
            }
            else if(strm(argv[i], "-rep")){
                i++;
            }
            else if(strm(argv[i], "-mpi")){
                
            }
            else if(strm(argv[i], "-%tree")){
                i++;
                if(i<argc && argv[i][0] != '-'){
                    P->random_percentage = atoi(argv[i]);
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -'%%'tree: No Value\n");
                    exit(EXIT_FAILURE);
                }
            }
            else if(strm(argv[i], "-nattempt")){
                if(rep_activated==1 && P->rep == 0 && P->mpi_mode == 0){
                    i++;
                    if(i<argc && argv[i][0] != '-'){
                        P->nattemps = atoi(argv[i]);
                    }
                    else {
                        fprintf(stderr, "ERROR - Parameter -nattempt: No Value\n");
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -nattemp: No repeated option is enabled or MPI version is enabled\n");
                    exit(EXIT_FAILURE);
                }
            }
            else if(strm(argv[i], "-treefile")){ //To pass a tree id file
                if(rep_activated==1 && P->rep == 0 && P->mpi_mode == 0){
                    i++;
                    if(i<argc && argv[i][0] != '-'){
                        P->tree_id_activated = 1;
                        sprintf(P->tree_id_file, "%s", argv[i]);
                    }
                    else {
                        fprintf(stderr, "ERROR - Parameter -treefile: No Value\n");
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -treefile: No repeated option is enabled or MPI version is enabled\n");
                    exit(EXIT_FAILURE);
                }
            }
            else if(strm(argv[i], "-treelist")){ //To pass a tree id file
                i++;
                if(i<argc && argv[i][0] != '-'){
                        P->treelist = 1;
                        sprintf(P->tree_list, "%s", argv[i]);
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -treefile: No Value\n");
                    exit(EXIT_FAILURE);
                }
               
            }
            else if(strm(argv[i], "-output")){
                i++;
                if(i<argc && argv[i][0] != '-'){
                    sprintf(P->outdir, "%s/", argv[i]);
                    sprintf(P->workdir, "%swork/", P->outdir);
                }
                else {
                    exit(EXIT_FAILURE);
                }
            }
            else if(strm(argv[i], "-matrix") || strm(argv[i], "-gop") || strm(argv[i], "-gep") || strm(argv[i], "-egap") || strm(argv[i], "-pdb_strike")|| strm(argv[i], "-pdb_irmsd")){
                i++;
            }
            else {
                fprintf(stderr, "ERROR - Parameters: Wrong Parameters\n");
                exit(EXIT_FAILURE);
            }
        }        
    }
    
    return P;
}

/*Read the score parameters like substitution matrices, gap penalties, templates*/

Parameters *score_parameters(Parameters *P, int argc, char** argv){
    
    int i;
    
    if(strm(P->score_method, "sp")){
        for(i=1; i<argc; i++){
            if(strm(argv[i], "-matrix")){
                i++;
                if(i<argc && argv[i][0] != '-'){
                    if (strm(argv[i], "pam250mt"))sprintf(P->mat, "pam250mt");
                    else if (strm(argv[i], "idmat"))sprintf(P->mat, "idmat");
                    else if (strm(argv[i], "dna_idmat"))sprintf(P->mat, "idmat");
                    else if (strm(argv[i], "est_idmat"))sprintf(P->mat, "est_idmat");
                    else if (strm(argv[i], "md_350mt"))sprintf(P->mat, "md_350mt");
                    else if (strm(argv[i], "md_250mt"))sprintf(P->mat, "md_250mt");
                    else if (strm(argv[i], "md_120mt"))sprintf(P->mat, "md_120mt");
                    else if (strm(argv[i], "md_40mt" ))sprintf(P->mat, " md_40mt");
                    else if (strm(argv[i], "pam350mt" ))sprintf(P->mat, "pam350mt");
                    else if (strm(argv[i], "pam160mt" ))sprintf(P->mat, "pam160mt");
                    else if (strm(argv[i], "pam120mt" ))sprintf(P->mat, "pam120mt");

                    else if (strm(argv[i], "blosum80mt" ))sprintf(P->mat, "blosum80mt");
                    else if (strm(argv[i], "blosum62mt" ))sprintf(P->mat, "blosum62mt");
                    else if (strm(argv[i], "exon2mt" ))sprintf(P->mat, "blosum62mt");
                    else if (strm(argv[i], "blosum62mt3" ))sprintf(P->mat, "blosum62mt3");

                    else if (strm(argv[i], "blosum62mt2" ))sprintf(P->mat, "blosum62mt2");
                    else if (strm(argv[i], "blosum55mt" ))sprintf(P->mat, "blosum55mt");
                    else if (strm(argv[i], "blosum50mt" ))sprintf(P->mat, "blosum50mt");
                    else if (strm(argv[i], "blosum45mt" ))sprintf(P->mat, "blosum45mt");

                    else if (strm(argv[i], "blosum40mt" ))sprintf(P->mat, "blosum40mt");
                    else if (strm(argv[i], "blosum30mt" ))sprintf(P->mat, "blosum30mt");
                    else if (strm(argv[i], "beta_mat" ))sprintf(P->mat, "beta_mat");
                    else if (strm(argv[i], "alpha_mat" ))sprintf(P->mat, "alpha_mat");
                    else if (strm(argv[i], "coil_mat" ))sprintf(P->mat, "coil_mat");

                    else if (strm(argv[i], "rblosum80mt" ))sprintf(P->mat, "rblosum80mt");
                    else if (strm(argv[i], "rblosum62mt" ))sprintf(P->mat, "rblosum62mt");
                    else if (strm(argv[i], "rblosum30mt" ))sprintf(P->mat, "rblosum30mt");

                    else if (strm(argv[i], "rpam250mt" ))sprintf(P->mat, "rpam250mt");
                    else if (strm(argv[i], "rpam350mt" ))sprintf(P->mat, "rpam350mt");
                    else if (strm(argv[i], "rpam160mt" ))sprintf(P->mat, "rpam160mt");
                    else if (strm(argv[i], "rpam120mt" ))sprintf(P->mat, "rpam120mt");

                    else if (strm(argv[i], "tmpam250mt" ))sprintf(P->mat, "tmpam250mt");
                    else if (strm(argv[i], "rtmpam250mt" ))sprintf(P->mat, "rtmpam250mt");

                    else if (strm(argv[i], "rbeta_mat" ))sprintf(P->mat, "rbeta_mat");
                    else if (strm(argv[i], "ralpha_mat" ))sprintf(P->mat, "ralpha_mat");
                    else if (strm(argv[i], "rcoil_mat" ))sprintf(P->mat, "rcoil_mat");
                    else if (strm (argv[i], "jtttm250mt"))sprintf(P->mat, "jtttm250mt");

                    else if (strm (argv[i], "promoter_tf1"))sprintf(P->mat, "promoter_tf1");
                    else {
                        fprintf(stderr, "ERROR - Parameter -mat: Wrong matrix\n");
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -mat: No Value\n");
                    exit(EXIT_FAILURE);
                }

            }
            else if (strm(argv[i], "-gop")){
                i++;
                if(i<argc){
                    P->gop = atof(argv[i]);
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -gop: No Value\n");
                    exit(EXIT_FAILURE);
                }
                
            }
            else if (strm(argv[i], "-gep")){
                i++;
                if(i<argc){
                    P->gep = atof(argv[i]);
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -gep: No Value\n");
                    exit(EXIT_FAILURE);
                }              
            }
        }
    }
    else if(strm(P->score_method, "normd")){
        
        P->gop = 0.0;
        P->gep = 0.1;
        P->egap = 0.0;
        sprintf(P->mat, "blosum62.bla");
        
        for(i=1; i<argc; i++){
            if(strm(argv[i], "-matrix")){
                i++;
                if(i<argc && argv[i][0] != '-'){
                    if (strm(argv[i], "blosum62"))sprintf(P->mat, "blosum62.bla");
                    else if (strm(argv[i], "gon250"))sprintf(P->mat, "gon250.bla");
                    else {
                        fprintf(stderr, "ERROR - Parameter -mat: Wrong matrix\n");
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -mat: No Value\n");
                    exit(EXIT_FAILURE);
                }    
            }
            else if (strm(argv[i], "-gop")){
                i++;
                if(i<argc){
                    P->gop = atof(argv[i]);
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -gop: No Value\n");
                    exit(EXIT_FAILURE);
                }
                
            }
            else if (strm(argv[i], "-gep")){
                i++;
                if(i<argc){
                    P->gep = atof(argv[i]);
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -gep: No Value\n");
                    exit(EXIT_FAILURE);
                }
            }
            else if (strm(argv[i], "-egap")){
                i++;
                if(i<argc){
                    P->egap = atof(argv[i]);
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -egap: No Value\n");
                    exit(EXIT_FAILURE);
                }
            }
           
        }
    }
    else if(strm(P->score_method, "strike")){
        int file_true=0;
        for(i=1; i<argc; i++){
            if(strm(argv[i], "-pdb_strike")){
                file_true = 1;
                i++;
                if(i<argc && argv[i][0] != '-'){
                    sprintf(P->str_file, "%s", argv[i]);
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -pdb_strike: No Value\n");
                    exit(EXIT_FAILURE);
                }
            }
        }
        if(file_true == 0 && (strm(P->score_method, "strike"))){
            fprintf(stderr, "ERROR - No connection or template file, necessary to strike, irmsd scores\n");
            exit(EXIT_FAILURE);
        }
        
    }
    else if(strm(P->score_method, "irmsd")){
        int file_true=0;
        for(i=1; i<argc; i++){
            if(strm(argv[i], "-pdb_irmsd")){
                file_true = 1;
                i++;
                if(i<argc && argv[i][0] != '-'){
                    sprintf(P->irmsd_file, "%s", argv[i]);
                }
                else {
                    fprintf(stderr, "ERROR - Parameter -pdb_irmsd: No Value\n");
                    exit(EXIT_FAILURE);
                }
            } 
        }
        if(file_true == 1 && strm(P->score_method, "irmsd")){
            sprintf(P->score_method, "irmsdt");
        }
    }
        
    return P;   
}

int search_mode(int argc, char** argv, char *mode){
    
    int i=0;
    if(argc < 2){
        fprintf(stderr, "ERROR - Number of parameters wrong\n\n");
        printhelp();
        exit(EXIT_FAILURE);
    }
    
    for(i=1; i<argc; i++){
        if(strm(argv[i], mode)){
            return 1;
        } 
    }
    
    return 0;
}

/*REDOOOOOOO*/
void printhelp(){
    fprintf(stdout, "\n******************HELP MENU***************************************************\n\n");
    /*fprintf(stdout, "- USAGE: $Multiple_Trees Sequence_file.fasta (Parameters)\n");
    fprintf(stdout, "\t-seq (Sequence file): The number of guide trees (required).\n");
    fprintf(stdout, "\t-ntree (N): The number of guide trees (optional - default=10).\n");
    fprintf(stdout, "\t-dm (method): The distance matrix method (ktup) (optional - default=ktup)\n");
    fprintf(stdout, "\t-msa (method): The aligner method (t_coffee, clustaw, clustalo) (optional - default=t_coffee)\n");
    fprintf(stdout, "\t-score (method): The score method (sp, normd, coffee, triplet, strike, irmsd) (optional - defualt=sp)\n");
    fprintf(stdout, "\t-output: To pass the name and the path of the output file (optional).");
    fprintf(stdout, "\t-rep (N): 0 - Repeated trees not acceptred. 1 - Repeated trees accepted (optional - default=0).\n");
    fprintf(stdout, "\t- If rep=0:\n");
    fprintf(stdout, "\t\t-nattempt (N): Number of attempts to find a no repeated tree (optional - default=500)\n");
    fprintf(stdout, "\t\t-treefile (Trees_id list file): To pass a file with a list of tree_ids that will be considered repeated trees (optional)\n");
    fprintf(stdout, "\n- OUTPUT:\n");
    fprintf(stdout, "\t- The best newick tree\n");
    fprintf(stdout, "\t- The best fasta alignment\n");
    fprintf(stdout, "\n- EXAMPLES:\n");
    fprintf(stdout, "\t- $Multiple_trees -seq sample.fasta\n");
    fprintf(stdout, "\t- $Multiple_trees -seq sample.fasta -ntree 100 -msa t_coffee -score normd -output dir/\n");
    fprintf(stdout, "\t- $Multiple_trees -seq sample.fasta -ntree 100 -rep 0 -nattempt 1000 -output dir/\n");
    fprintf(stdout, "\n**************************************************************************\n\n");*/
}
