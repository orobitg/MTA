
#include "distance_matrix_utils.h"

static int tsize;

/*Declare and free Distance matrix*/

Distance_matrix *declare_distance_matrix(int nseq){

    Distance_matrix *DM;

    DM=vcalloc ( 1, sizeof (Distance_matrix));
    DM->nseq=nseq;
    DM->similarity_matrix = declare_int(nseq, nseq);
    DM->distance_matrix = declare_int (nseq, nseq);
    DM->score_similarity_matrix=declare_int (nseq, nseq);

    return DM;
}

void free_distance_matrix(Distance_matrix *DM){

    free_int(DM->similarity_matrix, -1);
    free_int(DM->distance_matrix, -1);
    free_int(DM->score_similarity_matrix, -1);

    vfree(DM);
}

/*Create the distance matrix*/

Distance_matrix *make_distance_matrix(Sequence *S, char *mode){

    Distance_matrix *DM;
    int **sim_table=NULL;
    int ktup=6;
    int a=0, b=0;
    float score=0;
    int  id_score;
    int *ns;
    int **l_s;
    float id;

    fprintf(stdout, "\n---> Generationg the distance matrix\n");
    DM = declare_distance_matrix(S->nseq);
    //printf("MODE: %s\n", mode);
    if(strm(mode, "ktup")){
        sim_table = ktup_dist_mat(S->seq, S->nseq, ktup, S->type);
    }
    else if(strm(mode, "random")){
        sim_table = random_dist_mat(S->seq, S->nseq);
    }
    else {
        fprintf(stderr, "ERROR - Bad distance matrix mode\n");
        exit(EXIT_FAILURE);
    }
    //printf("SIM table fet\n");
    ns=vcalloc ( 2, sizeof(int));
    l_s=declare_int ( 2, 1);
    ns[0]=ns[1]=1;
    l_s[0][0]=0;
    l_s[1][0]=1;		

    for(a=0; a<S->nseq; a++){
        for(b=a; b<S->nseq; b++){
            if(b == a){
                DM->similarity_matrix[a][b]=MAXID;
            }
            else{
                l_s[0][0]=a;
                l_s[1][0]=b;
                id=sim_table[a][b];
                score=id*SCORE_K;

                /*Sim mat*/
                DM->similarity_matrix[a][b]=DM->similarity_matrix[b][a]=(int)(id);
                /*Dist mat*/
                DM->distance_matrix[a][b]=DM->distance_matrix[b][a]=MAXID-(int)(id);
                /*Score mat*/

                DM->score_similarity_matrix[a][b]=DM->score_similarity_matrix[b][a]=(int)score;
                id_score=id;

            }
        }
    }
    
    vfree (ns);
    free_int(l_s, -1);
    free_int(sim_table, -1);
    fprintf(stdout, "---> DONE\n");
    return DM;
}

/*Calculate random distances - Test matrix.*/

int **random_dist_mat(char **seq, int nseq){
    int **pscore;
    int i=0, j=0;

    pscore=declare_int ( nseq, nseq);
    srand((unsigned)time(0));

    for(i=0; i<nseq; i++){
        for(j=i; j<nseq; j++){
            if(i == j){
                pscore[i][j] = pscore[j][i] = 0;
            }
            else {
               //pscore[i][j] = pscore[j][i] = lowest + int(range*rand()/(RAND_MAX +1.0));
                pscore[i][j] = pscore[j][i] = (rand()%100)+1;
            }

        }
    }

    return pscore;
}

/*Calculate k-tup distances*/

int **ktup_dist_mat(char **seq, int nseq, int ktup, char *type){
  //Adapted from MAFFT 5: fast ktup
    int **pointt,*code=NULL, **pscore;
    int i, l, j, minl;
    double **mtx, score0;

    if(!seq || nseq==0)return NULL;
    for(minl=strlen(seq[0]),l=0,i=0;i<nseq; i++){
        int len;
        len=strlen (seq[i]);
        minl=MIN(minl, len);
        l=MAX(l,len);
    }
    ktup=MIN(minl, ktup);
    pointt=declare_int (nseq, l+1);
    mtx=declare_double (nseq, nseq);
    pscore=declare_int ( nseq, nseq);

    for(i=0; i<nseq; i++){
        makepointtable( pointt[i], code=code_seq (seq[i], type),ktup);
    }
    tsize=(int)pow(code[-1], ktup);

    for(i=0; i<nseq; i++){
        short *table1;
        table1=vcalloc ( tsize,sizeof (short));
        makecompositiontable( table1, pointt[i]);
        for (j=i; j<nseq; j++){
            mtx[i][j] = commonsextet(table1, pointt[j] );
        }
        vfree (table1);
    }
    for(i=0; i<nseq; i++ ){
        score0 = mtx[i][i];
        for(j=0; j<nseq; j++){
            pscore[i][j] = (int)( ( score0 - mtx[MIN(i,j)][MAX(i,j)] ) / score0 * 3 * 10.0 + 0.5 );
        }
    }
    for(i=0; i<nseq-1; i++)
    for(j=i+1; j<nseq; j++ ){
        pscore[i][j] = pscore[j][i]=100-MIN( pscore[i][j], pscore[j][i] );
    }
    
    return pscore;
}

int *makepointtable( int *pointt, int *n, int ktup ){

    int point, a, ng;
    register int *p;
    static int *prod;

    ng=n[-1];

    if(!prod){
        prod=vcalloc ( ktup, sizeof (int));
        for(a=0; a<ktup; a++){
            prod[ktup-a-1]=(int)pow(n[-1],a);
	}
    }
    p = n;

    for(point=0,a=0; a<ktup; a++){
        point+= *n++ *prod[a];
    }

    *pointt++ = point;

    while(*n != END_ARRAY){
        point -= *p++ * prod[0];
        point *= ng;
        point += *n++;
        *pointt++ = point;
    }
    *pointt = END_ARRAY;
    
    return pointt;
}

int *code_seq (char *seq, char *type){

    static int *code;
    static int *aa, ng;
    int a, b, l;

    if (!aa){
        char **gl;
        if(strm (type, "DNA") || strm (type, "RNA")){
            gl=declare_char (4,5);
            sprintf ( gl[ng++], "Aa");
            sprintf ( gl[ng++], "Gg");
            sprintf ( gl[ng++], "TtUu");
            sprintf ( gl[ng++], "Cc");
	}
        else{
            gl=make_group_aa_mafft(&ng, "mafft");
	}
        aa=vcalloc(256, sizeof (int));
        for(a=0; a<ng; a++){
            for(b=0; b< strlen (gl[a]); b++){
                aa[(int)gl[a][b]]=a;
	    }
	}
        free_char (gl, -1);
    }

    l=strlen (seq);

    if(code) code--;
    if(!code || read_array_size (code, sizeof (int))<(l+2)){
        vfree (code);
        code=vcalloc (l+2, sizeof (int));
    }
    code[0]=ng;
    code++;
    for(a=0; a<l; a++){
        code[a]=aa[(int)seq[a]];
    }
    code[a]=END_ARRAY;

    return code;
}

void makecompositiontable(short *table, int *pointt){

    int point;

    while((point = *pointt++) != END_ARRAY){
        table[point]++;
    }
}

int commonsextet(short *table, int *pointt){

    int value = 0;
    short tmp;
    int point;
    static short *memo = NULL;
    static int *ct = NULL;
    static int *cp;

    if(!memo){
        memo = vcalloc( tsize+1, sizeof( short ) );
        ct = vcalloc( tsize+1, sizeof( int ) );
    }

    cp = ct;
    while((point = *pointt++ )!= END_ARRAY){
        tmp = memo[point]++;
        if(tmp < table[point]){
            value++;
        }
        if(tmp == 0){
          *cp++ = point;
        }
    }
    *cp = END_ARRAY;
    cp =  ct;
    while(*cp != END_ARRAY){
            memo[*cp++] = 0;
    }

    return( value );
}

char** make_group_aa_mafft(int *ngroup, char *mode){
/*mode:         indicates which matrix will be used for the grouping*/
/*n_group:      pointer to the number of groups                     */
/*return value: an array of strings containing the AA of each group */
 
    char **group_list;
    char *matrix_name;
    int extend=0;
    matrix_name=vcalloc(100, sizeof (char));

    if(ngroup[0]==-1) extend=1;

    ngroup[0]=0;
    group_list=declare_char(100, 27);

    if(extend){
        sprintf ( group_list[ngroup[0]++], "gG");
        sprintf ( group_list[ngroup[0]++], "pP");
        sprintf ( group_list[ngroup[0]++], "aA");
        sprintf ( group_list[ngroup[0]++], "cC");
        sprintf ( group_list[ngroup[0]++], "dD");
        sprintf ( group_list[ngroup[0]++], "eE");

        sprintf ( group_list[ngroup[0]++], "fF");
        sprintf ( group_list[ngroup[0]++], "hH");
        sprintf ( group_list[ngroup[0]++], "iI");
        sprintf ( group_list[ngroup[0]++], "kK");
        sprintf ( group_list[ngroup[0]++], "lL");
        sprintf ( group_list[ngroup[0]++], "mM");
        sprintf ( group_list[ngroup[0]++], "nN");
        sprintf ( group_list[ngroup[0]++], "qQ");
        sprintf ( group_list[ngroup[0]++], "rR");

        sprintf ( group_list[ngroup[0]++], "sS");
        sprintf ( group_list[ngroup[0]++], "tT");
        sprintf ( group_list[ngroup[0]++], "vV");
        sprintf ( group_list[ngroup[0]++], "wW");
        sprintf ( group_list[ngroup[0]++], "*");
    }

    if(mode && mode[0]=='_'){ mode++; sprintf(matrix_name, "%s", mode); }

    if(strm(mode, "mafft")){
        sprintf ( group_list[ngroup[0]++],"agjopstAGJOPST");
        sprintf ( group_list[ngroup[0]++],"ilmvILMV");
        sprintf ( group_list[ngroup[0]++],"bdenqzBDENQZ");
        sprintf ( group_list[ngroup[0]++],"hkrHKR");
        sprintf ( group_list[ngroup[0]++],"fwyFWY");
        sprintf ( group_list[ngroup[0]++],"cC");
        vfree (matrix_name);
        return group_list;
    }
    else{
        fprintf(stderr, "ERROR - Incorrect Mode\n");
        exit(EXIT_FAILURE);
    }
}
