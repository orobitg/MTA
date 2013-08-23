
#include <stdio.h>
#include "evaluation_utils.h"
#include "utils.h"

/*Evaluate an alignment with a chosen score*/
double calculate_score(Parameters *P, char *alnname, int ntree, double *sc_list){

    FILE *fp;
    int ret;
    char *score_method, *scorefile, *libname;
    
    sc_list[0]=0.0;
    sc_list[1]=0.0;
    
    score_method = (char *) vcalloc(ALLPATH, sizeof(char));
    scorefile = (char *) vcalloc(FILENAMELEN, sizeof(char));
       
    if(strm(P->score_method, "sp")){
        if((P->PB)->tcoffee_bin != NULL){
            sprintf(scorefile, "%s/%s_%d.sp", TMP, (P->F)->name, ntree);
            sprintf(score_method, "%s -other_pg fastal -i %s --eval_aln -g %f -e %f -a --mat %s | grep Score: | cut -d' ' -f2 > %s", (P->PB)->tcoffee_bin, alnname, P->gop, P->gep, P->mat, scorefile);
        }
        else {
            fprintf(stderr, "Error - SP Score not found. Wrong T-Coffee path\n");
            exit(EXIT_FAILURE);
        }
    }
    else if(strm(P->score_method, "normd")){
        if((P->PB)->normd_bin != NULL){
            sprintf(scorefile, "%s/%s_%d.normd", TMP, (P->F)->name, ntree);
            sprintf(score_method, "%s %s > %s", (P->PB)->normd_bin, alnname, scorefile);
        }
        else {
            fprintf(stderr, "Error - NorMD Score not found. Wrong NorMD path\n");
            exit(EXIT_FAILURE);
        }
    }
    else if(strm(P->score_method, "triplet")){ //ARREGLAR
        if((P->PB)->tc_score_bin != NULL){
            sprintf(scorefile, "%s/%s_%d.triplet", TMP, (P->F)->name, ntree);
            libname = (char *) vcalloc(FILENAMELEN, sizeof(char));
/*
        if(P->mpi_mode == 0){    
            sprintf(libname, "%s%s_ev.lib", P->outdir, (P->F)->name);
            if(ntree == 0){
                sprintf(score_method, "%s%s -infile=%s -output=html -out_lib %s -outfile=%s/%s_%d.html > /dev/null 2>&1", MT_HOME, TCOFFEE_BIN, alnname, libname, TMP, (P->F)->name, ntree);
                ret = system(score_method);
                sprintf(score_method, "%s%s -infile=%s -output=html -lib %s -evaluate_mode=triplet -score -outfile=%s/%s_%d.html > /dev/null 2>&1", MT_HOME, TCOFFEE_BIN, alnname, libname, TMP, (P->F)->name, ntree);
            }
            else {
                sprintf(score_method, "%s%s -infile=%s -output=html -lib %s -evaluate_mode=triplet -score -outfile=%s/%s_%d.html > /dev/null 2>&1", MT_HOME, TCOFFEE_BIN, alnname, libname, TMP, (P->F)->name, ntree);
            }
        }
*/
       // else {
        
                sprintf(score_method, "%s -infile=%s -output=html -evaluate_mode=triplet -n_core=1 -score -outfile=%s/%s_%d.html 2> /dev/null | grep 'TRIPLET=' | cut -d'=' -f2 > %s",  (P->PB)->tc_score_bin, alnname, TMP, (P->F)->name, ntree, scorefile);
       // }
        }
        else {
            fprintf(stderr, "Error - Triple Score not found. Wrong MTA_HOME path\n");
            exit(EXIT_FAILURE);
        }
  
    }
    else if(strm(P->score_method, "coffee")){ //ARREGLAR
        if((P->PB)->tc_score_bin != NULL){
            sprintf(scorefile, "%s/%s_%d.coffee", TMP, (P->F)->name, ntree);
            libname = (char *) vcalloc(FILENAMELEN, sizeof(char));

/*
        if(P->mpi_mode == 0){    
            sprintf(libname, "%s%s_ev.lib", P->outdir, (P->F)->name);
            if(ntree == 0){
                sprintf(score_method, "%s%s -infile=%s -output=html -out_lib %s -outfile=%s/%s_%d.html > /dev/null 2>&1", MT_HOME, TCOFFEE_BIN, alnname, libname, TMP, (P->F)->name, ntree);
                ret = system(score_method);
                sprintf(score_method, "%s%s -infile=%s -output=html -lib %s -evaluate_mode=fast -score -outfile=%s/%s_%d.html > /dev/null 2>&1", MT_HOME, TCOFFEE_BIN, alnname, libname, TMP, (P->F)->name, ntree);
            }
            else {
                sprintf(score_method, "%s%s -infile=%s -output=html -lib %s -evaluate_mode=fast -score -outfile=%s/%s_%d.html > /dev/null 2>&1", MT_HOME, TCOFFEE_BIN, alnname, libname, TMP, (P->F)->name, ntree);
            }
        }
        else {
*/

                sprintf(score_method, "%s -infile=%s -output=html -evaluate_mode=fast -n_core=1 -score -outfile=%s/%s_%d.html 2> /dev/null | grep 'COFFEE=' | cut -d'=' -f2 > %s", (P->PB)->tc_score_bin, alnname, TMP, (P->F)->name, ntree, scorefile);
      //  }
        }
        else {
            fprintf(stderr, "Error - Coffee Score not found. Wrong MTA_HOME path\n");
            exit(EXIT_FAILURE);
        }
  
    }
    else if(strm(P->score_method, "strike")){
        if((P->PB)->strike_bin != NULL){
                sprintf(scorefile, "%s/%s_%d.strike", TMP, (P->F)->name, ntree);
                sprintf(score_method, "%s -a %s -c %s | tail -1 > %s", (P->PB)->strike_bin, alnname, P->str_file, scorefile);
        }
        else {
            fprintf(stderr, "Error - STRIKE Score not found. Correct STRIKE path\n");
            exit(EXIT_FAILURE);
        }

    }
    else if(strm(P->score_method, "irmsdt")){
        if((P->PB)->tcoffee_bin != NULL){
            sprintf(scorefile, "%s/%s_%d.irmsd", TMP, (P->F)->name, ntree);        
            sprintf(score_method, "%s -other_pg irmsd -aln %s -template_file %s 2>/dev/null | grep 'TOTAL' | grep 'NiRMSD' | cut -d':' -f2 | cut -d'A' -f1 | sed 's/ //g' > %s", (P->PB)->tcoffee_bin, alnname, P->irmsd_file, scorefile);
        }
        else {
            fprintf(stderr, "Error - iRMSDT Score not found. Wrong T-Coffee path\n");
            exit(EXIT_FAILURE);
        }
    }
    else if(strm(P->score_method, "irmsd")){
        if((P->PB)->tcoffee_bin != NULL){
            sprintf(scorefile, "%s/%s_%d.irmsd", TMP, (P->F)->name, ntree);        
            sprintf(score_method, "%s -other_pg irmsd -aln %s 2>/dev/null | grep 'TOTAL' | grep 'NiRMSD' | cut -d':' -f2 | cut -d'A' -f1 | sed 's/ //g' > %s", (P->PB)->tcoffee_bin, alnname, scorefile);
        //arreglar lo dels espais sed
        }
        else {
            fprintf(stderr, "Error - iRMSDT Score not found. Wrong T-Coffee path\n");
            exit(EXIT_FAILURE);
        }
    }
    else if(strm(P->score_method, "mcc") || strm(P->score_method, "wsc") || strm(P->score_method, "mcc_irmsd") || strm(P->score_method, "wsc_irmsd")){
        
        sprintf(score_method, "%s", P->score_method);
        
        sprintf(P->score_method, "sp");
        sprintf(P->mat, "blosum62mt");
        P->gop = -11;
        P->gep = -1;
        sc_list[2] = calculate_score(P, alnname, ntree, sc_list);

        sprintf(P->mat, "blosum62mt");
        P->gop = -6.6;
        P->gep = -0.9;
        sc_list[3] = calculate_score(P, alnname, ntree, sc_list);
        
        sprintf(P->mat, "pam250mt"); 
        P->gop = -13.8;
        P->gep = -0.2;
        sc_list[4] = calculate_score(P, alnname, ntree, sc_list);    
        
        sprintf(P->score_method, "normd");  
        sc_list[5] = calculate_score(P, alnname, ntree, sc_list);
        
        sprintf(P->score_method, "strike");     
        sc_list[6] = calculate_score(P, alnname, ntree, sc_list);
        
        sprintf(P->score_method, "triplet");  
        sc_list[7]=calculate_score(P, alnname, ntree, sc_list);
        if(strm(P->score_method, "mcc") || strm(P->score_method, "wsc")){
            sprintf(P->score_method, "coffee"); 
            sc_list[8]=calculate_score(P, alnname, ntree, sc_list);
        }
        else if(strm(P->score_method, "mcc_irmsd") || strm(P->score_method, "wsc_irmsd")){
            sprintf(P->score_method, "irmsdt");//2 parametres per tipus de fitxer
            sc_list[8]=calculate_score(P, alnname, ntree, sc_list);
        }

        sprintf(P->score_method, "%s", score_method);
    }
    else {
        fprintf(stderr, "ERROR - BAD SCORE METHOD\n");
    }

    if(!strm(P->score_method, "mcc") && !strm(P->score_method, "wsc") || !strm(P->score_method, "mcc_irmsd") || !strm(P->score_method, "wsc_irmsd")){
        ret = system(score_method);
        fp=openfile(scorefile, "r");
        ret=fscanf(fp, "%lf", &sc_list[1]);
        if(strm(P->score_method, "irmsdt") || strm(P->score_method, "irmsd")){
            if(sc_list[1]!=0){
                sc_list[1] = 1/sc_list[1];
            }
            else {
                sc_list[1]=0.0;
            }
        }
        remove_file(scorefile);
        fclose(fp);
    }

    if(strm(P->score_method, "triplet") || strm(P->score_method, "coffee")){
        sprintf(score_method, "%s/%s_%d.html", TMP, (P->F)->name, ntree);
        remove_file(score_method);
/*
        if(ntree == P->ntree-1 && P->mpi_mode == 0){
            sprintf(score_method, "rm %s", libname);
            ret = system(score_method);
        }
*/
        vfree(libname);
    }
    
    sc_list[0] = (double) ntree;
    //
    vfree(score_method);
    vfree(scorefile);

    return sc_list[1];
}

/*Score Genètics*/
void ga_completing_scores(Parameters *P, double **score_table, double *score_list){
    
    int i=0, j=0;
    double temp=0.0, interval=0.0;
    double *maxscore, *minscore;
    
    maxscore = (double *) vcalloc(NSCORE, sizeof(double));
    minscore = (double *) vcalloc(NSCORE, sizeof(double));

    for(i=0; i<NSCORE; i++){
        maxscore[i] = -999999.99999;
        minscore[i] = 999999.99999;
    }
        
    for(i=0; i<NSCORE; i++){
        for(j=0; j<P->ntree; j++){
            if(score_table[j][i] != 0.0){
                if(maxscore[i] < score_table[j][i]){
                    maxscore[i] = score_table[j][i];
                }
                if(minscore[i] > score_table[j][i]){
                    minscore[i] = score_table[j][i];
                }
            }
        }
    }
    /*Normalizar*/
    for(j=0; j<P->ntree; j++){
        for(i=0; i<NSCORE; i++){
            temp = score_table[j][i];
            interval = maxscore[i] - minscore[i];
            if(temp!=0){
                score_table[j][i] = (temp - minscore[i]) / interval;
            }
            else {
                score_table[j][i]=0;
            }
        }
    }
    
    for(i=0; i<NSCORE; i++){
        maxscore[i] = -1.0;
        minscore[i] = 2.0;
    }

    for(i=0; i<NSCORE; i++){      
        for(j=0; j<P->ntree; j++){
            score_table[P->ntree][i] = score_table[P->ntree][i] + score_table[j][i];
               // printf("%f ", data[i].score_list[ntree][k]);
                if(maxscore[i] < score_table[j][i]){
                    maxscore[i] = score_table[j][i];
                }
                if(minscore[i] > score_table[j][i]){
                    minscore[i] = score_table[j][i];
                }
        }
        score_table[P->ntree+3][0] = score_table[P->ntree+3][0] + score_table[P->ntree][i];
    }
              
    double maxgscore=-1.0, mingscore=1.0;
    for(i=0; i<NSCORE; i++){      
       //printf("%f ", data[i].score_list[ntree][k]);
        score_table[P->ntree][i] = score_table[P->ntree][i] / (float) P->ntree;
        score_table[P->ntree+1][i] = maxscore[i];
        score_table[P->ntree+2][i] = minscore[i];  

        if(maxscore[i] > maxgscore){
            maxgscore = maxscore[i];
        }
        if(minscore[i] < mingscore){
            mingscore = minscore[i];
        }
    }
    score_table[P->ntree+3][0] = score_table[P->ntree+3][0] / ((float) P->ntree * (NSCORE));
    score_table[P->ntree+3][1] = maxgscore;
    score_table[P->ntree+3][2] = mingscore;   
    
    if(strm(P->score_method, "wsc")){
        for(j=0; j<P->ntree; j++){
            score_list[j]=ga_wsc(score_table[j]);
        }
    }
    else if(P->score_method, "mcc"){
        for(j=0; j<P->ntree; j++){
            score_list[j]=ga_mcc2(score_table, j, P->ntree);
        }
    }
  
    vfree(maxscore);
    vfree(minscore);
}

double ga_wsc(double *score_table){
    
    double score;
    score = (0.14 * score_table[0]) + (0.07 * score_table[1]) + (0.07 * score_table[2]) + (0 * score_table[3]) + (0.64 * score_table[4]) + (0.07 * score_table[5]) + (0 * score_table[6]);
    return  score;    
}

double ga_wsc_clw_irmsd(double *score_table){
    
    double score;
    score = (0.01 * score_table[0]) + (0.14 * score_table[1]) + (0.01 * score_table[2]) + (0 * score_table[3]) + (0.37 * score_table[4]) + (0.01 * score_table[5]) + (98 * score_table[6]);
    return  score;    
}

double ga_mcc1(double **score_table, int tree, int ntree){
 
    double score;
    
    if(score_table[tree][3] > 54 && score_table[ntree+3][0] > 50 || score_table[tree][4] != 149 && score_table[tree][1] > 50 || 10 > score_table[ntree+1][0])
            score += 100;
    if(score_table[tree][3] != 42 /*&& getStatisticNth(exp, 6, 1) < 168*/ || score_table[tree][5] != score_table[tree][3] || score_table[tree][5] < score_table[tree][1] || score_table[tree][3] < 16)
            score -= score_table[tree][2];
    if(score_table[tree][5] < score_table[ntree+3][0] || score_table[tree][0] != score_table[ntree+3][0])
            score += score_table[tree][4] + score_table[tree][0] ;
    if(score_table[tree][3] == 2 && score_table[tree][1] != score_table[tree][3])
            score += 91 / 2 ;
    if(93 != score_table[tree][5])
            score += score_table[tree][3];
    if(score_table[tree][0] > score_table[tree][4])
            score -= score_table[tree][4];
    if(score_table[tree][4] > score_table[ntree][4])
            score *= score_table[tree][4];
    if(117 > score_table[ntree][1])
            score += score_table[tree][1];
       
    return score;
}

double ga_mcc2(double **score_table, int tree, int ntree){
 
    double score;
    
    if(120 > score_table[ntree+3][2] && score_table[tree][6] < score_table[tree][3] || score_table[tree][0] > score_table[tree][2] || score_table[ntree+2][1] != score_table[ntree+3][2] && 220 > score_table[tree][1])
        score += 154 ;
    if( score_table[tree][3] > 89 || 56 != score_table[tree][6] && score_table[tree][1] < score_table[tree][4] || score_table[tree][3] != score_table[tree][6] || score_table[tree][4] < score_table[ntree][5])
        score -= 100 ;
    if( score_table[ntree+1][5] > score_table[tree][0] || score_table[ntree][5] < score_table[ntree][3])
        score -= score_table[ntree][4] - 6 ;
    if( score_table[ntree][6] == 74 && score_table[tree][4] < score_table[tree][0])
        score += score_table[tree][3] * score_table[ntree][4] ;
    if( score_table[ntree][6] > 4)
        score += score_table[tree][1] ;
    if( 189 > score_table[tree][3])
        score += score_table[tree][4] ;
    if( score_table[tree][1] > score_table[tree][3])
        score += score_table[tree][5] ;
    if( score_table[tree][0] != score_table[ntree][0])
        score -= score_table[tree][2] ;
       
    return score;
}


