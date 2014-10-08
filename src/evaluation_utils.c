
#include <stdio.h>
#include "evaluation_utils.h"
#include "utils.h"

/*Evaluate an alignment with a chosen score*/
double calculate_score(Parameters *P, char *alnname, int ntree, double *sc_list){

    FILE *fp;
    int ret;
    char *score_method, *scorefile;
    
    sc_list[0]=0.0;
    sc_list[1]=0.0;
    
    score_method = (char *) vcalloc(ALLPATH, sizeof(char));
    scorefile = (char *) vcalloc(FILENAMELEN, sizeof(char));
    sprintf(scorefile, "%s/%s_%d.score", P->workdir, (P->F)->name, ntree);
    
    if(strm(P->score_method, "sp")){
        if((P->PB)->tcoffee_bin != NULL){      
            sprintf(score_method, "%s -other_pg fastal -i %s --eval_aln -g %f -e %f -a --mat %s | grep Score: | cut -d' ' -f2 > %s", (P->PB)->tcoffee_bin, alnname, P->gop, P->gep, P->mat, scorefile);
        }
        else {
            fprintf(stderr, "Error - SP Score not found. Wrong T-Coffee path\n");
            exit(EXIT_FAILURE);
        }
    }
    else if(strm(P->score_method, "normd")){
        if((P->PB)->normd_bin != NULL){
            sprintf(score_method, "%s %s > %s", (P->PB)->normd_bin, alnname, scorefile);
        }
        else {
            fprintf(stderr, "Error - NorMD Score not found. Wrong NorMD path\n");
            exit(EXIT_FAILURE);
        }
    }
    /*else if(strm(P->score_method, "tcs")){ //ARREGLAR
	   if((P->PB)->tc_score_bin != NULL){
		   sprintf(score_method, "%s -infile=%s -method=proba_pair -output=score_ascii -n_core=1 -outfile=%s/%s_%d.ascii 2> /dev/null",  (P->PB)->tc_score_bin, alnname, P->workdir, (P->F)->name, ntree);
	   }
	   else {
		   fprintf(stderr, "Error - Triple Score not found. Wrong MTA_HOME path\n");
		   exit(EXIT_FAILURE);
	   }

    }*/
    /*else if(strm(P->score_method, "triplet")){ //ARREGLAR
        if((P->PB)->tc_score_bin != NULL){
            sprintf(score_method, "%s -infile=%s -output=html -evaluate_mode=triplet -n_core=1 -score -outfile=%s/%s_%d.html 2> /dev/null | grep 'TRIPLET=' | cut -d'=' -f2 > %s",  (P->PB)->tc_score_bin, alnname, TMP, (P->F)->name, ntree, scorefile);
        }
        else {
            fprintf(stderr, "Error - Triple Score not found. Wrong MTA_HOME path\n");
            exit(EXIT_FAILURE);
        }
  
    }
    else if(strm(P->score_method, "coffee")){ //ARREGLAR
        if((P->PB)->tc_score_bin != NULL){
            sprintf(score_method, "%s -infile=%s -output=html -evaluate_mode=fast -n_core=1 -score -outfile=%s/%s_%d.html 2> /dev/null | grep 'COFFEE=' | cut -d'=' -f2 > %s", (P->PB)->tc_score_bin, alnname, TMP, (P->F)->name, ntree, scorefile);
        }
        else {
            fprintf(stderr, "Error - Coffee Score not found. Wrong MTA_HOME path\n");
            exit(EXIT_FAILURE);
        }
  
    }
    else if(strm(P->score_method, "strike")){
        if((P->PB)->strike_bin != NULL){
                sprintf(score_method, "%s -a %s -c %s | tail -1 > %s", (P->PB)->strike_bin, alnname, P->str_file, scorefile);
        }
        else {
            fprintf(stderr, "Error - STRIKE Score not found. Correct STRIKE path\n");
            exit(EXIT_FAILURE);
        }

    }
    else if(strm(P->score_method, "irmsdt")){
        if((P->PB)->tcoffee_bin != NULL){       
            sprintf(score_method, "%s -other_pg irmsd -aln %s -template_file %s 2>/dev/null | grep 'TOTAL' | grep 'NiRMSD' | cut -d':' -f2 | cut -d'A' -f1 | sed 's/ //g' > %s", (P->PB)->tcoffee_bin, alnname, P->irmsd_file, scorefile);
        }
        else {
            fprintf(stderr, "Error - iRMSD Score not found. Wrong T-Coffee path\n");
            exit(EXIT_FAILURE);
        }
    }
    else if(strm(P->score_method, "irmsd")){
        if((P->PB)->tcoffee_bin != NULL){      
            sprintf(score_method, "%s -other_pg irmsd -aln %s 2>/dev/null | grep 'TOTAL' | grep 'NiRMSD' | cut -d':' -f2 | cut -d'A' -f1 | sed 's/ //g' > %s", (P->PB)->tcoffee_bin, alnname, scorefile);
        }
        else {
            fprintf(stderr, "Error - iRMSD Score not found. Wrong T-Coffee path\n");
            exit(EXIT_FAILURE);
        }
    }*/

    ret = system(score_method);
    sc_list[0] = (double) ntree;
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


    if(strm(P->score_method, "triplet") || strm(P->score_method, "coffee")){
        sprintf(score_method, "%s/%s_%d.html", P->workdir, (P->F)->name, ntree);
        remove_file(score_method);
    }
    
    vfree(score_method);
    vfree(scorefile);

    return sc_list[1];
}
