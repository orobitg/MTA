#include "sequence_utils.h"

Sequence *declare_sequence(int nseq, int max_len, int min_len){
    
    Sequence *LS;

    LS=vcalloc(1, sizeof(Sequence));
    LS->seq=declare_char(nseq, max_len+1);
    LS->seq_names=declare_char(nseq,MAXNAMES+1);
    LS->type=vcalloc(30, sizeof (char));
    LS->seq_len=vcalloc(nseq, sizeof (int));
    LS->seq_id=vcalloc(nseq, sizeof (int));

    LS->max_len=max_len;
    LS->min_len=min_len;
    LS->nseq=nseq;

    return LS;
}

void free_sequence(Sequence *LS){

    if (!LS) return;
    
    free_char(LS->seq, -1);
    free_char(LS->seq_names,-1);
    vfree(LS->seq_len);
    vfree(LS->type);
    vfree(LS->seq_id);
    vfree(LS);
}

Sequence *read_fasta_sequences(char *fname){

    FILE *fp;
    Sequence *S;
    int nseq=0, max_len_seq=0, min_len_seq=0, max=0, coor=0, current=0;
    int l=0, a=0, p=0, i=0, clen=0;
    int c, nul;
    char *sub, *name;

    fprintf(stdout, "\n---> Reading the input sequences (Fasta format)\n");
    nseq=count_nseqs(fname, '>');
    fprintf(stdout, "\tN Seqs: %d\n", nseq);

    if (nseq==0){
        return NULL;
    }

    min_len_seq=max=count_n_char_in_file(fname);
    sub = (char *) vcalloc(max+1, sizeof (char));
    name = (char *) vcalloc(10000, sizeof(char));
    fp=openfile(fname, "r");

    c=fgetc(fp);
    while (c!=EOF){
        if (c=='>'){
            nul=fscanf_seq_name (fp, name);
            while((c=fgetc(fp))!='\n' && c!=EOF);
            while((c=fgetc(fp))!='>' && c!=EOF){
                if(isalnum(c) || is_gap(c)){
                    sub[clen++]=c;
                }
            }
            max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
            min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
            clen=0;
        }
        else {
            c=fgetc (fp);
        }
    }
    fclose (fp);
    vfree(name);
    vfree(sub);

    S = declare_sequence(nseq, max_len_seq, min_len_seq);

    fp=openfile(fname,"r");
    c=fgetc(fp);
    coor++;

    while (c!=EOF){
        if (c=='>'){
            coor+=fscanf_seq_name(fp, S->seq_names[current]);
            l=strlen (S->seq_names[current]);

            if (S->seq_names[current][l-1]==','|| S->seq_names[current][l-1]==';'){
                S->seq_names[current][l-1]='\0';
            }
            //seq_names[current]=translate_name(seq_names[current]);
            a=0;
            while ((c=fgetc(fp))!='\n' && c!=EOF && a<(COMMENT_SIZE-1)){
                coor++;
            }
            coor++;
            p=0;
            while ((c=fgetc(fp))!='>' && c!=EOF){
                coor++;
                if (!isspace(c)){
                    S->seq[current][p++]=c;
                }
            }
            coor++;
            S->seq[current][p]='\0';
            S->seq_len[current]=p;
            S->seq_id[current]=current;
            current++;
        }
        else {
            c=fgetc(fp);
            coor++;
        }

    }
    
    S=get_sequence_type(S);
    fprintf(stdout, "\tSequences: (%s)\n", S->type);
    for(i=0; i < S->nseq; i++){
        fprintf(stdout, "\t\tSeq_id: %d - %s\n", S->seq_id[i], S->seq_names[i]);
    }
    fprintf(stdout, "---> DONE\n");

    return S;
}

int count_nseqs(char *fname, char x){
    FILE *fp;
    int n, c;

    n=0;
    fp=openfile(fname, "r");
    while((c=fgetc(fp))!=EOF){
        n+=(c==x);
    }
    fclose (fp);
    return n;
}

int count_n_char_in_file(char *name){
      int  c, n;
      FILE *fp;

      n=0;
      fp=openfile(name, "r");
      while((c=fgetc(fp))!=EOF){
          n++;
      }
      fclose (fp);
      return n;
  }

int fscanf_seq_name ( FILE *fp, char *sname){

    static char *name;
    int r, nul;

    if (!sname){
        return 0;
    }
    if (!name){
        name = (char *) calloc(10000, sizeof (char));
    }
    nul=fscanf (fp, "%s", name);
    r=strlen (name);
    if(strlen(name)>MAXNAMES){
        fprintf(stderr, "\nWARNING: Seq Name Too long: [%s]. Truncated to %d", name, MAXNAMES);
    }
    name[MAXNAMES]='\0';
    sprintf(sname, "%s", name);

    return r;
}



Sequence *get_sequence_type(Sequence *S){

    if(!S){
        return NULL;
    }
    else{
        sprintf(S->type, "%s", get_array_type (S->nseq, S->seq));
    }

    return S;
}

char * get_array_type (int n, char **seq){

    char *buf, *buf2;
    int a, tot=0;
    buf2=vcalloc ( 100, sizeof (char));

    for(tot=0,a=0; a<n; a++){
        tot+=strlen (seq[a]);
    }
    buf=vcalloc (tot+1, sizeof (char));
    for ( a=0; a<n; a++){
        strcat (buf, seq[a]);
    }
    sprintf(buf2, "%s", get_string_type(buf));

    vfree (buf);
    return buf2;
}

char *get_string_type(char *S){

    int a, l;
    int protein=0,  dna=0,rna=0, tot=0;
    char *type;
    static char *ltype;
    static int warning;

    if(!ltype){
      declare_name(ltype);
    }

    declare_name(type);
    l=(int)strlen (S);

    if (l==0){
        sprintf ( type, "UNKNOWN");
        return type;
    }

    for(a=0; a<l; a++){
        if(!is_gap(S[a])){
            protein+=(is_aa (S[a]) && !is_dna(S[a]));
            dna+=(is_dna(S[a]));
            rna+=(is_rna(S[a]));
            tot++;
        }
    }

    protein=(protein*100)/tot;
    dna=(dna*100)/tot;
    rna=(rna*100)/tot;

    if(l<20 && warning==0){
        /*add_warning ( stderr, "WARNING: short sequences, use -type=DNA or PROTEIN");*/
        warning=1;
    }
    if(l<20 && ltype && ltype[0]){
        if(dna==100){
            sprintf ( type, "%s", ltype);
            ltype[0]='\0';
        }
        else
            sprintf ( type, "PROTEIN");
        }
    else if(dna>98 && rna>0){
        sprintf ( type, "RNA");
    }
    else if(dna>98){
        sprintf ( type, "DNA");
    }
    else{
        sprintf ( type, "PROTEIN");
    }
    sprintf ( ltype, "%s", type);

    return type;
}

int is_aa(char x){
    return (is_in_set (x, AA_ALPHABET) && !is_in_set (x, GAP_LIST));
}

int is_rna(char x){
    return (is_in_set (x, RNAONLY_ALPHABET)&& !is_in_set (x, GAP_LIST));
}

int is_dna(char x){
    return (is_in_set (x, DNA_ALPHABET)&& !is_in_set (x, GAP_LIST));
}

int is_gap(char x){
    return (is_in_set( x, GAP_LIST));
}

int is_in_set(char r, char *list){

    char s[2];
    s[0]=r;

    if (strstr(list, s)!=NULL){
        return 1;
    }
    else {
        return 0;
    }
}