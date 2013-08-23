/* 
 * File:   sequence_utils.h
 * Author: oro
 *
 * Created on 30 / desembre / 2010, 16:10
 */

#ifndef SEQUENCE_UTILS_H
#define	SEQUENCE_UTILS_H

#include "utils.h"

struct Sequence {
    char **seq;
    char **seq_names;
    int *seq_len;
    int *seq_id;
    int nseq;
    int min_len;
    int max_len;
    char *type;
};
typedef struct Sequence Sequence;

Sequence *declare_sequence(int nseq, int max_len, int min_len);
void free_sequence(Sequence *LS);

Sequence *read_fasta_sequences(char *fname);
int count_nseqs(char *fname, char x);
int count_n_char_in_file(char *name);
int fscanf_seq_name ( FILE *fp, char *sname);

Sequence *get_sequence_type(Sequence *S);
char * get_array_type (int n, char **seq);
char *get_string_type(char *S);

int is_aa(char x);
int is_rna(char x);
int is_dna(char x);
int is_gap( char x);
int is_in_set ( char r, char *list);
#endif	/* SEQUENCE_UTILS_H */

