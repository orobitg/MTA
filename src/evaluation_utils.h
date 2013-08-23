/* 
 * File:   evaluation_utils.h
 * Author: oro
 *
 * Created on 14 / febrer / 2013, 17:10
 */

#ifndef EVALUATION_UTILS_H
#define	EVALUATION_UTILS_H

#include "utils.h"
#include "parameters_utils.h"

double calculate_score(Parameters *P, char *alnname, int ntree, double *sc_list);
void ga_completing_scores(Parameters *P, double **score_table, double *score_list);
double ga_wsc(double *score_table);
double ga_mcc1(double **score_table, int tree, int ntree);
double ga_mcc2(double **score_table, int tree, int ntree);

#endif	/* EVALUATION_UTILS_H */

