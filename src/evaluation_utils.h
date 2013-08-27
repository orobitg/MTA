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

#endif	/* EVALUATION_UTILS_H */

