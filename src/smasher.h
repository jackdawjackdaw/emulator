#ifndef __INC_SMASHER__
#define __INC_SMASHER__

#include "stdio.h"
#include "multifit.h"


int smasher(gsl_matrix **split_ranges, int* nsplits, eopts* toplevel, int min_points, gsl_rng* random_number);
int compare_regions(const region* a, const region* b);
void extend_regions(gsl_matrix* split_ranges, int number_splits, eopts* toplevel, int dim);
#endif
