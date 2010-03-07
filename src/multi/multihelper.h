#ifndef __INC_MULTIHELPER__
#define __INC_MULTIHELPER__

#include "multifit.h"

void dump_result(emuResult *res, FILE *fptr);
void alloc_emuRes(emuResult *thing, eopts *options);
void free_eopts(eopts* options);
void free_emuRes(emuResult *thing);
void print_splits(gsl_matrix* splits, int n);

#endif	

