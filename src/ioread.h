#ifndef __INC_IOREAD_
#define __INC_IOREAD_

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "useful.h"


char** unconstrained_read(char* filename, int* line_count_final);
void copy_char_arrays(char** dest, char** src, int lx, int ly);
void free_char_array(char** array, int ly);


#endif
