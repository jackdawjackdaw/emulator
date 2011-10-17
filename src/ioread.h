#ifndef __INC_IOREAD_
#define __INC_IOREAD_

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"


/**
 * @file ioread.h
 * \brief fns for reading simple tab/space separated input from stdin or a file 
 * 
 * leave this alone
 */

char** unconstrained_read(char* filename, int* line_count_final);
void copy_char_arrays(char** dest, char** src, int lx, int ly);
void free_char_array(char** array, int ly);


#endif
