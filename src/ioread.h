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

void scan_int(FILE * fp, int * ip);
void scan_double(FILE * fp, double * dp);
#define print_int(FILE_ptr, an_int) (fprintf((FILE_ptr),"%d\n",(an_int)))
#define print_double(FILE_ptr, a_double) (fprintf((FILE_ptr), "%.17f\n",(a_double)))

#endif
