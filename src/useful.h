/* 11-1-09
 * CCS
 * useful.h, some basic functions one always uses in scientific code
 */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "gsl/gsl_rng.h"
#include "math.h"
#include "sys/time.h"

//*******************
// memory management
//*******************
double* alloc_2d_double_array(int size_x, int size_y, double init_val);
int* alloc_2d_int_array(int size_x, int size_y, int init_val);
void print_double_array(double* the_array, int n_x, int n_y);
void print_int_array(int* the_array, int n_x, int n_y);

//******************
// tries to get a seed from /dev/random or otherwise uses the system time
//*****************
unsigned long int get_seed(void);

//***************
// binary io
// some low level stuff, snarfed from treeio.c
//**************
void in_int(FILE *fptr, int *);
// read in a single double
void in_double(FILE *fptr, double *);
// read in a vector of doubles
void in_vector(FILE *fptr, double *);
// read in N doubles
void in_blob(FILE *fptr, int n_things, double *vec);

void out_int(FILE *fptr, int);
void out_double(FILE *fptr, double);
void out_vector(FILE *fptr, double*);

// this will change the way that out vec and in vec work
#define NDIM 3
