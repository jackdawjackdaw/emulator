#ifndef __INC_PERSIST__
#define __INC_PERSIST__

#include "stdio.h"
#include "stdlib.h"
#include "emulator.h"
#include "estimator.h"
#include "maximise.h"
#include "multifit.h"
#include "useful.h"

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






#endif
