#ifndef __INC_LBFGS_
#define __INC_LBFGS_

#include "stdio.h"
#include "unistd.h"
#include "math.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_rng.h"
#include "assert.h"
#include "../useful.h"

double fSphere(double *x, int nparams, void* args);
double fPowers( double *x, int nparams, void* args);
double fRosenbrock( double* x, int nparams, void* args);
void copy_gslvec_vec(gsl_vector* gsl_vec, double* vec, int n);
void copy_vec_gslvec(double* vec, gsl_vector* gsl_vec, int n);

void set_bounds(double *vec, gsl_matrix *ranges, int index, int nparams);
void set_nbd(int *vec, int nparams, int nmax);
void set_zero( double *vec, int size);
void print_vec(gsl_vector* x, int n);

void doBoundedBFGS( double(*fn)(double*, int, void*),													
										void(*gradientFn)(double*, double*, int, void*),
										gsl_matrix* ranges, 
										gsl_vector *xkInit, gsl_vector* xkFinal, int nparams, int nsteps, void* args);

extern int setulb_(int* n, int*m, double* x, double* l, double*u, int*nbd, double*f, double* g, \
										double* factr, double*pgtol, double* wa, int*iwa, char* task, int* iprint, char* csave,
										long int *lsave, int* isave, double* dsave, int nstr1, int nstr2);


void getGradientNumericLBFGS(double(*fn)(double*, int, void*), double* xk, double* gradient, int nparams, void* args);

#endif
