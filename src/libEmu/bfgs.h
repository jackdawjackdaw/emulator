#ifndef __INC_BFGS__
#define __INC_BFGS__

//! A C implementation of the BFGS algorithm as sketched out in mathematica
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "assert.h"

/** 
 * This should provide an additional maximisation method for libEmu
 */

/**
 * function pointers
 * float fn (float a, float b)
 * 
 * void foo( float (*ptr)(float, float)){
 * result = ptr(a, b);
 * }
 *
 * foo(&fn);
 */


/* could use these macros to provide the functiontypes
 * but it feels like this would rather obscure the methods 
 *
 * #define FNTYPE double (*fn)(gsl_vector*, int)
 * #define GRADTYPE void (*gradientFn)( double (*fn)(gsl_vector*, int), gsl_vector*, gsl_vector* ,int)
 */


/* test functions */
double fRosenbrock( gsl_vector* x, int nparams); 
double fSphere(gsl_vector *x, int nparams);
double fPowers( gsl_vector*x, int nparams);

/* actual bfgs functions */
void getGradientNumeric(double(*fn)(gsl_vector*, int), gsl_vector* xk, gsl_vector* gradient, int nparams);

void obtainStep(gsl_vector* step, double (*fn)(gsl_vector*, int),
								void (*gradientFn)( double (*fn)(gsl_vector*, int), gsl_vector*, gsl_vector* ,int),	\
								gsl_matrix* ranges,\
								gsl_vector* xk, int nparams, gsl_matrix* bkInv);

// can use macros void obtainStep(gsl_vector* step, FNTYPE, GRADTYPE,gsl_vector* xk, int nparams, gsl_matrix* bkInv);


void getYk(gsl_vector* yk, double(*fn)(gsl_vector*, int),
					 void (*gradientFn)( double (*fn)(gsl_vector*, int), gsl_vector*, gsl_vector* ,int), 
					 gsl_vector* xk, gsl_vector* xk1, int nparams);


void getNewBInverse(gsl_matrix *bprev,  gsl_vector *s, gsl_vector *y, int nparams);
void getNewB(gsl_matrix* b, gsl_vector* s, gsl_vector* y, int nparams);

double lineSearch(double(*fn)(gsl_vector*, int),\
					 void (*gradientFn)( double (*fn)(gsl_vector*, int), gsl_vector*, gsl_vector*, int), \
									gsl_vector* direction, gsl_vector* position, int nparams);


int armGold(double(*fn)(gsl_vector*, int), 
						void (*gradientFn)( double (*fn)(gsl_vector*, int), gsl_vector*, gsl_vector* ,int), 
						gsl_vector* direction, gsl_vector* position, double a, int nparams);

void doSimpleBFGS( double(*fn)(gsl_vector*, int),
									 void (*gradientFn)( double (*fn)(gsl_vector*, int), gsl_vector*, gsl_vector* ,int), 
									 gsl_matrix* ranges,
									 gsl_vector* xkInit, gsl_vector* xkFinal, gsl_matrix* Binit, int nparams, int nsteps);


int test_range_vector( gsl_vector *x, gsl_matrix *ranges, int nparams);

#endif
