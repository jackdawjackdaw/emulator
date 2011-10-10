#ifndef __EMULATOR_INC_
#define __EMULATOR_INC_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h> 
#include <gsl/gsl_sf.h>

#include "regression.h"
//#include "../optstruct.h"
//#include "../modelstruct.h"


/**
 * ccs, added a void* argument to the derivative functions, we need to pass more information in the 
 * case that we're varying alpha as well 
 */

/**
 * the fn ptr to the covariance function, this is the most called function in libEmu
 * you can change this when you setup the optstruct.
 */
double (*covariance_fn)(gsl_vector*, gsl_vector*, gsl_vector*, int, int);

void print_matrix(gsl_matrix* m, int nx, int ny);

double covariance_fn_gaussian(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams);

void derivative_l_gauss(gsl_matrix *dCdTheta, gsl_matrix* xmodel, 
												double thetaLength, int index, int nmodel_points, int nparams, void* args);


// same as covariance_fn_gaussian but with a variable power too
double covariance_fn_gaussian_alpha(gsl_vector *xm, gsl_vector *xn, gsl_vector *thetas, int nthetas, int nparams);

// this function needs the void* args to compute the gradient correctly
void derivative_l_gauss_alpha(gsl_matrix *dCdTheta, gsl_matrix* xmodel, 
															double thetaLength, int index, int nmodel_points, int nparams, void* args);


// a helper, not to be directly called
double powExpGradAlpha(double rtemp, double alpha,  int thetaLength);

double covariance_fn_matern_three(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams);

void derivative_l_matern_three(gsl_matrix *dCdTheta, gsl_matrix* xmodel, double thetaLength,
															 int index, int nmodel_points, int nparams, void* args);

double covariance_fn_matern_five(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams);

void derivative_l_matern_five(gsl_matrix *dCdTheta, gsl_matrix* xmodel, double thetaLength,
															int index, int nmodel_points, int nparams, void* args);


void makeKVector(gsl_vector* kvector, gsl_matrix *xmodel, gsl_vector *xnew, gsl_vector* thetas, int nmodel_points, int nthetas, int nparams);


void makeCovMatrix(gsl_matrix *cov_matrix, gsl_matrix *xmodel, gsl_vector* thetas, int nmodel_points, int nthetas, int nparams);

double makeEmulatedMean(gsl_matrix *inverse_cov_matrix, gsl_vector *training_vector, gsl_vector *kplus_vector, gsl_vector* h_vector, gsl_matrix* h_matrix, gsl_vector* beta_vector,  int nmodel_points);

double makeEmulatedVariance(gsl_matrix *inverse_cov_matrix, gsl_vector *kplus_vector, gsl_vector *h_vector, gsl_matrix *h_matrix, double kappa, int nmodel_points, int nregression_fns);



void initialise_new_x(gsl_matrix* new_x, int nparams, int nemulate_points, double emulate_min, double emulate_max);

/* deprecated */
/* double covariance_fn_gaussian_nondiag(gsl_vector* xm, gsl_vector*xn, gsl_vector*thetas, int nthetas, int nparams); */
/* double covariance_fn_matern(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams); */

#endif
