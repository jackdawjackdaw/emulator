#ifndef __EMULATOR_INC_
#define __EMULATOR_INC_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h> 
#include <gsl/gsl_sf.h>

#include "regression.h"

/**
 * the fn ptr to the covariance function, this is the most called function in libEmu
 * you can change this when you setup the optstruct.
 */
double (*covariance_fn)(gsl_vector*, gsl_vector*, gsl_vector*, int, int);

void print_matrix(gsl_matrix* m, int nx, int ny);

double covariance_fn_gaussian(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams);

void derivative_l_gauss(gsl_matrix *dCdTheta, gsl_matrix* xmodel, 
												double thetaLength, int index, int nmodel_points, int nparams);

// same as above but without clamping on small values
double covariance_fn_gaussian_exact(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams);

double covariance_fn_matern_three(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams);

void derivative_l_matern_three(gsl_matrix *dCdTheta, gsl_matrix* xmodel, double thetaLength,
															  int index, int nmodel_points, int nparams);

double covariance_fn_matern_five(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams);

void derivative_l_matern_five(gsl_matrix *dCdTheta, gsl_matrix* xmodel, double thetaLength,
															int index, int nmodel_points, int nparams);


void makeKVector(gsl_vector* kvector, gsl_matrix *xmodel, gsl_vector *xnew, gsl_vector* thetas, int nmodel_points, int nthetas, int nparams);

void makeKVector_fnptr(gsl_vector* kvector, gsl_matrix *xmodel, gsl_vector *xnew, gsl_vector* thetas, int nmodel_points, int nthetas, int nparams,
											 double (*covariance_fn_ptr)(gsl_vector*, gsl_vector*, gsl_vector*, int, int));


void makeCovMatrix(gsl_matrix *cov_matrix, gsl_matrix *xmodel, gsl_vector* thetas, int nmodel_points, int nthetas, int nparams);

void makeCovMatrix_fnptr(gsl_matrix *cov_matrix, gsl_matrix *xmodel, gsl_vector* thetas, int nmodel_points, int nthetas, int nparams, 
												 double (*covariance_fn_ptr)(gsl_vector*, gsl_vector*, gsl_vector*, int, int));


double makeEmulatedMean(gsl_matrix *inverse_cov_matrix, gsl_vector *training_vector, gsl_vector *kplus_vector, gsl_vector* h_vector, gsl_matrix* h_matrix, gsl_vector* beta_vector,  int nmodel_points);

double makeEmulatedVariance(gsl_matrix *inverse_cov_matrix, gsl_vector *kplus_vector, gsl_vector *h_vector, gsl_matrix *h_matrix, double kappa, int nmodel_points, int nregression_fns);



void initialise_new_x(gsl_matrix* new_x, int nparams, int nemulate_points, double emulate_min, double emulate_max);

/* deprecated */
/* double covariance_fn_gaussian_nondiag(gsl_vector* xm, gsl_vector*xn, gsl_vector*thetas, int nthetas, int nparams); */
/* double covariance_fn_matern(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams); */

#endif
