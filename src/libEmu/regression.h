#ifndef _INC_REGRESSION_
#define _INC_REGRESSION_
/**
 * contains the functions needed to handle the linear regression
 * this entails, creating hvectors (for given x vectors in the parameter space)
 * estimatingBeta, this needs the inverse covariance matrix and the training data
 * creating hmatrix, a done once kind of job, slap together hvectors for each of the design points
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h> 
#include <gsl/gsl_sf.h>

void makeHVector(gsl_vector *h_vector, gsl_vector *x_location, int nparams);
void makeHMatrix(gsl_matrix *h_matrix, gsl_matrix *xmodel, int nmodel_points, int nparams, int nregresion_fns);
void estimateBeta(gsl_vector *beta_vector, gsl_matrix *h_matrix, gsl_matrix* cinverse, gsl_vector *trainingvector, int nmodel_points, int nregression_fns);

void makeHVector_linear(	gsl_vector *h_vector, gsl_vector *x_location, int nparams);
void makeHVector_quadratic( gsl_vector *h_vector, gsl_vector *x_location, int nparams);
void makeHVector_trivial(gsl_vector *h_vector, gsl_vector *x_location, int nparams);
void makeHVector_cubic( gsl_vector *h_vector, gsl_vector *x_location, int nparams);
#endif
