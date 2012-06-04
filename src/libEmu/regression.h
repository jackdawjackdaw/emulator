#ifndef _INC_REGRESSION_
#define _INC_REGRESSION_
/**
 * contains the functions needed to handle the linear regression
 * this entails, creating hvectors (for given x vectors in the parameter space)
 * estimatingBeta, this needs the inverse covariance matrix and the training data
 * creating hmatrix, a done once kind of job, slap together hvectors for each of the design points
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h> 
#include <gsl/gsl_sf.h> 


/**
 * \brief creates a vector of regression coefficients at x_location
 * 
 * create a h_vector evaluated at given x_location., 
 * this currently only creates a set of linear basis fns, you could add another loop which 
 * would add squares etc 
 * 
 * on calling h_vector should be initliased to a vector of length nregression_fns which 
 * in the default case nregression_fns = nparams + 1 
 * 
 * the value of options->nregression_fns is going to be very important for this
 * 
 * this pointer is set to one of the fns in this file by optstruct.c:setup_regression
 */
void (*makeHVector)(gsl_vector *h_vector, gsl_vector *x_location, int nparams);


//void makeHVector(gsl_vector *h_vector, gsl_vector *x_location, int nparams);
void makeHMatrix(gsl_matrix *h_matrix, gsl_matrix *xmodel, int nmodel_points, int nparams, int nregresion_fns);

void makeHMatrix_fnptr(gsl_matrix *h_matrix, gsl_matrix *xmodel, int nmodel_points, int nparams, int nregression_fns,
											 void (*makeHVector_ptr)(gsl_vector *h_vector, gsl_vector *x_location, int nparams));

void estimateBeta(gsl_vector *beta_vector, gsl_matrix *h_matrix, gsl_matrix* cinverse, gsl_vector *trainingvector, int nmodel_points, int nregression_fns);

void makeHVector_linear(	gsl_vector *h_vector, gsl_vector *x_location, int nparams);
void makeHVector_quadratic( gsl_vector *h_vector, gsl_vector *x_location, int nparams);
void makeHVector_trivial(gsl_vector *h_vector, gsl_vector *x_location, int nparams);
void makeHVector_cubic( gsl_vector *h_vector, gsl_vector *x_location, int nparams);
#endif
