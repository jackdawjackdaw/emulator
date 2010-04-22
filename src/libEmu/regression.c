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

/** 
 * create a h_vector evaluated at given x_location., 
 * this currently only creates a set of linear basis fns, you could add another loop which 
 * would add squares etc 
 * 
 * on calling h_vector should be initliased to a vector of length nregression_fns which 
 * in the trivial case nregression_fns = nparams + 1 
 */
void makeHVector(gsl_vector *h_vector, gsl_vector *x_location, int nparams){
	int i;
	gsl_vector_set(h_vector, 0, 1); // the first element is always a constant
	for(i = 0; i < nparams; i++) 
		gsl_vector_set(h_vector, i, gsl_vector_get(x_location, i));
}

/**
 * create the h_matrix, which is the a matrix of h_vectors evaluated at each 
 * design point from xmodel
 * where h_matrix is defined to be n_model_points x nregression_fns
 */
void makeHMatrix(gsl_matrix *h_matrix, gsl_matrix *xmodel, int nmodel_points, int nparams, int nregresion_fns){
	int i; 
	gsl_vector *h_vec = gsl_vector_alloc(nregresion_fns);
	gsl_vector_view xmodel_row_i;
	for(i = 0; i < nmodel_points; i++){
		xmodel_row_i = gsl_matrix_row(xmodel, i);
		makeHVector(h_vec, &xmodel_row_i.vector, nparams);
		gsl_matrix_set_row(h_matrix, i, h_vec);
	}
	gsl_vector_free(h_vec);
}

/**
 * estimate the values of the coefficients beta for a given training vector and inverse covariance matrix
 * 
 * beta = (h_matrix^{T}.cinverse.hmatrix)^{-1} . (h_matrix^{T}.cinverse.yvector)
 */
void estimateBeta(gsl_vector *beta_vector, gsl_matrix *h_matrix, gsl_matrix* cinverse, gsl_vector *trainingvector, int nmodel_points, int nregression_fns){
	int i;
	gsl_matrix *htrans_cinverse = gsl_matrix_alloc(nregression_fns, nmodel_points);
	gsl_matrix *temp_denominator = gsl_matrix_alloc(nregresion_fns, nregression_fns);
	gsl_vector *temp_numerator = gsl_matrix_alloc(nregression_fns);

	// first calculate h_matrix^{T}.cinverse (since we'll use this twice)
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, h_matrix, inverse_cov_matrix, htrans_cinverse);
	// now calculate htrans_cinverse.hmatrix
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, htrans_cinverse, h_matrix, temp_denominator);

	gsl_linalg_cholesky_decomp(temp_denominator);
	gsl_linalg_cholesky_invert(temp_denominator); // temp_denominator = (h_matrix^{T}.cinverse.hmatrix)^{-1}
	
	// temp_numerator = htrans_cinverse.trainingvector
	gsl_blas_dgemv(CblasNoTrans, 1.0, htrans_cinverse, trainingvector, temp_numerator);
	
	// and finally we set the value of the beta vector
	// beta = temp_denominator.temp_numerator
	gsl_blas_dgemv(CblasNoTrans, 1.0, temp_denominator, temp_numerator, beta_vector);
	
	gsl_vector_free(temp_numerator);
	gsl_matrix_free(temp_denominator);
	gsl_matrix_free(htrans_cinverse);
}
