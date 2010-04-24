#include "regression.h"

/** 
 * create a h_vector evaluated at given x_location., 
 * this currently only creates a set of linear basis fns, you could add another loop which 
 * would add squares etc 
 * 
 * on calling h_vector should be initliased to a vector of length nregression_fns which 
 * in the trivial case nregression_fns = nparams + 1 
 * 
 */
void makeHVector(gsl_vector *h_vector, gsl_vector *x_location, int nparams){
	makeHVector_linear(h_vector, x_location, nparams);
}

void makeHVector_linear(	gsl_vector *h_vector, gsl_vector *x_location, int nparams){
	int i;
	int offset = nparams + 1;
	double temp_val;
	gsl_vector_set(h_vector, 0, 1); // the first element is always a constant
	for(i = 0; i < nparams; i++) {
		temp_val = gsl_vector_get(x_location, i);
		gsl_vector_set(h_vector, i+1, temp_val);
	}
}

/**
 * an example of a basis of more complicated regression functions,
 * fiddling around with this might be useful if you *really* need prediction and 
 * you have some idea of what the regression spectrum should be like, i.e up to a certain
 * order in polys etc, 
 * 
 * sticking with the linear model is probably fine
 */
void makeHVector_quadratic( gsl_vector *h_vector, gsl_vector *x_location, int nparams){
	int i;
	int offset = nparams + 1;
	double temp_val;
	gsl_vector_set(h_vector, 0, 1); // the first element is always a constant
	for(i = 0; i < nparams; i++) {
		temp_val = gsl_vector_get(x_location, i);
		gsl_vector_set(h_vector, i+1, temp_val);
	}
	for(i = 0; i < nparams; i++){
		temp_val = gsl_vector_get(x_location, i);
		temp_val = temp_val * temp_val;
		gsl_vector_set(h_vector, nparams+i+1, temp_val);
	}
}	


/**
 * create the h_matrix, which is the a matrix of h_vectors evaluated at each 
 * design point from xmodel
 * where h_matrix is defined to be n_model_points x nregression_fns
 */
void makeHMatrix(gsl_matrix *h_matrix, gsl_matrix *xmodel, int nmodel_points, int nparams, int nregression_fns){
	int i,j; 
	gsl_vector *h_vec = gsl_vector_alloc(nregression_fns);
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
	gsl_matrix *temp_denominator = gsl_matrix_alloc(nregression_fns, nregression_fns);
	gsl_vector *temp_numerator = gsl_vector_alloc(nregression_fns);

	// first calculate h_matrix^{T}.cinverse (since we'll use this twice)
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, h_matrix, cinverse, 0.0, htrans_cinverse);
	// now calculate htrans_cinverse.hmatrix
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, htrans_cinverse, h_matrix, 0.0 ,temp_denominator);
	
	gsl_linalg_cholesky_decomp(temp_denominator);
	gsl_linalg_cholesky_invert(temp_denominator); // temp_denominator = (h_matrix^{T}.cinverse.hmatrix)^{-1}
	
	// temp_numerator = htrans_cinverse.trainingvector
	gsl_blas_dgemv(CblasNoTrans, 1.0, htrans_cinverse, trainingvector,0.0, temp_numerator);
	
	// and finally we set the value of the beta vector
	// beta = temp_denominator.temp_numerator
	gsl_blas_dgemv(CblasNoTrans, 1.0, temp_denominator, temp_numerator, 0.0, beta_vector);
	
	gsl_vector_free(temp_numerator);
	gsl_matrix_free(temp_denominator);
	gsl_matrix_free(htrans_cinverse);
}
