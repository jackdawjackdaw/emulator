#include "estimator.h"
#include "emulator.h"
#include "gsl/gsl_linalg.h"

double model_function(double x){
	return(x*exp(-x));
}

// test the gradient and log likelyhoods against the MM functions
int main (void){
	int i, j;
	int nthetas = 4;
	int nparams = 1;
	int nmodel_points = 11;
	int lu_signum = 0;
	int theta_index = 0;
	gsl_matrix *model_test = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix *covariance_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix *cinverse_cholesky = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix *temp_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix *cinverse = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_vector *training_points = gsl_vector_alloc(nmodel_points);
	gsl_vector *thetas = gsl_vector_alloc(nthetas);
	gsl_permutation *c_LU_permutation = gsl_permutation_alloc(nmodel_points);
	
	double likelyhood_test = 0.0;
	double gradient_test = 0.0;
	double det_cinverse = 0.0;
	
	// do the model and training vecs
		for(i = 0; i < nmodel_points; i++){
		for(j = 0; j < nparams; j++){
			gsl_matrix_set(model_test, i, j, 0.1*((double)i));
			gsl_vector_set(training_points, i, model_function(gsl_matrix_get(model_test, i, j)));
		}
	}
		
 	printf("modelpoints:\n");
	print_matrix(model_test, nmodel_points, nparams);
	printf("\n");

	printf("trainingvec:\n");
	for(i = 0; i < nmodel_points; i++){
		printf("%g ", gsl_vector_get(training_points, i));
	}
	printf("\n");
	
	// these test values used in the MM file grad-test.nb
	gsl_vector_set(thetas, 0, 0.04);
	gsl_vector_set(thetas, 1, 0.03);
	gsl_vector_set(thetas, 2, 0.001);
	gsl_vector_set(thetas, 3, 0.03);
	
	// make the covariance matrix
	makeCovMatrix(covariance_matrix, model_test, thetas, nmodel_points, nthetas, nparams);
	
	print_matrix(covariance_matrix, nmodel_points, nmodel_points);

	printf("\n");

	// be very careful with the arguments here, or you won't get a good copy
	gsl_matrix_memcpy( temp_matrix, covariance_matrix);
	gsl_matrix_memcpy( cinverse_cholesky, covariance_matrix);
	//print_matrix(temp_matrix, nmodel_points, nmodel_points);
	// now we have to make cinverse, and the determinant
	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
	// now invert
	gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse);
	
	printf("cinverse -> LU\n");
	print_matrix(cinverse, nmodel_points, nmodel_points);

	// now cinverse is how we want it, but have to decomp this to get the determinant (sigh)
	gsl_matrix_memcpy(temp_matrix, cinverse);
	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
	det_cinverse = gsl_linalg_LU_det(temp_matrix, lu_signum);

	//gsl_linalg_cholesky_decomp(cinverse_cholesky);
	//gsl_linalg_cholesky_invert(cinverse_cholesky); // documented but not in the lib?
	
	//printf("\n");
	//printf("cinverse -> cholesky\n");
	//print_matrix(cinverse_cholesky, nmodel_points, nmodel_points);

	// the gradient in the 0 direction should be -75.0812
	gradient_test = getGradient(cinverse, model_test, training_points, thetas, theta_index, nmodel_points, nthetas, nparams);
	likelyhood_test = getLogLikelyhood(cinverse, det_cinverse, model_test, training_points, thetas, nmodel_points, nthetas, nparams);
	
	printf("likelyhood = %g\n", likelyhood_test);
	printf("det_cinverse = %g\n", det_cinverse);
	printf("theta_index = %d\tgradient = %g\n", theta_index, gradient_test);
	
	gsl_matrix_free(model_test);
	gsl_matrix_free(covariance_matrix);
	gsl_matrix_free(temp_matrix);
	gsl_matrix_free(cinverse);
	gsl_vector_free(training_points);
	gsl_vector_free(thetas);
	gsl_permutation_free(c_LU_permutation);
	return(0);
}
	 
	
