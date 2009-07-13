#include "emulator.h"
#include "math.h"


double model_function(double x){
	return(x*exp(-x));
}

int main (void){
	int i,j;
	int nparams = 1;
	int nthetas = 3 + nparams;
	int nmodel_points = 10;
	double test_covariance = 0.0;
	gsl_vector_view xm_view;
	gsl_vector_view xn_view;
	

	gsl_matrix *xmodel_test = gsl_matrix_alloc(nmodel_points, nparams); 
	gsl_matrix *covariance_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_vector *kvec_1 = gsl_vector_alloc(nmodel_points);
	gsl_vector *xnew = gsl_vector_alloc(nparams);	
	gsl_vector *training_points = gsl_vector_alloc(nmodel_points);
	gsl_vector *thetas = gsl_vector_alloc(nthetas);

	for(i = 0; i < nparams; i++) {gsl_vector_set(xnew, i, 0.35);}
	for(i = 0; i < nthetas; i++) {gsl_vector_set(thetas, i, 0.2);}


	for(i = 0; i < nmodel_points; i++){
		for(j = 0; j < nparams; j++){
			gsl_matrix_set(xmodel_test, i, j, 0.1*((double)i));
			gsl_vector_set(training_points, i, model_function(gsl_matrix_get(xmodel_test, i, j)));
		}
	}
	
	printf("modelpoints:\n");
	print_matrix(xmodel_test, nmodel_points, nparams);
	printf("\n");

	printf("trainingvec:\n");
	for(i = 0; i < nmodel_points; i++){
		printf("%g ", gsl_vector_get(training_points, i));
	}
	printf("\n");

	xm_view  = gsl_matrix_row(xmodel_test, 3);
	xn_view = gsl_matrix_row(xmodel_test, 2);
	test_covariance = covariance_fn(&xm_view.vector, &xn_view.vector, thetas, nthetas, nparams);

	makeKVector(kvec_1, xmodel_test, xnew, thetas, nmodel_points, nthetas, nparams);
	makeCovMatrix(covariance_matrix, xmodel_test, thetas, nmodel_points, nthetas, nparams); 
	printf("covariance of (3,2) = %g\n", test_covariance);

	
	printf("kvector:\n");	
	for (i = 0; i < nmodel_points; i++){
		printf("%g ", gsl_vector_get(kvec_1, i));
	}
	printf("\n");

	printf("covariance matrix: \n");

	print_matrix(covariance_matrix, nmodel_points, nmodel_points);
	printf("\n");

	// tidy up
	gsl_matrix_free(xmodel_test);
	gsl_matrix_free(covariance_matrix);
	gsl_vector_free(kvec_1);
	gsl_vector_free(training_points);
	gsl_vector_free(xnew);
		
	return(0);
}
