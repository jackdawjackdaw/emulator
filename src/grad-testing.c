#include "maximise.h"

double model_function(double x){
	return(x*exp(-x));
}

int main (void){
	// testing the gradDesc method
	int i,j;
	int nparams = 1;
	int nthetas = 3 + nparams;
	int nmodel_points = 10;
	double test_covariance = 0.0;
	const gsl_rng_type *T;
	gsl_rng *random_number;
	gsl_vector_view xm_view;
	gsl_vector_view xn_view;

	gsl_matrix *xmodel_test = gsl_matrix_alloc(nmodel_points, nparams); 
	gsl_vector *xnew = gsl_vector_alloc(nparams);	
	gsl_vector *training_points = gsl_vector_alloc(nmodel_points);
	gsl_vector *thetas = gsl_vector_alloc(nthetas);


	T = gsl_rng_default;
	random_number = gsl_rng_alloc(T);


	for(i = 0; i < nparams; i++){gsl_vector_set(thetas, i, 0.1)}; // whatever don't care initially
	
	for(i = 0; i < nmodel_points; i++){
		for(j = 0; j < nparams; j++){
			gsl_matrix_set(xmodel_test, i, j, 0.1*((double)i));
			gsl_vector_set(training_points, i, model_function(gsl_matrix_get(xmodel_test, i, j)));
		}
	}
	

	return(0);
}
