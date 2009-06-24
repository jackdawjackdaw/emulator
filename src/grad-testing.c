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
	
	/*gsl_vector_view xm_view;
		gsl_vector_view xn_view;*/

	gsl_matrix *xmodel_test = gsl_matrix_alloc(nmodel_points, nparams); 
	gsl_vector *xnew = gsl_vector_alloc(nparams);	
	gsl_vector *training_points = gsl_vector_alloc(nmodel_points);
	gsl_vector *thetas = gsl_vector_alloc(nthetas);

	gsl_vector *nelder_thetas = gsl_vector_alloc(nthetas);

	// things for the gradient descent
	int max_tries = 1;
	int number_steps = 10;
	double gamma = 1.0e-6;
	gsl_matrix *grad_ranges = gsl_matrix_alloc(nthetas, 2);

	T = gsl_rng_default;
	random_number = gsl_rng_alloc(T);


	for(i = 0; i < nparams; i++){gsl_vector_set(thetas, i, 0.1);}; // whatever don't care initially
	
	// set the training points and the model
	for(i = 0; i < nmodel_points; i++){
		for(j = 0; j < nparams; j++){
			gsl_matrix_set(xmodel_test, i, j, 0.1*((double)i));
			gsl_vector_set(training_points, i, model_function(gsl_matrix_get(xmodel_test, i, j)));
		}
	}

	
	// set the ranges for doing the gradient descent
	// just done 0..1 for each param right now
	/*for(i = 0; i < nthetas; i++){
		gsl_matrix_set(grad_ranges, i, 0, 0.0);
		gsl_matrix_set(grad_ranges, i, 1, 1.0);
		}*/
	gsl_matrix_set(grad_ranges, 0, 0, 0.0);
	gsl_matrix_set(grad_ranges, 0, 1, 0.1);
	gsl_matrix_set(grad_ranges, 1, 0, 0.0);
	gsl_matrix_set(grad_ranges, 1, 1, 0.1);
	gsl_matrix_set(grad_ranges, 2, 0, 0.0);
	gsl_matrix_set(grad_ranges, 2, 1, 0.1);
	gsl_matrix_set(grad_ranges, 3, 0, 0.0);
	gsl_matrix_set(grad_ranges, 3, 1, 0.4);

	/*
	// possibly worlds most gigantic function call!
	gradDesc(random_number, max_tries, number_steps, gamma, grad_ranges, xmodel_test, training_points, thetas, nmodel_points, nthetas, nparams);
	
	
	printf("best thetas");
	for(i = 0; i < nthetas; i++){
		printf("%g\n", gsl_vector_get(thetas, i));
	}
	*/


	printf("NELDERMEAD METHOD\n");
	// now do it again, but this time with neldermead
	nelderMead(random_number, max_tries, number_steps, nelder_thetas, grad_ranges, xmodel_test, training_points, thetas, nmodel_points, nthetas, nparams);
	
	printf("best thetas");
	for(i = 0; i < nthetas; i++){
		printf("%g\n", gsl_vector_get(nelder_thetas, i));
	}
	
	gsl_vector_free(xnew);
	gsl_vector_free(training_points);
	gsl_vector_free(thetas);
	gsl_matrix_free(xmodel_test);
 	gsl_rng_free(random_number);
	return(0);
}
