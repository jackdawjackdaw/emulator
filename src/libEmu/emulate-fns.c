#include "emulate-fns.h"


#include "emulator.h"
#include "regression.h"


/**
 * emulates the model at through emulate_min -> emulate_max in each dimension
 * @param results -> malloc'd and filled by testPts.txt
 * @return results has the calculted new_x, emulated_mean, emulated_var
 */
void emulate_model_results(modelstruct *the_model, optstruct* options, resultstruct* results){

	int i;
	double determinant_c = 0.0;

	gsl_matrix *c_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_matrix *cinverse = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);

	gsl_vector *beta_vector = gsl_vector_alloc(options->nregression_fns);
	gsl_matrix *h_matrix = gsl_matrix_alloc(options->nmodel_points, options->nregression_fns);

	gsl_matrix *temp_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);

	
	// create the covariance matrix and save a copy in temp_matrix
	makeCovMatrix(c_matrix, the_model->xmodel, the_model->thetas,options->nmodel_points, options->nthetas, options->nparams);
	gsl_matrix_memcpy(temp_matrix, c_matrix);

	//print_matrix(temp_matrix, options->nmodel_points, options->nmodel_points);
	
	chol_inverse_cov_matrix(options, temp_matrix, cinverse, &determinant_c);

	// regression cpts
	makeHMatrix(h_matrix, the_model->xmodel, options->nmodel_points, options->nparams, options->nregression_fns);
	estimateBeta(beta_vector, h_matrix, cinverse, the_model->training_vector, options->nmodel_points, options->nregression_fns);

	fprintf(stderr, "regression cpts: ");
	for(i = 0; i < options->nregression_fns; i++)
		fprintf(stderr, "%g ", gsl_vector_get(beta_vector, i));
	fprintf(stderr, "\n");
	
	// set the new_x values
	// NO longer, now we will use the values given in the resultstruct
	/* initialise_new_x(results->new_x, options->nparams, options->nemulate_points, options->emulate_min, options->emulate_max); */

	for(i = 0; i < options->nemulate_points; i++){
		printf("%g\n", gsl_matrix_get(results->new_x, i,0));
	}
	
	for(i = 0; i < options->nemulate_points; i++){		
		emulate_ith_location(the_model, options, results, i, h_matrix, cinverse, beta_vector);
	}
	
	gsl_matrix_free(c_matrix);
	gsl_matrix_free(cinverse);
	gsl_matrix_free(h_matrix);
	gsl_matrix_free(temp_matrix);
	gsl_vector_free(beta_vector);
}


/** 
 * this is the same as calling emulateAtPoint a bunch of times but you 
 * save on the setup costs of init'ing all the memory, should be quite
 * a lot faster for a large number of points
 * 
 * emulate the model at a list of points
 * @return *the_variance is set to the variance of each location in the list
 * @return *the_mean is set to the mean at each loc
 */
void emulateAtPointList(modelstruct *the_model, gsl_matrix* point_list, optstruct* options,
												double* the_mean, double* the_variance){
	int i;
	double determinant_c = 0.0;
	
	resultstruct results;

	gsl_matrix *c_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_matrix *cinverse = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);

	gsl_vector *beta_vector = gsl_vector_alloc(options->nregression_fns);
	gsl_matrix *h_matrix = gsl_matrix_alloc(options->nmodel_points, options->nregression_fns);

	gsl_matrix *temp_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);

	

	// allocate the resultsstruct
	alloc_resultstruct(&results, options);
	//results.new_x = gsl_matrix_alloc(options->nemulate_points, options->nparams);

	//fprintf(stderr, "%d %d\n", point_list->size1, point_list->size2);
	// copy in the x-values we want
	gsl_matrix_memcpy(results.new_x, point_list);

	/* for(i = 0; i < options->nparams; i++){ */
	/* 	printf("%g", gsl_matrix_get(results.new_x,0,i)); */
	/* } */
	/* printf("\n"); */

	
	// create the covariance matrix and save a copy in temp_matrix
	makeCovMatrix(c_matrix, the_model->xmodel, the_model->thetas,options->nmodel_points, options->nthetas, options->nparams);
	gsl_matrix_memcpy(temp_matrix, c_matrix);
	
	chol_inverse_cov_matrix(options, temp_matrix, cinverse, &determinant_c);

	// regression cpts
	makeHMatrix(h_matrix, the_model->xmodel, options->nmodel_points, options->nparams, options->nregression_fns);
	estimateBeta(beta_vector, h_matrix, cinverse, the_model->training_vector, options->nmodel_points, options->nregression_fns);

	
	for(i = 0; i < options->nemulate_points; i++){
		emulate_ith_location(the_model, options, &results, i, h_matrix, cinverse, beta_vector);
		the_mean[i] = gsl_vector_get(results.emulated_mean, i);
		the_variance[i] = gsl_vector_get(results.emulated_var, i);
	}


	
	gsl_matrix_free(c_matrix);
	gsl_matrix_free(cinverse);
	gsl_matrix_free(h_matrix);
	gsl_matrix_free(temp_matrix);
	gsl_vector_free(beta_vector);

	free_resultstruct(&results);
}


/**
 * emulate the model at the location of the_point
 * @return *the_variance is set to the variance at this location
 * @return *the_mean is set to the emulated mean at this loc
 */
void emulateAtPoint(modelstruct *the_model, gsl_vector* the_point, optstruct* options, double* the_mean, double* the_variance){

	int i;
	double determinant_c = 0.0;

	gsl_matrix *c_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_matrix *cinverse = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_vector *beta_vector = gsl_vector_alloc(options->nregression_fns);
	gsl_matrix *h_matrix = gsl_matrix_alloc(options->nmodel_points, options->nregression_fns);
	gsl_matrix *temp_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	double kappa;
	double temp_mean, temp_var;
	gsl_vector_view new_x_row;
	gsl_vector *kplus = gsl_vector_alloc(options->nmodel_points);
	gsl_vector *h_vector = gsl_vector_alloc(options->nregression_fns);

	
	// create the covariance matrix and save a copy in temp_matrix
	makeCovMatrix(c_matrix, the_model->xmodel, the_model->thetas,options->nmodel_points, options->nthetas, options->nparams);
	gsl_matrix_memcpy(temp_matrix, c_matrix);
	
	chol_inverse_cov_matrix(options, temp_matrix, cinverse, &determinant_c);

	// regression cpts
	makeHMatrix(h_matrix, the_model->xmodel, options->nmodel_points, options->nparams, options->nregression_fns);
	estimateBeta(beta_vector, h_matrix, cinverse, the_model->training_vector, options->nmodel_points, options->nregression_fns);
	
	makeKVector(kplus, the_model->xmodel, the_point, the_model->thetas, options->nmodel_points, options->nthetas, options->nparams);
	makeHVector(h_vector, the_point, options->nparams);
	
	temp_mean = makeEmulatedMean(cinverse, the_model->training_vector, kplus, h_vector, h_matrix, beta_vector, options->nmodel_points);
	kappa = covariance_fn(the_point, the_point, the_model->thetas, options->nthetas, options->nparams);
	temp_var = makeEmulatedVariance(cinverse, kplus, h_vector, h_matrix, kappa, options->nmodel_points, options->nregression_fns);

	if(temp_mean > 100){
		printf("EF:mean:%lf\tvar:%lf\n", temp_mean, temp_var);
		print_vector_quiet(the_point,options->nparams);
		printf("kplus:\n");
		print_vector_quiet(kplus, options->nmodel_points);
		printf("hvector:\n");
		print_vector_quiet(h_vector, options->nregression_fns);
		printf("beta:\n");
		print_vector_quiet(beta_vector, options->nregression_fns);
		//exit(1);
	}
	

	*the_mean = temp_mean;
	*the_variance = temp_var;

	gsl_vector_free(kplus);
	gsl_vector_free(h_vector);
	gsl_matrix_free(c_matrix);
	gsl_matrix_free(cinverse);
	gsl_matrix_free(h_matrix);
	gsl_matrix_free(temp_matrix);
	gsl_vector_free(beta_vector);
	
}

/**
 * emulate the model at the location the_point using the precomputed inverse cov_matrix, beta_vector, h_matrix and 
 * so forth
 */

void emulateQuick( modelstruct *the_model, gsl_vector* the_point, optstruct* options, 
									 double* mean_out, double* var_out, gsl_matrix *h_matrix, gsl_matrix *cinverse, gsl_vector* beta_vector ){

	double kappa;
	double temp_mean, temp_var;
	gsl_vector_view new_x_row;
	gsl_vector *kplus = gsl_vector_alloc(options->nmodel_points);
	gsl_vector *h_vector = gsl_vector_alloc(options->nregression_fns);

	
	makeKVector(kplus, the_model->xmodel, the_point, the_model->thetas, options->nmodel_points, options->nthetas, options->nparams);

	makeHVector(h_vector, the_point, options->nparams);
	
	temp_mean = makeEmulatedMean(cinverse, the_model->training_vector, kplus, h_vector, h_matrix, beta_vector, options->nmodel_points);

	kappa = covariance_fn(the_point, the_point, the_model->thetas, options->nthetas, options->nparams);

	temp_var = makeEmulatedVariance(cinverse, kplus, h_vector, h_matrix, kappa, options->nmodel_points, options->nregression_fns);

	//fprintf(stderr, "temp_mean %lf\ttemp_var %lf\n", temp_mean, temp_var);
	// save the results
	*mean_out = temp_mean;
	*var_out = temp_var;

	gsl_vector_free(kplus);
	gsl_vector_free(h_vector);

}


/**
 * emulate the model at the ith entry in results->new_x
 */
void emulate_ith_location(modelstruct *the_model, optstruct *options, resultstruct *results,int i, gsl_matrix* h_matrix, gsl_matrix* cinverse, gsl_vector *beta_vector){
	double kappa;
	double temp_mean, temp_var;
	gsl_vector_view new_x_row;
	gsl_vector *kplus = gsl_vector_alloc(options->nmodel_points);
	gsl_vector *h_vector = gsl_vector_alloc(options->nregression_fns);

	// read the new x location 
	new_x_row = gsl_matrix_row(results->new_x, i);
	//fprintf(stderr, "i(%d) new_x_row: ", i);
	//print_vector_quiet(&new_x_row.vector, options->nparams);

	
	makeKVector(kplus, the_model->xmodel, &new_x_row.vector, the_model->thetas, options->nmodel_points, options->nthetas, options->nparams);

	makeHVector(h_vector, &new_x_row.vector, options->nparams);
	
	temp_mean = makeEmulatedMean(cinverse, the_model->training_vector, kplus, h_vector, h_matrix, beta_vector, options->nmodel_points);


	kappa = covariance_fn(&new_x_row.vector, &new_x_row.vector, the_model->thetas, options->nthetas, options->nparams);

	temp_var = makeEmulatedVariance(cinverse, kplus, h_vector, h_matrix, kappa, options->nmodel_points, options->nregression_fns);

	//fprintf(stderr, "temp_mean %lf\ttemp_var %lf\n", temp_mean, temp_var);

	gsl_vector_set(results->emulated_mean, i, temp_mean);
	gsl_vector_set(results->emulated_var, i, temp_var);

	gsl_vector_free(kplus);
	gsl_vector_free(h_vector);
}
	
/**
 * cholesky decompose temp_matrix, invert it and copy the result into result_matrix
 * the determinant of the cholesky decomp is also calculated
 * 
 * it might be worth adding a little noise here, or catching the gsl error and trying to
 * recover?
 */

void chol_inverse_cov_matrix(optstruct* options, gsl_matrix* temp_matrix, gsl_matrix* result_matrix, double* final_determinant_c){
	int cholesky_test, i;
	double determinant_c = 0.0;
	gsl_error_handler_t *temp_handler; // disable the default handler
	// do a cholesky decomp of the cov matrix, LU is not stable for ill conditioned matrices
	temp_handler = gsl_set_error_handler_off();
	cholesky_test = gsl_linalg_cholesky_decomp(temp_matrix);
	if(cholesky_test == GSL_EDOM){
		fprintf(stderr, "trying to cholesky a non postive def matrix, in emulate-fns.c sorry...\n");
		exit(1);
	}
	gsl_set_error_handler(temp_handler);

	// find the determinant and then invert 
	// the determinant is just the trace squared
	determinant_c = 1.0;
	for(i = 0; i < options->nmodel_points; i++)
		determinant_c *= gsl_matrix_get(temp_matrix, i, i);
	determinant_c = determinant_c * determinant_c;

	//printf("det CHOL:%g\n", determinant_c);	
	gsl_linalg_cholesky_invert(temp_matrix);
	gsl_matrix_memcpy(result_matrix, temp_matrix);
	*final_determinant_c = determinant_c;
}
