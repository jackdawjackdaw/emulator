#include "rbind.h"

void callEmulator(double* xmodel_in, int *nparams_in, double* training_in, int *nmodelpts, int *nthetas, double* final_x, int* nemupts, \
									double* finaly, double* finalvar, double* rangemin, double* rangemax){
	double *thetas_vec;
	thetas_vec = MallocChecked(sizeof(double)*(*nthetas));
	// this will run the estimator and fill in the thetas vector with the best ones
	callEstimate(xmodel_in, nparams_in, training_in, nmodelpts, nthetas, thetas_vec);
	callEmulate(xmodel_in, nparams_in, training_in, nmodelpts, thetas_vec, nthetas, final_x, nemupts, finaly, finalvar, rangemin, rangemax);
	free(thetas_vec);
}


/**
 * calculate the thetas for a model
 * 
 * post-refactoring this is wonderfully simple
 */
void callEstimate(double* xmodel_in, int* nparams_in, double* training_in, int *nmodelpts, int *nthetas_in, double* final_thetas){
	int i;
	optstruct options;
	modelstruct the_model;

	
	// setup the optstruct
	options.nmodel_points = *nmodelpts;
	options.nparams = *nparams_in;
	options.nthetas = *nthetas_in;
	options.emulate_min = EMULATEMINDEFAULT;
	options.emulate_max = EMULATEMAXDEFAULT;
	options.grad_ranges = gsl_matrix_alloc(options.nthetas, 2);
	options.nregression_fns = 1;		// simple constant regression
	setup_cov_fn(&options);
	setup_optimization_ranges(&options);

	alloc_modelstruct(&the_model, &options);

	// fill in xmodel 
	convertDoubleToMatrix(the_model.xmodel, xmodel_in, options.nparams, options.nmodel_points);
	// fill in the training vec
	convertDoubleToVector(the_model.training_vector, training_in, options.nmodel_points);

	// actually do the estimation using libEmu
	estimate_thetas_threaded(&the_model, &options);

	// set the final results
	for(i = 0; i < options.nthetas; i++)
		final_thetas[i] = gsl_vector_get(the_model.thetas, i);

	// tidy up
	free_modelstruct(&the_model);
	free_optstruct(&options);

}

/**
 * run the emulator against a given set of data and given hyperparams theta
 */
void callEmulate(double* xmodel_in, int* nparams_in, double* training_in, int* nmodelpts, double* thetas_in, int* nthetas_in, double* final_emulated_x, int *nemupts_in, double* final_emulated_y, double* final_emulated_variance, double* range_min_in, double*range_max_in){
	
	optstruct options;
	modelstruct the_model;
	resultstruct results;

	int i, j;
	options.nmodel_points = *nmodelpts;
	options.nparams = *nparams_in;
	options.nemulate_points = *nemupts_in;
	options.nthetas = *nthetas_in;
	// fuck this needs to be set automatically :(
	options.nregression_fns = 1;
	options.emulate_min = *range_min_in;
	options.emulate_max = *range_max_in;
	setup_cov_fn(&options);
	setup_optimization_ranges(&options);	// this strictly isn't needed for emulator

	alloc_modelstruct(&the_model, &options);

	alloc_resultstruct(&results, &options);
	
	// fill in thetas
	convertDoubleToVector(the_model.thetas, thetas_in, options.nthetas);
	// fill in xmodel 
	convertDoubleToMatrix(the_model.xmodel, xmodel_in, options.nparams, options.nmodel_points);
	// fill in the training vec
	convertDoubleToVector(the_model.training_vector, training_in, options.nmodel_points);

	/* for(i = 0; i < options.nmodel_points; i++){ */
	/* 	printf("%lf\n", gsl_vector_get(the_model.training_vector, i)); */
	/* } */

	/* printf("nparams:%d\n", options.nparams); */
	/* for(i = 0; i < options.nparams; i++){ */
	/* 	for(j = 0; j < options.nmodel_points; j++){ */
	/* 		printf("%lf\t", gsl_matrix_get(the_model.xmodel, j, i)); */
	/* 	} */
	/* 	printf("\n"); */
	/* } */
						 

	// and call out to libEmu
	emulate_model_results(&the_model, &options, &results);
	
	// Fill in emulated_y, emulated_variance
	for(i = 0; i < options.nemulate_points; i++){
		final_emulated_y[i] = gsl_vector_get(results.emulated_mean, i);
		final_emulated_variance[i] = gsl_vector_get(results.emulated_var, i);
	}

	// fill in final emulated_x
	for(j = 0; j < options.nparams; j++){
		for(i = 0; i < options.nemulate_points; i++){
			final_emulated_x[i+j*options.nmodel_points] = gsl_matrix_get(results.new_x, i, j);
		}
	}

	printf("hello from rbind-emulator\n");

	for(i = 0; i < options.nthetas; i++)
		printf("%g\t", gsl_vector_get(the_model.thetas, i));
	printf("\n");
					 
	printf("nemulate pts = %d\n", options.nemulate_points);

	// check that the results struct agrees with final_emulated etc
	for(i = 0; i < options.nemulate_points; i++){
		for(j = 0; j < options.nparams; j++){
			printf("%d:%g\t%g\t", i+j*options.nemulate_points, gsl_matrix_get(results.new_x, i, j), final_emulated_x[i+j*options.nemulate_points]);
		}
		
		printf("%g\t%g\t", gsl_vector_get(results.emulated_mean, i), final_emulated_y[i]);
		printf("%g\t%g\n", gsl_vector_get(results.emulated_var, i), final_emulated_variance[i]);
	}
	

	// tidy up
	free_modelstruct(&the_model);
	free_optstruct(&options);
	free_resultstruct(&results);
}



/**
 * this is a binding to call the function libEmu/maximise.c:evalLikelyhood
 * double evalLikelyhood(gsl_vector *vertex, gsl_matrix *xmodel, gsl_vector *trainingvector, 
 * int nmodel_points, int nthetas, int nparams) 
 * 
 * this is mostly copied from above
 */

void callEvalLikelyhood(double * xmodel_in, int* nparams_in, double* training_in, \
													int *nmodelpts_in, int* nthetas_in, double* thetas_in, \
													double* answer){
	optstruct options;
	modelstruct the_model;
	double the_likelyhood = 0.0;
	struct estimate_thetas_params params;
	const gsl_rng_type *T;
	
	T = gsl_rng_default;

	params.the_model = MallocChecked(sizeof(modelstruct));
	params.options = MallocChecked(sizeof(optstruct));

	
	options.nmodel_points = *nmodelpts_in;
	options.nparams = *nparams_in;
	options.nthetas = *nthetas_in;
	options.nregression_fns = options.nparams+1;
	setup_cov_fn(&options);
	setup_optimization_ranges(&options);
	
	alloc_modelstruct(&the_model, &options);
	convertDoubleToMatrix(the_model.xmodel, xmodel_in, options.nparams, options.nmodel_points);
	convertDoubleToVector(the_model.training_vector, training_in, options.nmodel_points);
	convertDoubleToVector(the_model.thetas, thetas_in, options.nthetas);

	// copy in the structures we just created
	copy_optstruct(params.options, &options);
	alloc_modelstruct(params.the_model, &options);
	copy_modelstruct(params.the_model, &the_model);
	
	params.random_number = gsl_rng_alloc(T);
	gsl_rng_set(params.random_number, get_seed_noblock());

	// this calls the log likelyhood
	the_likelyhood = evalLikelyhoodLBFGS_struct(&params);

	*answer = the_likelyhood;
	
	// tidy up
	gsl_rng_free(params.random_number);
	free_modelstruct(params.the_model);
	free_modelstruct(&the_model);
	gsl_matrix_free(options.grad_ranges);
	gsl_matrix_free(params.options->grad_ranges);
	gsl_matrix_free(params.h_matrix);
	free(params.the_model);
	free(params.options);
	 
	

}

//! creates the coeffs for a lagrange poly interpolation
/**
 *  @param xin are the points at which we require our polynomial to pass through
 *  @param valin are the values of the function we wish to interpolate at xin
 *  @param npts how many of xin there are, this also determines the order of the L poly
 *  @param where we want this to be evaluated
 *  @return desired_point_in is set to the correct value
 *
 * this is clearly possible in R but i just can't quite get it right
 */
void lagrange_interp(double* xin, double* valin, int* npts_in, double* desired_point_in){
	int i, j;
	double retVal = 0;
	double weight = 0.0;
	double npts = *npts_in;
	double desired_point = *desired_point_in;
	
	for(i=0; i<npts; ++i){
		weight = 1.0;
		for(j=0; j < npts; ++j){
			if( j != i){
				weight *= (desired_point - xin[j]) / (xin[i] - xin[j]);
			}
		}
		retVal += weight*valin[i];
	}
	// push the answer back into the desired_point_in
	*desired_point_in = retVal;
	
}




/*
 * driving to atlanta: fixed bug in this
 * 
 * takes an ALLOCATED gsl_matrix and copies the input vector into it, 
 * doesn't check anything
 */
/** dimensions are interleaved! */
void convertDoubleToMatrix(gsl_matrix* the_matrix, double* input, int nx, int ny){
	int i, j;
	for(j = 0; j < ny; j++){
		for(i =0; i < nx; i++){
			gsl_matrix_set(the_matrix, j, i, input[nx*j+i]);
		}
	}
}

										 
/*
 * takes an ALLOCATED gsl_vector and copies the input vector into it
 * non checking again
 */
void convertDoubleToVector(gsl_vector* the_vec, double* input, int nx){
	int i;
	for(i =0; i < nx; i++){
		gsl_vector_set(the_vec, i, input[i]);
	}
}

