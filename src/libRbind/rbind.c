#include "rbind.h"



/**
 * estimate the optimum covariance length scales to reproduce the training data
 * 
 * the user can force the nugget to be fixed, in this case its value will not be optimized
 * 
 * @param xmodel_in: an r vector (flattened from a row-major array) of length nmodelpts_in*nparams_in which is converted into a col-major array
 * by a call to convertDoubleToMatrix. This contains the locations in the design space where the training runs of the model were evaluated 
 * @param nparams_in: dimensionality of the parameter (design) space
 * @param training_in: a flattened r vector of length nmodelpts_in, the output of runnign the model at the xmodel_in design locations
 * @param nmodelpts: number of locations where the model was trained
 * @param nthetas_in: the number of hyperparameters to be used, this needs to be correct for the specified cov fn, otherwise we'll have 
 * some kind of memory fault when we try and store the results in final_thetas
 * @param final_thetas: an nthetas_in blank vector
 * @param use_fixed_nugget: 0 -> nugget is set to be some relatively small value. 
 *           | x > 0 -> nugget is constrained to be within ~ 10% of x
 * @param cov_fn_index_in takes values: 1 -> POWEREXPCOVFN, 2-> MATERN32, 3->MATERN53 anything else 
 * will result in a failure
 * @param regression_order_in fixes the order of the linear regression model applied to the mean 
 * of the process, values of {0,1,2,3} are acceptable, other values result in using a trivial (order 0 model)
 * 
 * @return final_thetas is filled with the best set of thetas from the estimation process
 */
void callEstimate(double* xmodel_in, int* nparams_in, double* training_in, int *nmodelpts, int *nthetas_in, double* final_thetas,
									int* use_fixed_nugget, 
									double* fixed_nugget_in, int* cov_fn_index_in, int* regression_order_in){
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
	options.cov_fn_index = *cov_fn_index_in;
	options.regression_order = *regression_order_in;
	
	options.use_data_scales = 1; // use scales set by the data (this is a good idea)

	if(*use_fixed_nugget == 1){
		options.fixed_nugget_mode = 1;
		options.fixed_nugget = *fixed_nugget_in;
	} else {
		options.fixed_nugget_mode = 0;
		options.fixed_nugget = 0;
	}

	setup_cov_fn(&options);
	setup_regression(&options);

	alloc_modelstruct(&the_model, &options);

	// fill in xmodel, this is the right way!
	// note that xmodel has to be transposd in EmuRbind.R or else disaster
	printf("nx = %d, ny = %d\n", options.nparams, options.nmodel_points);
	convertDoubleToMatrix(the_model.xmodel, xmodel_in, options.nparams, options.nmodel_points);
	// fill in the training vec
	convertDoubleToVector(the_model.training_vector, training_in, options.nmodel_points);

	fill_sample_scales(&the_model, &options);
	setup_optimization_ranges(&options, &the_model);

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
 * compute the mean and variance at a set of points, for a certain emulator
 * the emulator is implicitly specified in terms of the training set (xmodel_in, training_in)
 * the hyper-parameters thetas_in and the choice of cov fn (cov_fn_index) and regression model (regression_order)
 * 
 * @param xmodel_in: an r vector (flattened from a row-major array) of length nmodelpts_in*nparams_in which is converted into a col-major array
 * by a call to convertDoubleToMatrix. This contains the locations in the design space where the training runs of the model were evaluated 
 * @param nparams_in: dimensionality of the parameter (design) space
 * @param points_in is a double vector (flattend from an R array of nparams x nemupoints), to be inflated back to life by convertDoubleToMatrix
 * @param nemupoints: is a list of how many points to emulate
 * @param training_in: a flattened r vector of length nmodelpts_in, the output of runnign the model at the xmodel_in design locations
 * @param nomodelpts: number of training locations in design space 
 * @param thetas_in: the hyperparameters which partially spec the emulator we're evaluating
 * @param nthetas_in: the number of hyperparameters to be used, this needs to be correct for the specified cov fn, otherwise we'll have 
 * some kind of memory fault when we try and store the results in final_thetas
 * @param final_emulated_y: (hopefully) blank vector of length nemupoints
 * @param final_emulated_variance: (hopefully) blank vector of length nemupoints
 * @param cov_fn_index_in takes values: 1 -> POWEREXPCOVFN, 2-> MATERN32, 3->MATERN53 anything else 
 * will result in a failure
 * @param regression_order_in fixes the order of the linear regression model applied to the mean 
 * of the process, values of {0,1,2,3} are acceptable, other values result in using a trivial (order 0 model)

 * 
 * @return final_emulated_y is a vector of the mean at each point
 * @return final_emulated_variance is a vector of the variance at each point
 */
void callEmulateAtList(double *xmodel_in, int *nparams_in, double* points_in, int *nemupoints, double* training_in,
											 int *nmodelpts, double* thetas_in, int *nthetas_in, double* final_emulated_y, 
											 double* final_emulated_variance, int* cov_fn_index_in, int* regression_order_in){
	optstruct options;
	modelstruct the_model; 
	//double the_mean, the_variance;
	gsl_matrix *the_point_array;
	int i, j;

	options.nmodel_points = *nmodelpts;
	options.nparams = *nparams_in;
	options.nemulate_points = *nemupoints;
	options.nthetas = *nthetas_in;
	options.cov_fn_index = *cov_fn_index_in;
	options.regression_order = *regression_order_in;

	/****************************************/
	/* defaults not used in this process    */
	options.use_data_scales = 1; 
	options.emulate_min = 0; 
	options.emulate_max = 1; 
	/****************************************/

	setup_cov_fn(&options);
	setup_regression(&options);


	//fprintf(stderr, "#callEmulate at list, nparams %d, nemulate_points %d\n", *nparams_in, *nemupoints);

	alloc_modelstruct(&the_model, &options);
	
	fill_sample_scales(&the_model, &options);
	setup_optimization_ranges(&options, &the_model); // not used

	//the_point_array = gsl_matrix_alloc(options.nparams, options.nemulate_points);
	the_point_array = gsl_matrix_alloc(options.nemulate_points, options.nparams);


	// fill in the point array
	convertDoubleToMatrix(the_point_array, points_in, options.nparams, options.nemulate_points);
	//convertDoubleToMatrix(the_point_array, points_in, options.nemulate_points, options.nparams);
	

	// fill in thetas
	convertDoubleToVector(the_model.thetas, thetas_in, options.nthetas);
	// fill in xmodel 
	convertDoubleToMatrix(the_model.xmodel, xmodel_in, options.nparams, options.nmodel_points);
	// fill in the training vec
	convertDoubleToVector(the_model.training_vector, training_in, options.nmodel_points);
	

	emulateAtPointList(&the_model, the_point_array, &options, final_emulated_y, final_emulated_variance);	
	
	gsl_matrix_free(the_point_array);
	// tidy up
	free_modelstruct(&the_model);
	free_optstruct(&options);

}

/**
 * computes the mean and variance at point_in, similar to callEmulateAtList
 * 
 * @param xmodel_in: an r vector (flattened from a row-major array) of length nmodelpts_in*nparams_in which is converted into a col-major array
 * by a call to convertDoubleToMatrix. This contains the locations in the design space where the training runs of the model were evaluated 
 * @param nparams_in: dimensionality of the parameter (design) space
 * @param point_in is a double array of nparams length, the point that we wish to eval the emulator at 
 * @param nemupoints: is a list of how many points to emulate
 * @param training_in: a flattened r vector of length nmodelpts_in, the output of runnign the model at the xmodel_in design locations
 * @param nmodelpts: number of training locations in design space 
 * @param thetas_in: the hyperparameters which partially spec the emulator we're evaluating
 * @param nthetas_in: the number of hyperparameters to be used, this needs to be correct for the specified cov fn, otherwise we'll have 
 * some kind of memory fault when we try and store the results in final_thetas
 * @param final_emulated_y: (hopefully) blank vector of length nemupoints
 * @param final_emulated_variance: (hopefully) blank vector of length nemupoints
 * @param cov_fn_index_in takes values: 1 -> POWEREXPCOVFN, 2-> MATERN32, 3->MATERN53 anything else 
 * will result in a failure
 * @param regression_order_in fixes the order of the linear regression model applied to the mean 
 * of the process, values of {0,1,2,3} are acceptable, other values result in using a trivial (order 0 model)
 *
 * @return final_emulated_y is the mean at point
 * @return final_emulated_variance is the variance at point
 *
 */
void callEmulateAtPt(double* xmodel_in, int* nparams_in, double* point_in, double* training_in, 
										 int* nmodelpts, double* thetas_in, int* nthetas_in, double* final_emulated_y, 
										 double* final_emulated_variance, int*cov_fn_index_in, int*regression_order_in){
	optstruct options;
	modelstruct the_model;
	double the_mean, the_variance;
	gsl_vector *the_point;
	int i;
	
	options.nmodel_points = *nmodelpts;
	options.nparams = *nparams_in;
	options.nemulate_points = 1;
	options.nthetas = *nthetas_in;

	options.cov_fn_index = *cov_fn_index_in;
	options.regression_order = *regression_order_in;
	/****************************************/
	/* defaults not used in this process    */
	options.use_data_scales = 1; 
	options.emulate_min = 0; 
	options.emulate_max = 1;
	/****************************************/
	setup_cov_fn(&options);
	setup_regression(&options);

	alloc_modelstruct(&the_model, &options);
	
	the_point = gsl_vector_alloc(options.nparams);

	convertDoubleToVector(the_point, point_in, options.nparams);

	// fill in thetas
	convertDoubleToVector(the_model.thetas, thetas_in, options.nthetas);
	// fill in xmodel 
	convertDoubleToMatrix(the_model.xmodel, xmodel_in, options.nparams, options.nmodel_points);
	// fill in the training vec
	convertDoubleToVector(the_model.training_vector, training_in, options.nmodel_points);

	fill_sample_scales(&the_model, &options);
	setup_optimization_ranges(&options, &the_model);	// this strictly isn't needed for emulator

	emulateAtPoint(&the_model, the_point, &options, &the_mean, &the_variance);

	//print_vector_quiet(the_point,options.nparams);
	fprintf(stderr, "#mean:%lf\tvar:%lf\n", the_mean, the_variance);
	*final_emulated_y  = the_mean;
	*final_emulated_variance = the_variance;

	gsl_vector_free(the_point);
	// tidy up
	free_modelstruct(&the_model);
	free_optstruct(&options);
}

/**
 * compute the log likelihood of a given configuration in theta space for the specified 
 * potential emulator, useful for validating/exploring what callEstimate gives
 *
 * @param xmodel_in: an r vector (flattened from a row-major array) of length nmodelpts_in*nparams_in which is converted into a col-major array
 * by a call to convertDoubleToMatrix. This contains the locations in the design space where the training runs of the model were evaluated 
 * @param nparams_in: dimensionality of the parameter (design) space
 * @param pointList_in is a double array of nevalPoints_in*nparams length, the points that we wish to eval the lhood at 
 * @param nevalPoints_in: how many points to evaluate
 * @param training_in: a flattened r vector of length nmodelpts_in, the output of runnign the model at the xmodel_in design locations
 * @param nmodelpts_in: number of training locations in design space 
 * @param nthetas_in: the number of hyperparameters to be used, this needs to be correct for the specified cov fn, otherwise we'll have 
 * some kind of memory fault when we try and store the results in final_thetas
 * @param answer: (hopefully) blank vector of length nevalPoints_in
 *
 * @param cov_fn_index_in takes values: 1 -> POWEREXPCOVFN, 2-> MATERN32, 3->MATERN53 anything else 
 * will result in a failure
 * @param regression_order_in fixes the order of the linear regression model applied to the mean 
 * of the process, values of {0,1,2,3} are acceptable, other values result in using a trivial (order 0 model)
 * 
 * @return answer: a list of nevalPoints_in lhood values
 */

void callEvalLhoodList(double *xmodel_in, int *nparams_in, double *pointList_in,
											 int *nevalPoints_in, double *training_in, int *nmodelPoints_in,
											 int *nthetas_in, double *answer,
											 int* cov_fn_index_in, int *regression_order_in)
{
	
	optstruct options;
	modelstruct the_model;
	double the_likelyhood = 0.0;
	int nevalPts  = *nevalPoints_in;
	gsl_matrix *the_point_array;
	double *xinput;
	int i,j;
	struct estimate_thetas_params params;

	const gsl_rng_type *T;
	T = gsl_rng_default;

	params.the_model = MallocChecked(sizeof(modelstruct));
	params.options = MallocChecked(sizeof(optstruct));

	/* printf("# callEvalLhoodList called:\nnparams %d\nnevalPoints %d\nnmodelPoints %d\nnthetas %d\n", */
	/* 			 *nparams_in, *nevalPoints_in, *nmodelPoints_in, *nthetas_in); */

	options.nmodel_points = *nmodelPoints_in;
	options.nparams = *nparams_in;
	options.nthetas = *nthetas_in;
	options.nemulate_points =1;
	options.cov_fn_index = *cov_fn_index_in;
	options.regression_order = *regression_order_in;
	/****************************************/
	/* defaults not used in this process    */
	options.use_data_scales = 1; 
	options.emulate_min = 0; 
	options.emulate_max = 1;
	/****************************************/
	setup_cov_fn(&options);
	setup_regression(&options);

	alloc_modelstruct(&the_model, &options);
	

	the_point_array = gsl_matrix_alloc(nevalPts, options.nthetas);
	convertDoubleToMatrix(the_point_array, pointList_in, options.nthetas, nevalPts);

	xinput = MallocChecked(sizeof(double)*options.nthetas);

	convertDoubleToMatrix(the_model.xmodel, xmodel_in, options.nparams, options.nmodel_points);
	convertDoubleToVector(the_model.training_vector, training_in, options.nmodel_points);


	fill_sample_scales(&the_model, &options);
	setup_optimization_ranges(&options, &the_model); // not actually used



	for(i = 0; i < options.nthetas; i++){
		gsl_vector_set(the_model.thetas, i, 0.0);
	}

	// copy the structures we just created into params
	copy_optstruct(params.options, &options);
	alloc_modelstruct(params.the_model, &options);
	copy_modelstruct(params.the_model, &the_model);

	// setup the rng in params
	params.random_number = gsl_rng_alloc(T);
	gsl_rng_set(params.random_number, get_seed_noblock());


	// setup the h matrix
	params.h_matrix = gsl_matrix_alloc(options.nmodel_points, options.nregression_fns);

	makeHMatrix(params.h_matrix, the_model.xmodel,options.nmodel_points, options.nparams, options.nregression_fns);


	for(i = 0; i < nevalPts; i++){
		for(j = 0; j < options.nthetas; j++){
			xinput[j] = gsl_matrix_get(the_point_array, i, j);
		}
		// we can get away with changing xinput each time since the params thetas are not used
		the_likelyhood = evalFnLBFGS(xinput, options.nthetas, &params);
		answer[i] = the_likelyhood;
	}
	
	// free params
	gsl_rng_free(params.random_number);
	free_modelstruct(params.the_model);
	free_optstruct(params.options);
	gsl_matrix_free(params.h_matrix);

	free_modelstruct(&the_model);
	free_optstruct(&options);
	gsl_matrix_free(the_point_array);
	free(xinput);
}

void fill_sample_scales(modelstruct* the_model, optstruct* options)
{
	int i, j; 
	double average_value, min_value;
	gsl_vector* differences = gsl_vector_alloc(options->nmodel_points-1);
	
	// compute the average separations
	for(i = 0; i < options->nparams; i++){
		average_value = 0;
		for(j = 0; j < (options->nmodel_points-1); j++){
			gsl_vector_set(differences, j, fabs(gsl_matrix_get(the_model->xmodel, j+1, i) - 
																					gsl_matrix_get(the_model->xmodel, j, i))) ;
			average_value += gsl_vector_get(differences, j);
		}
		// compute the min difference
		min_value = gsl_vector_min(differences);
		// compute the average difference
		average_value /= (options->nmodel_points-1);
		if(min_value < 1E-5){
			min_value = 0.00001;
		}
		gsl_vector_set(the_model->sample_scales, i, min_value);
		//fprintf(stderr, "# param %d min-value %lf average %lf\n", i, min_value, average_value);
	}
	
	gsl_vector_free(differences);

}



/*
 * takes an ALLOCATED gsl_matrix and copies the input vector into it, 
 * doesn't check anything
 * 
 * matrices in R are stored like this
 *      [,1] [,2] [,3] [,4]
 * [1,]    1    4    7    0
 * [2,]    2    5    8    0
 * [3,]    3    6    9    0
 * 
 * the input vector will look like this:
 * nx = 4	ny= 3
 * IN:0 = 1.000000
 * IN:1 = 2.000000
 * IN:2 = 3.000000
 * IN:3 = 4.000000
 * IN:4 = 5.000000
 * IN:5 = 6.000000
 * IN:6 = 7.000000
 * IN:7 = 8.000000
 * IN:8 = 9.000000
 * IN:9 = 0.000000
 * IN:10 = 0.000000
 * IN:11 = 0.000000
 *
 * as far as c is concerned matrices are going to be indexed by y and then x, i.e
 * a faithful copy of the input matrix should have:
 * 
 * m[1,1] = 1, m[2,1] = , m[3,2] =6 etc
 * does this work?
 */
void convertDoubleToMatrix(gsl_matrix* the_matrix, double* input, int nx, int ny){
	int i, j;
	/* for(i = 0; i < nx*ny; i++){ */
	/* 	printf("IN:%d = %lf\n", i, input[i]); */
	/* } */
	/* printf("\n"); */
	/* printf("nx = %d\tny=%d\n", nx, ny); */
	for(i = 0; i < nx; i++){
		for(j =0; j < ny; j++){
			gsl_matrix_set(the_matrix, j, i, input[j+ny*i]);
			/* printf("j+ny*i = %d\n", j+ny*i); */
			/* printf("MAT(%d,%d)=%lf\n", j+1, i+1, gsl_matrix_get(the_matrix,j,i)); */
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

void testConvert(double* matrix, int *nx_in, int *ny_in){
	int nx = *nx_in;
	int ny = *ny_in;
	printf("nx = %d\tny= %d\n", nx, ny);
	gsl_matrix *test = gsl_matrix_alloc(ny, nx);
	convertDoubleToMatrix(test, matrix, nx, ny);
	gsl_matrix_free(test);
}


/// graveyard ////


/** 
 * this is deprecated 
 * use callEstimate to get the thetas you need and then callEmulateAtPt or callEmulateAtList to
 * get the predicted mean(s) and variance(s) you desire.
 */
/* void callEmulator(double* xmodel_in, int *nparams_in, double* training_in, int *nmodelpts, int *nthetas, double* final_x, int* nemupts, \ */
/* 									double* finaly, double* finalvar, double* rangemin, double* rangemax){ */
/* 	double *thetas_vec; */
/* 	int no_nugget  = 0; */
/* 	thetas_vec = MallocChecked(sizeof(double)*(*nthetas)); */
/* 	// this will run the estimator and fill in the thetas vector with the best ones */
/* 	callEstimate(xmodel_in, nparams_in, training_in, nmodelpts, nthetas, thetas_vec, &no_nugget, NULL); */
/* 	callEmulate(xmodel_in, nparams_in, training_in, nmodelpts, thetas_vec, nthetas, final_x, nemupts, finaly, finalvar, rangemin, rangemax); */
/* 	free(thetas_vec); */
/* } */



/**
 * run the emulator against a given set of data and given hyperparams theta
 * this is deprecated, use callEmulateatPt or callEmulateatList
 */
/* void callEmulate(double* xmodel_in, int* nparams_in, double* training_in, int* nmodelpts, double* thetas_in, int* nthetas_in, double* final_emulated_x, int *nemupts_in, double* final_emulated_y, double* final_emulated_variance, double* range_min_in, double*range_max_in){ */
	
/* 	optstruct options; */
/* 	modelstruct the_model; */
/* 	resultstruct results; */

/* 	int i, j; */
/* 	options.nmodel_points = *nmodelpts; */
/* 	options.nparams = *nparams_in; */
/* 	options.nemulate_points = *nemupts_in; */
/* 	options.nthetas = *nthetas_in; */
/* 	// fuck this needs to be set automatically :( */
/* 	options.nregression_fns = 1 + options.nparams; */
/* 	options.emulate_min = *range_min_in; */
/* 	options.emulate_max = *range_max_in; */
/* 	setup_cov_fn(&options); */


/* 	alloc_modelstruct(&the_model, &options); */

/* 	alloc_resultstruct(&results, &options); */
	
/* 	// fill in thetas */
/* 	convertDoubleToVector(the_model.thetas, thetas_in, options.nthetas); */
/* 	// fill in xmodel  */
/* 	convertDoubleToMatrix(the_model.xmodel, xmodel_in, options.nparams, options.nmodel_points); */
/* 	// fill in the training vec */
/* 	convertDoubleToVector(the_model.training_vector, training_in, options.nmodel_points); */

/* 	// first we have to fill the sample ranges */
/* 	fill_sample_scales(&the_model, &options); */
/* 	setup_optimization_ranges(&options, &the_model);	// this strictly isn't needed for emulator */


/* 	// and call out to libEmu */
/* 	emulate_model_results(&the_model, &options, &results); */
	
/* 	// Fill in emulated_y, emulated_variance */
/* 	for(i = 0; i < options.nemulate_points; i++){ */
/* 		final_emulated_y[i] = gsl_vector_get(results.emulated_mean, i); */
/* 		final_emulated_variance[i] = gsl_vector_get(results.emulated_var, i); */
/* 	} */

/* 	// fill in final emulated_x */
/* 	for(j = 0; j < options.nparams; j++){ */
/* 		for(i = 0; i < options.nemulate_points; i++){ */
/* 			final_emulated_x[i+j*options.nmodel_points] = gsl_matrix_get(results.new_x, i, j); */
/* 		} */
/* 	} */

/* 	printf("hello from rbind-emulator\n"); */

/* 	for(i = 0; i < options.nthetas; i++) */
/* 		printf("%g\t", gsl_vector_get(the_model.thetas, i)); */
/* 	printf("\n"); */
					 
/* 	printf("nemulate pts = %d\n", options.nemulate_points); */

/* 	// check that the results struct agrees with final_emulated etc */
/* 	for(i = 0; i < options.nemulate_points; i++){ */
/* 		for(j = 0; j < options.nparams; j++){ */
/* 			printf("%d:%g\t%g\t", i+j*options.nemulate_points, gsl_matrix_get(results.new_x, i, j), final_emulated_x[i+j*options.nemulate_points]); */
/* 		} */
		
/* 		printf("%g\t%g\t", gsl_vector_get(results.emulated_mean, i), final_emulated_y[i]); */
/* 		printf("%g\t%g\n", gsl_vector_get(results.emulated_var, i), final_emulated_variance[i]); */
/* 	} */
	

/* 	// tidy up */
/* 	free_modelstruct(&the_model); */
/* 	free_optstruct(&options); */
/* 	free_resultstruct(&results); */
/* } */



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
/* void lagrange_interp(double* xin, double* valin, int* npts_in, double* desired_point_in){ */
/* 	int i, j; */
/* 	double retVal = 0; */
/* 	double weight = 0.0; */
/* 	double npts = *npts_in; */
/* 	double desired_point = *desired_point_in; */
	
/* 	for(i=0; i<npts; ++i){ */
/* 		weight = 1.0; */
/* 		for(j=0; j < npts; ++j){ */
/* 			if( j != i){ */
/* 				weight *= (desired_point - xin[j]) / (xin[i] - xin[j]); */
/* 			} */
/* 		} */
/* 		retVal += weight*valin[i]; */
/* 	} */
/* 	// push the answer back into the desired_point_in */
/* 	*desired_point_in = retVal; */
	
/* } */
