#include "rbind.h"

/**
 * \todo make similar or wrapper functions which would allow any C compatible ffi 
 * without R's weirdness about vectors etc to setup and use the emulator
 * 
 * intended uses: python, paraview, scott? 
 */


/**
 * estimate the optimum covariance length scales to reproduce the training data
 * do this first
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
 * compute the mean and variance at a set of points (the list), for a certain emulator. 
 * this is much faster than repeated calls to callEmulateAtPoint
 * 
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
	options.fixed_nugget_mode = 0;
	options.fixed_nugget = 0;

	/****************************************/

	setup_cov_fn(&options);
	setup_regression(&options);


	//fprintf(stderr, "#callEmulate at list, nparams %d, nemulate_points %d\n", *nparams_in, *nemupoints);

	alloc_modelstruct(&the_model, &options);
	
	fill_sample_scales(&the_model, &options);
	setup_optimization_ranges(&options, &the_model); // not used

	//the_point_array = gsl_matrix_alloc(options.nparams, options.nemulate_points);
	the_point_array = gsl_matrix_alloc(options.nemulate_points, options.nparams);


	// fill in the_point_array from points_in
	// for(i, ...)
	// for(j, ...)
	// gsl_matrix_set(the_point_array, j, i, points_in[j+ny*i]);
	// so the points are evaluated in the natural order of the matrix
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
	options.fixed_nugget_mode = 0;
	options.fixed_nugget = 0;
	options.use_data_scales = 1; 
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
 * sometimes we will want to make MANY calls to emulate at a set of at the time of calling this 
 * function,  unknown points, 
 * for instance if we're doing some kind of monte-carlo sampling of the posterior of our emulator.
 * 
 * this function will setup an optstruct, modelstruct and some other useful bits of memory
 * which can then be used  by callEmulateMC to make fast calls to the emulator without all the setup
 * and teardown.
 * 
 * finally you'll want to call destroyEmulateMC to stop some kind of horrible leakyness
 * 
 * @param xmodel_in: an r vector (flattened from a row-major array) of length nmodelpts_in*nparams_in which is converted into a col-major array
 * by a call to convertDoubleToMatrix. This contains the locations in the design space where the training runs of the model were evaluated 
 * @param nparams_in: dimensionality of the parameter (design) space
 * @param nemupoints: is a list of how many points to emulate
 * @param training_in: a flattened r vector of length nmodelpts_in, the output of runnign the model at the xmodel_in design locations
 * @param nmodelpts: number of training locations in design space 
 * @param thetas_in: the hyperparameters which partially spec the emulator we're evaluating
 * @param nthetas_in: the number of hyperparameters to be used, this needs to be correct for the specified cov fn, otherwise we'll have 
 * some kind of memory fault when we try and store the results in final_thetas
 *
 * @param cov_fn_index_in takes values: 1 -> POWEREXPCOVFN, 2-> MATERN32, 3->MATERN53 anything else 
 * will result in a failure
 * @param regression_order_in fixes the order of the linear regression model applied to the mean 
 * of the process, values of {0,1,2,3} are acceptable, other values result in using a trivial (order 0 model)

 */
void setupEmulateMC(double* xmodel_in, int* nparams_in,  double* training_in, 
										 int* nmodelpts, double* thetas_in, int* nthetas_in, 
										int *cov_fn_index_in, int*regression_order_in){

	// call the helper with the mono-variate emuMCData objet
	setupEmulateMCHelper(&emuMCData, xmodel_in, nparams_in, training_in, 
									 nmodelpts, thetas_in, nthetas_in, cov_fn_index_in,
									 regression_order_in);
}

/** 
 * internal: called by setupEmulateMC and setupEmulateMCMulti
 */
void setupEmulateMCHelper(struct emulateMCData* emuMCData, double* xmodel_in, 
													int* nparams_in,  double* training_in, 
										 int* nmodelpts, double* thetas_in, int* nthetas_in, 
										int *cov_fn_index_in, int*regression_order_in){
	// we'll setup these structures first and then copy them to the
	// ones at the scope of rbind.h
	optstruct* options;
	modelstruct* the_model;
	double determinant_c;
	gsl_matrix *c_matrix ;
	gsl_matrix *cinverse;
	gsl_vector *beta_vector;
	gsl_matrix *h_matrix; 
	gsl_matrix *temp_matrix;

	emuMCData->options = MallocChecked(sizeof(optstruct));
	emuMCData->model = MallocChecked(sizeof(modelstruct));

	options = emuMCData->options;
	the_model = emuMCData->model;

	
	options->nmodel_points = *nmodelpts;
	options->nparams = *nparams_in;
	options->nemulate_points = 1; // not used
	options->nthetas = *nthetas_in;

	options->cov_fn_index = *cov_fn_index_in;
	options->regression_order = *regression_order_in;

	/****************************************/
	/* defaults not used in this process    */
	options->fixed_nugget_mode = 0;
	options->fixed_nugget = 0;
	options->use_data_scales = 1; 
	/****************************************/
	setup_cov_fn(options);
	setup_regression(options);

	alloc_modelstruct(the_model, options);
	
	temp_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);

	// fill in thetas
	convertDoubleToVector(the_model->thetas, thetas_in, options->nthetas);
	// fill in xmodel 
	convertDoubleToMatrix(the_model->xmodel, xmodel_in, options->nparams, options->nmodel_points);
	// fill in the training vec
	convertDoubleToVector(the_model->training_vector, training_in, options->nmodel_points);

	fill_sample_scales(the_model, options);
	setup_optimization_ranges(options, the_model);	// this strictly isn't needed for emulator



	// now we can also generate the covariance matrix we'll use all the time when computing the 
	// posterior mean and variance.
	// specifically: cov_matrix, inverse_cov_matrix, h_matrix, beta_vector
	emuMCData->cov_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	emuMCData->cov_matrix_inverse = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	emuMCData->beta_vector = gsl_vector_alloc(options->nregression_fns);
	emuMCData->h_matrix = gsl_matrix_alloc(options->nmodel_points, options->nregression_fns);

	c_matrix = emuMCData->cov_matrix;
	cinverse = emuMCData->cov_matrix_inverse;
	beta_vector = emuMCData->beta_vector;
	h_matrix = emuMCData->h_matrix;

	
	// compute the cov matrix
	makeCovMatrix(c_matrix, the_model->xmodel, the_model->thetas,options->nmodel_points, options->nthetas, options->nparams);
	gsl_matrix_memcpy(temp_matrix, c_matrix);
	// now invert the cov_matrix
	chol_inverse_cov_matrix(options, temp_matrix, cinverse, &determinant_c);
	// do the regression stuff
	makeHMatrix(h_matrix, the_model->xmodel, options->nmodel_points, options->nparams, options->nregression_fns);
	estimateBeta(beta_vector, h_matrix, cinverse, the_model->training_vector, options->nmodel_points, options->nregression_fns);
	
}

/**
 * rapidly sample the posterior distribution for the emluator
 *  at the point_in, we need to have called setupEmulateMC first
 * 
 * @param point_in: a nparams array representing the location in parameter space at which we 
 * wish to sample the emulator. 
 * @return mean_out: the estimated mean of the emulator at this point
 * @return var_out: the estimated variance of the emulator at this point
 */
void callEmulateMC(double* point_in, double* mean_out, double* var_out){
	double temp_mean;
	double temp_var;

	gsl_vector* the_point;

	optstruct* options = emuMCData.options;
	modelstruct* the_model = emuMCData.model;
	gsl_matrix* cmat = emuMCData.cov_matrix;
	gsl_matrix* cmatrix_inv = emuMCData.cov_matrix_inverse;
	gsl_vector* beta_vector = emuMCData.beta_vector;
	gsl_matrix* h_matrix = emuMCData.h_matrix;

	assert(cmat != NULL);
	assert(cmatrix_inv != NULL);
	assert(beta_vector != NULL);
	assert(h_matrix != NULL);

	
	the_point = gsl_vector_alloc(options->nparams);
	convertDoubleToVector(the_point, point_in, options->nparams);

	emulateQuick(the_model, the_point, options, &temp_mean, &temp_var, 
							 h_matrix, cmatrix_inv, beta_vector);
	
	*mean_out = temp_mean;
	*var_out = temp_var;

	// actually, let's not do this
	//fprintf(stderr, "# mean:%lf\tvar:%lf\n", temp_mean, temp_var);
	gsl_vector_free(the_point);
}

/**
 * free all the memory we allocated in setupEmulateMC
 */
void freeEmulateMC(void){
	free_optstruct(emuMCData.options);
	free_modelstruct(emuMCData.model);
	free(emuMCData.options);
	free(emuMCData.model);
	gsl_matrix_free(emuMCData.cov_matrix);
	gsl_matrix_free(emuMCData.cov_matrix_inverse);
	gsl_vector_free(emuMCData.beta_vector);
	gsl_matrix_free(emuMCData.h_matrix);
}


/**
 ******************************************************************************************
 * multivariate MC calls.
 ******************************************************************************************/



/**
 * essentially the same as setupEmulateMC but for multidimensional data 
 *
 * this function will setup an optstruct, modelstruct and some other useful bits of memory
 * which can then be used  by callEmulateMCMulti to make fast calls to the emulator without all the setup
 * and teardown.
 * 
 * finally you'll want to call destroyEmulateMCMulti to stop some kind of horrible leakyness
 * 
 * @param xmodel_in: an r vector (flattened from a row-major array) of length nmodelpts_in*nparams_in which is converted into a col-major array
 * by a call to convertDoubleToMatrix. This contains the locations in the design space where the training runs of the model were evaluated 
 * @param nparams_in: dimensionality of the parameter (design) space
 * @param nemupoints: is a list of how many points to emulate
 * @param training_in: flat r matrix of nrows=nmodelpts_in ncols=nydims_in, the output of running
 * @param nmodelpts: number of training locations in design space 
 * @param thetas_in: the hyperparameters which partially spec the emulator we're evaluating,  a flattened r matrix of nrows=nydims ncols=nthetas
 * @param nthetas_in: the number of hyperparameters to be used, this needs to be correct for the specified cov fn, otherwise we'll have 
 * some kind of memory fault when we try and store the results in final_thetas
 *
 * @param cov_fn_index_in takes values: 1 -> POWEREXPCOVFN, 2-> MATERN32, 3->MATERN53 anything else 
 * will result in a failure
 * @param regression_order_in fixes the order of the linear regression model applied to the mean 
 * of the process, values of {0,1,2,3} are acceptable, other values result in using a trivial (order 0 model)

 */
void setupEmulateMCMulti(double* xmodel_in, int* nparams_in,  
												 double* training_in, int* nydims_in,
												 int* nmodelpts_in, double* thetas_in, int* nthetas_in, 
												 int *cov_fn_index_in, int*regression_order_in){

	int nydims = *nydims_in;
	int nparams = *nparams_in;
	int nmodelpts = *nmodelpts_in;
	int nthetas = *nthetas_in;
	int index;
	int i,j;
	double *training_double_vec = MallocChecked(sizeof(double)*nmodelpts);
	double *theta_double_vec = MallocChecked(sizeof(double)*nthetas);

	// we want to allocate nydims copies of a emulateMCData structure
	// i'm worried we haven't done that below
	emuMCDataMulti = MallocChecked(sizeof(struct emulateMCData)*nydims);
	
	// coerce the flattened R data back to the matrix shapes we care about
	gsl_matrix *training_matrix = gsl_matrix_alloc(nmodelpts, nydims);
	gsl_matrix *thetas_matrix = gsl_matrix_alloc(nydims, nthetas);

	// annoyingly there's a mix up between nx and ny...
	convertDoubleToMatrix(training_matrix, training_in, nydims, nmodelpts);
	convertDoubleToMatrix(thetas_matrix, thetas_in, nthetas, nydims);

	// now we setup each emuMCDataMulti structure
	for (index = 0; index < nydims; ++index){
		
		for (i = 0; i < nthetas; ++i)
			theta_double_vec[i] = gsl_matrix_get(thetas_matrix, index, i);

		for(i = 0; i < nmodelpts; ++i)
			training_double_vec[i] = gsl_matrix_get(training_matrix, i, index);
		
		assert(&(emuMCDataMulti[index]) != NULL);
		setupEmulateMCHelper(&(emuMCDataMulti[index]), xmodel_in, nparams_in,
												 training_double_vec, nmodelpts_in, theta_double_vec, 
												 nthetas_in, cov_fn_index_in, regression_order_in);
	}			
		
	gsl_matrix_free(training_matrix);
	gsl_matrix_free(thetas_matrix);
}

/**
 * rapidly emulates the a multivariate model  at the location point_in
 * @param point_in: a nparams array representing the location in parameter space at which we 
 * wish to sample the emulator. 
 * @param nydims_in: the dim of our multivar model
 * @return final_mean: vec of the estimated mean of the emulator at this point 
 * @return var_out: vec of the estimated variance of the emulator at this point
 */
void callEmulateMCMulti(double* point_in, int* nydims_in, double* final_mean, double* final_var){
	int nydims = *nydims_in;
	int nparams, index;
	gsl_vector* the_point;
	double *temp_mean = MallocChecked(sizeof(double)*nydims);
	double *temp_var = MallocChecked(sizeof(double)*nydims);

	optstruct* options = emuMCDataMulti[0].options; //this is the same for each dimension
	modelstruct* the_model;
	gsl_matrix* cmat;
	gsl_matrix* cmatrix_inv;
	gsl_vector* beta_vector;
	gsl_matrix* h_matrix;
		
	nparams = options->nparams;

	the_point = gsl_vector_alloc(nparams);
	convertDoubleToVector(the_point, point_in, nparams);

	for(index = 0; index <nydims; index++){
		options = emuMCDataMulti[index].options;
		the_model = emuMCDataMulti[index].model;
		cmat = emuMCDataMulti[index].cov_matrix;
		cmatrix_inv = emuMCDataMulti[index].cov_matrix_inverse;
		beta_vector = emuMCDataMulti[index].beta_vector;
		h_matrix = emuMCDataMulti[index].h_matrix;

		assert(cmat != NULL);
		assert(cmatrix_inv != NULL);
		assert(beta_vector != NULL);
		assert(h_matrix != NULL);

		emulateQuick(the_model, the_point, options, &(temp_mean[index]), &(temp_var[index]), 
								 h_matrix, cmatrix_inv, beta_vector);
		
		final_mean[index] = temp_mean[index];
		final_var[index] = temp_var[index];
	}

	free(temp_mean);
	free(temp_var);
	gsl_vector_free(the_point);
} 

/* free all the memory we alloccd 
 */
void freeEmulateMCMulti(int *nydims_in){
	int nydims = *nydims_in;
	int index;
	for(index = 0; index < nydims; index++){
		free_optstruct(emuMCDataMulti[index].options);
		free_modelstruct(emuMCDataMulti[index].model);
		free(emuMCDataMulti[index].options);
		free(emuMCDataMulti[index].model);
		gsl_matrix_free(emuMCDataMulti[index].cov_matrix);
		gsl_matrix_free(emuMCDataMulti[index].cov_matrix_inverse);
		gsl_vector_free(emuMCDataMulti[index].beta_vector);
		gsl_matrix_free(emuMCDataMulti[index].h_matrix);
	}
}


/**
 ******************************************************************************************
 * likelihood evaluation
 ******************************************************************************************/


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
	gsl_vector_view point_vec_view;

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
	options.fixed_nugget_mode = 0;
	options.fixed_nugget = 0;
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
		/* for(j = 0; j < options.nthetas; j++){ */
		/* 	xinput[j] = gsl_matrix_get(the_point_array, i, j); */
		/* } */
		point_vec_view = gsl_matrix_row(the_point_array, i);
		// we can get away with changing xinput each time since the params thetas are not used
		the_likelyhood = evalFnMulti(&(point_vec_view.vector),  &params);
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


/**
 * print out the version number and some usage info to stdout
 * to be used in R
 */
void callInfo(void){
	printf("# libRbind build version: VERSION_NUMBER\n");
	printf("# covFn's:\n");
	printf("# 1 POWEREXP\n");
	printf("# 2 MATERN 32\n");
	printf("# 3 MATERN 52\n");
	printf("# regression order's supported: 0-3\n");

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


