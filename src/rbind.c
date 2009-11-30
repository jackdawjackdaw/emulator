#include "libEmu/emulator.h"
#include "libEmu/estimator.h"
#include "libEmu/maximise.h"
#include "estimate_threaded.h"
#include "multifit.h"
#include "useful.h"
#include "stdio.h"

void convertDoubleToMatrix(gsl_matrix* the_matrix, double* input, int ny, int nx);
void convertDoubleToVector(gsl_vector* the_vec, double* input, int nx);
void setEmulatorOptions(int* useGauss_in, double* alpha_in);
// this lives in libEmu/emulator.c it's important!
extern emulator_opts the_emulator_options;

/**
 * just enough setup and teardown to call the emulator directly from R
 */

//! change the covariance function options
/**
 * set the emulator options
 * @param useGauss, if this is one we'll use the gaussian covariance function
 * otherwise use the matern function
 * @param alpha, if we're using the gaussian covariance function this 
 * will set the power in the exponent, you should set it somewhere between 1 and 2 
 * if you're using the matern one and you set alpha to either 1/2, 3/2, 5/2 
 * we get fast versions of besselK and the emulator will be faster (but 
 * it will have one fixed param)
 * \bold the matern behaviour is not implemented yet 
 */
void setEmulatorOptions(int* useGauss_in, double* alpha_in){
	int useGauss = *useGauss_in;
	double alpha = *alpha_in;
	printf("%d\n", useGauss);
	printf("%g\n", alpha);


	if(useGauss == 1.0){
		the_emulator_options.usematern=0;
		if(alpha > 0 && alpha < 2){
			the_emulator_options.alpha = alpha;
		} else {
			the_emulator_options.alpha = 1.9;
		}
		fprintf(stderr, "using Gaussian Cov, with alpha = %g\n", the_emulator_options.alpha);
	}
	else if(useGauss==0.0){
		the_emulator_options.usematern=1;
		fprintf(stderr, "using Matern Cov, no alpha support yet\n");
	} else {
		fprintf(stderr, " bad call to setEmulatorOptions\n");
		exit(1);
	}
}



//! callEmulator from R
/**
 * so R can't handle passing anything other than arrays, 
 * as such anything which is canonically a 2d array (such as xmodel) will be passed 
 * in and out as a flat vector.
 */
void callEmulator(double* xmodel_in, int* nparams_in,  double* training_in, int *nmodelpts, int* nthetas_in, double* final_emulated_x, int *nemupts_in, double* final_emulated_y, double* final_emulated_variance , double* range_min_in, double* range_max_in ){
	int i, j;
	int nmodel_points = *nmodelpts;
	int nparams = *nparams_in;
	int nemupts = *nemupts_in;
	int nthetas = *nthetas_in;
	gsl_matrix *xmodel = gsl_matrix_alloc(nmodel_points, nparams);
	gsl_matrix *new_x = gsl_matrix_alloc(nemupts, nparams);
	gsl_vector *training_vec = gsl_vector_alloc(nmodel_points);
	gsl_vector *emulated_variance = gsl_vector_alloc(nemupts);
	gsl_vector *emulated_y = gsl_vector_alloc(nemupts);
	gsl_vector *thetas = gsl_vector_alloc(nthetas);
	//\todo sethis up
	gsl_rng* random; 
	const gsl_rng_type *rng_type;

	// need these to call the routines which actually do the emulation
	emuResult theResult;
	eopts theOptions;

	// reset the common vecs, this might be messing things up?
	for(i = 0; i < nemupts; i++){
		final_emulated_y[i] = 0.0;
		final_emulated_variance[i] = 0.0;
	}
	for(i =0; i < nemupts*nparams;i++){
		final_emulated_x[i] = 0.0;
	}
	

	fprintf(stderr, "%d\t%d\t%d\t%d\n", nmodel_points, nparams, nemupts, nthetas);

	// setup the rng
	rng_type = gsl_rng_default;
	random = gsl_rng_alloc(rng_type);
	gsl_rng_set(random, get_seed_noblock());

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// this is key
	// fills in a structure in libEmu which 
	// sets gaussian or matern cov fn and 
	// the alpha option for the gaussian
	//set_emulator_defaults(&the_emulator_options);
	// show the default options in the lib
	//print_emulator_options(&the_emulator_options);
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	// fill in xmodel 
	//! \todo fix this for nparams >1 !
	convertDoubleToMatrix(xmodel, xmodel_in,nparams, nmodel_points);

	// fill in the training vec
	convertDoubleToVector(training_vec, training_in, nmodel_points);

	
	print_matrix(xmodel, nmodel_points, nparams);

	// fill in the options
	theOptions.nmodel_points = nmodel_points;
	theOptions.nemu_points = nemupts;
	theOptions.nparams = nparams;
	theOptions.range_min = *range_min_in;
	theOptions.range_max = *range_max_in;
	theOptions.xmodel = xmodel;
	theOptions.nthetas = nthetas;
	theOptions.training = training_vec;
	theOptions.thetas = thetas;
	sprintf(theOptions.filename,"NOTHING");

	// alloc the result bits
	theResult.nemu_points = nemupts;
	theResult.nparams = nparams;
	theResult.new_x = new_x;
	theResult.new_mean = emulated_y; 
	theResult.new_var = emulated_variance;

	// actually do the emulation, this function 
	// can be found in multifit.c
	evaluate_region_threaded(&theResult, &theOptions, random);
	//evaluate_region(&theResult, &theOptions, random);

	fprintf(stderr, "back from evaluate region\n");

	// Fill in emulated_y, emulated_variance
	for(i = 0; i < nemupts; i++){
		final_emulated_y[i] = gsl_vector_get(emulated_y, i);
		final_emulated_variance[i] = gsl_vector_get(emulated_variance, i);
	}

	// fill in final emulated_x
	// this doesn't seem to work in >1d
	// should be interleaved
	for(j = 0; j < nparams; j++){
		for(i = 0; i < nemupts; i++){
			final_emulated_x[i+j*nmodel_points] = gsl_matrix_get(new_x, i, j);
		}
	}

	
	
	// tidy up
	gsl_matrix_free(xmodel);
	gsl_matrix_free(new_x);
	gsl_vector_free(emulated_variance);
	gsl_vector_free(emulated_y);
	gsl_vector_free(thetas);
	gsl_vector_free(training_vec);
	gsl_rng_free(random);
}



/**
 * just calculate the thetas for a model
 */
void callEstimate(double* xmodel_in, int* nparams_in, double* training_in, int *nmodelpts, int *nthetas_in, double* final_thetas){
	int i;
	int nmodel_points = *nmodelpts;
	int nparams = *nparams_in;
	int nthetas = *nthetas_in;
	gsl_matrix *xmodel = gsl_matrix_alloc(nmodel_points, nparams);
	gsl_vector *training_vec = gsl_vector_alloc(nmodel_points);
	gsl_vector *thetas = gsl_vector_alloc(nthetas);

	gsl_rng* random; 
	const gsl_rng_type *rng_type;
	
	eopts theOptions;
	// setup the rng
	rng_type = gsl_rng_default;
	random = gsl_rng_alloc(rng_type);
	gsl_rng_set(random, get_seed());

	// fill in xmodel 
	convertDoubleToMatrix(xmodel, xmodel_in, nparams, nmodel_points);
	// fill in the training vec
	convertDoubleToVector(training_vec, training_in, nmodel_points);

	// fill in the options
	theOptions.nmodel_points = nmodel_points;
	theOptions.nemu_points = 200;
	theOptions.nparams = nparams;
	// these two don't do anything, right?
	theOptions.range_min = 0;
	theOptions.range_max = 1;

	theOptions.xmodel = xmodel;
	theOptions.nthetas = nthetas;
	theOptions.training = training_vec;
	theOptions.thetas = thetas;

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// this is key
	// fills in a structure in libEmu which 
	// sets gaussian or matern cov fn and 
	// the alpha option for the gaussian
	set_emulator_defaults(&the_emulator_options);
	// show the default options in the lib
	print_emulator_options(&the_emulator_options);
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	// call out to the estimation func
	estimate_region_threaded(&theOptions);


	// set the final results
	for(i = 0; i < nthetas; i++)
		final_thetas[i] = gsl_vector_get(thetas, i);


	// tidy up
	gsl_matrix_free(xmodel);
	gsl_vector_free(thetas);
	gsl_vector_free(training_vec);
	gsl_rng_free(random);

}

/**
 * run the emulator against a given set of data and given hyperparams theta
 */
void callEmulate(double* xmodel_in, int* nparams_in, double* training_in, int* nmodelpts, double* thetas_in, int* nthetas_in, double* final_emulated_x, int *nemupts_in, double* final_emulated_y, double* final_emulated_variance, double* range_min_in, double*range_max_in){
	
	int i, j;
	int nmodel_points = *nmodelpts;
	int nparams = *nparams_in;
	int nemupts = *nemupts_in;
	int nthetas = *nthetas_in;
	fprintf(stderr,"nthetas = %d\n", nthetas);
	fprintf(stderr, "nparams = %d\n", nparams);
	fprintf(stderr, "nemupts = %d\n", nemupts);
	gsl_matrix *xmodel = gsl_matrix_alloc(nmodel_points, nparams);
	gsl_matrix *new_x = gsl_matrix_alloc(nemupts, nparams);
	gsl_vector *training_vec = gsl_vector_alloc(nmodel_points);
	gsl_vector *emulated_variance = gsl_vector_alloc(nemupts);
	gsl_vector *emulated_y = gsl_vector_alloc(nemupts);
	gsl_vector *thetas = gsl_vector_alloc(nthetas);
	//\todo sethis up
	gsl_rng* random; 
	const gsl_rng_type *rng_type;

	// need these to call the routines which actually do the emulation
	emuResult theResult;
	eopts theOptions;


	// setup the rng
	rng_type = gsl_rng_default;
	random = gsl_rng_alloc(rng_type);
	// don't wait for a truly random seeed :(
	gsl_rng_set(random, get_seed_noblock());


	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// this is key
	// fills in a structure in libEmu which 
	// sets gaussian or matern cov fn and 
	// the alpha option for the gaussian
	set_emulator_defaults(&the_emulator_options);
	// show the default options in the lib
	print_emulator_options(&the_emulator_options);
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	// fill in thetas
	convertDoubleToVector(thetas, thetas_in, nthetas);
	// fill in xmodel 
	convertDoubleToMatrix(xmodel, xmodel_in, nparams, nmodel_points);
	// fill in the training vec
	convertDoubleToVector(training_vec, training_in, nmodel_points);

	// fill in the options
	theOptions.nmodel_points = nmodel_points;
	theOptions.nemu_points = nemupts;
	theOptions.nparams = nparams;
	theOptions.range_min = *range_min_in;
	theOptions.range_max = *range_max_in;
	theOptions.xmodel = xmodel;
	theOptions.nthetas = nthetas;
	theOptions.training = training_vec;
	theOptions.thetas = thetas;

	// alloc the result bits
	theResult.nemu_points = nemupts;
	theResult.nparams = nparams;
	theResult.new_x = new_x;
	theResult.new_mean = emulated_y; 
	theResult.new_var = emulated_variance;

	// call out to multifit.c:emulate_region
	emulate_region(theResult.new_x, theResult.new_mean, theResult.new_var, &theOptions);
	
	// Fill in emulated_y, emulated_variance
	for(i = 0; i < nemupts; i++){
		final_emulated_y[i] = gsl_vector_get(emulated_y, i);
		final_emulated_variance[i] = gsl_vector_get(emulated_variance, i);
	}

	// fill in final emulated_x
	for(j = 0; j < nparams; j++){
		for(i = 0; i < nemupts; i++){
			final_emulated_x[i+j*nmodel_points] = gsl_matrix_get(new_x, i, j);
		}
	}

	
	
	// tidy up
	gsl_matrix_free(xmodel);
	gsl_matrix_free(new_x);
	gsl_vector_free(emulated_variance);
	gsl_vector_free(emulated_y);
	gsl_vector_free(thetas);
	gsl_vector_free(training_vec);
	gsl_rng_free(random);
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

	int nmodel_points = *nmodelpts_in;
	int nparams = *nparams_in;
	int nthetas = *nthetas_in;
	gsl_matrix *xmodel = gsl_matrix_alloc(nmodel_points, nparams);
	gsl_vector *training_vec = gsl_vector_alloc(nmodel_points);
	gsl_vector *thetas = gsl_vector_alloc(nthetas);
	double the_likelyhood = 0.0;

	convertDoubleToMatrix(xmodel, xmodel_in, nparams, nmodel_points);
	convertDoubleToVector(training_vec, training_in, nmodel_points);
	convertDoubleToVector(thetas, thetas_in, nthetas);

	// this calls the log likelyhood
	the_likelyhood = evalLikelyhood(thetas, xmodel, training_vec, nmodel_points, nthetas, nparams);

	*answer = the_likelyhood;

	gsl_matrix_free(xmodel);
	gsl_vector_free(training_vec);
	gsl_vector_free(thetas);

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
