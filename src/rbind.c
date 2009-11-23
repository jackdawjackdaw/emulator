#include "libEmu/emulator.h"
#include "libEmu/estimator.h"
#include "libEmu/maximise.h"
#include "multifit.h"
#include "useful.h"

void convertDoubleToMatrix(gsl_matrix* the_matrix, double* input, int ny, int nx);
void convertDoubleToVector(gsl_vector* the_vec, double* input, int nx);

/**
 * just enough setup and teardown to call the emulator directly from R
 */

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

	//fprintf(stderr, "%d\t%d\t%d\t%d\n", nmodel_points, nparams, nemupts, nthetas);

	// setup the rng
	rng_type = gsl_rng_default;
	random = gsl_rng_alloc(rng_type);
	gsl_rng_set(random, get_seed());

	// fill in xmodel 
	convertDoubleToMatrix(xmodel, xmodel_in, nparams, nmodel_points);
	// fill in the training vec
	convertDoubleToVector(training_vec, training_in, nmodel_points);


/* 	for(j=0; j < nparams; j++){ */
/* 		for(i = 0; i < nmodel_points; i++){ */
/* 			fprintf(stderr, "%g\n" , gsl_matrix_get(xmodel, i, j)); */
/* 		} */
/* 	} */
	
/* 	// fill in the training vec */
/* 	for(i = 0; i < nmodel_points; i++){ */
/* 		fprintf(stderr, "%g\n", gsl_vector_get(training_vec, i)); */
/* 	} */
	


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

	// actually do the emulation, this function 
	// can be found in multifit.c
	evaluate_region_threaded(&theResult, &theOptions, random);

	fprintf(stderr, "back from evaluate region\n");

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

/*void callEmulator(double* xmodel_in, int* nparams_in,  double* training_in, int *nmodelpts, 
int* nthetas_in, double* final_emulated_x, int *nemupts_in, double* final_emulated_y, 
double* final_emulated_variance , double* range_min_in, double* range_max_in ){ */

/*
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
 * takes an ALLOCATED gsl_matrix and copies the input vector into it, 
 * doesn't check anything
 */

void convertDoubleToMatrix(gsl_matrix* the_matrix, double* input, int ny, int nx){
	int i, j;
	for(j = 0; j < ny; j++){
		for(i =0; i < nx; i++){
			gsl_matrix_set(the_matrix, i, j, input[i+j*nx]);
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
