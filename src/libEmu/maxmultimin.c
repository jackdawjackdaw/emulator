#include "maxmultimin.h"

/**
 * @file
 * @author C.Coleman-Smith, cec24@phy.duke.edu
 * @date oct-11
 * @version 0.1
 * @section DESCRIPTION
 * 
 * contains routines to use the multimin maximisation methods
 * these are provided by the gsl, and include simplex (nelder-mead) and
 * lbfgs methods. 
 * 
 * This is currently experimental work but it is desirable as it removes dependences upon
 * fortran!
 * 
 * we want to interface with estimate_threaded, so we expose the function:
 * maxWithMultimin(struct estimate_thetas_params *params)
 * 
 * when optimization has finished params->the_model->thetas is set to the winning 
 * set and params->best_value should be set to the loglikelihood of this set
 * 
 * 
 */

#define SCREWUPVALUE -2000

/**
 * call this to max the given model (params) with multimin, this is copied 
 * from maxWithLBFGS
 */
void maxWithMultiMin(struct estimate_thetas_params *params){
	int tries = 0;
	double likelyHood = 0.0;
	double bestLikelyHood = SCREWUPVALUE;

	gsl_vector *xInit = gsl_vector_alloc(params->options->nthetas);
	gsl_vector *xFinal = gsl_vector_alloc(params->options->nthetas);
	gsl_vector *xBest = gsl_vector_alloc(params->options->nthetas);
	double *tempVec = malloc(sizeof(double)*params->options->nthetas);
	pthread_t self = pthread_self(); // this is your thread id

	// init our vectors
	gsl_vector_set_zero(xBest);
	gsl_vector_set_zero(xFinal);
	set_random_initial_value(params->random_number, xInit, params->options->grad_ranges, params->options->nthetas);

	// generate the regression details for our params object
	params->h_matrix = gsl_matrix_alloc(params->options->nmodel_points, params->options->nregression_fns);
	makeHMatrix(params->h_matrix, params->the_model->xmodel,params->options->nmodel_points, params->options->nparams, params->options->nregression_fns);

	/* printf("max_tries = %d\n", params->max_tries); */
	while(tries < params->max_tries) {
		// do the actual optimization using MultiMin
		doOptimizeMultiMiN(&evalFnMulti, &getGradientExactMulti, params->options->grad_ranges, xInit, xFinal, params->options->nthetas, 500, (void*)params);
		
		copy_gslvec_vec(xFinal, tempVec, params->options->nthetas);
		likelihood = evalFnMulti(tempVec, params->options->nthetas, (void*)params);
		
		/*annoying!
		 *fprintPt(stdout, self);
		 *printf(":L = %g:try = %d\n", likelihood, tries);
		 */
		
		if(likelihood > bestLikelyHood && (isnan(likelihood) == 0 && isinf(likelihood) == 0)){
			bestLikelyHood = likelihood;
			gsl_vector_memcpy(xBest, xFinal);

			/* fprintPt(stdout, self); */
			/* printf(":best = %g\n", bestLikelyHood); */
		}
		tries++;
		set_random_initial_value(params->random_number, xInit, params->options->grad_ranges, params->options->nthetas);
		gsl_vector_set_zero(xFinal);
	}

	//printf("Final Best = %g\n", bestLikleyHood);
	if(bestLikelyHood == SCREWUPVALUE) {
		fprintf(stderr, "maximisation didn't work at all, relax your ranges\n");
	}

	gsl_vector_memcpy(params->the_model->thetas, xBest);
	
	gsl_vector_free(xInit);
	gsl_vector_free(xFinal);
	gsl_vector_free(xBest);
	free(tempVec);

}

/**
 * to specify a system for multimin we should supply:
 * the evaluation function itself: double (*f)(gsl_vector* x, void *params)
 * the grad fn: void (*df)(gsl_vector* x, void* params, gsl_vector*g) - sets the grad vector
 * a combo fn: void (*fdf)(gsl_vector* x, void* params, double* f, gsl_vector*g) - sets f and g
 * 
 */

// the actual eval fn
double evalFnMulti(gsl_vector* x, void* params){
}

double gradFnMulti(gsl_vector* x, void* params, 

// actually carry out the multimin steps
void doBoundedMultiMin( double(*fn)(gsl_vector*, void*),													\
												void(*gradientFn)(gsl_vector*,void*, double*, gsl_vector*),
										gsl_matrix* ranges, 
										gsl_vector *xkInit, gsl_vector* xkFinal, int nparams, int nsteps, void* args){

}





/**
 * init the vector x to a set of random values which are sampled from a 
 * uniform dist (from rand) generated on the ranges given by the 
 * matrix ranges
 */
void set_random_initial_value(gsl_rng* rand, gsl_vector* x, gsl_matrix* ranges,int  nthetas){
	int i;
	double range_min; 
	double range_max;
	double the_value;
	
	for(i = 0; i < nthetas; i++){
		range_min = gsl_matrix_get(ranges, i, 0);
		range_max = gsl_matrix_get(ranges, i, 1);
		// set the input vector to a random value in the range
		the_value = gsl_rng_uniform(rand) * (range_max - range_min) + range_min;
		//printf("theta %d set to %g\n", i, the_value);
		//gsl_vector_set(x, gsl_rng_uniform(rand)*(range_max-range_min)+range_min, i);
		gsl_vector_set(x, i, the_value);
	}
}



