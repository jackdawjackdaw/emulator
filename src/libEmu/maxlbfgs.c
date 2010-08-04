#include "maxlbfgs.h"

#define SCREWUPVALUE -20000

/**
 * @file 
 * @author Chris Coleman-Smith cec24@phy.duke.edu
 * @version 1.1
 * @section DESCRIPTION
 * 
 * contains routines to use the lbfgs maximisation method
 * this is a constrained modified-newtons method, see the 
 * file lbfgs/routines.f and the files lbfgs/drivers[1-3].f for 
 * more information about how it works and what it does 
 *
 * the maxWithBLAH interface of rand, max_tries etc has become kind of standard, see the neldermead, and 
 * plain bfgs methods.
 *  
 * April '10 updated with support for linear regression model
 */
void maxWithLBFGS(struct estimate_thetas_params *params){
									
	/* 
	 * this is the main routine for driving the lbfgs method
	 */

	int tries = 0;
	double likelyHood = 0.0;
	double bestLikelyHood = SCREWUPVALUE;
	gsl_vector *xInit = gsl_vector_alloc(params->options->nthetas);
	gsl_vector *xFinal = gsl_vector_alloc(params->options->nthetas);
	gsl_vector *xBest = gsl_vector_alloc(params->options->nthetas);
	double *tempVec = malloc(sizeof(double)*params->options->nthetas);

	params->h_matrix = gsl_matrix_alloc(params->options->nmodel_points, params->options->nregression_fns);

	makeHMatrix(params->h_matrix, params->the_model->xmodel,params->options->nmodel_points, params->options->nparams, params->options->nregression_fns);
	
	gsl_vector_set_zero(xBest);
	gsl_vector_set_zero(xFinal);
	set_random_initial_value(params->random_number, xInit, params->options->grad_ranges, params->options->nthetas);
	
	printf("max_tries = %d\n", params->max_tries);
	while(tries < params->max_tries) {
		doBoundedBFGS(&evalFnLBFGS, &getGradientNumericLBFGS, params->options->grad_ranges, xInit, xFinal, params->options->nthetas, 500, (void*)params);
		
		copy_gslvec_vec(xFinal, tempVec, params->options->nthetas);
		likelyHood = -1*evalFnLBFGS(tempVec, params->options->nthetas, (void*)params);
		printf("%lu:L = %g\n", (unsigned long)pthread_self(), likelyHood);
		printf("try = %d\n", tries);
		if(likelyHood > bestLikelyHood && (isnan(likelyHood) == 0 && isinf(likelyHood) == 0)){
			bestLikelyHood = likelyHood;
			gsl_vector_memcpy(xBest, xFinal);

			printf("%lu:best = %g\n", (unsigned long)pthread_self(), bestLikelyHood);
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
 * returns the loglikleyhood for the final step of estimate_thread_function
 *
 * used in estimate_thread_function to set "my_theta_val"
 * this is a wrapper for evalFnLBFGS
 * 
 * see estimate_threaded.h for a spec of estimate_thetas_params
 */
double evalLikelyhoodLBFGS_struct(struct estimate_thetas_params *params){
	int i;
	double *xinput = MallocChecked(sizeof(double)*params->options->nthetas);
	double likelihood = 0.0;

	for(i = 0; i < params->options->nthetas; i++) 	/* setup xinput*/
		xinput[i] = gsl_vector_get(params->the_model->thetas, i);
	
	likelihood = evalFnLBFGS(xinput, params->options->nthetas, params);;

	free(xinput);
	return(likelihood);
}


	
/**
 * Calculates the loglikelihood for a given set of thetas
 * 
 * @params xinput -> a flat double vector of the position in parameter space to be evaluated
 * @params nthetas -> length of xinuput, the number of hyperparams in the statistical model
 * @params args -> a voided (struct evalFnLBFGSArgs)
 * @return the loglikleyhood of xinput
 */
double evalFnLBFGS(double *xinput, int nthetas, void* args){
	struct estimate_thetas_params *params = (struct estimate_thetas_params*) args;
	
	gsl_matrix* covariance_matrix = gsl_matrix_alloc(params->options->nmodel_points, params->options->nmodel_points);
	gsl_matrix* cinverse = gsl_matrix_alloc(params->options->nmodel_points, params->options->nmodel_points);
	gsl_matrix* temp_matrix = gsl_matrix_alloc(params->options->nmodel_points, params->options->nmodel_points);

	gsl_vector *xk = gsl_vector_alloc(nthetas);
	gsl_permutation *c_LU_permutation = gsl_permutation_alloc(params->options->nmodel_points);	
	double determinant_c = 0.0;
	double temp_val = 0.0;
	int lu_signum = 0, i;
	int cholesky_test = 0;
		
	// copy the given double vector into the gsl vector
	copy_vec_gslvec(xinput, xk, params->options->nthetas );

	// make the covariance matrix 
	// using the random initial conditions! (xold not thetas)
	makeCovMatrix(covariance_matrix, params->the_model->xmodel, xk, params->options->nmodel_points, nthetas, params->options->nparams, params->options->covariance_fn);
	gsl_matrix_memcpy(temp_matrix, covariance_matrix);
	//print_matrix(temp_matrix, params->options->nmodel_points, params->options->nmodel_points);

	// do a cholesky decomp of the cov matrix, LU is not stable for ill conditioned matrices
	cholesky_test = gsl_linalg_cholesky_decomp(temp_matrix);
	if(cholesky_test == GSL_EDOM){
		fprintf(stderr, "trying to cholesky a non postive def matrix, sorry...\n");
		exit(1);
	}
	// find the determinant and then invert 
	// the determinant is just the trace squared
	determinant_c = 1.0;
	for(i = 0; i < params->options->nmodel_points; i++)
		determinant_c *= gsl_matrix_get(temp_matrix, i, i);
	determinant_c = determinant_c * determinant_c;

	//printf("det CHOL:%g\n", determinant_c);	
	gsl_linalg_cholesky_invert(temp_matrix);
	gsl_matrix_memcpy(cinverse, temp_matrix);

	// temp_val is now the likelyhood for this answer
	temp_val = getLogLikelyhood(cinverse, determinant_c, params->the_model->xmodel, params->the_model->training_vector, xk, params->h_matrix, params->options->nmodel_points, params->options->nthetas, params->options->nparams, params->options->nregression_fns);
		
	/* fprintf(stderr,"L:%f\n", temp_val);					 */
	/* print_vector_quiet(xk, nthetas); */
		
	gsl_matrix_free(covariance_matrix);
	gsl_matrix_free(cinverse);
	gsl_matrix_free(temp_matrix);
	gsl_permutation_free(c_LU_permutation);
	gsl_vector_free(xk);
	return(-1*temp_val);
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
