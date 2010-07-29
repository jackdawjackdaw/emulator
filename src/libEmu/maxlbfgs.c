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

void maxWithLBFGS(gsl_rng *rand, int max_tries, int nsteps, gsl_matrix *ranges, gsl_matrix *xmodel,
									gsl_vector *trainingvector, gsl_vector *thetas, int nmodel_points, int nthetas, int nparams, int nregression_fns){
	/* 
	 * this is the main routine for driving the lbfgs method
	 */

	int tries = 0;
	double likelyHood = 0.0;
	double bestLikelyHood = SCREWUPVALUE;
	gsl_vector *xInit = gsl_vector_alloc(nthetas);
	gsl_vector *xFinal = gsl_vector_alloc(nthetas);
	gsl_vector *xBest = gsl_vector_alloc(nthetas);
	gsl_matrix *h_matrix = gsl_matrix_alloc(nmodel_points, nregression_fns);
	double *tempVec = malloc(sizeof(double)*nthetas);

	makeHMatrix(h_matrix, xmodel,nmodel_points, nparams, nregression_fns);

	/* setup the arguments which will be used in the eval and gradient functions */
	struct evalFnLBFGSArgs eval_fn_args;
	eval_fn_args.nparams = nparams;
	eval_fn_args.nmodel_points = nmodel_points;
	eval_fn_args.xmodel = gsl_matrix_alloc(nmodel_points, nparams);
	eval_fn_args.training_vector = gsl_vector_alloc(nmodel_points);
	eval_fn_args.h_matrix = gsl_matrix_alloc(nmodel_points, nregression_fns);
	eval_fn_args.nregression_fns = nregression_fns;
	gsl_matrix_memcpy(eval_fn_args.xmodel, xmodel);
	gsl_matrix_memcpy(eval_fn_args.h_matrix, h_matrix);
	gsl_vector_memcpy(eval_fn_args.training_vector, trainingvector);


	gsl_vector_set_zero(xBest);
	gsl_vector_set_zero(xFinal);
	set_random_initial_value(rand, xInit, ranges, nthetas);
	
	printf("max_tries = %d\n", max_tries);
	while(tries < max_tries) {
		doBoundedBFGS(&evalFnLBFGS, &getGradientNumericLBFGS, ranges, xInit, xFinal, nthetas, 500, (void*)&eval_fn_args);
		
		copy_gslvec_vec(xFinal, tempVec, nthetas);
		likelyHood = -1*evalFnLBFGS(tempVec, nthetas, (void*)&eval_fn_args);
		printf("%lu:L = %g\n", (unsigned long)pthread_self(), likelyHood);
		printf("try = %d\n", tries);
		if(likelyHood > bestLikelyHood && (isnan(likelyHood) == 0 && isinf(likelyHood) == 0)){
			bestLikelyHood = likelyHood;
			gsl_vector_memcpy(xBest, xFinal);

			printf("%lu:best = %g\n", (unsigned long)pthread_self(), bestLikelyHood);
		}
		tries++;
		set_random_initial_value(rand, xInit, ranges, nthetas);
		gsl_vector_set_zero(xFinal);
	}
	
	//printf("Final Best = %g\n", bestLikleyHood);
	if(bestLikelyHood == SCREWUPVALUE) {
		fprintf(stderr, "maximisation didn't work at all, relax your ranges\n");
	}

	gsl_vector_memcpy(thetas, xBest);
	
	gsl_matrix_free(eval_fn_args.xmodel);
	gsl_matrix_free(eval_fn_args.h_matrix);
	gsl_vector_free(eval_fn_args.training_vector);
	gsl_matrix_free(h_matrix);
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
	double *xinput = MallocChecked(sizeof(double)*params->nthetas);
	double likelihood = 0.0;
	struct evalFnLBFGSArgs arguments;

	for(i = 0; i < params->nthetas; i++) 	/* setup xinput*/
		xinput[i] = gsl_vector_get(params->thetas, i);
	
	/* setup arguments */
	setup_evalFnLBFGSArgs(&arguments, params);

	likelihood = evalFnLBFGS(xinput, params->nthetas, &arguments);

	gsl_matrix_free(arguments.xmodel);
	gsl_matrix_free(arguments.h_matrix);
	gsl_vector_free(arguments.training_vector);
	free(xinput);
	return(likelihood);
}
	

/**
 * inits a evalFnLBFGSArgs structure from an estimate_thetas_params structure
 * 
 * struct evalFnLBFGSArgs{
 * 	int nparams;
 *	int nmodel_points;
 * 	int nregression_fns;
 *	gsl_matrix* xmodel;
 * 	gsl_vector* training_vector;
 *	gsl_matrix* h_matrix;
 * } evalFnLBFGSArgs;
 *
 */
void setup_evalFnLBFGSArgs(struct evalFnLBFGSArgs *arguments, struct estimate_thetas_params *params){
	/* setup the gsl vecs/matrices  */
	arguments->xmodel = gsl_matrix_alloc(params->nmodel_points, params->nparams);
	arguments->training_vector = gsl_vector_alloc(params->nmodel_points);
	arguments->h_matrix = gsl_matrix_alloc(params->nmodel_points, params->nregression_fns);

	arguments->nparams = params->nparams;
	arguments->nmodel_points = params->nmodel_points;
	arguments->nregression_fns = params->nregression_fns;
	
	gsl_matrix_memcpy(arguments->xmodel, params->model_input);
	gsl_vector_memcpy(arguments->training_vector, params->training_vector);

	/* and finally setup the hmatrix (for the regression)*/
	makeHMatrix(arguments->h_matrix, arguments->xmodel, arguments->nmodel_points, arguments->nparams, arguments->nregression_fns);
		 

	
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
	struct evalFnLBFGSArgs *params = (struct evalFnLBFGSArgs*) args;
	
	gsl_matrix* covariance_matrix = gsl_matrix_alloc(params->nmodel_points, params->nmodel_points);
	gsl_matrix* cinverse = gsl_matrix_alloc(params->nmodel_points, params->nmodel_points);
	gsl_matrix* temp_matrix = gsl_matrix_alloc(params->nmodel_points, params->nmodel_points);

	gsl_vector *xk = gsl_vector_alloc(nthetas);
	gsl_permutation *c_LU_permutation = gsl_permutation_alloc(params->nmodel_points);	
	double determinant_c = 0.0;
	double temp_val = 0.0;
	int lu_signum = 0, i;
	int cholesky_test = 0;
		
	// copy the given double vector into the gsl vector
	copy_vec_gslvec(xinput, xk, nthetas );

	// make the covariance matrix 
	// using the random initial conditions! (xold not thetas)
	makeCovMatrix(covariance_matrix, params->xmodel, xk, params->nmodel_points, nthetas, params->nparams);
	gsl_matrix_memcpy(temp_matrix, covariance_matrix);

#define _CHOLDECOMP
#ifndef _CHOLDECOMP
	// this is not stable it seems
	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
	determinant_c = gsl_linalg_LU_decomp(temp_matrix, lu_signum);
	gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse); // now we have the inverse

#else 
	// do cholesky (should be twice as fast)
	// do the decomp and then run along
	cholesky_test = gsl_linalg_cholesky_decomp(temp_matrix);
	if(cholesky_test == GSL_EDOM){
		fprintf(stderr, "trying to cholesky a non postive def matrix, sorry...\n");
		exit(1);
	}
	// find the determinant and then invert 
	// the determinant is just the trace squared
	determinant_c = 1.0;
	for(i = 0; i < params->nmodel_points; i++)
		determinant_c *= gsl_matrix_get(temp_matrix, i, i);
	determinant_c = determinant_c * determinant_c;

	//printf("det CHOL:%g\n", determinant_c);	
	gsl_linalg_cholesky_invert(temp_matrix);
	gsl_matrix_memcpy(cinverse, temp_matrix);

#endif
		
	// temp_val is now the likelyhood for this answer
	temp_val = getLogLikelyhood(cinverse, determinant_c, params->xmodel, params->training_vector, xk, params->h_matrix, params->nmodel_points, nthetas, params->nparams, params->nregression_fns);

		
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
