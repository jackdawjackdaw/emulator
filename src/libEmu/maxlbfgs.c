#include "maxlbfgs.h"

#define SCREWUPVALUE -100

/**
 * @file 
 * @author Chris Coleman-Smith cec24@phy.duke.edu
 * @version 1
 * @section DESCRIPTION
 * 
 * contains routines to use the lbfgs maximisation method
 * this is a constrained modified-newtons method, see the 
 * file lbfgs/routines.f and the files lbfgs/drivers[1-3].f for 
 * more information about how it works and what it does 
 *
 * the maxWithBLAH interface of rand, max_tries etc has become kind of standard, see the neldermead, and 
 * plain bfgs methods.
 */

void maxWithLBFGS(gsl_rng *rand, int max_tries, int nsteps, gsl_matrix *ranges, gsl_matrix *xmodel,
									gsl_vector *trainingvector, gsl_vector *thetas, int nmodel_points, int nthetas, int nparams){
	/* 
	 * this is the main routine for driving the lbfgs method
	 */

	int tries = 0;
	double likelyHood = 0.0;
	double bestLikelyHood = SCREWUPVALUE;
	gsl_vector *xInit = gsl_vector_alloc(nthetas);
	gsl_vector *xFinal = gsl_vector_alloc(nthetas);
	gsl_vector *xBest = gsl_vector_alloc(nthetas);
	double *tempVec = malloc(sizeof(double)*nthetas);
	

	/* setup the arguments which will be used in the eval and gradient functions */
	struct evalFnLBFGSArgs eval_fn_args;
	eval_fn_args.nparams = nparams;
	eval_fn_args.nmodel_points = nmodel_points;
	eval_fn_args.xmodel = gsl_matrix_alloc(nmodel_points, nparams);
	eval_fn_args.training_vector = gsl_vector_alloc(nmodel_points);
	gsl_matrix_memcpy(eval_fn_args.xmodel, xmodel);
	gsl_vector_memcpy(eval_fn_args.training_vector, trainingvector);

	gsl_vector_set_zero(xBest);
	gsl_vector_set_zero(xFinal);
	set_random_initial_value(rand, xInit, ranges, nthetas);
	
	while(tries < max_tries) {
		doBoundedBFGS(&evalFnLBFGS, &getGradientNumericLBFGS, ranges, xInit, xFinal, nthetas, nsteps, (void*)&eval_fn_args);
		
		copy_gslvec_vec(xFinal, tempVec, nthetas);
		likelyHood = -1*evalFnLBFGS(tempVec, nthetas, (void*)&eval_fn_args);
		printf("%lu:L = %g\n", pthread_self(), likelyHood);
		if(likelyHood > bestLikelyHood){
			bestLikelyHood = likelyHood;
			gsl_vector_memcpy(xBest, xFinal);

			printf("%lu:best = %g\n", pthread_self(), bestLikelyHood);
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
	gsl_vector_free(eval_fn_args.training_vector);
	gsl_vector_free(xInit);
	gsl_vector_free(xFinal);
	gsl_vector_free(xBest);
	free(tempVec);
}


	

 
double evalFnLBFGS(double *xinput, int nthetas, void* args){
	struct evalFnLBFGSArgs *params = (struct evalFnLBFGSArgs*) args;
	
	gsl_matrix* covariance_matrix = gsl_matrix_alloc(params->nmodel_points, params->nmodel_points);
	gsl_matrix* cinverse = gsl_matrix_alloc(params->nmodel_points, params->nmodel_points);
	gsl_matrix* temp_matrix = gsl_matrix_alloc(params->nmodel_points, params->nmodel_points);
	gsl_vector *xk = gsl_vector_alloc(nthetas);
	gsl_permutation *c_LU_permutation = gsl_permutation_alloc(params->nmodel_points);
	double cinverse_det = 0.0;
	double temp_val = 0.0;
	int lu_signum = 0, i;
	int cholesky_test = 0;
		
	// copy the given double vector into the gsl vector
	copy_vec_gslvec(xinput, xk, nthetas );

	// make the covariance matrix 
	// using the random initial conditions! (xold not thetas)
	makeCovMatrix(covariance_matrix, params->xmodel, xk, params->nmodel_points, nthetas, params->nparams);
	gsl_matrix_memcpy(temp_matrix, covariance_matrix);


#ifndef _CHOLDECOMP
	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
	gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse); // now we have the inverse
		
	// now we want to calc the logLikelyhood and see how it is
	gsl_matrix_memcpy(temp_matrix, cinverse);
	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
	cinverse_det = gsl_linalg_LU_det(temp_matrix, lu_signum);
#else 
	cholesky_test = gsl_linalg_cholesky_decomp(temp_matrix);
	if(cholesky_test == GSL_EDOM){
		fprintf(stderr, "trying to cholesky a non postive def matrix, sorry...\n");
		exit(1);
	}
	gsl_linalg_cholesky_invert(temp_matrix);
	gsl_matrix_memcpy(cinverse, temp_matrix);
	
	// now get the determinant of the inverse
	// for a cholesky decomp matrix L.L^T = A the determinant of A is 
	// the square of the product of the diagonal elements of L	
	cholesky_test = gsl_linalg_cholesky_decomp(temp_matrix);
	if(cholesky_test == GSL_EDOM){
		fprintf(stderr, "trying to cholesky a non postive def matrix, sorry...\n");
		exit(1);
	}
	
	cinverse_det = 1.0;
	for(i = 0; i < params->nmodel_points; i++){
		cinverse_det *= gsl_matrix_get(temp_matrix, i,i);
	}
	cinverse_det = cinverse_det * cinverse_det;
#endif
		
	// temp_val is now the likelyhood for this answer
	temp_val = getLogLikelyhood(cinverse, cinverse_det, params->xmodel, params->training_vector, xk, params->nmodel_points, nthetas, params->nparams);

		
	//fprintf(stderr,"L:%f\n", temp_val);
		
	gsl_matrix_free(covariance_matrix);
	gsl_matrix_free(cinverse);
	gsl_matrix_free(temp_matrix);
	gsl_permutation_free(c_LU_permutation);
	gsl_vector_free(xk);
	return(-1*temp_val);
}