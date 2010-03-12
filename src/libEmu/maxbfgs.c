#include "maxbfgs.h"
#include "pthread.h"

#define _CHOLDECOMP

/** 
 * @file 
 * @author Chris Coleman-Smith cec24@phy.duke.edu
 * @version 0.1
 * @section DESCRIPTION
 * 
 * contains the routines to attempt the maximum likelyhood estimation
 * using the bfgs gradient based method.
1 * At present we need the gradient of the covariance function to be explicitally 
 * provided as a function pointer, but this could really be derived by automatic-differentiation
 */


void maxWithBFGS(gsl_rng *rand, int max_tries, int nsteps, gsl_matrix *ranges, gsl_matrix* xmodel,
								 gsl_vector *trainingvector, gsl_vector* thetas, int nmodel_points, int nthetas, int nparams){

	// a wrapper to evaluate, relies on inherited scope!
	double evalFn(gsl_vector* xk, int nt){
		gsl_matrix* covariance_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points);
		gsl_matrix* cinverse = gsl_matrix_alloc(nmodel_points, nmodel_points);
		gsl_matrix* temp_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points);
		gsl_permutation *c_LU_permutation = gsl_permutation_alloc(nmodel_points);
		double cinverse_det = 0.0;
		double temp_val = 0.0;
		int lu_signum = 0, i;
		int cholesky_test = 0;
		
		// make the covariance matrix 
		// using the random initial conditions! (xold not thetas)
		makeCovMatrix(covariance_matrix, xmodel, xk, nmodel_points, nthetas, nparams);
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
		for(i = 0; i < nmodel_points; i++){
			cinverse_det *= gsl_matrix_get(temp_matrix, i,i);
		}
		cinverse_det = cinverse_det * cinverse_det;
		#endif
		
		// temp_val is now the likelyhood for this answer
		temp_val = getLogLikelyhood(cinverse, cinverse_det, xmodel, trainingvector, xk, nmodel_points, nthetas, nt);

		
		//fprintf(stderr,"L:%f\n", temp_val);
		
		gsl_matrix_free(covariance_matrix);
		gsl_matrix_free(cinverse);
		gsl_matrix_free(temp_matrix);
		gsl_permutation_free(c_LU_permutation);
		return(-1*temp_val);
	}
	
	/* don't use this if you're not using the power-exp covfn, since estimator.c:getGradient is HARD CODED
	 * to only return the grad of the power-exp 
	 * you are WARNED
	 */
	/* // a wrapper to calculate the gradient */
	/* void gradFn( double(*fn)(gsl_vector*, int), gsl_vector* xk, gsl_vector* gradient, int nt){ */
	/* 	// now the function estimator.c:getGradient in its wisdom only returns the gradient along one  */
	/* 	// dimension so we have to call it nparams times, joy */
	/* 	gsl_matrix* covariance_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points); */
	/* 	gsl_matrix* cinverse = gsl_matrix_alloc(nmodel_points, nmodel_points); */
	/* 	gsl_matrix* temp_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points); */
	/* 	gsl_permutation *c_LU_permutation = gsl_permutation_alloc(nmodel_points); */
	/* 	int cholesky_test = 0; */
	/* 	double cinverse_det = 0.0; */
	/* 	double gradTemp = 0.0; */
	/* 	int lu_signum = 0; */
	/* 	int i; */
		
	/* 	// make the covariance matrix  */
	/* 	// using the random initial conditions! (xold not thetas) */
	/* 	makeCovMatrix(covariance_matrix, xmodel, xk, nmodel_points, nthetas, nparams); */
	/* 	gsl_matrix_memcpy(temp_matrix, covariance_matrix); */
	/* 	#ifndef _CHOLDECOMP */
	/* 	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum); */
	/* 	gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse); // now we have the inverse */
	/* 	#else */
	/* 	cholesky_test = gsl_linalg_cholesky_decomp(temp_matrix); */
	/* 	if(cholesky_test == GSL_EDOM){ */
	/* 	fprintf(stderr, "trying to cholesky a non postive def matrix, sorry...\n"); */
	/* 	exit(1); */
	/* 	} */
	/* 	gsl_linalg_cholesky_invert(temp_matrix); */
	/* 	gsl_matrix_memcpy(cinverse, temp_matrix); */
	/* 	#endif */

	/* 		/\* oh man how could you even begin to think this is ok if you're not using the gaussian? */
	/* 			 fucking hell */
	/* 		*\/ */
	/* 		for(i = 0; i < nparams; i++){ */
	/* 			gradTemp = getGradient(cinverse, xmodel, trainingvector, xk, i, nmodel_points,  nt, nparams); */
	/* 			gsl_vector_set(gradient, i, gradTemp); */
	/* 		} */
			
	/* 	gsl_matrix_free(covariance_matrix); */
	/* 	gsl_matrix_free(cinverse); */
	/* 	gsl_matrix_free(temp_matrix); */
	/* 	gsl_permutation_free(c_LU_permutation); */
	/* } */
	
	int tries  = 0;
	double likelyHood = 0.0;
	double bestLikleyHood = -100;
	gsl_vector *xInit = gsl_vector_alloc(nthetas);
	gsl_vector *xFinal = gsl_vector_alloc(nthetas);
	gsl_vector *xBest = gsl_vector_alloc(nthetas);
	gsl_matrix *Binit = gsl_matrix_alloc(nthetas, nthetas);

	gsl_vector_set_zero(xBest);
	gsl_vector_set_zero(xFinal);
	
	set_random_initial_value(rand, xInit, ranges, nthetas);
	
	gsl_matrix_set_identity(Binit);

	while(tries < max_tries){
		doSimpleBFGS(&evalFn, &getGradientNumeric, ranges,  xInit, xFinal, Binit, nthetas, nsteps);
		

		// now if the function fell out of range it'll come back here with some stupid
		// answer, so we need to test the range again and then do something
		if(test_range_vector(xFinal, ranges, nthetas) == 0){
			// it's in range, so we're ok
			//vector_print(xFinal, nthetas);
			likelyHood = -1*evalFn(xFinal, nthetas);
			printf("%lu:L = %g\n", pthread_self(),likelyHood);
			if(likelyHood > bestLikleyHood){
				bestLikleyHood = likelyHood;
				gsl_vector_memcpy(xBest, xFinal);

				printf("%lu:best = %g\n", pthread_self(),bestLikleyHood);
			}
		} 
		tries++;
		
		// tidy up to run again
		set_random_initial_value(rand, xInit, ranges, nthetas);
		gsl_matrix_set_identity(Binit);
		gsl_vector_set_zero(xFinal);
	}

	printf(" Final best = %g\n", bestLikleyHood);

	if(bestLikleyHood == -100){
		fprintf(stderr, "maximisation didn't work at all, relax your ranges\n");
	}

	// should now save the best values
	gsl_vector_memcpy(thetas, xBest);

	gsl_matrix_free(Binit);
	gsl_vector_free(xFinal);
	gsl_vector_free(xInit);
}

