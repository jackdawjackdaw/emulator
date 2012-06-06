#include "maxmultimin.h"

/**
 * @file
 * @author C.Coleman-Smith, cec24@phy.duke.edu
 * @date oct-11
 * @version 0.
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
 * there's a lot of duplication here...
 */

/**
 * to specify a system for multimin we should supply:
 * the evaluation function itself: double (*f)(gsl_vector* x, void *params)
 * the grad fn: void (*df)(gsl_vector* x, void* params, gsl_vector*g) - sets the grad vector
 * a combo fn: void (*fdf)(gsl_vector* x, void* params, double* f, gsl_vector*g) - sets f and g
 * 
 */



#define SCREWUPVALUE -2000

/**
 * call this to max the given model (params) with multimin, this is copied 
 * from maxWithLBFGS
 *
 * we will only numerically estimate thetas_1..Nparams, the zeroth theta will instead be set 
 * by the data, using the function estimateSigma
 */
void maxWithMultiMin(struct estimate_thetas_params *params){
	int tries = 0;
	// make local bindings of commonly used values in params
	int nthetas = params->options->nthetas;
	int nmodel_points = params->options->nmodel_points;
	int nparams = params->options->nparams;
	int nregression_fns = params->options->nregression_fns;
	int i;

	double likelihood = 0.0;
	double bestLHood = SCREWUPVALUE;

	gsl_vector *xInit = gsl_vector_alloc(nthetas);
	gsl_vector *xFinal = gsl_vector_alloc(nthetas);
	gsl_vector *xBest = gsl_vector_alloc(nthetas);
	gsl_vector *xTest = gsl_vector_alloc(nthetas-1);
	double *tempVec = malloc(sizeof(double)*nthetas);
	double sigma = 0.0;
	pthread_t self = pthread_self(); // this is your thread id

	// init our vectors
	gsl_vector_set_zero(xBest);
	gsl_vector_set_zero(xFinal);
	//set_random_init_value(params->random_number, xInit, params->options->grad_ranges, nthetas);

	// generate the regression details for our params object
	params->h_matrix = gsl_matrix_alloc(nmodel_points, nregression_fns);
	makeHMatrix(params->h_matrix, params->the_model->xmodel,nmodel_points, nparams, nregression_fns);

	/* printf("max_tries = %d\n", params->max_tries); */
	while(tries < params->max_tries) {
		set_random_init_value(params->random_number, xInit, params->options->grad_ranges, nthetas);
		// do the actual optimization using MultiMin
		
		doOptimizeMultiMin(&evalFnMulti, // computes the log-likelihood 
											 &gradFnMulti,  // computes the gradient
											 &evalFnGradMulti, // does both at once
											 xInit, xFinal, (void*)params);
		
		/** 
		 * need to send a theta_vector without the amplitude if you want to 
		 * get the lhood from the evalFn directly
		 */
		for(i = 0; i < nthetas-1; i++)
			gsl_vector_set(xTest, i, gsl_vector_get(xFinal, i+1));

		likelihood = -1*evalFnMulti(xTest, (void*)params);
		
		/*annoying!
		fprintPt(stdout, self);
		printf(":L = %g:try = %d\n", likelihood, tries);
		*/
		 
		if(likelihood > bestLHood && (isnan(likelihood) == 0 && isinf(likelihood) == 0)){
			bestLHood = likelihood;
			gsl_vector_memcpy(xBest, xFinal);
			
			/* fprintPt(stdout, self); */
			/* printf(":best = %g\n", bestLikelyHood); */
		}
		tries++;

		gsl_vector_set_zero(xFinal);
	}

	//printf("Final Best = %g\n", bestLikleyHood);
	if(bestLHood == SCREWUPVALUE) {
		fprintf(stderr, "maximisation didn't work at all, relax your ranges\n");
	}

	gsl_vector_memcpy(params->the_model->thetas, xBest);
	params->lhood_current = bestLHood;
	
	gsl_vector_free(xInit);
	gsl_vector_free(xFinal);
	gsl_vector_free(xBest);
	gsl_vector_free(xTest);
	free(tempVec);

}

/**
 * estimate the training sample covariance, like estimateSigma but here the cov matrix is directly 
 * inverted. 
 * this takes a full length theta
 * 
 * @return sigma^2 or GSL_NAN if inversion fails, this is NOT log scaled
 * @param thetas, full vector of thetas {amp, nug, length_1, length_2,...} (including amplitude field)
 * @param params_in a correctly setup estimate_thetas_params structure
 * 
 */
double estimateSigmaFull(gsl_vector *thetas, void* params_in){
	double sigmaValue = 0.0;

	struct estimate_thetas_params *params = (struct estimate_thetas_params*) params_in;
	
	int nmodel_points = params->options->nmodel_points;
	int nthetas = params->options->nthetas;
	int nparams = params->options->nparams;
	int nthetas_opt = nthetas - 1;
	int i,j;
	int cholesky_test;
	gsl_matrix* xmodel = params->the_model->xmodel;
	gsl_matrix* cmatrix = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_error_handler_t *temp_handler;


	gsl_vector* theta_local = gsl_vector_alloc(nthetas);
	gsl_vector_set_zero(theta_local); // set the amp to zero
	for(i = 1; i < nthetas; ++i)
		gsl_vector_set(theta_local, i, gsl_vector_get(thetas, i-1));

	makeCovMatrix(cmatrix, xmodel, theta_local, nmodel_points, nthetas, nparams);

	temp_handler = gsl_set_error_handler_off();
	cholesky_test = gsl_linalg_cholesky_decomp(cmatrix);
	if(cholesky_test == GSL_EDOM){
		FILE *fptr;
		fprintf(stderr, "estSigmaFull\n");
		fprintf(stderr, "trying to cholesky a non postive def matrix, sorry...\n");
		fprintf(stderr, "matrix dumped to chol-err.dat\n");
		fptr = fopen("chol-err.dat", "w");
		for (i = 0; i < nmodel_points; ++i){
			for (j = 0; j < nmodel_points; ++j){
				fprintf(fptr, "%lf", gsl_matrix_get(cmatrix, i, j));				
				}
			fprintf(fptr, "\n");
		}
		
		// clean up and return
		gsl_matrix_free(cmatrix);
		gsl_vector_free(theta_local);
		// if the fn cannot be evaluated we return GSL_NAN?
		return(GSL_NAN); // umm, should maybe jump out instead?
	}
	gsl_set_error_handler(temp_handler);
	
	gsl_linalg_cholesky_invert(cmatrix);
	// now compute sigma
	sigmaValue = estimateSigma(cmatrix, params_in);

	gsl_matrix_free(cmatrix);
	gsl_vector_free(theta_local);
	return(sigmaValue);
}


/** 
 * estimate the process variance
 *
 * sigma^2 = 1/Nmodel_points (yt - \hat{mu})^t %*% cinverse %*% (yt - \hat{mu})
 * where yt is the training vector and \hat{mu} is the estimated mean 
 * this is not log scaled.
 * 
 * @return sigma^2
 * @param cinverse inverted covariance matrix, computed with the cov-fn amplitude set to 1!
 * @param params_in a correctly initialized estimate_thetas_params structure
 */
double estimateSigma(gsl_matrix* cinverse, void* params_in){
	struct estimate_thetas_params *params = (struct estimate_thetas_params*) params_in;
	
	int nmodel_points = params->options->nmodel_points;
	int nparams = params->options->nparams;
	int nregression_fns = params->options->nregression_fns;
	double estimated_mean_val;

	gsl_vector *yt_local = gsl_vector_alloc(nmodel_points);
	gsl_vector* vec_temp = gsl_vector_alloc(nmodel_points);

	gsl_vector *h_vector = gsl_vector_alloc(nregression_fns);
	gsl_vector *beta_vector = gsl_vector_alloc(nregression_fns);
	gsl_vector *estimated_mean = gsl_vector_alloc(nmodel_points);
	gsl_vector *train_sub_mean = gsl_vector_alloc(nmodel_points);

	
	gsl_vector_memcpy(yt_local, params->the_model->training_vector);

	// ptrs to useful things
	gsl_matrix* h_matrix = params->h_matrix;
	gsl_matrix* xmodel = params->the_model->xmodel;
	gsl_vector_view xmodel_row;	

	int i;
	double sigma_est = 0.0;

	/* we need to calculate the mean vector for this set of thetas 
	 * estMean[i] = hvector(training[i]).betavector
	 */
	estimateBeta(beta_vector, h_matrix, cinverse,  yt_local, nmodel_points, nregression_fns);
	for(i = 0; i < nmodel_points; i++){
		xmodel_row = gsl_matrix_row(xmodel, i);
		makeHVector(h_vector, &xmodel_row.vector, nparams);
		//print_vector_quiet(h_vector, nregression_fns);
		gsl_blas_ddot(beta_vector, h_vector, &estimated_mean_val);
		gsl_vector_set(estimated_mean, i, estimated_mean_val);
	}

	/* train_sub_mean = trainingvector- estimated_mean */
	gsl_vector_memcpy(train_sub_mean, yt_local);
	gsl_vector_sub(train_sub_mean, estimated_mean);
	
	// vec_temp = cinverse %*% (yt - est_mean)
	gsl_blas_dgemv(CblasNoTrans, 1.0, cinverse, train_sub_mean, 0.0, vec_temp);
	// sigma_est = yt_local %*% vec_temp
	gsl_blas_ddot(yt_local, vec_temp, &sigma_est);
		
	sigma_est /= (double)nmodel_points; // scale
	

	gsl_vector_free(yt_local);
	gsl_vector_free(vec_temp);
	gsl_vector_free(h_vector);
	gsl_vector_free(beta_vector);
	gsl_vector_free(estimated_mean);
	gsl_vector_free(train_sub_mean);
	return(sigma_est);
}


/**
 * 
 * the eval fn, returns the loglikelihood for a given set of vectors at the location theta_vec
 * for the model specified by the params struct
 * 
 * when we fix sigma then theta_vec starts from the nugget, as far as the max routine cares
 * so we need to be very careful to copy it into entries 1..nthetas in the 
 * local vectors
 * 
 * @return loglikelihood for theta_vec_less_amp or GSL_NAN
 * @param theta_vec_less_amp -> length nthetas-1 , {nug, theta_1, theta_2,...} (no ampltitude term)
 */
double evalFnMulti(const gsl_vector* theta_vec_less_amp, void* params_in){
	struct estimate_thetas_params *params = (struct estimate_thetas_params*) params_in;
	
	int nmodel_points = params->options->nmodel_points;
	int nthetas = params->options->nthetas;
	int nparams = params->options->nparams;
	
	// need to disable the gsl error handler
	gsl_error_handler_t *temp_handler;
	gsl_matrix* covariance_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix* cinverse = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix* temp_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points);

	gsl_vector *theta_local = gsl_vector_alloc(nthetas);
	
	double determinant_c = 0.0;
	double temp_val = 0.0;
	double sigma_est = 0.0;
	int i, j;
	int cholesky_test = 0;

	// if we're fixing sigma we want to first compute the cov matrix 
	// without the amplitude
	gsl_vector_set(theta_local, 0, 0.0); 
	for(i = 1; i < nthetas; i++) // copy in the rest of the thetas
		gsl_vector_set(theta_local, i, gsl_vector_get(theta_vec_less_amp, i-1));
	
		
	// make the covariance matrix 
	makeCovMatrix(covariance_matrix, params->the_model->xmodel, theta_local, nmodel_points, nthetas, nparams);
	gsl_matrix_memcpy(temp_matrix, covariance_matrix);

	//print_matrix(temp_matrix, params->options->nmodel_points, params->options->nmodel_points);

	temp_handler = gsl_set_error_handler_off();
	cholesky_test = gsl_linalg_cholesky_decomp(temp_matrix);

	if(cholesky_test == GSL_EDOM){
		FILE *fptr;
		fprintf(stderr, "evalFnMulti\n");
		fprintf(stderr, "trying to cholesky a non postive def matrix, sorry...\n");		fprintf(stderr, "matrix dumped to chol-err.dat\n");
		fptr = fopen("chol-err.dat", "w");
		fprintf(fptr, "#thetas: ");
		for(i = 0; i < nthetas; i++)
			fprintf(fptr, "%lf\t", gsl_vector_get(theta_local, i));
		fprintf(fptr, "\n");
					
		for(i = 0; i < nmodel_points; ++i){
			for(j = 0; j < nmodel_points; ++j){
				fprintf(fptr, "%lf ", gsl_matrix_get(temp_matrix, i, j));				
			}
			fprintf(fptr, "\n");
		}
		fclose(fptr);
		
		gsl_matrix_free(covariance_matrix);
		gsl_matrix_free(cinverse);
		gsl_matrix_free(temp_matrix);
		gsl_vector_free(theta_local);
		return(GSL_NAN);
	}
	gsl_set_error_handler(temp_handler);

	// find the determinant and then invert 
	// the determinant is just the trace squared
	determinant_c = 1.0;
	for(i = 0; i < nmodel_points; i++)
		determinant_c *= gsl_matrix_get(temp_matrix, i, i);
	determinant_c = determinant_c * determinant_c;

	//printf("det CHOL:%g\n", determinant_c);	
	gsl_linalg_cholesky_invert(temp_matrix);
	gsl_matrix_memcpy(cinverse, temp_matrix);

	// now we can estimate the process scale (sigma) and pass this through into the ll
	// note that we have to log our estimate
	sigma_est = log(estimateSigma(cinverse, params_in));
	/* fprintf(stderr, "#sigma_est(eval): %g\n", sigma_est); */
	/* print_vector_quiet(theta_local, nthetas); */
	// now we store the estimated sigma, so that we get the correct likelihood
	gsl_vector_set(theta_local, 0, sigma_est);


	// compute the likelyhood for these thetas and amplitude
	temp_val = getLogLikelyhood(cinverse, determinant_c, 
															params->the_model->xmodel, 
															params->the_model->training_vector, 
															theta_local, 
															params->h_matrix, nmodel_points, 
															nthetas, nparams, params->options->nregression_fns);
		
	
	/* printf("evalfn: %g\t", -1*temp_val); */
	/* for(i = 0; i < nthetas; i++) */
	/* 	printf("%g ", gsl_vector_get(theta_local, i)); */
	/* printf("\n"); */


	gsl_matrix_free(covariance_matrix);
	gsl_matrix_free(cinverse);
	gsl_matrix_free(temp_matrix);
	gsl_vector_free(theta_local);
	return(-1*temp_val);
}

/**
 * this gets the exact gradient, for a theta_vec
 *
 * the gradient is a vector in each of the theta directions
 * theta_1 (the nugget)
 * theta_2..theta_n (the directions of the cov-fn)
 *
 * the gradient in the length directions (theta_3 etc) is computed exactly, by functions 
 * named derivative_l_<cov_fn_name> in emulator.c. The function pointer  maxlbfgs.h:makeGradMatLength is 
 * set to the correct derivative function by a call to optstruct.c:setup_cov_fn
 * 
 * 
 * dL/dTheta = -1/2 tr (Cn^{-1} dCn/dtheta) + 1/2 tn^{t}Cn^{-1} dCn/dtheta Cn^{-1} tn
 *  where t_n are the n training points
 * 
 * this relies upon makeGradMatLength which is a fn ptr
 * 
 * note->we set theta_0 by estimation, as such we dont need to comppute the gradient in this 
 * direction
 */
void gradFnMulti(const gsl_vector* theta_vec_less_amp, void* params_in, gsl_vector * grad_vec){
	int nmpoints, nthetas, i, j;
	int nparams;
	double amp = 0.0, nug = 0.0; 
	double sigma_est = 0.0;
	double grad_temp = 0.0;
	int nthetas_opt;


	struct estimate_thetas_params *params = (struct estimate_thetas_params*)params_in;
	nmpoints = params->options->nmodel_points;
	nthetas = params->options->nthetas;
	nparams = params->options->nparams;
	nthetas_opt = nthetas - 1; // we optimize one less than we actually have

	int cholesky_test;
	gsl_error_handler_t *temp_handler;
	gsl_vector* theta_local = gsl_vector_alloc(nthetas);
	gsl_matrix* covariance_matrix = gsl_matrix_alloc(nmpoints, nmpoints);
	gsl_matrix* cinverse = gsl_matrix_alloc(nmpoints, nmpoints);
	gsl_matrix* temp_matrix = gsl_matrix_alloc(nmpoints, nmpoints);
	
	//gsl_vector_memcpy(theta_local, theta_vec);
	//we need to estimate the sigma for our thetas, start by setting it to 0
	// recall that it is defined exponentiated
	gsl_vector_set(theta_local, 0, 0.0); 
	for(i = 1; i < nthetas; i++) // copy in the rest of the thetas
		gsl_vector_set(theta_local, i, gsl_vector_get(theta_vec_less_amp, i-1));
	
	// using the random initial conditions! (xold not thetas)
	makeCovMatrix(covariance_matrix, params->the_model->xmodel, theta_local, nmpoints, nthetas, nparams);
	gsl_matrix_memcpy(temp_matrix, covariance_matrix);

	/* printf("#thetas: "); */
	/* for(i = 0; i < nthetas; ++i) */
	/* 	printf("%g ", gsl_vector_get(theta_local, i)); */
	/* printf("\n"); */

	//print_matrix(temp_matrix, nmpoints, nmpoints);

	// do a cholesky decomp of the cov matrix
	temp_handler = gsl_set_error_handler_off();
	cholesky_test = gsl_linalg_cholesky_decomp(temp_matrix);
	if(cholesky_test == GSL_EDOM){
		FILE *fptr;
		fprintf(stderr, "gradFnMulti\n");
		fprintf(stderr, "trying to cholesky a non postive def matrix, sorry...\n");
		fprintf(stderr, "matrix dumped to chol-err.dat\n");
		fptr = fopen("chol-err.dat", "w");
		for (i = 0; i < nmpoints; ++i){
			for (j = 0; j < nmpoints; ++j){
				fprintf(fptr, "%lf", gsl_matrix_get(temp_matrix, i, j));				
				}
			fprintf(fptr, "\n");
		}
		
		// lets clean up and jump out of here
		gsl_matrix_free(covariance_matrix);
		gsl_matrix_free(cinverse);
		gsl_matrix_free(temp_matrix);
		gsl_vector_free(theta_local);

		// if the fn cannot be evaluated we return GSL_NAN?
		//return(GSL_NAN); // umm, should maybe jump out instead?
		// FIXME: do somthing to indicate failure.
		return; //void function.
	}
	gsl_set_error_handler(temp_handler);

	gsl_linalg_cholesky_invert(temp_matrix);
	gsl_matrix_memcpy(cinverse, temp_matrix);
	
	// now we can estimate sigma again  (note that we have to log it)
	sigma_est = log(estimateSigma(cinverse, params_in));
	/* fprintf(stderr, "#sigma_est(grad): %g\n", sigma_est); */
	// now we store the estimated sigma, so that we get the correct likelihood
	gsl_vector_set(theta_local, 0, sigma_est);

	// C(x,y) = c(x,y) + theta_1
	// the first component, dC/dN_0 = 1_delta(ij) (easssy)
	// the rest: dC/dN_{h} = amp * dC/d\theta

	
	// first, get the scale parameter and the nugget
	amp = exp(gsl_vector_get(theta_local, 0)); 
	nug = exp(gsl_vector_get(theta_local, 1)); 
	gsl_matrix_memcpy(temp_matrix, covariance_matrix);


	// first the nugget term, if we don't log scale we set this to the identity
	gsl_matrix_set_identity(temp_matrix);
	// since we're log scaling the nugget we need to set the cov-matrix accordingly
	gsl_matrix_scale(temp_matrix, nug);

 	/* printf("#dC/dtheta_1\n"); */
	/* print_matrix(temp_matrix, nmpoints, nmpoints); */

	grad_temp= -1.0*getGradientCn(temp_matrix, cinverse, params->the_model->training_vector, nmpoints,  nthetas);
	gsl_vector_set(grad_vec, 0, grad_temp);
	
	for(i = 2; i < nthetas; i++){
		makeGradMatLength(temp_matrix, params->the_model->xmodel, gsl_vector_get(theta_local, i) , i, nmpoints, nparams);
		gsl_matrix_scale(temp_matrix, amp);
		/* printf("#dC/dtheta_2\n"); */
		/* print_matrix(temp_matrix, nmpoints, nmpoints); */
		grad_temp = -1.0*getGradientCn(temp_matrix, cinverse, params->the_model->training_vector, nmpoints,  nthetas);
		// we have to offset the gradient vector locations by -1
		gsl_vector_set(grad_vec, i-1, grad_temp);
	}


	/* printf("# gradient: "); */
	/* for(i = 0; i < nthetas_opt; ++i) */
	/* 	printf("%g ", gsl_vector_get(grad_vec,i)); */
	/* printf("\n"); */

	gsl_vector_free(theta_local);
	gsl_matrix_free(covariance_matrix);
	gsl_matrix_free(cinverse);
	gsl_matrix_free(temp_matrix);
}

/**
 * 
 * computes the gradient for a given hyper-parameter direction, the important arguments here are:
 * we rely on the CONVENTION that all cov-fns are laid out as follows:
 * C(x,y) = theta_0 * c(x,y, theta_2,theta_3,...) + theta_1 
 * i.e each cov fn has a scale given by the first hyper param and a nugget term given by the 
 * second. 
 * this is important when we compute the gradient
 * 
 * 
 * @param dCdtheta a matrix of the covariance function at each point differentiated wrt to the theta 
 * that we want the gradient in.
 * @param cinverse the inverted cov matrix
 * @return the value of the gradient in this direction
 * 
 * grad = -1/2*trace( cinverse %*% dCdtheta ) + 1/2 tn^{t} %*% cinverse %*% dCdtheta %*% cinverse %*% tn 
 *
 * see Rasumussen eqn 5.9 for a little simplification of the algebra here
 */
double getGradientCn(gsl_matrix * dCdtheta, gsl_matrix *cinverse,  gsl_vector* training_vector,int nmodel_points, int nthetas){
	double grad = 0;
	double trace = 0;
	int i;
	// this is going to be Cn^{-1}*dCn/dTheta
	gsl_matrix *temp = gsl_matrix_calloc(nmodel_points, nmodel_points);
	// this is Cn{-1} * training
	gsl_vector *v = gsl_vector_alloc(nmodel_points);
	gsl_vector *w = gsl_vector_alloc(nmodel_points);
	/* — Function: int gsl_blas_dsymm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C) */
	//gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, cinverse, dCdtheta, 0.0, temp);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, cinverse, dCdtheta, 0.0, temp);

	
	for(i = 0; i < nmodel_points; i++)
		trace += gsl_matrix_get(temp, i, i);
	trace *= -(0.5);
	
	//printf("trace = %lf\n", trace);

	/* — Function: int gsl_blas_dsymv (CBLAS_UPLO_t Uplo, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y) */
	//gsl_blas_dsymv(CblasUpper, 0.5, cinverse, training_vector, 0.0, v); // v = cinverse %*% training_vector
	gsl_blas_dgemv(CblasNoTrans, 0.5, cinverse, training_vector, 0.0, v);
	//gsl_blas_dsymv(CblasUpper, 1.0, temp, v, 0.0, w); // w = (cinverse %*% dCdthtea) %*% v
	gsl_blas_dgemv(CblasNoTrans, 1.0, temp, v, 0.0, w); // w = (cinverse %*% dCdthtea) %*% v
	
	gsl_blas_ddot(training_vector, w, &grad);

	//printf("grad = %lf\n", grad);
	
	grad += trace;

	gsl_matrix_free(temp);
	gsl_vector_free(v);
	gsl_vector_free(w);
	return(grad);
}


/**
 * compute the function value and the gradient vector in one call
 * required by the multimin lib
 */
void evalFnGradMulti(const gsl_vector* theta_vec, void* params, double* fnval, gsl_vector * grad_vec){
		*fnval = evalFnMulti(theta_vec, params);
	gradFnMulti(theta_vec, params, grad_vec);
}



/**
 * minimize the supplied functions, taking thetaInit as a starting point.
 * 
 * @return thetaFinal is set to the full ntheta {amp, nug, theta_1, theta_2,...} vector.
 * @param double(*fn)(...) the evaluation function
 * @param void(*gradientFn)(...) computes the gradient
 * @param void(*fnGradFn)(...) computes the grad and the fn value
 * @param thetaInit the initial location in the optimization space, nthetas long (includes amp term)
 * @param thetaFinal the final location, nthetas long (includes estimated amp term)
 */
void doOptimizeMultiMin( double(*fn)(const gsl_vector*, void*),													\
												void(*gradientFn)(const gsl_vector*,void*, gsl_vector*),
												void(*fnGradFn)(const gsl_vector*, void*, double*, gsl_vector*),
												gsl_vector *thetaInit, gsl_vector* thetaFinal, void* args){
	
	struct estimate_thetas_params *params = (struct estimate_thetas_params*)args;

	int nparams = params->options->nparams;
	int nmpoints = params->options->nmodel_points;
	int nthetas = params->options->nthetas;
	int nthetas_opt = nthetas - 1; 
	int status;
	int stepcount = 0;
	int stepmax = 30; // tune this if a huge number of steps seem to be happening
	int i;

	double sigma_final = 0.0;

	double stepSizeInit = 0.1;
	double tolerance = 0.05; // sets the accuracy of the line search. 
	double fnValue = 0.0;
	double norm = 0.0; // norm of the gradient
	/*
	 * we stop when |g| < epsAbs, not sure how to set this yet
	 */
	double epsAbs = 0.02; 

	/* fprintf(stderr, "#doOptimizeMultiMin: nthetas %d\tnmpoints %d\tnparams %d\n",  */
	/* 				nthetas, nmpoints, nparams); */
	
	gsl_vector *tempTest = gsl_vector_alloc(nthetas); // get values out of the min to check on progress
	/*
	 * because we're fixing the scale with the data we're actually optimizing one less
	 * theta than in the unconstrained case
	 * this means we need to carefully copy in the original theta values which 
	 * correspond to the nugget and length scales while avoiding copying in the original amp
	 */
	gsl_vector *thetaTemp = gsl_vector_alloc(nthetas_opt);
	for(i = 0; i < nthetas_opt; i ++){
		//gsl_vector_memcpy(thetaTemp, thetaInit);
		gsl_vector_set(thetaTemp, i, gsl_vector_get(thetaInit, i+1));
	}

	// our multimin fn
	gsl_multimin_function_fdf multiminFn;
	multiminFn.n = nthetas_opt; // the number of dimensions
	multiminFn.f = fn;
	multiminFn.df = gradientFn;
	multiminFn.fdf = fnGradFn;
	multiminFn.params = args;

	// the multiminimizer
	const gsl_multimin_fdfminimizer_type *min_type = gsl_multimin_fdfminimizer_vector_bfgs2; //set this properly
	//const gsl_multimin_fdfminimizer_type *min_type = gsl_multimin_fdfminimizer_conjugate_fr; //set this properly

	gsl_multimin_fdfminimizer *fdfmin = gsl_multimin_fdfminimizer_alloc(min_type, nthetas_opt);
	
	/* printf ("#multimin: using a '%s' minimizer with %d dimensions\n", */
	/* 				gsl_multimin_fdfminimizer_name (fdfmin), nthetas_opt); */

	status = gsl_multimin_fdfminimizer_set(fdfmin, &multiminFn, thetaTemp, stepSizeInit, tolerance);
	/* printf("set_status: %d\n", status); */

	// iterate our minimizer
	do{
		status = gsl_multimin_fdfminimizer_iterate(fdfmin);
		/* printf("status =%d\n", status); */
		
		if(status == GSL_ENOPROG && stepcount > 0){
			// no progress on this step
			//fprintf(stderr, "#mutimin: no progress. nsteps: %d\n", stepcount);
			break;
		}

		fnValue = gsl_multimin_fdfminimizer_minimum(fdfmin);
		tempTest = gsl_multimin_fdfminimizer_gradient(fdfmin);
		/* // output current info */
		/* fprintf(stderr, "#(%d) f: %lf x:", stepcount, fnValue); */
		/* for (i = 0; i < nthetas_opt; ++i) */
		/* 	fprintf(stderr, "%lf ", gsl_vector_get(fdfmin->x, i)); */
		/* fprintf(stderr, "#grad: "); */
		/* for (i = 0; i < nthetas_opt; ++i) */
		/* 	fprintf(stderr, "%lf ", gsl_vector_get(fdfmin->gradient, i)); */
		
		/* norm = 0.0; */
		/* for(i = 0; i < nthetas_opt; ++i){ */
		/* 	norm+= gsl_vector_get(fdfmin->gradient,i)*gsl_vector_get(fdfmin->gradient,i); */
		/* } */
		/* fprintf(stderr, "norm: %g\n", sqrt(norm)); */
		
		
		// this will return GSL_SUCCESS if we are done
		status = gsl_multimin_test_gradient(fdfmin->gradient, epsAbs);

		stepcount++;
	} while(status == GSL_CONTINUE && stepcount < stepmax);

	if(stepcount == stepmax)
		fprintf(stderr, "#multimin: no converge at stepmax %d\n", stepmax);
	

	fnValue = gsl_multimin_fdfminimizer_minimum(fdfmin);
	/* fprintf(stderr, "#multimin: best_value %lf\n", fnValue); */


	// get the min location
	tempTest = gsl_multimin_fdfminimizer_x(fdfmin);

	/* for (i = 0; i < nthetas_opt; ++i) */
	/* 	fprintf(stderr, "%lf ", gsl_vector_get(tempTest, i)); */
	/* fprintf(stderr, "\n"); */
	
	// now estimate sigma one last time
	sigma_final = log(estimateSigmaFull(tempTest, args));
	if(sigma_final == GSL_NAN){
		fprintf(stderr, "# cannot estimate sigma for these x\n");
		gsl_vector_free(thetaTemp);
		gsl_vector_free(tempTest);
		gsl_multimin_fdfminimizer_free(fdfmin);
		exit(1);
	}
	/* fprintf(stderr, "#FINAL estimated sigma: %g\n", sigma_final); */

	gsl_vector_set(thetaFinal, 0, sigma_final);
	for(i = 0; i < nthetas_opt; i++)
		gsl_vector_set(thetaFinal, i+1, gsl_vector_get(tempTest, i));
	

	// clear up
	gsl_multimin_fdfminimizer_free(fdfmin);
	gsl_vector_free(thetaTemp);
	//gsl_vector_free(tempTest);


}





/**
 * init the vector x to a set of random values which are sampled from a 
 * uniform dist (from rand) generated on the ranges given by the 
 * matrix ranges
 */
void set_random_init_value(gsl_rng* rand, gsl_vector* x, gsl_matrix* ranges,int  nthetas){
	int i;
	double range_min; 
	double range_max;
	double the_value;
	
	for(i = 0; i < nthetas; i++){
		range_min = gsl_matrix_get(ranges, i, 0);
		range_max = gsl_matrix_get(ranges, i, 1);
		// set the input vector to a random value in the range
		the_value = gsl_rng_uniform(rand) * (range_max - range_min) + range_min;
		/* printf("theta(%d) init-val:  %g\n", i, the_value); */
		//gsl_vector_set(x, gsl_rng_uniform(rand)*(range_max-range_min)+range_min, i);
		gsl_vector_set(x, i, the_value);
	}
}



