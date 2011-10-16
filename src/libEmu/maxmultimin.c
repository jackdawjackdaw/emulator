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
 * currently this roughly works, but the routine seems to terminate very quickly, well
 * before the gradient is satisfactory
 * 
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
 */
void maxWithMultiMin(struct estimate_thetas_params *params){
	int tries = 0;
	// make local bindings of commonly used values in params
	int nthetas = params->options->nthetas;
	int nmodel_points = params->options->nmodel_points;
	int nparams = params->options->nparams;
	int nregression_fns = params->options->nregression_fns;

	double likelihood = 0.0;
	double bestLHood = SCREWUPVALUE;

	gsl_vector *xInit = gsl_vector_alloc(nthetas);
	gsl_vector *xFinal = gsl_vector_alloc(nthetas);
	gsl_vector *xBest = gsl_vector_alloc(nthetas);
	double *tempVec = malloc(sizeof(double)*nthetas);
	pthread_t self = pthread_self(); // this is your thread id

	// init our vectors
	gsl_vector_set_zero(xBest);
	gsl_vector_set_zero(xFinal);
	set_random_init_value(params->random_number, xInit, params->options->grad_ranges, nthetas);

	// generate the regression details for our params object
	params->h_matrix = gsl_matrix_alloc(nmodel_points, nregression_fns);
	makeHMatrix(params->h_matrix, params->the_model->xmodel,nmodel_points, nparams, nregression_fns);

	/* printf("max_tries = %d\n", params->max_tries); */
	while(tries < params->max_tries) {
		// do the actual optimization using MultiMin
		
		doOptimizeMultiMin(&evalFnMulti, // computes the log-likelihood 
											 &gradFnMulti,  // computes the gradient
											 &evalFnGradMulti, // does both at once
											 xInit, xFinal, (void*)params);
		
		likelihood = -1*evalFnMulti(xFinal, (void*)params);
		
		/*annoying!
		**/
		fprintPt(stdout, self);
		printf(":L = %g:try = %d\n", likelihood, tries);
		 
		if(likelihood > bestLHood && (isnan(likelihood) == 0 && isinf(likelihood) == 0)){
			bestLHood = likelihood;
			gsl_vector_memcpy(xBest, xFinal);
			
			/* fprintPt(stdout, self); */
			/* printf(":best = %g\n", bestLikelyHood); */
		}
		tries++;
		set_random_init_value(params->random_number, xInit, params->options->grad_ranges, nthetas);
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
	free(tempVec);

}


/*
 * the actual eval fn, returns the loglikelihood for a given set of vectors at the location theta_vec
 * for the model specified by the params struct
 */
double evalFnMulti(const gsl_vector* theta_vec, void* params_in){
	struct estimate_thetas_params *params = (struct estimate_thetas_params*) params_in;
	
	int nmodel_points = params->options->nmodel_points;
	int nthetas = params->options->nthetas;
	int nparams = params->options->nparams;
	
	// need to disable the gsl error handler
	gsl_error_handler_t *temp_handler;
	gsl_matrix* covariance_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix* cinverse = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix* temp_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points);

	// make a copy of thetas
	gsl_vector *theta_local = gsl_vector_alloc(nthetas);
	gsl_vector_memcpy(theta_local, theta_vec);
	
	double determinant_c = 0.0;
	double temp_val = 0.0;
	int i, j;
	int cholesky_test = 0;
		
	// copy the given double vector into the gsl vector
	//copy_vec_gslvec(xinput, xk, params->nthetas );

	// make the covariance matrix 
	// using the random initial conditions! (xold not thetas)
	makeCovMatrix(covariance_matrix, params->the_model->xmodel, theta_local, nmodel_points, nthetas, nparams);
	gsl_matrix_memcpy(temp_matrix, covariance_matrix);
	
	//print_matrix(temp_matrix, params->options->nmodel_points, params->options->nmodel_points);

	// do a cholesky decomp of the cov matrix, LU is not stable for ill conditioned matrices
	temp_handler = gsl_set_error_handler_off();
	cholesky_test = gsl_linalg_cholesky_decomp(temp_matrix);

	if(cholesky_test == GSL_EDOM){
		FILE *fptr;
		fprintf(stderr, "evalFnMulti\n");
		fprintf(stderr, "trying to cholesky a non postive def matrix, sorry...\n");
		fprintf(stderr, "matrix dumped to chol-err.dat\n");
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
		exit(1);
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

	// temp_val is now the likelyhood for this answer
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
 * theta_0 (the scale of the cov fn)
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
 */
void gradFnMulti(const gsl_vector* theta_vec, void* params_in, gsl_vector * grad_vec){
	int nmpoints, nthetas, i, j;
	int nparams;
	double amp, nug; 
	double grad_temp = 0.0;

	struct estimate_thetas_params *params = (struct estimate_thetas_params*)params_in;
	nmpoints = params->options->nmodel_points;
	nthetas = params->options->nthetas;
	nparams = params->options->nparams;
	int cholesky_test;
	gsl_error_handler_t *temp_handler;
	gsl_vector* theta_local = gsl_vector_alloc(nthetas);
	gsl_matrix* covariance_matrix = gsl_matrix_alloc(nmpoints, nmpoints);
	gsl_matrix* cinverse = gsl_matrix_alloc(nmpoints, nmpoints);
	gsl_matrix* temp_matrix = gsl_matrix_alloc(nmpoints, nmpoints);
	
	gsl_vector_memcpy(theta_local, theta_vec);

	// using the random initial conditions! (xold not thetas)
	makeCovMatrix(covariance_matrix, params->the_model->xmodel, theta_local, nmpoints, nthetas, nparams);
	gsl_matrix_memcpy(temp_matrix, covariance_matrix);

	#ifdef DEBUG1
	printf("#thetas: ");
	for(i = 0; i < nthetas; ++i)
		printf("%g ", gsl_vector_get(theta_local, i));
	printf("\n");
	#endif
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

		exit(1);
	}
	gsl_set_error_handler(temp_handler);

	gsl_linalg_cholesky_invert(temp_matrix);
	gsl_matrix_memcpy(cinverse, temp_matrix);

	// C(x,y) = theta_0 c(x,y) + theta_1
	// make the first component, here dC/dN_0 = 1/theta0 * (C  - theta_1)
	// the second component, dC/dN_1 = 1_delta(ij) (easssy)
	// the rest: dC/dN_{h} = amp * dC/d\theta
	//
	// but since we log scaled theta_0 and theta_1 we need to do 
	// dC/dN_0 = 1/exp(theta0) * c * d(exp(theta_0), theta_0) = (C - theta_1)
	
	// first, get the scale parameter and the nugget
	amp = exp(gsl_vector_get(theta_local, 0)); 
	nug = exp(gsl_vector_get(theta_local, 1)); 
	gsl_matrix_memcpy(temp_matrix, covariance_matrix);
	// now we need to sub out the contribution from the nugget, in the bad way...
	for(i = 0; i < nmpoints; i++)
		gsl_matrix_set(temp_matrix, i, i, gsl_matrix_get(temp_matrix,i,i) - nug);

	/* printf("#dC/dtheta_0\n"); */
	/* print_matrix(temp_matrix, nmpoints, nmpoints); */

	grad_temp = -1.0*getGradientCn(temp_matrix, cinverse, params->the_model->training_vector, nmpoints,  nthetas);
	gsl_vector_set(grad_vec, 0, grad_temp);

	// next, the nugget, if we don't log scale we set this to the identity
	gsl_matrix_set_identity(temp_matrix);
	// since we're log scaling the nugget we need to set the cov-matrix accordingly
	gsl_matrix_scale(temp_matrix, nug);

 	/* printf("#dC/dtheta_1\n"); */
	/* print_matrix(temp_matrix, nmpoints, nmpoints); */

	grad_temp= -1.0*getGradientCn(temp_matrix, cinverse, params->the_model->training_vector, nmpoints,  nthetas);
	gsl_vector_set(grad_vec, 1, grad_temp);
	
	for(i = 2; i < nthetas; i++){
		makeGradMatLength(temp_matrix, params->the_model->xmodel, gsl_vector_get(theta_local, i) , i, nmpoints, nparams);
		gsl_matrix_scale(temp_matrix, amp);
		/* printf("#dC/dtheta_2\n"); */
		/* print_matrix(temp_matrix, nmpoints, nmpoints); */
		grad_temp = -1.0*getGradientCn(temp_matrix, cinverse, params->the_model->training_vector, nmpoints,  nthetas);
		gsl_vector_set(grad_vec, i, grad_temp);
	}

	#ifdef DEBUG1
	printf("# gradient: ");
	for(i = 0; i < nthetas; ++i)
		printf("%g ", gsl_vector_get(grad_vec,i));
	printf("\n");
	#endif

	
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


// compute the function value and the gradient vector in one call
void evalFnGradMulti(const gsl_vector* theta_vec, void* params, double* fnval, gsl_vector * grad_vec){
	printf("# called fdf\n");
	*fnval = evalFnMulti(theta_vec, params);
	gradFnMulti(theta_vec, params, grad_vec);
}


// compute the loglikelihood just using params
double evalLikelihood_struct(struct estimate_thetas_params *params){
	int i;
	double likelihood = 0.0;
	int nthetas = params->options->nthetas;
	gsl_vector *thetas_local = gsl_vector_alloc(nthetas);

	gsl_vector_memcpy(thetas_local, params->the_model->thetas);

	printf("struct thetas: ");
	for(i=0; i< nthetas;i++)
		printf("%g ", gsl_vector_get(thetas_local, i));
	printf("\n");

	likelihood = -1*evalFnMulti(thetas_local, params);
	printf("lhood: %g\n", likelihood);

	gsl_vector_free(thetas_local);
	return(likelihood);
}



// actually carry out the multimin steps
void doOptimizeMultiMin( double(*fn)(const gsl_vector*, void*),													\
												void(*gradientFn)(const gsl_vector*,void*, gsl_vector*),
												void(*fnGradFn)(const gsl_vector*, void*, double*, gsl_vector*),
												gsl_vector *thetaInit, gsl_vector* thetaFinal, void* args){
	
	struct estimate_thetas_params *params = (struct estimate_thetas_params*)args;
	
	int nparams = params->options->nparams;
	int nmpoints = params->options->nmodel_points;
	int nthetas = params->options->nthetas;
	int status;
	int stepcount = 0;
	int stepmax = 30;
	int i;

	double stepSizeInit = 0.1;
	double tolerance = 0.2; // sets the accuracy of the line search.
	double fnValue = 0.0;
	double norm = 0.0; // norm of the gradient
	/*
	 * we stop when |g| < epsAbs, not sure how to set this yet
	 */
	double epsAbs = 0.05; 

	fprintf(stderr, "#doOptimizeMultiMin: nthetas %d\tnmpoints %d\tnparams %d\n", 
					nthetas, nmpoints, nparams);
	
	gsl_vector *tempTest = gsl_vector_alloc(nthetas); // get values out of the min to check on progress
	
	gsl_vector *thetaTemp = gsl_vector_alloc(nthetas);
	gsl_vector_memcpy(thetaTemp, thetaInit);
	//gsl_vector_set_zero(thetaTemp);

	// our multimin fn
	gsl_multimin_function_fdf multiminFn;
	multiminFn.n = nthetas; // the number of dimensions
	/* multiminFn.f = fn; */
	/* multiminFn.df = gradientFn; */
	/* multiminFn.fdf = fnGradFn; */
	multiminFn.f = &evalFnMulti;
	multiminFn.df = &gradFnMulti;
	multiminFn.fdf = &evalFnGradMulti;
	
	multiminFn.params = args;

	// the multiminimizer, this is not working very well currently...
	const gsl_multimin_fdfminimizer_type *min_type = gsl_multimin_fdfminimizer_vector_bfgs2; //set this properly
	//const gsl_multimin_fdfminimizer_type *min_type = gsl_multimin_fdfminimizer_conjugate_fr; //set this properly

	gsl_multimin_fdfminimizer *fdfmin = gsl_multimin_fdfminimizer_alloc(min_type, nthetas);
	
	
	
	printf ("#multimin: using a '%s' minimizer with %d dimensions\n",
					gsl_multimin_fdfminimizer_name (fdfmin), nthetas);

	status = gsl_multimin_fdfminimizer_set(fdfmin, &multiminFn, thetaTemp, stepSizeInit, tolerance);
	printf("set_status: %d\n", status);

	// iterate our minimizer
	do{
		status = gsl_multimin_fdfminimizer_iterate(fdfmin);
		printf("status =%d\n", status);
		
		if(status == GSL_ENOPROG && stepcount > 0){
			// no progress on this step
			fprintf(stderr, "#mutimin: no progress. nsteps: %d\n", stepcount);
			break;
		}

		fnValue = gsl_multimin_fdfminimizer_minimum(fdfmin);
		tempTest = gsl_multimin_fdfminimizer_gradient(fdfmin);
		// output current info
		fprintf(stderr, "#(%d) f: %lf x:", stepcount, fnValue);
		for (i = 0; i < nthetas; ++i)
			fprintf(stderr, "%lf ", gsl_vector_get(fdfmin->x, i));
		fprintf(stderr, "#grad: ");
		for (i = 0; i < nthetas; ++i)
			fprintf(stderr, "%lf ", gsl_vector_get(fdfmin->gradient, i));
		
		norm = 0.0;
		for(i = 0; i < nthetas; ++i){
			norm+= gsl_vector_get(fdfmin->gradient,i)*gsl_vector_get(fdfmin->gradient,i);
		}
		fprintf(stderr, "norm: %g\n", sqrt(norm));
		
		
		
		// this will return GSL_SUCCESS if we are done
		status = gsl_multimin_test_gradient(fdfmin->gradient, epsAbs);

		stepcount++;
	} while(status == GSL_CONTINUE && stepcount < stepmax);

	if(stepcount == stepmax)
		fprintf(stderr, "#multimin: no converge at stepmax %d\n", stepmax);
		

	fnValue = gsl_multimin_fdfminimizer_minimum(fdfmin);
	fprintf(stderr, "#multimin: best_value %lf\n", fnValue);

	for (i = 0; i < nthetas; ++i)
		fprintf(stderr, "%lf ", gsl_vector_get(fdfmin->x, i));
	fprintf(stderr, "\n\n");
	
	gsl_vector_memcpy(thetaFinal, fdfmin->x);

	

	// clear up
	gsl_vector_free(thetaTemp);
	gsl_multimin_fdfminimizer_free(fdfmin);

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
		printf("theta(%d) init-val:  %g\n", i, the_value);
		//gsl_vector_set(x, gsl_rng_uniform(rand)*(range_max-range_min)+range_min, i);
		gsl_vector_set(x, i, the_value);
	}
}



