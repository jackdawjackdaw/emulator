#include "estimator.h"

/** 
 * @file 
 * @author Chris Coleman-Smith cec24@phy.duke.edu
 * @version 1.0
 * 
 */ 



//! calc the log likelyhood for a given cinverse and hyperparameters
/**
 * calculate the log likelyhood for the system, given the inverse covariance matrix and it's determinant.
 * If you change the covariance function around a bit this should probably still work without too much 
 * trouble. You are advised to test this however
 *  
 * somehow log(det_cinverse) is > 0 and judging from the MM code it should be negative,
 * have hacked it so that it is, but this is abit weird
 * 
 * L = (-1/2)*Log[Det[cinverse]]  - (1/2)*(trainingvector-est_mean).cinverse.(trainingvector-est_mean) - (nmodel_points/2)*Log[2*Pi]
 * 
 * @return the log likelyhood for a given set of hyperparams theta (and cinverse) 
 * @param cinverse -> the inverse covariance matrix for this evaluation,
 * calculated through LU decomp somewhere else, since we want to use it in getGradient too, without re-calcing it
 * @param det_cinverse -> the determinant of the inverse covariance matrix, calculated at the same time as the inverse
 * @param xmodel -> matrix of model evaluation points 
 * @param trainingvector -> vector of model output values
 * @param thetas  -> vector of hyperparameters, these may not be needed here, but they are what we are actually hoping to optimise 
 * @param h_matrix -> matrix of the regression vector evaluated at all the design points (xmodel)
 * @param nregression_fns the number of regression fns used (length of hvector)
 */
double getLogLikelyhood(gsl_matrix *cinverse, double det_cmatrix,  gsl_matrix *xmodel, gsl_vector *trainingvector, gsl_vector *thetas, gsl_matrix *h_matrix, int nmodel_points, int nthetas, int nparams, int nregression_fns){
	int i;
	double the_likelyhood = 0.0;
	double vector_matrix_vector_product = 0.0;
	double log_2_pi = 1.83788;
	gsl_vector_view xmodel_row;
	gsl_vector *result_holder = gsl_vector_alloc(nmodel_points);
	gsl_vector *h_vector = gsl_vector_alloc(nregression_fns);
	gsl_vector *beta_vector = gsl_vector_alloc(nregression_fns);
	gsl_vector *estimated_mean = gsl_vector_alloc(nmodel_points);
	gsl_vector *train_sub_mean = gsl_vector_alloc(nmodel_points);
	double estimated_mean_val = 0.0;
	// but this should be determinant of C not inverse of C!
	double log_det_c = log(det_cmatrix);

	

	// not sure why we want the fabs here
	//log_det_c = fabs(log_det_c);

	/* we need to calculate the mean vector for this set of thetas 
	 * estMean[i] = hvector(training[i]).betavector
	 */
	estimateBeta(beta_vector, h_matrix, cinverse,  trainingvector, nmodel_points, nregression_fns);


	for(i = 0; i < nmodel_points; i++){
		xmodel_row = gsl_matrix_row(xmodel, i);
		makeHVector(h_vector, &xmodel_row.vector, nparams);
		//print_vector_quiet(h_vector, nregression_fns);
		gsl_blas_ddot(beta_vector, h_vector, &estimated_mean_val);
		gsl_vector_set(estimated_mean, i, estimated_mean_val);
	}

	/* train_sub_mean = trainingvector- estimated_mean */
	gsl_vector_memcpy(train_sub_mean, trainingvector);
	gsl_vector_sub(train_sub_mean, estimated_mean);

	//print_vector_quiet(beta_vector, nregression_fns); <- verified
	//print_vector_quiet(trainingvector, nmodel_points); <- verified 
	//print_vector_quiet(estimated_mean, nmodel_points); 
	/* print_vector_quiet(train_sub_mean, nmodel_points); */
	// DEBUG 

	// the  log likelyhood is a given by
	// L = (-1/2)*Log[Det[cinverse]]  - (1/2)*trainingvector.cinverse.trainingvector - (nmodel_points/2)*Log[2*Pi]
	the_likelyhood = -(1.0/2.0)*log_det_c -  (nmodel_points/2.0)*log_2_pi;

	gsl_blas_dgemv(CblasNoTrans, 1.0, cinverse, train_sub_mean, 0.0, result_holder);
	gsl_blas_ddot(train_sub_mean, result_holder, &vector_matrix_vector_product);
	
	//printf("%g\n", (log_2_pi)*(nmodel_points/2.0));
	/* printf("det_cmatrix = %g\n", det_cmatrix); */
	/* printf("parta:(-1/2)*log_det_c = %g\n", (-0.5)*log_det_c); */
	/* printf("partb:(training-mean).cinverse.(training-mean) = %g\n", (-0.5)*vector_matrix_vector_product); */
	/* printf("partc:(-nmodel_points/2.0)*log_2_pi = %g\n", -(nmodel_points/2.0)*log_2_pi); */

	the_likelyhood += vector_matrix_vector_product*(-1.0/2.0);

	gsl_vector_free(result_holder);
	gsl_vector_free(h_vector);
	gsl_vector_free(beta_vector);
	gsl_vector_free(estimated_mean);
	gsl_vector_free(train_sub_mean);
	return(the_likelyhood);
}



/* ***************************** DEPRECATED **************************************** 
 * 
 * There should be no calls for analytic gradients, this is too confusingly implemented 
 * right now to be useful.
 * 
 * use a getGradientNumeric for this , like in lbfgs.c 
 * 
 */
