#include "estimator.h"

/** 
 * @file 
 * @author Chris Coleman-Smith cec24@phy.duke.edu
 * @version 0.1
 * @section DESCRIPTION
 * 
 * these are the functions which we will put together to form our evalution function
 * for a single step of the maximum likelyhood estimator, we will need to create 
 * the inverse of the covariance matrix and then return the gradient vector or 
 * the loglikelyhood depending on the method used to maximise it (the log likelyhood). 
 *
 * as such:
 * these functions are VERY important since they entirely define the resulting shape 
 * of your MLE estimate, if you want to change the emulation process you may need to change
 * the getGradient function. if you use a gradient method it is clear that you WILL
 * have to change these
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
 * L = (-1/2)*Log[Det[cinverse]]  - (1/2)*trainingvector.cinverse.trainingvector - (nmodel_points/2)*Log[2*Pi]
 * 
 * @return the log likelyhood for a given set of hyperparams theta (and cinverse) 
 * @param cinverse -> the inverse covariance matrix for this evaluation,
 * calculated through LU decomp somewhere else, since we want to use it in getGradient too, without re-calcing it
 * @param det_cinverse -> the determinant of the inverse covariance matrix, calculated at the same time as the inverse
 * @param xmodel -> matrix of model evaluation points 
 * @param trainingvector -> vector of model output values
 * @param thetas  -> vector of hyperparameters, these may not be needed here, but they are what we are actually hoping to optimise 
 */
double getLogLikelyhood(gsl_matrix *cinverse, double det_cinverse,  gsl_matrix *xmodel, gsl_vector *trainingvector, gsl_vector *thetas, int nmodel_points, int nthetas, int nparams){
	double the_likelyhood = 0.0;
	double vector_matrix_vector_product = 0.0;
	double log_2_pi = 1.83788;
	gsl_vector *result_holder = gsl_vector_alloc(nmodel_points);
	double log_det_c = log(det_cinverse);

	if(log_det_c > 0.0){
		// force postive
		log_det_c = log_det_c *(-1);
	}

	//printf("%g\n", log_det_c);
	// the  log likelyhood is a given by
	// L = (-1/2)*Log[Det[cinverse]]  - (1/2)*trainingvector.cinverse.trainingvector - (nmodel_points/2)*Log[2*Pi]
	the_likelyhood += -(1.0/2.0)*log_det_c -  (nmodel_points/2.0)*log_2_pi;
	
	gsl_blas_dgemv(CblasNoTrans, 1.0, cinverse, trainingvector, 0.0, result_holder);
	gsl_blas_ddot(trainingvector, result_holder, &vector_matrix_vector_product);
	
	//printf("%g\n", (log_2_pi)*(nmodel_points/2.0));

	the_likelyhood += vector_matrix_vector_product*(-1.0/2.0);

	gsl_vector_free(result_holder);
	return(the_likelyhood);
}

//! get the gradient of the log likelyhood in a specific direction (given by index)
/** 
 * calculate the gradient of the log likleyhood in a given direction, 
 * this is the analytical gradient, worked out for the specific form of the covariance function 
 * (here t -> thetas)
 * C = t1 * exp(-(x1-x2)^2/t4) + t2 + delta(x1,s2)*t3
 * so if you change the covariance function s.t the derivatives wrt to the various tvalues changes
 * you will have to adjust this function. 
 * 
 * the gradient is: (-1/2)*Tr(cinverse.dCdT) + (1/2)*trainingvector.cinverse.dCdT.cinverse.trainingvector
 *
 *
 * @return the derivative of the log likelyhood, evaluated at the hyperparams theta, in the direction given by index
 * @param cinverse -> the inverted covariance matrix
 * @param xmodel -> the model evaluation points
 * @param trainingvector -> vector of the model output
 * @param index -> which direction to calculate the gradient in, i.e wrt to which theta do we differentiate the covariance_fn
 */
double getGradient(gsl_matrix *cinverse, gsl_matrix *xmodel, gsl_vector *trainingvector, gsl_vector *thetas, int index, int nmodel_points, int nthetas, int nparams){
	int i;
	gsl_vector *result_holder = gsl_vector_alloc(nmodel_points);
	gsl_matrix *dcdt = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix *cinverse_dcdt = gsl_matrix_alloc(nmodel_points, nmodel_points);
	double big_matrix_product = 0.0;
	double the_gradient = 0.0;
	double trace = 0.0;
	getdCdt(dcdt, xmodel, thetas, index, nmodel_points, nthetas, nparams);


	// the gradient is: (-1/2)*Tr(cinverse.dCdT) + (1/2)*trainingvector.cinverse.dCdT.cinverse.trainingvector
	
	// do the trace part,
	// first we calc cinverse_dcdt
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, cinverse, dcdt, 0.0, cinverse_dcdt);
	// now get the trace
	for(i = 0; i < nmodel_points; i++){
		trace += gsl_matrix_get(cinverse_dcdt, i, i);
	}
	the_gradient += trace*(-1.0/2.0);

	// now do the other bit, we can recyle some of the matricies, so this will probably look a bit confusing
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, cinverse_dcdt, cinverse, 0.0, dcdt); // now dcdt holds cinverse.dCdt.cinverse
	gsl_blas_dgemv(CblasNoTrans, 1.0, dcdt, trainingvector, 0.0, result_holder);
	gsl_blas_ddot(trainingvector, result_holder,  &big_matrix_product);
	the_gradient += big_matrix_product*(1.0/2.0);
	
	gsl_vector_free(result_holder);
	gsl_matrix_free(dcdt);
	gsl_matrix_free(cinverse_dcdt);
	return(the_gradient);

}


//! get differentiated covariance matrix for a direction given by index
/** 
 * calculate the covariance matrix in the case where we have differentiated the covariance function by 
 * the theta parameter labelled by index. 
 * This is incredibly covariance_fn dependant! You must look at this and adjust it if you adjust the covariance_fn!
 * you will also need to adjust detCovFn_1 and detCovFn_higher also
 * 
 * of course, if you use a MLE method which does not need the gradient you can safely ignore getGradient and this function
 * 
 * @param dcdt -> this is set to the result of diff(C, index)
 * @param xmodel -> the model evaluation points
 * @param thetas -> the hyperparameters, important here as they will determine the final values of dcdt
 * @return dcdt is set to the given derivative 
 */
void getdCdt(gsl_matrix* dcdt, gsl_matrix* xmodel, gsl_vector* thetas, int index, int nmodel_points, int nthetas, int nparams){
	assert(index <= nthetas-1);
	// declare so that they can only be called from this scope
	double detCovFn_1(gsl_vector *xm, gsl_vector *xn, gsl_vector* thetas, int nthetas, int nparams);
	double detCovFn_higher(gsl_vector *xm, gsl_vector *xn, gsl_vector* thetas, int index, int nthetas, int nparams);
	// this is very covariance specific, the other functions here arn't so this is where you need to pay attention if 
	// you are changing things, the idea is the return the covariance matrix where the covariance function 
	// has been differentiated by theta[index]
	// 
	// if this becomes unpossible, approximating the gradient through finite-diffs might be ok
	// but i don't know what happens to the errors there, probably nothing healthy
	int i, j;
	gsl_vector_view xmodel_row_i;
	gsl_vector_view xmodel_row_j;
	for(i = 0; i < nmodel_points; i++){
		for(j=0; j < nmodel_points; j++){
			xmodel_row_i = gsl_matrix_row(xmodel, i);
			xmodel_row_j = gsl_matrix_row(xmodel, j);				
			if(index == 0){ // the amplitude
				gsl_matrix_set(dcdt, i, j, detCovFn_1(&xmodel_row_i.vector, &xmodel_row_j.vector, thetas, nthetas, nparams));
			} else if (index == 1){ // the additional constant
				gsl_matrix_set(dcdt, i, j, 1.0);
			} else if (index == 2){ // the kernel term
				if(i == j){
					gsl_matrix_set(dcdt, i, j, 1.0);
				} else {
					gsl_matrix_set(dcdt, i, j, 0.0);
				}
			} else if( index >= 3){
				gsl_matrix_set(dcdt, i, j, detCovFn_higher(&xmodel_row_i.vector, &xmodel_row_j.vector, thetas, index, nthetas, nparams));
			}											 
		}
	}		
}



//! return the derivative of the covariance function with respect to the first constant theta (the amplitude) 
/**
 * 
 * return diff(covariance_fn, theta0), the derivative of the covaraince fn wrt the amplitude
 * 
 * @return the derivative as a double
 * @param xm,xn -> vectors of model params to be compared
 * @param thetas -> vector of the model hyper parameters
 * @param nthetas -> length of thetas
 * @param nparams -> number of input parameters (length of the xm,xn vectors)
 */
double detCovFn_1(gsl_vector *xm, gsl_vector *xn, gsl_vector* thetas, int nthetas, int nparams){
	int i;
	double covariance = 0.0;
	double xm_temp = 0.0;
	double xn_temp = 0.0;
	double r_temp = 0.0;
	for (i = 0; i < nparams; i++){
		xm_temp = gsl_vector_get(xm, i);
		xn_temp = gsl_vector_get(xn, i);
		r_temp = gsl_vector_get(thetas, i+3);
		r_temp = r_temp*r_temp;
		covariance += exp((-1.0/2.0)*((xm_temp-xn_temp)*(xm_temp-xn_temp))/(r_temp));
	}
	// this is the derivative, don't actually use any of the basic thetas
	return(covariance);
}


//! return the derivative of the covariance function with respect to the higher thetas
/**
 * return the deriv of covariance_fn wrt to the higher order thetas, these 
 *  are the ones which set the correlation length etc, very important!
 * @return the derivative as a double
 * @param xm,xn -> vectors of model params to be compared
 * @param thetas -> vector of the model hyper parameters
 * @param index -> which of the higher thetas we are looking at. 
 * @param nthetas -> length of thetas
 * @param nparams -> number of input parameters (length of the xm,xn vectors)
 */
double detCovFn_higher(gsl_vector *xm, gsl_vector *xn, gsl_vector* thetas, int index, int nthetas, int nparams){
	int i;
	double covariance = 0.0;
	double xm_temp = 0.0;
	double xn_temp = 0.0;
	double r_temp = 0.0;
	// this part can be the same as before
	for (i = 0; i < nparams; i++){
		xm_temp = gsl_vector_get(xm, i);
		xn_temp = gsl_vector_get(xn, i);
		r_temp = gsl_vector_get(thetas, i+3);
		r_temp = r_temp*r_temp;
		covariance += exp((-1.0/2.0)*((xm_temp-xn_temp)*(xm_temp-xn_temp))/(r_temp));
	}
	
	// get the right things from the model parameters
	// i.e if we have a 2d x, we will have a theta -> (0,1,2,3,4)
	// where thetas 4 and 5 are the length scale indicies
	// thefore the fn will be called with index 4 or 3
	xm_temp = gsl_vector_get(xm, index-3);
	xn_temp = gsl_vector_get(xn, index-3);
	covariance *= pow(gsl_vector_get(thetas, index), -3.0)*(gsl_vector_get(thetas, 0))*(xm_temp-xn_temp)*(xm_temp-xn_temp);
	return(covariance);
}
