#include "estimator.h"

// these are the functions which we will put together to form our evalution function
// for a single step of the maximum likelyhood estimator, we will need to create 
// the inverse of the covariance matrix and then return the gradient vector or 
// the loglikelyhood depending on the method used to maximise it (the log likelyhood). 
//
// as such:
// these functions are VERY important since they entirely define the resulting shape 
// of your MLE estimate, if you want to change the emulation process you may need to change
// the getGradient function. if you use a gradient method it is clear that you WILL
// have to change these
// 



// i think i'm passing more arguments than i probably need to
// oh well
// the functions we'll need for the MLE estimation


// @param cinverse the inverse covariance matrix for this evaluation,
// calculated through LU decomp somewhere else, since we want to use it in getGradient too, without re-calcing it
// 
// @param det_cinverse the determinant, calculated at the same time as the inverse
double getLogLikelyhood(gsl_matrix *cinverse, double det_cinverse,  gsl_matrix *xmodel, gsl_vector *trainingvector, gsl_vector *thetas, int nmodel_points, int nthetas, int nparams){
	double the_likelyhood = 0.0;
	double vector_matrix_vector_product = 0.0;
	double log_2_pi = 0.798179868;
	gsl_vector *result_holder = gsl_vector_alloc(nmodel_points);
	// the  log likelyhood is a given by
	// L = (-1/2)*Log[Det[cinverse]]  - (1/2)*trainingvector.cinverse.trainingvector - (nmodel_points/2)*Log[2*Pi]
	the_likelyhood += -(1.0/2.0)*Log(det_cinverse) -  (nnmodel_points/2.0)*log_2_pi;
	
	gsl_blas_dgemv(CblasNoTrans, 1.0, inverse_cov_matrix, trainingvector, 0.0, result_holder);
	gsl_blas_ddot(trainingvector, result_holder, &vector_matrix_vector_product);
	
	the_likelyhood += vector_matrix_vector_product*(-1.0/2.0);

	gsl_vector_free(result_holder);
	return(the_likelyhood);
}

// get the gradient in a specific direction (given by index)
double getGradient(gsl_matrix *cinverse, gsl_matrix *xmodel, gsl_vector *trainingvector, gsl_vector *thetas, int index, int nmodel_points, int nthetas, int nparams){
	int i;
	gsl_vector *result_holder = gsl_vector_alloc(nmodel_points);
	gsl_matrix *dcdt = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix *cinverse_dcdt = gsl_matrix_alloc(nmodel_points, nmodel_points);
	double big_matrix_product = 0.0;
	double the_gradient = 0.0;
	double trace = 0.0
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
	gsl_blas_dgemv(CblasNoTrans, 1.0, dct, training_vector, 0.0, result_holder);
	gsl_blas_ddot(training_vector, result_holder,  &big_matrix_product);
	the_gradient += big_matrix_product*(1.0/2.0);
	
	gsl_vector_free(result_holder);
	gsl_vector_free(dcdt);
	gsl_vector_free(cinvese_dcdt);
	return(the_gradient);

}


double getdCdt(gsl_matrix* dcdt, gsl_matrix* xmodel, gsl_vector* thetas, int index, int nmodel_points, int nthetas, int nparams){
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
				gsl_matrix_set(dcdt, i, j, detCovFn_higher(&xmodel_row_i.vector, &xmodel_row_j.vector, thetats, thetas, nparams));
			}											 
		}
	}		
}



//! return the derivative of the covariance function with respect to the first constant theta (the amplitude) 
//@return the derivative as a double
//@param xm,xn -> vectors of model params to be compared
//@param thetas -> vector of the model hyper parameters
//@param nthetas -> length of thetas
//@param nparams -> number of input parameters (length of the xm,xn vectors)
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


//! return the derivative of the covariance function with respect to the higher thetas, these 
//! are the ones which set the correlation length etc, very important!
//@return the derivative as a double
//@param xm,xn -> vectors of model params to be compared
//@param thetas -> vector of the model hyper parameters
//@param index -> which of the higher thetas we are looking at. 
//@param nthetas -> length of thetas
//@param nparams -> number of input parameters (length of the xm,xn vectors)
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
	xm_temp = gsl_vector_get(xn, index-3);
	covariance *= (1/gsl_vector_get(thetas, index))*(gsl_vector_get(thetas, 0))*(xm-xn)*(xm-xn);
	return(covariance);
}
