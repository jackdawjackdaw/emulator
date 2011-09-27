#include "emulator.h"
/** 
 * @file 
 * @author Chris Coleman-Smith cec24@phy.duke.edu
 * @version 1.1
 * @section DESCRIPTION
 * 
 * This file contains the main functions needed to run the emulator, given a set of 
 * correct hyperparameters (theta) then the makeEmulatedMean and makeEmulatedVariance
 * functions can be used  
 *
 * the following convention for the hyperparameters (thetas) must be observed for 
 * any new covariance function added to the code-base
 * C(x,y) = theta_0 * c(x,y, theta_2,theta_3,...) + theta_1 
 * i.e each cov fn has a scale given by the first hyper param and a nugget term given by the 
 * second. 
 * if you don't do this then the gradient cannot simply be computed by dropping in a new dCn_dL matrix
 * as created by fns derivative_l_<name> for the different directions, this will result in woe.
 * 
 * sep-1 2011
 * updated cov fns to remove spurious double argument and allow for selection of 
 * covaraince fn from optstruct
 * 
 * sep-12 2011
 * 1)each cov fn is now joined by a derivative_l_<name> fn which takes a matix 
 * and returns a matrix of the derivative of the cov fn w.r.t one of the length scales
 * evaluated at all the model_points
 * these MUST BE UPDATED when you change the cov-fns
 * 2) remove exponential theta scaling within cov-fns, this is fucking confusing
 */

//! print the given matrix to the stdout
/**
 * prints the matrix to the stdout, assumes elements are doubles
 * @param nx -> width
 * @param ny -> height
 */
void print_matrix(gsl_matrix* m, int nx, int ny){
	int i,j;
	for(i = 0; i < nx; i++){
		for(j= 0; j < ny; j++){
			//printf("%g ", gsl_matrix_get(m, i, j));
			fprintf(stderr, "%g ", gsl_matrix_get(m,i,j));
		}
		//printf("\n");
		fprintf(stderr, "\n");
 	}
}


//! calculate the covariance between a set of input points
/** 
 * calculate the covariance between the two given vectors (xm, xn) 
 * where xm and xn are vectors of length nparams and their elements represent the various
 * input parameters used to evaluate the model at this point. 
 * The usual way to get xm and xn is to take a row slice from the model_points matrix representing a 
 * single evaluation of the model with the parameters. 
 * 
 *
 * the hyperparameters are used in this function, 
 * theta0 -> amplitude of the covariance gaussian
 * theta1 -> the nugget (only for diagonal terms) 
 * theta2 (and higher) -> scale parameter for the gaussian
 *
 * @param thetas -> hyperparameters 
 * @param nparams -> the legnth of xm and xn
 * @return the covariance
 * @param thetas -> vector of the hyperparameters of the process. 
 * @param nthetas -> length of thetas
 *
 * Updated 21 July 09
 * Changed the covariance function from being a pure gaussian to exp(-(x-y)^alpha)  
 * where alpha should be something on the range [1..2] and is set by the above #define
 * In lit (o'hagan) it's suggested that one might actually optimise for alpha also.
 * 
 * WELL i did, but then it sucked...
 * Also removed what was theta2, the offset term. This is suspect (wolpert)
 * 
 * Using the neldermead optimisation, this function is all you need to change
 * to use the grad-desc type methods which rely on the derivative of C we need 
 * to adjust some of the calls there. 
 * 
 * updated: march 23rd 11
 * optimisations by fixing alpha = 2.0 allowing us to remove calls to pow 
 * exp calls are reduced by clamping the range over which we actually compute the
 * exponential from 0..CLAMPVALUE, the rest are taken to be zero.
 * 
 * updated: sep 12 11,
 * 1) why are the thetas exponentiated? from now on all thetas 
 * are not to be scaled by cov fns
 * 
 *
 */

/** CLAMPVALUE is the cutoff for calculating the exponential covariance, anythign less than this
 * is just a waste of precision and computation
 * 
 * -10 is too large
 * -25 is about right, but you can fiddle around with this a bit to pull a tiny smidgin more speed out of the code
 * -100 is too small
 */
#define CLAMPVALUE -15
double covariance_fn_gaussian(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams){
	// calc the covariance for a given set of input points
	int i, truecount  = 0;
	double covariance = 0.0;
	double exponent = 0.0;
	double xm_temp = 0.0;
	double xn_temp = 0.0;
	double r_temp = 0.0;
	double dist_temp = 0.0;
	//double amp = exp(gsl_vector_get(thetas, 0));
	double amp = gsl_vector_get(thetas, 0);
	double nug = gsl_vector_get(thetas, 1);
	
	for(i = 0; i < nparams; i++){
		xm_temp = gsl_vector_get(xm, i);  // get the elements from the gsl vector, just makes things a little clearer
		xn_temp = gsl_vector_get(xn, i);

		//r_temp = gsl_vector_get(thetas, i+2);
    r_temp = exp(gsl_vector_get(thetas, i+2));

		// fix alpha = 2.0
		//r_temp = pow(r_temp , alpha); 
		r_temp = r_temp * r_temp;
		// gaussian term				

		// change from pow to explicit multiplication, for alpha = 2.0
		//exponent += (-1.0/2.0)*pow(fabs(xm_temp-xn_temp), alpha)/(r_temp);
		dist_temp = fabs(xm_temp-xn_temp);
		exponent += (-1.0/2.0)*dist_temp*dist_temp/(r_temp);

		//DEBUGprintf("%g\n", covariance);
		if (dist_temp < 0.0000000001){
			truecount++; 		
		}
	}
	
	// we're going to clamp exp(exponent) == 0 if exponent < -100
	if(exponent < CLAMPVALUE){
		covariance = 0.0;
	} else {
		covariance = exp(exponent)*amp;
	}

	/** 
	 * the nugget is only added to the diagonal covariance terms,
	 * the parts where xm == xn (vectorwise)
	 */
	if(truecount == nparams) {
		// i.e the two vectors are hopefully the same
		covariance += gsl_vector_get(thetas,1);
	}


	return(covariance);
}

/**
 * same as covariance_fn_gaussian but without any clamping of the exponential
 * 
 * we can use the same gradient here
 */
double covariance_fn_gaussian_exact(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams){
	// calc the covariance for a given set of input points
	int i, truecount  = 0;
	double covariance = 0.0;
	double exponent = 0.0;
	double xm_temp = 0.0;
	double xn_temp = 0.0;
	double r_temp = 0.0;
	double dist_temp = 0.0;
	//double amp = exp(gsl_vector_get(thetas, 0));
	double amp = gsl_vector_get(thetas, 0);
	double nug = gsl_vector_get(thetas, 1);
	
	for(i = 0; i < nparams; i++){
		xm_temp = gsl_vector_get(xm, i);  // get the elements from the gsl vector, just makes things a little clearer
		xn_temp = gsl_vector_get(xn, i);

		//r_temp = gsl_vector_get(thetas, i+2);
    r_temp = exp(gsl_vector_get(thetas, i+2));

		// fix alpha = 2.0
		//r_temp = pow(r_temp , alpha); 
		r_temp = r_temp * r_temp;
		// gaussian term				

		// change from pow to explicit multiplication, for alpha = 2.0
		//exponent += (-1.0/2.0)*pow(fabs(xm_temp-xn_temp), alpha)/(r_temp);
		dist_temp = fabs(xm_temp-xn_temp);
		exponent += (-1.0/2.0)*dist_temp*dist_temp/(r_temp);

		//DEBUGprintf("%g\n", covariance);
		if (dist_temp < 0.0000000001){
			truecount++; 		
		}
	}

	covariance = exp(exponent)*amp;

	/** 
	 * the nugget is only added to the diagonal covariance terms,
	 * the parts where xm == xn (vectorwise)
	 */
	if(truecount == nparams) {
		// i.e the two vectors are hopefully the same
		covariance += gsl_vector_get(thetas,1);
	}


	return(covariance);
}



/* 
 * compute a matrix of dC/dtheta-length / theta_0 for 
 * the gaussian covariance function c(x,y) = theta_0 * exp( - 1/2 * (x-y)^2 / theta_L^2 ) + theta_1
 * where theta_L = exp(theta_2) is the scale along one of the nparams dimensions 
 *  
 * if r^2 = (x-y)^2 then
 *
 * dC/dtheta_l / theta_0  = \frac{r^2 } {\theta_L^3}  exp(-1/2 * r^2 / theta_L^2 ) 
 * 
 * but we really want dC/dtheta
 * 
 * dC/dtheta / exp(theta_0) = exp(-0.5 * exp(-2*theta_L) * r^2  - 2*theta_L) * r^2
 * 
 * if alpha is ever modified above, this should also be changed
 */
void derivative_l_gauss(gsl_matrix *dCdTheta, gsl_matrix* xmodel, 
												double thetaLength, int index, int nmodel_points, int nparams){ 
	int i, j;
	double scale;
	double rtemp;
	const int nthetasConstant = 2;
	int indexScaled = index - nthetasConstant;
	double thetaLCubed;
	double partialCov = 0.0;
	double expTheta = exp(-2*thetaLength);

	//thetaLCubed = thetaLength * thetaLength * thetaLength;

	for(i = 0; i < nmodel_points; i++){
		for(j = 0; j < nmodel_points; j++){

			// is this the correct indexing into xmodel?
			rtemp = gsl_matrix_get(xmodel, i, indexScaled) - gsl_matrix_get(xmodel, j, indexScaled); 
			/*
			 * this is the correct form for the un logged thetas 
			 * dC/dtheta_l
			 */
			
			/* scale = (rtemp * rtemp) / (thetaLCubed); */
			/* partialCov = exp(-(0.5*rtemp*rtemp)/ thetaLength*thetaLength); */
			/* gsl_matrix_set(dCdTheta, i, j, scale*partialCov); */

			/*
			 * if we have defined our thetas on a log scale then
			 * we need to use the dC/dtheta expression
			 */
			partialCov = exp(-0.5*expTheta *rtemp*rtemp - 2*thetaLength)*rtemp*rtemp;
			gsl_matrix_set(dCdTheta, i, j, partialCov);
		}
	}
	
}




//! the matern cov fn but with nu set to 3/2
/**
 * uses 3 thetas only, 0 and 1 are for the vert scale and the nugget and 
 * 3 is the actual length scale 
 * 
 * this is a radially symmetric covariance fn so we only care about the distance
 * r between point pairs
 */
double covariance_fn_matern_three(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams){
	double covariance = 0.0;
	int i, truecount = 0;
	double xm_temp = 0.0;
	double xn_temp = 0.0;
	double distance = 0.0;
	double temp_dist = 0.0;
	
	assert(nthetas >= 3); // can throw away the upper ones no problem

	// map the thetas onto some local variables so the formula is more transparent
	double amp = gsl_vector_get(thetas, 0);
	double nugget = gsl_vector_get(thetas,1); 
	double rho = exp(gsl_vector_get(thetas, 2));

	double tempx = 0.0;
	double root3 = 1.732050808;

	// calculate the euclidean distance between the two points;
	for(i = 0; i < nparams; i++){
		xm_temp = gsl_vector_get(xm, i);
		xn_temp = gsl_vector_get(xn, i);
		// this is currently the distance squared
		temp_dist = fabs(xm_temp-xn_temp);
		distance += temp_dist * temp_dist;
		if(temp_dist < 0.0000000000000001){
			truecount++;
		}			
	}
	// reduce back to the right dimensions
	distance = sqrt(distance);

	if(distance > 0.0){
		covariance = amp*(1 + root3*(distance/rho))*exp(-root3*(distance/rho));
	} else {
		covariance = amp;
	}

	// this means we're on a diagonal term, golly but i write bad code :(
	if(truecount == nparams){
		covariance += nugget;
	}
	return(covariance);
}

/**
 * compute a matrix of dC/dtheta-length/theta_0 for 
 * the matern 3/2 covariance function c(x,y) = theta_0 * ( 1 + \sqrt{3}r/theta_2) * exp(-\sqrt{3}*r/theta_2) + theta_1
 * where theta_2 is the scale along one of the nparams dimensions 
 * 
 * dc(r)/dtheta_2/theta_0 = (  exp(-\sqrt{3}*r/theta_2) * (3*r^2 / theta_3^3) )
 * 
 * @param covsub = c(x,y) - theta_1
 * @param thetaLength <- value of param we're evaluating gradient at
 * @param index <- must be 2 (since matern 3/2 only has 3 indices)
 * @return dCdTheta <- matrix of dc(r)/dtheta_2
 * 
 */
void derivative_l_matern_three(gsl_matrix *dCdTheta, gsl_matrix* xmodel, double thetaLength, 
															 int index, int nmodel_points, int nparams){ 
	assert(index == 2);
	int i, j, paramIndex;
	double derivVal = 0.0;
	double x_temp = 0.0, y_temp = 0.0;
	gsl_vector_view xmodel_row_i;
	gsl_vector_view xmodel_row_j;
	double root3 = 1.732050808;
	double rtemp = 0.0;
	const int nthetasConstant = 2;
	double thetaLCubed = thetaLength * thetaLength * thetaLength;

	for(i = 0; i < nmodel_points; i++){
		for(j = 0; j < nmodel_points; j++){
			/* extract the nparams long rows of the xmodel matrix */
			xmodel_row_i = gsl_matrix_row(xmodel, i);
			xmodel_row_j = gsl_matrix_row(xmodel, j);
			
			// compute the euclidean distance between the two points
			for(paramIndex = 0; paramIndex < nparams; paramIndex++){
				x_temp = gsl_vector_get(&xmodel_row_i.vector, paramIndex);
				y_temp = gsl_vector_get(&xmodel_row_j.vector, paramIndex);
				rtemp += (x_temp - y_temp) * (x_temp-y_temp);
			}
			rtemp=sqrt(rtemp); 

			derivVal = 3.0*exp(-root3*rtemp/thetaLength)*(rtemp*rtemp / thetaLCubed);

			gsl_matrix_set(dCdTheta, i, j, derivVal);
		}
	}
	
}


	
//! the matern cov fn but with nu set to 5/2
double covariance_fn_matern_five(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams){
	double covariance = 0.0;
	int i, truecount = 0;
	double xm_temp = 0.0;
	double xn_temp = 0.0;
	double distance = 0.0;
	
	assert(nthetas >= 3); // can throw away the upper ones no problem

	// map the thetas onto some local variables so the formula is more transparent
	double amp = gsl_vector_get(thetas, 0);
	double nugget = gsl_vector_get(thetas,1);
	double rho = exp(gsl_vector_get(thetas, 2));

	double tempx = 0.0;
	double root5 = 2.236067978;

	double d_over_r = 0.0;
		
	// calculate the euclidean distance between the two points;
	for(i = 0; i < nparams; i++){
		xm_temp = gsl_vector_get(xm, i);
		xn_temp = gsl_vector_get(xn, i);
		// this is currently the distance squared
		distance += pow(fabs(xm_temp - xn_temp), 2.0);
		if(fabs(xm_temp - xn_temp) < 0.0000000000000001){
			truecount++;
		}			
	}
	// reduce back to the right dimensions
	distance = sqrt(distance);
	d_over_r = distance / rho;
	if(distance > 0.0){
		covariance = amp*(1+root5*(d_over_r)+(5.0/3.0)*(d_over_r)*(d_over_r))*exp(-root5*(d_over_r));
	} else if(distance == 0){
		covariance = amp;
	}
	
	// this means we're on a diagonal term, golly but i write bad code :(
	if(truecount == nparams){
		covariance += nugget;
	}
	return(covariance);
}


/** 
 * compute a matrix of dC/dtheta-length/theta_0 for 
 * the matern 5/2 covariance function 
 * c(x,y) = theta_0 * (1+root5 *(r/theta_2) + (5.0/3.0)*(r/theta_2)^2 ) * exp(-root5 * (r/theta_2))
 * where theta_2 is the scale along one of the only radial dimension
 * 
 * dc / dtheta_2 / theta_0 = (r^2/theta_2^3)*exp(-\sqrt{5}*r/theta_2) * ( 3.72678 *r  + 1.66667 theta_2)
 * 
 * @param thetaLength <- value of param we're evaluating gradient at
 * @param index <- must be 2 (since matern 3/2 only has 3 indices)
 * @param nparams <- the number of design space dimensions we have in our xmodel
 * @return dCdTheta <- matrix of dc(r)/dtheta_2
 * 
 */
void derivative_l_matern_five(gsl_matrix *dCdTheta, gsl_matrix* xmodel, double thetaLength, 
															 int index, int nmodel_points, int nparams){ 
	assert(index == 2);
	int i, j, paramIndex;
	double derivVal = 0.0;
	double x_temp = 0.0, y_temp = 0.0;
	gsl_vector_view xmodel_row_i;
	gsl_vector_view xmodel_row_j;
	double root5 = 2.2360680;
	double rtemp = 0.0;
	double rsq = 0.0;
	const int nthetasConstant = 2;
	double thetaLCubed = thetaLength * thetaLength * thetaLength;

	for(i = 0; i < nmodel_points; i++){
		for(j = 0; j < nmodel_points; j++){
			/* extract the nparams long rows of the xmodel matrix */
			xmodel_row_i = gsl_matrix_row(xmodel, i);
			xmodel_row_j = gsl_matrix_row(xmodel, j);
			
			// compute the euclidean distance between the two points
			for(paramIndex = 0; paramIndex < nparams; paramIndex++){
				x_temp = gsl_vector_get(&xmodel_row_i.vector, paramIndex);
				y_temp = gsl_vector_get(&xmodel_row_j.vector, paramIndex);
				rtemp += (x_temp - y_temp) * (x_temp-y_temp);
			}
			rsq = rtemp;
			rtemp=sqrt(rtemp); 

			derivVal = (rsq / (thetaLCubed)) * exp(-root5 * rtemp / thetaLength) *
				(3.72768*rtemp + 1.66667*thetaLength);

			gsl_matrix_set(dCdTheta, i, j, derivVal);
		}
	}
	
}



//! calculate a vector of the covariances of a given point against all of the model points
/**
 * calculate a vector of (c(xnew, xmodel[1]), c(xnew, xmodel[2]) ....) 
 * this can be thought of as the final row or column of the covariance matrix if you were
 * to augment it with the additional point xnew. 
 *
 * @param xmodel -> the matrix representing all of the model evaluations (nmodel_points x nparams)
 * @param xnew -> the new point at which the kvector is to be calculated
 * @param kvector -> on return this is the result
 * @param thetas -> the hyperparams of the gp used in the estimation.
 * @return kvector is set to the result
 */
void makeKVector(gsl_vector* kvector, gsl_matrix *xmodel, gsl_vector *xnew, gsl_vector* thetas, int nmodel_points, int nthetas, int nparams){
	int i = 0; 
	gsl_vector_view  xmodel_row;
	double cov;
	for (i = 0; i < nmodel_points; i++){
		xmodel_row = gsl_matrix_row(xmodel, i);
		// send the rows from the xmodel matrix to the kvector, these have nparams width
		cov = covariance_fn(&xmodel_row.vector, xnew, thetas, nthetas, nparams);
		if(cov < 1E-10){
			cov = 0.0;
		}
		gsl_vector_set(kvector, i, cov);
	}
}


//! create the covariance matrix
/**
 * calculate the covariance matrix for a given set of model points and hyperparameters theta
 * since xmodel is taken to have have nparams columns and rows with the values of these params for each of the training points
 * we pass slices to the covariance function to create the matrix
 * 
 * @param cov_matrix overwritten by the calculated covariances.
 * @param xmodel matrix of the model points (nmodel_points x nparams)
 * @param thetas hyperparameters
 */
void makeCovMatrix(gsl_matrix *cov_matrix, gsl_matrix *xmodel, gsl_vector* thetas, int nmodel_points, int nthetas, int nparams){
	int i,j;
	double covariance; 
	gsl_vector_view xmodel_row_i;
	gsl_vector_view xmodel_row_j;

	for(i = 0; i < nmodel_points; i++){
		for(j = 0; j < nmodel_points; j++){
			xmodel_row_i = gsl_matrix_row(xmodel, i);
			xmodel_row_j = gsl_matrix_row(xmodel, j);
			covariance = covariance_fn(&xmodel_row_i.vector, &xmodel_row_j.vector, thetas, nthetas,nparams);
			//printf("(%d,%d) cov: %g\n", i, j, covariance);
			gsl_matrix_set(cov_matrix, i,j, covariance);
		}
	}
	//print_matrix(cov_matrix, nmodel_points, nmodel_points);
}
	

//! calculate the emulated mean 
/**
 * calculate the mean at t+1 using the current training vec, the inverse matrix (calculated elsewhere!) and the
 * kplus vector
 * @return the mean at the point used to caclulate kplus_vector
 * @param inverse_cov_matrix the inverse of the covariance matrix
 * @param training_vector -> the scalar output of the model
 * @param kplus_vector -> a K vector, from makeKVector, evaluated at the t+1 point
 * @param h_vector -> the vector of linear regression functions evaluated at t+1
 * @param h_matrix -> matrix of linear regression functions evaluated at all of the design points
 * @param beta_vector -> vector of linear regression coeffs, this should be constant once you have a set of thetas
 * @param nmodel_points -> the length of training_vector and also the length of inverse_cov_matrix
 * 
 * so now emulated_mean = h_vector.beta + kplus_vector.inverse_cov_matrix.(training_vector - h_matrix.beta_vector)
 */
double makeEmulatedMean(gsl_matrix *inverse_cov_matrix, gsl_vector *training_vector, gsl_vector *kplus_vector, gsl_vector* h_vector, gsl_matrix* h_matrix, gsl_vector* beta_vector,  int nmodel_points){
	gsl_vector *result_1 = gsl_vector_alloc(nmodel_points);
	gsl_vector *result_2 = gsl_vector_alloc(nmodel_points);
	double emulated_mean = 0.0;
	double regression_cpt = 0.0; // transpose(h).beta 
	double residual_cpt=0.0;  // -k_vector.inverse_cov_matrix.Hmatrix.beta
	//— Function: int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
	// result_holder = inverse_cov_matrix . training_vector (matrix, vec multiply)
	gsl_blas_dgemv(CblasNoTrans, 1.0, inverse_cov_matrix, training_vector, 0.0, result_1);
	//— Function: int gsl_blas_sdsdot (float alpha, const gsl_vector_float * x, const gsl_vector_float * y, float * result)
	// emulatedMean = kplusvector . result_1
	gsl_blas_ddot(kplus_vector, result_1, &emulated_mean) ;
	// regression_cpt = transpose(h).beta
	gsl_blas_ddot(h_vector, beta_vector, &regression_cpt);

	gsl_vector_set_zero(result_1);
	
	// residual_cpt = -k_vector.inverse_cov_matrix.Hmatrix.beta
	// do result_1 = Hmatrix.beta first (this comes out as a nmodel_points long vector)
	gsl_blas_dgemv(CblasNoTrans, 1.0, h_matrix, beta_vector, 0.0, result_1); 
	// result_2 = inverse_cov_matrix.result_1
	gsl_blas_dgemv(CblasNoTrans, 1.0, inverse_cov_matrix, result_1, 0.0, result_2);
	// -residual_cpt = kplusvector . result_2
	gsl_blas_ddot(kplus_vector, result_2, &residual_cpt);

	// now we have the correct result
	// m(x) = h(x)^{T}\hat\beta + t(x)^{T}A^{-1}(y - H\hat\beta)
	//emulated_mean += regression_cpt - residual_cpt;

	gsl_vector_free(result_1);
	gsl_vector_free(result_2);
	return(regression_cpt + emulated_mean -residual_cpt);
}


//! create the emulated variance
/**
 * calculate the variance at t+1 using the kplus vector (for that point), the inverse matrix, and kappa the 
 * constant variance for the process.
 * @return -> the new variance
 * @param kappa -> the constant variance for the process.
 * @param h_vector -> regression functions evaluated at t+1
 * @param h_matrix -> matrix of the regression fns evaluated at the design points
 * 
 * the regression cpt to this is a bit more complicated than for the mean we have
 * c*(x,x') = c(x,x') - t(x)^{T}A^{-1}t(x') + (h(x)^{T} - t(x)^{T}A^{-1}H)(H^{T}A^{-1}H)^{-1}(h(x')^{T} - t(x')^{T}A^{-1}H)
 *
 */
double makeEmulatedVariance(gsl_matrix *inverse_cov_matrix, gsl_vector *kplus_vector, gsl_vector *h_vector, gsl_matrix *h_matrix, double kappa, int nmodel_points, int nregression_fns){
	double emulated_variance;
	double regression_cpt;
	gsl_vector *result_nreg = gsl_vector_alloc(nregression_fns);
	gsl_vector *result_nreg_2 = gsl_vector_alloc(nregression_fns);
	gsl_vector *result_holder = gsl_vector_alloc(nmodel_points);
	gsl_matrix *result_minverse_dot_h = gsl_matrix_alloc(nmodel_points, nregression_fns);
	gsl_matrix *result_inverse_h_minverse_h = gsl_matrix_alloc(nregression_fns, nregression_fns);	
	gsl_error_handler_t *temp_handler;
	int i,j, cholesky_test;
	
	// result_nreg = (h(x)^T - t(x)^T A^(-1).H) 
	// where t -> kplus_
	//int gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, inverse_cov_matrix, h_matrix, 0.0, result_minverse_dot_h);
	// now result_nreg = t(x)^(T).result_minverse_dot_h
	// but since we can't do vector.matrix we do
	// result_nreg = result_minverse_dot_h^(T).t(x)
	gsl_blas_dgemv(CblasTrans, -1.0, result_minverse_dot_h, kplus_vector, 0.0, result_nreg);
	// note that we scaled the previous matrix multiply by -1.0
	// gsl_vector_add(a,b) := a <- a + b so
	gsl_vector_add(result_nreg, h_vector);
	
	// result_inverse_h_minverse_h := (H^{t} .A^{-1} . H)^{-1}
	//                              = H^{t} . result_minverse.h
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, h_matrix, result_minverse_dot_h, 0.0, result_inverse_h_minverse_h);
	
	// now result_inverse_h_minverse_h is symmetric and so we can consider it as triangular etc and then
	// do linear solves etc to get the inverse and all that malarky, but since it's only nreg * nreg its 
	// probably not worth the pain in the blass (arf)
	temp_handler = gsl_set_error_handler_off();
	cholesky_test = gsl_linalg_cholesky_decomp(result_inverse_h_minverse_h);
	if(cholesky_test == GSL_EDOM){
		FILE *fptr;		
		fprintf(stderr, "trying to cholesky a non postive def matrix, sorry...\n");
		fprintf(stderr, "matrix dumped to chol-err.dat\n");
		fptr = fopen("chol-err.dat", "w");
		for(i = 0; i < nmodel_points; ++i){
			for(j = 0; j < nmodel_points; ++j){
				fprintf(fptr, "%lf", gsl_matrix_get(h_matrix, i, j));				
			}
			fprintf(fptr, "\n");
		}
		fclose(fptr);
		exit(1);
	}
	gsl_set_error_handler(temp_handler);


	gsl_linalg_cholesky_invert(result_inverse_h_minverse_h); // now it really is the inverse
	
	// now compute regression_cpt = result_nreg . (result_inverse_h_minverse_h)^{-1} . result_nreg
	// result_nreg_2 = result_inverse_h_minverse_h^{-1}.result_nreg
	gsl_blas_dgemv(CblasNoTrans, 1.0, result_inverse_h_minverse_h, result_nreg, 0.0, result_nreg_2);
	// regression_cpt = result_nreg. result_nreg_2
	gsl_blas_ddot(result_nreg, result_nreg_2, &regression_cpt);

	gsl_blas_dgemv(CblasNoTrans, 1.0, inverse_cov_matrix, kplus_vector, 0.0, result_holder);
	gsl_blas_ddot(kplus_vector, result_holder, &emulated_variance);
	gsl_vector_free(result_holder);
	gsl_vector_free(result_nreg);
	gsl_vector_free(result_nreg_2);
	gsl_matrix_free(result_minverse_dot_h);
	gsl_matrix_free(result_inverse_h_minverse_h);
	return(kappa-emulated_variance + regression_cpt);
}


//! fill in the new_x matrix for 1 and 2d
/**
 * Given an empty matrix of nparams x nemulate_points this function inits the points into a square lattice
 * however so far this only works for nparams = 1, 2
 */
void initialise_new_x(gsl_matrix* new_x, int nparams, int nemulate_points, double emulate_min, double emulate_max){
	int i, j;
	int n_side;
	double step_size;
	
	if (nparams == 1){
		step_size = (emulate_max - emulate_min) / ((double)nemulate_points);	
		for(i = 0; i < nemulate_points;i++){
			gsl_matrix_set(new_x, i, 0, step_size*((double)i)+emulate_min);
		}
	} else if(nparams == 2){
		n_side = floor(sqrt(nemulate_points));
		step_size = (emulate_max - emulate_min) / ((double)n_side);	
		for(i = 0; i < n_side; i++){
			for(j = 0; j < n_side; j++){
				gsl_matrix_set(new_x, i*n_side+j,0, step_size*((double)(i))+emulate_min);
				gsl_matrix_set(new_x, i*n_side+j, 1, step_size*((double)(j))+emulate_min);
			}
		}
	} else{
		fprintf(stderr, "oops there's no support for %d'd problems yet!\n", nparams);
	}
	//print_matrix(new_x, n_emu_points, nparams);
	
}	


/**
 * deprecated functions
 * /

/* //! calculate the covariance between xm and xn but with a full matrix of covariance parameters between */
/* //! xm and xn, not just diagonal entries */
/* /\**  */
/*  * this needs nthetas to go lke nparams ^2 + 2 */
/*  * doesn't seem to work because this covariance matrix is not symmetric. */
/*  * doh */
/*  *\/ */
/* double covariance_fn_gaussian_nondiag(gsl_vector* xm, gsl_vector*xn, gsl_vector*thetas, int nthetas, int nparams ){ */
/* 	double the_covariance = 0.0;   */
/* 	double nugget = gsl_vector_get(thetas, 1); */
/* 	double amplitude = gsl_vector_get(thetas,0); */
/* 	double sigma_temp = 0.0; */
/* 	double distance = 0.0; */
/* 	double r_temp = 0.0; */
/* 	double small_no = 1E-10; */
/* 	int i, j, diagcount = 0;  */

/* 	double alpha = 1.90; */

/* 	assert((int)thetas->size == (nparams*nparams) + 2); */

/* 	for(i = 0; i < nparams; i++){ */
/* 		for(j = 0; j < nparams; j++){ */
/* 			// the vairance matrix elements are index-offest by 2, since we have the nugget and amplitude stored in the same vector */
/* 			sigma_temp = gsl_vector_get(thetas, i+2); */
/* 			r_temp = gsl_vector_get(xm,i) - gsl_vector_get(xn, j); */
/* 			if(sigma_temp  > small_no){ */
/* 				distance += pow(fabs(r_temp), alpha)/sigma_temp; */
/* 			}  */
/* 		} */
/* 		if(fabs(gsl_vector_get(xm,i) - gsl_vector_get(xn,i)) < 1E-10) diagcount++; */
/* 	} */
	
/* 	the_covariance = amplitude * exp(-distance); */
/* 	if(diagcount == nparams) the_covariance += nugget; */
/* 	//fprintf(stderr, "%g\n", the_covariance); */
/* 	return(the_covariance); */
/* } */
	



/* // this only works with exactly 4 thetas! no matter how many params there are */
/* //! calculates the covariance function using the matern metric. http://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function */
/* /\** */
/*  * you can switch which covariance function is called in covariance_fn */
/*  *  */
/*  * not well tested yet */
/*  *\/ */
/* double covariance_fn_matern(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams){ */
/* 	double covariance = 0.0; */
/* 	int i, truecount = 0; */
/* 	double xm_temp = 0.0; */
/* 	double xn_temp = 0.0; */
/* 	double distance = 0.0; */
	
/* 	assert(nthetas >= 4); // can throw away the upper ones no problem */
	
/* 	// map the thetas onto some local variables so the formula is more transparent */
/* 	double sigsquared = gsl_vector_get(thetas, 0); */
/* 	double rho = gsl_vector_get(thetas, 1); */
/* 	double nu = gsl_vector_get(thetas, 2); */
/* 	double nugget = gsl_vector_get(thetas,3); */
/* 	double tempx = 0.0; */

/* 	// calculate the euclidean distance between the two points; */
/* 	for(i = 0; i < nparams; i++){ */
/* 		xm_temp = gsl_vector_get(xm, i); */
/* 		xn_temp = gsl_vector_get(xn, i); */
/* 		// this is currently the distance squared */
/* 		distance += pow(fabs(xm_temp - xn_temp), 2.0); */
/* 		if(fabs(xm_temp - xn_temp) < 0.0000000000000001){ */
/* 			truecount++; */
/* 		}			 */
/* 	} */
/* 	// reduce back to the right dimensions */
/* 	distance = sqrt(distance); */
	
/* 	//fprintf(stderr, "nu = %g, x = %g\n", nu, 2*sqrt(nu)*(distance/rho)); */

/* 	tempx = (2.0*sqrt(nu)*distance/rho); */

/* 	if(distance > 0.0){ */
/* 		// debug */
/* 		if(nu < 0){ */
/* 			fprintf(stderr, "nu = %g, x = %g\n", nu, tempx); */
/* 		} */
/* 		assert(nu > 0); // otherwise it'll crash anyway */
/* 		//fprintf(stderr, "tempx = %g\n", tempx); */
/* 		/\* */
/* 		 * besselK -> 0 really fast so we just cut of at 300 and replace that part with zero */
/* 		 * this might not be close enough but it should stop the overflow */
/* 		 *\/ */
/* 		if(tempx < 300){  */
/* 			covariance = sigsquared * (1.0/(gsl_sf_gamma(nu)*pow(2.0, nu - 1.0))); */
/* 			covariance *= pow(tempx, nu)*(gsl_sf_bessel_Knu(nu, tempx)); */
/* 		} else { */
/* 			covariance = 0.0; */
/* 		} */
/* 		/\*  */
/* 		 * in the limit that they're at the same place just */
/* 		 * apply the regular process variance  */
/* 		 *\/ */
/* 	} else if (distance == 0){  */
/* 		covariance = sigsquared;  */
/* 	} else { */
/* 		fprintf(stderr, "distance is negative in covariance_fn_matern!\n"); */
/* 		fprintf(stderr, "nu = %g, x = %g\n", nu, 2*sqrt(nu)*(distance/rho)); */
/* 		exit(1); */
/* 	} */
/* 	// this means it's diagonal, i.e the distance is less than zero */
/* 	if(truecount == nparams){ */
/* 		covariance += nugget; */
/* 	} */

/* 	return(covariance); */
/* } */
