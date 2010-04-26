#include "emulator.h"
/** 
 * @file 
 * @author Chris Coleman-Smith cec24@phy.duke.edu
 * @version 0.1
 * @section DESCRIPTION
 * 
 * This file contains the main functions needed to run the emulator, given a set of 
 * correct hyperparameters (theta) then the makeEmulatedMean and makeEmulatedVariance
 * functions can be used  
 */

//! a GLOBAL structure which will be queried by the functions to set params
emulator_opts the_emulator_options;


//! this MUST BE CALLED FIRST before any emulatoring
void set_emulator_defaults(emulator_opts* x){
	// set some default behaviour
	x->alpha = 1.9;
	// don't use any of the matern functions
	x->usematern = 0;
	x->usematern_three = 0;
	x->usematern_five = 0;
}
	

void print_emulator_options(emulator_opts* x){
	if(x->usematern ==0 && x->usematern_three == 0 && x->usematern_five==0){
		fprintf(stderr, "using a power-exp covariance function\n");
		fprintf(stderr, "alpha = %g\n", x->alpha);				
	} else if (x->usematern == 1){
		fprintf(stderr, "using a *FULL* matern  covariance function\n");
	} else if (x->usematern_three ==1){
		fprintf(stderr, "using 3/2 matern  covariance function\n");
	} else {
		fprintf(stderr, "using 5/2 matern  covariance function\n");
	}
}



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


//! wrapper fn, calls the appropriate covariance. Change this by hand!
/** 
 * this is now just a wrapper which calls the appropriate covariance function
 * it seems that the matern version doesn't work very well compared to the 
 * gaussian covariance function, at least not on the ising model.
 */
inline double covariance_fn(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams){
	if(the_emulator_options.usematern == 1){
		return(covariance_fn_matern(xm, xn, thetas, nthetas, nparams));
	} else if(the_emulator_options.usematern_three ==1){
		return(covariance_fn_matern_three(xm,xn, thetas, nthetas, nparams));
	} else if(the_emulator_options.usematern_five ==1){
		return(covariance_fn_matern_five(xm,xn, thetas, nthetas, nparams));
	} else {
		// this is the default option
		return(covariance_fn_gaussian(xm, xn , thetas, nthetas, nparams, the_emulator_options.alpha));
	}
}

/** 
 * and you'll kneel in the market place and draw your vast mandala
 *
 */
 


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
 * theta2 -> a confusing offset without which this doens't work at all
 * theta3 (and higher) -> scale parameter for the gaussian
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
 */
double covariance_fn_gaussian(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams, double alpha){
	// calc the covariance for a given set of input points
	int i, truecount  = 0;
	double covariance = 0.0;
	double exponent = 0.0;
	double xm_temp = 0.0;
	double xn_temp = 0.0;
	double r_temp = 0.0;
	double amp = exp(gsl_vector_get(thetas, 0));
	double nug = gsl_vector_get(thetas, 1);
	
	for(i = 0; i < nparams; i++){
		xm_temp = gsl_vector_get(xm, i);  // get the elements from the gsl vector, just makes things a little clearer
		xn_temp = gsl_vector_get(xn, i);
		r_temp = exp(gsl_vector_get(thetas, i+2));
		r_temp = pow(r_temp , alpha); 
		// gaussian term				
		exponent += (-1.0/2.0)*pow(xm_temp-xn_temp, alpha)/(r_temp);
		//DEBUGprintf("%g\n", covariance);
		if (fabs(xm_temp - xn_temp) < 0.0000000000000001){
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
		// add the nugget
		covariance += gsl_vector_get(thetas,1);
	}


	return(covariance);
}


// this only works with exactly 4 thetas! irregardless of how many params there are
//! calculates the covariance function using the matern metric. http://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
/**
 * you can switch which covariance function is called in covariance_fn
 * 
 * not well tested yet
 */
double covariance_fn_matern(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams){
	double covariance = 0.0;
	int i, truecount = 0;
	double xm_temp = 0.0;
	double xn_temp = 0.0;
	double distance = 0.0;
	
	assert(nthetas >= 4); // can throw away the upper ones no problem
	
	// map the thetas onto some local variables so the formula is more transparent
	double sigsquared = gsl_vector_get(thetas, 0);
	double rho = gsl_vector_get(thetas, 1);
	double nu = gsl_vector_get(thetas, 2);
	double nugget = gsl_vector_get(thetas,3);
	double tempx = 0.0;

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
	
	//fprintf(stderr, "nu = %g, x = %g\n", nu, 2*sqrt(nu)*(distance/rho));

	tempx = (2.0*sqrt(nu)*distance/rho);

	if(distance > 0.0){
		// debug
		if(nu < 0){
			fprintf(stderr, "nu = %g, x = %g\n", nu, tempx);
		}
		assert(nu > 0); // otherwise it'll crash anyway
		//fprintf(stderr, "tempx = %g\n", tempx);
		/*
		 * besselK -> 0 really fast so we just cut of at 300 and replace that part with zero
		 * this might not be close enough but it should stop the overflow
		 */
		if(tempx < 300){ 
			covariance = sigsquared * (1.0/(gsl_sf_gamma(nu)*pow(2.0, nu - 1.0)));
			covariance *= pow(tempx, nu)*(gsl_sf_bessel_Knu(nu, tempx));
		} else {
			covariance = 0.0;
		}
		/* 
		 * in the limit that they're at the same place just
		 * apply the regular process variance 
		 */
	} else if (distance == 0){ 
		covariance = sigsquared; 
	} else {
		fprintf(stderr, "distance is negative in covariance_fn_matern!\n");
		fprintf(stderr, "nu = %g, x = %g\n", nu, 2*sqrt(nu)*(distance/rho));
		exit(1);
	}
	// this means it's diagonal, i.e the distance is less than zero
	if(truecount == nparams){
		covariance += nugget;
	}

	return(covariance);
}

//! the matern cov fn but with nu set to 3/2
/**
 * uses 3 thetas only, 0 and 1 are for the vert scale and the length scale and 
 * 3 is the nugget 
 * 
 * this doesn't use quite the same normalisation as the other matern function
 * i don't think this will make much difference although you can't directly 
 * compare the "best" hyperparams
 */
double covariance_fn_matern_three(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams){
	double covariance = 0.0;
	int i, truecount = 0;
	double xm_temp = 0.0;
	double xn_temp = 0.0;
	double distance = 0.0;
	
	assert(nthetas >= 3); // can throw away the upper ones no problem

	// map the thetas onto some local variables so the formula is more transparent
	double sigsquared = gsl_vector_get(thetas, 0);
	double rho = gsl_vector_get(thetas, 1);
	double nugget = gsl_vector_get(thetas,2);
	double tempx = 0.0;
	double root3 = 1.732050808;

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

	if(distance > 0.0){
		covariance = sigsquared*(1 + root3*(distance/rho))*exp(-root3*(distance/rho));
	} else if(distance == 0){
		covariance = sigsquared;
	}
	
	// this means we're on a diagonal term, golly but i write bad code :(
	if(truecount == nparams){
		covariance += nugget;
	}
	return(covariance);
}

	
//! the matern cov fn but with nu set to 5/2
/**
 * uses 3 thetas only, 0 and 1 are for the vert scale and the length scale and 
 * 3 is the nugget 
 * 
 * 
 * same hyperparam normalisation as covariance_fn_matern_three,
 */
double covariance_fn_matern_five(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams){
	double covariance = 0.0;
	int i, truecount = 0;
	double xm_temp = 0.0;
	double xn_temp = 0.0;
	double distance = 0.0;
	
	assert(nthetas >= 3); // can throw away the upper ones no problem

	// map the thetas onto some local variables so the formula is more transparent
	double sigsquared = gsl_vector_get(thetas, 0);
	double rho = gsl_vector_get(thetas, 1);
	double nugget = gsl_vector_get(thetas,2);
	double tempx = 0.0;
	double root5 = 2.236067978;

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

	if(distance > 0.0){
		covariance = sigsquared*(1+root5*(distance/rho)+(5.0/3.0)*pow((distance/rho),2.0))*exp(-root5*(distance/rho));
	} else if(distance == 0){
		covariance = sigsquared;
	}
	
	// this means we're on a diagonal term, golly but i write bad code :(
	if(truecount == nparams){
		covariance += nugget;
	}
	return(covariance);
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
	for (i = 0; i < nmodel_points; i++){
		xmodel_row = gsl_matrix_row(xmodel, i);
		// send the rows from the xmodel matrix to the kvector, these have nparams width
		gsl_vector_set(kvector, i, covariance_fn(&xmodel_row.vector, xnew, thetas, nthetas, nparams));
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
			covariance = covariance_fn(&xmodel_row_i.vector, &xmodel_row_j.vector, thetas, nthetas, nparams);
			gsl_matrix_set(cov_matrix, i,j, covariance);
		}
	}
}
	

//! caculate the emulated mean 
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
	gsl_linalg_cholesky_decomp(result_inverse_h_minverse_h);
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
