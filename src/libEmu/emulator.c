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
 * the hyperparameters are used in this function, 
 * theta0 -> amplitude of the covariance gaussian
 * theta1 -> offset from 0 
 * theta2 -> the nugget (only for diagonal terms) 
 * theta3 (and higher) -> scale parameter for the gaussian
 *
 * @param thetas -> hyperparameters 
 * @param nparams -> the legnth of xm and xn
 * @return the covariance
 * @param thetas -> vector of the hyperparameters of the process. 
 * @param nthetas -> length of thetas
 */
double covariance_fn(gsl_vector *xm, gsl_vector* xn, gsl_vector* thetas, int nthetas, int nparams){
	// calc the covariance for a given set of input points
	int i, truecount  = 0;
	double covariance = 0.0;
	double xm_temp = 0.0;
	double xn_temp = 0.0;
	double r_temp = 0.0;
	for(i = 0; i < nparams; i++){
		xm_temp = gsl_vector_get(xm, i);  // get the elements from the gsl vector, just makes things a little clearer
		xn_temp = gsl_vector_get(xn, i);
		r_temp = gsl_vector_get(thetas, i+3);
		r_temp = r_temp * r_temp; 
		// gaussian term
		// the type coersion here is VERY bad, if you use a rational 1/2 then death
		covariance += exp((-1.0/2.0)*((xm_temp-xn_temp)*(xm_temp-xn_temp))/(r_temp));
		//DEBUGprintf("%g\n", covariance);
		/*
		 * this is slightly dangerous float comparison
		 */
		if (fabs(xm_temp - xn_temp) < 0.0000000000000001){
			truecount++; 		
		}
	}
	covariance = covariance * gsl_vector_get(thetas,0) + gsl_vector_get(thetas,1);

	/** 
	 * the nugget is only added to the diagonal covariance terms,
	 * the parts where xm == xn (vectorwise)
	 */
	if(truecount == nparams) {
		// i.e the two vectors are hopefully the same
		// add the nugget
		covariance += gsl_vector_get(thetas,2);
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
 * @param nmodel_points -> the length of training_vector and also the length of inverse_cov_matrix
 */
double makeEmulatedMean(gsl_matrix *inverse_cov_matrix, gsl_vector *training_vector, gsl_vector *kplus_vector, int nmodel_points){
	gsl_vector *result_holder = gsl_vector_alloc(nmodel_points);
	double emulated_mean;
	//— Function: int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
	// result_holder = inverse_cov_matrix . training_vector (matrix, vec multiply)
	gsl_blas_dgemv(CblasNoTrans, 1.0, inverse_cov_matrix, training_vector, 0.0, result_holder);
	//— Function: int gsl_blas_sdsdot (float alpha, const gsl_vector_float * x, const gsl_vector_float * y, float * result)
	// emulatedMean = kplusvector . result_holder
	gsl_blas_ddot(kplus_vector, result_holder, &emulated_mean);
	gsl_vector_free(result_holder);
	return(emulated_mean);
}


//! create the emulated variance
/**
 * calculate the variance at t+1 using the kplus vector (for that point), the inverse matrix, and kappa the 
 * constant variance for the process.
 * @return -> the new variance
 * @param kappa -> the constant variance for the process.
 */
double makeEmulatedVariance(gsl_matrix *inverse_cov_matrix, gsl_vector *kplus_vector, double kappa, int nmodel_points){
	double emulated_variance;
	gsl_vector *result_holder = gsl_vector_alloc(nmodel_points);
	gsl_blas_dgemv(CblasNoTrans, 1.0, inverse_cov_matrix, kplus_vector, 0.0, result_holder);
	gsl_blas_ddot(kplus_vector, result_holder, &emulated_variance);
	gsl_vector_free(result_holder);
	return(kappa-emulated_variance);
}





