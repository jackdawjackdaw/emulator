#include "emulator.h"

void print_matrix(gsl_matrix* m, int nx, int ny){
	int i,j;
	for(i = 0; i < nx; i++){
		for(j= 0; j < ny; j++){
			printf("%g ", gsl_matrix_get(m, i, j));
		}
		printf("\n");
 	}
}


//! calculate the covariance between a set of input points
// note that the xm is a vector view (from the 

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
		if (fabs(xm_temp - xn_temp) < 0.0000000000000001){
			truecount++; 		
		}
	}
	covariance = covariance * gsl_vector_get(thetas,0) + gsl_vector_get(thetas,1);

	if(truecount == nparams) {
		// i.e the two vectors are hopefully the same
		// add the nugget
		covariance += gsl_vector_get(thetas,2);
	}

	return(covariance);
}

// kvector is of length nmodel_points
// xmodel is a matrix of nmode_points * nparams
void makeKVector(gsl_vector* kvector, gsl_matrix *xmodel, gsl_vector *xnew, gsl_vector* thetas, int nmodel_points, int nthetas, int nparams){
	int i = 0; 
	gsl_vector_view  xmodel_row;
	for (i = 0; i < nmodel_points; i++){
		xmodel_row = gsl_matrix_row(xmodel, i);
		// send the rows from the xmodel matrix to the kvector, these have nparams width
		gsl_vector_set(kvector, i, covariance_fn(&xmodel_row.vector, xnew, thetas, nthetas, nparams));
	}
}


// create the covariance matrix
// note the use of vector views again. 
// since xmodel is taken to have have nparams columns and rows with the values of these params for each of the training points
// we pass slices to the covariance function to create the matrix
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
	


// calculate the mean at t+1 using the current training vec, the inverse matrix (calculated elsewhere!) and the
// kplus vector
// @return -> the mean
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


// calculate the variance at t+1 using the kplus vector (for that point), the inverse matrix, and kappa the 
// constant variance for the process.
// @return -> the new variance
double makeEmulatedVariance(gsl_matrix *inverse_cov_matrix, gsl_vector *kplus_vector, double kappa, int nmodel_points){
	double emulated_variance;
	gsl_vector *result_holder = gsl_vector_alloc(nmodel_points);
	gsl_blas_dgemv(CblasNoTrans, 1.0, inverse_cov_matrix, kplus_vector, 0.0, result_holder);
	gsl_blas_ddot(kplus_vector, result_holder, &emulated_variance);
	gsl_vector_free(result_holder);
	return(emulated_variance);
}





