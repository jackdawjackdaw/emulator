#include "maximise.h"

/** 
 * @file 
 * @author Chris Coleman-Smith cec24@phy.duke.edu
 * @version 0.1
 * @section DESCRIPTION
 * 
 * contains the routines to attempt the maximum likleyhood estimation
 */

//! attempts a *stupid* gradient descent (ascent) maximisation
/**
 * does gradeint descent on the system
 * 
 * @param ranges -> a matrix (nthetas * 2) where each row is the max and min allowable range to search for hyperparams in, use this to seed the 
 * initial search etc.
 */
void gradDesc(int max_tries, int nsteps, int gamma, gsl_matrix* ranges,  gsl_matrix *xmodel, gsl_vector *trainingvector, gsl_vector *thetas, int nmodel_points, int nthetas, int nparams){
	int i,j;
	int tries = 0;
	double best_value = 0.0;
	// vectors for the grad desc
	gsl_vector *gradient_vec = gsl_vector_alloc(nthetas);
	gsl_vector *xNew = gsl_vector_alloc(nthetas);
	gsl_vector *xOld = gsl_vector_alloc(nthetas);
	gsl_vector *best_vector = gsl_vector_alloc(nthetas);
	
	// matricies etc for the actual estimation
	gsl_matrix *covariance_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix *cinverse = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix *temp_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_permutation *c_LU_permutation = gsl_permutation_alloc(nmodel_points);
	double cinverse_det = 0.0;
	double the_gradient = 0.0;
	int lu_signum = 0;
	
	while(tries < max_tries){
		set_random_initial_value(xOld, ranges, nthetas);
		
		for(i = 0; i < nsteps; i++){
			
			// make the covariance matrix
			makeCovmatrix(covariance_matrix, xmodel, thetas, nmodel_points, nthetas, nparams);
			gsl_matrix_memcpy(temp_matrix, covariance_matrix);
			gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
			gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse); // now we have the inverse
			
			// now get the determinant of the inverse
			// actually don't need that here, need it when we're seeing to keep the estimate
			//		gsl_matrix_memcpy(temp_matrix, cinverse);
			//gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
			//cinverse_det = gsl_linalg_LU_det(temp_matrix, lu_signum);

			for(j = 0; j < nthetas; j++){ /* calc the gradient vector*/
				the_gradient =  getGradient(cinverse, xmodel, trainingvector, thetas, i, nmodel_points, nthetas, nparams);
				gsl_vector_set(gradient_vec, i, the_gradient);
			}
			


















