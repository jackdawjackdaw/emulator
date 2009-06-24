#include "maximise.h"

#define FAILVALUE -100

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
void gradDesc(gsl_rng *rand, int max_tries, int nsteps, double gamma, gsl_matrix* ranges,  gsl_matrix *xmodel, gsl_vector *trainingvector, gsl_vector *thetas, int nmodel_points, int nthetas, int nparams){
	int i,j;
	int tries = 0;
	double best_value = -50.0;
	double temp_val = 0.0;

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
		set_random_initial_value(rand, xOld, ranges, nthetas);

		for(i = 0; i < nsteps; i++){
			
			// make the covariance matrix 
			// using the random initial conditions! (xold not thetas)
			makeCovMatrix(covariance_matrix, xmodel, xOld, nmodel_points, nthetas, nparams);
			gsl_matrix_memcpy(temp_matrix, covariance_matrix);
			gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
			gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse); // now we have the inverse
			
			// now get the determinant of the inverse
			// actually don't need that here, need it when we're seeing to keep the estimate
			//gsl_matrix_memcpy(temp_matrix, cinverse);
			//gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
			//cinverse_det = gsl_linalg_LU_det(temp_matrix, lu_signum);

			for(j = 0; j < nthetas; j++){ /* calc the gradient vector*/
				the_gradient =  getGradient(cinverse, xmodel, trainingvector, xOld, j, nmodel_points, nthetas, nparams);
				gsl_vector_set(gradient_vec, j, the_gradient);
			}

			
			// xNew = xOld + gamma*the_gradient(xold)
			// there's probably a neat function for this but i can't be bothered
			for(j = 0; j < nthetas; j++){
				temp_val = gsl_vector_get(xOld, j) + gamma*gsl_vector_get(gradient_vec, j);
				gsl_vector_set(xNew, j, temp_val);
			}
			
			if( range_check(xNew, ranges, nthetas) == 1){
				gsl_vector_memcpy(xOld, xNew); // set the next xOld value
			} else {
				for(j = 0; j < nthetas; j++){
					gsl_vector_set(xNew, j, FAILVALUE);					
				}
				fprintf(stderr, "gradDesc fell out of range\n");
				break; // get out of this loop, it didn't work
			}
		}
		tries++;
		if ( vector_components_equal(xNew, FAILVALUE, nthetas) != 1){
			// the thing didn't get set to the failval
			// now we want to calc the logLikelyhood and see how it is
			gsl_matrix_memcpy(temp_matrix, cinverse);
			gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
			cinverse_det = gsl_linalg_LU_det(temp_matrix, lu_signum);
			// temp_val is now the likelyhood for this answer
			temp_val = getLogLikelyhood(cinverse, cinverse_det, xmodel, trainingvector, xNew, nmodel_points, nthetas, nparams);
			printf("log likelyhood = %g\n", temp_val);

			if(temp_val > best_value){
				//debug
				best_value = temp_val;
				printf("best: %g\n", best_value);
				gsl_vector_memcpy(best_vector, xNew);
			}			
		} else {
			fprintf(stderr, "caught a failed run, ignoring\n");
		}
	}
	
	printf("final-best = %g\n", best_value);

	// now we copy the best value into thetas
	gsl_vector_memcpy(thetas, best_vector);
	
	// clean up extra junk
	gsl_vector_free(best_vector);
	gsl_vector_free(xOld);
	gsl_vector_free(xNew);
	gsl_vector_free(gradient_vec);
	
	gsl_matrix_free(covariance_matrix);
	gsl_matrix_free(cinverse);
	gsl_matrix_free(temp_matrix);
	gsl_permutation_free(c_LU_permutation);
	
}

//! see if a given vector has it's components in the correct range
/**
 * @return returns 1 if true
 */
int range_check(gsl_vector* x, gsl_matrix* ranges, int nthetas){
	int i;
	double range_min;
	double range_max;
	int truecount = 0;
	for(i = 0; i < nthetas; i++){
		range_min = gsl_matrix_get(ranges, i, 0);
		range_max = gsl_matrix_get(ranges, i, 1);
		if  ((gsl_vector_get(x, i) <= range_max) && (gsl_vector_get(x,i) > range_min)){
			truecount++;
		}
	}
	if (truecount == nthetas){
		// all values of x are in range
		return(1);
	} else {
		return(0);
	}
}

//! see if a given vector has all of its components equal to the given value
/**
 * @return 1 if the test passes, zero otherwise 
 */
int vector_components_equal(gsl_vector *x, double test_value, int nthetas){
	int i;
	int truecount = 0;
	for(i = 0; i < nthetas; i++){
		// double comparision, but should be ok here
		if(gsl_vector_get(x, i) == test_value){
			truecount++;
		}
	}	
	if(truecount == nthetas){
		return(1);
	} else {
		return(0);
	}
}
 
	
//! set a random starting point in nthetas dimensions, sampled from the ranges
/**
 * @param rand a pre setup gsl_rng
 * @return x is set to the correct values
 */
void set_random_initial_value(gsl_rng* rand, gsl_vector* x, gsl_matrix* ranges,int  nthetas){
	int i;
	double range_min; 
	double range_max;
	double the_value;
	
	for(i = 0; i < nthetas; i++){
		range_min = gsl_matrix_get(ranges, i, 0);
		range_max = gsl_matrix_get(ranges, i, 1);
		// set the input vector to a random value in the range
		the_value = gsl_rng_uniform(rand) * (range_max - range_min) + range_min;
		//printf("theta %d set to %g\n", i, the_value);
		//gsl_vector_set(x, gsl_rng_uniform(rand)*(range_max-range_min)+range_min, i);
		gsl_vector_set(x, i, the_value);
	}
}
	








