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
	

int compare_evals(const evalunit *a, const evalunit *b){
	double temp = a->value - b->value;
	if(temp > 0){
		return 1;
	} else if(temp < 0){
		return -1;
	} else 
		return 0;
}

//! a simplex method for maximis
void nelderMead(gsl_rng *rand, int max_tries, int nsteps, gsl_vector* the_answer,  gsl_matrix* ranges, gsl_matrix* xmodel, gsl_vector *trainingvector, gsl_vector* thetas, int nmodel_points, int nthetas, int nparams){
	int i,j;
	int t = 0;
	int tries = 0;
	int nverticies  = nthetas+1;
	int nanflag = 0;
	double alpha = 1.0;
	double gamma = 2.0;
	double rho = 0.5;
	double sigma = 0.5;
	evalunit* evalList;
	gsl_vector_view vertex;
	gsl_vector *com = gsl_vector_alloc(nthetas); // the centre of mass of the simplex
	gsl_vector *temp_vector = gsl_vector_alloc(nthetas);
	gsl_vector *xreflected = gsl_vector_alloc(nthetas);
	gsl_vector *xcontracted = gsl_vector_alloc(nthetas);
	gsl_vector *xextended = gsl_vector_alloc(nthetas);
	// simplex is nthetas+1 points (each of which have dim nthetas)
	gsl_matrix *verticies = gsl_matrix_alloc(nverticies, nthetas);
	gsl_matrix *new_verticies = gsl_matrix_alloc(nverticies, nthetas);
	int last_vertex = nverticies -1; // since list starts from 0

	double best_value = -50.0;
	double answer_likelyhood;
	gsl_vector *best_vector = gsl_vector_alloc(nthetas);

	double the_likelyhood = 0.0;
	double evalReflected = 0.0;
	double evalContracted = 0.0;
	double evalExtended = 0.0;
	
	evalList = malloc(sizeof(evalunit)*nverticies);

	
	while(tries < max_tries){
		// assign initial values to the vertices
		for(i = 0; i < nverticies; i++){
			vertex = gsl_matrix_row(verticies, i);
			set_random_initial_value(rand, &vertex.vector, ranges, nthetas);
		}
		



		printf("INITIAL VERTICIES\n");
		print_matrix(verticies, nverticies, nthetas);
		


								
		while(t < nsteps){				
			// fill the evalList
			for(i = 0; i < nverticies; i++){
				vertex = gsl_matrix_row(verticies, i);
				the_likelyhood = evalLikelyhood(&vertex.vector, xmodel, trainingvector, nmodel_points, nthetas, nparams);
				evalList[i].index = i;
				evalList[i].value = the_likelyhood;
			}
									
			// now we have to sort the evallist by value
			qsort(evalList, nverticies, sizeof(evalunit), (void*)compare_evals);
			
			for(i = 0; i < nverticies; i++){ 
				if((isnan(evalList[i].value)) || isinf(evalList[i].value)){
					fprintf(stderr, "nan!\n");
					nanflag = 1;
				}
			}
			if(nanflag == 1){ break;};
				
			// now we use the sorted evalList to sort the vertex matrix
			sort_vertex_list(verticies, evalList, nverticies, nthetas);			
			// calculate the com of the simplex
			// but we have to ignore the last point
			calc_com(verticies, com, nverticies, nthetas);
			
			printf("VERTICIES\n");
			print_matrix(verticies, nverticies, nthetas);

			printf("COM IS\t");
			for(i = 0; i < nthetas; i++){ printf("%g\t", gsl_vector_get(com, i));}; printf("\n");
			
			// set xreflected
			// xreflected = com + alpha*(com - verticies[last])
			gsl_vector_memcpy(xreflected, com); 
			gsl_vector_memcpy(temp_vector, com);
			vertex = gsl_matrix_row(verticies, last_vertex); // get the last one from the list
			gsl_vector_sub(temp_vector, &vertex.vector);
			gsl_vector_scale(temp_vector, alpha);
			gsl_vector_add(xreflected, temp_vector);
			
			// now we have to evaluate the likelyhood of xreflected
			evalReflected = evalLikelyhood(xreflected, xmodel, trainingvector, nmodel_points, nthetas, nparams);
			
			// now the heuristics begin
			if(evalList[0].value <= evalReflected && evalReflected <= evalList[last_vertex -1].value){
				// so the reflected vertex is better than the worst, so we keep it
				gsl_matrix_set_row(verticies, last_vertex, xreflected); 
			} else if( evalReflected < evalList[0].value){
				// i.e evalReflected is the best we have so far
				// setup xExtended = com + gamma*(com-verticies[last])
				gsl_vector_memcpy(xextended, com);
				gsl_vector_memcpy(temp_vector, com);
				vertex = gsl_matrix_row(verticies, last_vertex);
				gsl_vector_sub(temp_vector, &vertex.vector);
				gsl_vector_scale(temp_vector, gamma);
				gsl_vector_add(xextended, temp_vector);
				
				evalExtended = evalLikelyhood(xextended, xmodel, trainingvector, nmodel_points, nthetas, nparams);
				
				if(evalExtended < evalReflected){
					gsl_matrix_set_row(verticies, last_vertex, xextended);
				} else {
					gsl_matrix_set_row(verticies, last_vertex, xreflected);
				}
			} else if(evalReflected > evalList[last_vertex-1].value){
				// setup xcontracted = last_vertex + rho*(com-last_vertex);
				vertex = gsl_matrix_row(verticies, last_vertex);
				gsl_vector_memcpy(xcontracted, &vertex.vector);
				gsl_vector_memcpy(temp_vector, com);
				gsl_vector_sub(temp_vector, &vertex.vector);
				gsl_vector_scale(temp_vector, rho);
				gsl_vector_add(xcontracted, temp_vector);
				
				evalContracted = evalLikelyhood(xcontracted, xmodel, trainingvector, nmodel_points, nthetas, nparams);
				
				if(evalContracted < evalList[last_vertex].value){
					gsl_matrix_set_row(verticies, last_vertex, xcontracted);
				} else {
					make_new_vlist(new_verticies, verticies, sigma, nverticies, nthetas);
				}
			}
				
			t++;
			
			printf("EVALS:\t");
			for(i = 0; i < nverticies; i++){ printf("%g\t", evalList[i].value);} 
			printf("\n");		 
		 
		}
		// now it's done, calc the com one more time and then run with it
		if(nanflag == 0){
			calc_com(verticies, com, nverticies, nthetas);
			answer_likelyhood = evalLikelyhood(com, xmodel, trainingvector, nmodel_points, nthetas, nparams);
			if(answer_likelyhood > best_value){
				best_value = answer_likelyhood;
				printf("best_value = %g\n", answer_likelyhood);
				gsl_vector_memcpy(best_vector, com);
			}
		} else {
			nanflag = 0;
		}
		tries++;
	}
	printf("final_answer = %g\n", answer_likelyhood);
	gsl_vector_memcpy(the_answer, best_vector);

	
	free(evalList);
	gsl_vector_free(best_vector);
	gsl_vector_free(com);
	gsl_vector_free(temp_vector);
	gsl_vector_free(xreflected);
	gsl_vector_free(xcontracted);
	gsl_vector_free(xextended);
	gsl_matrix_free(verticies);
	gsl_matrix_free(new_verticies);

}


void make_new_vlist(gsl_matrix* new_verticies, gsl_matrix* verticies, double sigma, int nverticies, int nthetas){
	int i;
	gsl_vector *temp_vector = gsl_vector_alloc(nthetas);
	gsl_vector *new_vector = gsl_vector_alloc(nthetas);
	gsl_vector *first_vertex = gsl_vector_alloc(nthetas);
	gsl_vector_view vertex;
	vertex = gsl_matrix_row(verticies, 0);
	gsl_vector_memcpy(first_vertex, &vertex.vector);
	
	// new_verticies = first_vertex + sigma*(vertex[i] - first_vertex)
	for(i = 1; i < nverticies; i++){
		vertex = gsl_matrix_row(verticies, i);
		gsl_vector_memcpy(temp_vector, &vertex.vector);
		gsl_vector_memcpy(new_vector, temp_vector);
		gsl_vector_sub(temp_vector, first_vertex);
		gsl_vector_scale(temp_vector, sigma);
		gsl_vector_add(new_vector, temp_vector);
		gsl_matrix_set_row(new_verticies, i,  new_vector);
	}
	
	gsl_vector_free(temp_vector);
	gsl_vector_free(new_vector);
	gsl_vector_free(first_vertex);
}


//! get the likelyhood for a given vertex
/**
 * returns the loglikelyhood for a given vertex
 * each time you do this you have to invert the matrix, so this method becomes expensive in high nthetas
 * 
 */
double evalLikelyhood(gsl_vector *vertex, gsl_matrix *xmodel, gsl_vector *trainingvector, int nmodel_points, int nthetas, int nparams){
	// evaluate the likelyhood at a given point
	gsl_matrix* covariance_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix* cinverse = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_matrix* temp_matrix = gsl_matrix_alloc(nmodel_points, nmodel_points);
	gsl_permutation *c_LU_permutation = gsl_permutation_alloc(nmodel_points);
	int lu_signum =0;
	double cinverse_det = 0.0;
	double the_likelyhood = 0.0;
	

	// make the covariance matrix 
	// using the random initial conditions! (xold not thetas)
	makeCovMatrix(covariance_matrix, xmodel, vertex, nmodel_points, nthetas, nparams);
	gsl_matrix_memcpy(temp_matrix, covariance_matrix);
	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
	gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse); // now we have the inverse
	// now get the determinant of the inverse			
	gsl_matrix_memcpy(temp_matrix, cinverse);
	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
	cinverse_det = gsl_linalg_LU_det(temp_matrix, lu_signum);
	
	the_likelyhood  = getLogLikelyhood(cinverse, cinverse_det, xmodel, trainingvector, vertex, nmodel_points, nthetas, nparams);
	

	gsl_matrix_free(covariance_matrix);
	gsl_matrix_free(cinverse);
	gsl_matrix_free(temp_matrix);
	gsl_permutation_free(c_LU_permutation);
	
	return(the_likelyhood);
}


//! sort the verticies by the sorted evallist
/** 
 * sort the vertex matrix using the information in the sorted evvallist,
 * if the evallist hasn't been sorted by value this is not going to be very useful
 *
 * @param verticies a matrix nverticies (cols) nparams rows
 * @param evallist a list of evalunit structs, length -> nverticies
 * @return verticies is sorted into the correct order
 */
void sort_vertex_list(gsl_matrix *verticies, evalunit* evalList, int nverticies, int nthetas){
	gsl_matrix * temp_matrix = gsl_matrix_alloc(nverticies, nthetas);
	gsl_vector_view vertex;
	int i;
	int the_index = 0;
	
	// use the indicies from the sorted evallist to 
	// copy the verticies, rowwise into the tempmatrix
	for(i = 0; i < nverticies; i++){
		the_index = evalList[i].index;
		vertex = gsl_matrix_row(verticies, the_index);
		gsl_matrix_set_row(temp_matrix, nthetas, &vertex.vector);
	}
	// copy the sorted matrix back into the vertex list
	gsl_matrix_memcpy(verticies, temp_matrix);
	gsl_matrix_free(temp_matrix);
}


//! calculate the centre of mass of the simplex
/** 
 * calculte the geometric centre of mass of the simplex
 * ignoring the last point
 * 
 * @return com is a vector set to the centre of mass
 */
void calc_com(gsl_matrix *verticies, gsl_vector *com, int nverticies, int nthetas){
	int i, j;
	double temp_val = 0.0;
	for(i = 0; i < nthetas; i++){
		for(j = 0; j < (nverticies-1); j++){
			temp_val += gsl_matrix_get(verticies, j,i);
		}
		gsl_vector_set(com, i, temp_val / ((double)nverticies));
		temp_val = 0.0;
	}
}

	
