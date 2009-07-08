
#include "multifit.h"

/* typedef struct eopts{ */
/* 	int nmodel_points; */
/* 	int nemu_points; */
/* 	int nparams; */
/* 	int nthetas; */
/* 	double range_min; */
/* 	double range_max; */
/* 	gsl_matrix *xmodel; */
/* 	gsl_vector *training; */
/* 	gsl_vector *thetas; */
/* } eopts; */

// this will only work in 1d for now
//! emulate the given region, returning the new_x, emulated_mean and emulated_variances
/**
 * emulate a given region, the eopts struct has to be set correctly and also 
 * the hyperparameters have to have been set
 * @return new_x are the emulated points, emulated_mean is the mean at these points
 * and emualted_variance is the variance at these points */
void emulate_region(gsl_matrix *new_x, gsl_vector* emulated_mean, gsl_vector* emulated_variance , eopts* options){
	int i, j;
	double kappa = 0.0;
	double step_size = (options->range_max - options->range_min) /((double)(options->nemu_points));
	double temp_mean = 0.0;
	double temp_var = 0.0;
	gsl_matrix *c_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_matrix *cinverse = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_vector *kplus = gsl_vector_alloc(options->nmodel_points);
	gsl_vector_view new_x_row;

	gsl_matrix *temp_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_permutation *c_LU_permutation = gsl_permutation_alloc(options->nmodel_points);
	int lu_signum = 0;

	makeCovMatrix(c_matrix, options->xmodel, options->thetas,options->nmodel_points, options->nthetas, options->nparams);
	gsl_matrix_memcpy(temp_matrix, c_matrix);
	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
	gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse);
	
	// set the new_x values
	for(i = 0; i < options->nemu_points; i++){
		// this doesn't make sense for many params!
		for(j = 0; j < options->nparams; j++){	
			gsl_matrix_set(new_x, i, j,step_size*((double)i)+options->range_min);			
		}
	}
	
	for(i = 0; i < options->nemu_points; i++){
		new_x_row = gsl_matrix_row(new_x, i);
		makeKVector(kplus, options->xmodel, &new_x_row.vector, options->thetas, options->nmodel_points, options->nthetas, options->nparams);
		temp_mean = makeEmulatedMean(cinverse, options->training, kplus, options->nmodel_points);
		kappa = covariance_fn(&new_x_row.vector, &new_x_row.vector, options->thetas, options->nthetas, options->nparams);
		temp_var = makeEmulatedVariance(cinverse, kplus, kappa, options->nmodel_points);
		gsl_vector_set(emulated_mean, i, temp_mean);
		gsl_vector_set(emulated_variance, i, temp_var);
	}

	gsl_matrix_free(c_matrix);
	gsl_matrix_free(cinverse);
	gsl_vector_free(kplus);
	gsl_matrix_free(temp_matrix);
	gsl_permutation_free(c_LU_permutation);
}

//! estimate the hyperparams for a region
/**
 * @param random is a gsl_rng which has already been setup and seeded
 * @return options->thetas is set to the best hyperparams
 */
void estimate_region(eopts* options, gsl_rng *random){
	int max_tries = 20;
	int i;
	int number_steps = 40;
	gsl_matrix *grad_ranges = gsl_matrix_alloc(options->nthetas,2);
  
	for(i = 0; i < options->nthetas; i++){
		gsl_matrix_set(grad_ranges, i, 0, 0.0);
		gsl_matrix_set(grad_ranges, i, 1, 1.1);
	}

	gsl_matrix_set(grad_ranges, 3, 0, 0.0);
	gsl_matrix_set(grad_ranges, 3,1, 1.1);
 		
	nelderMead(random, max_tries, number_steps, options->thetas, grad_ranges, options->xmodel, options->training, options->nmodel_points, options->nthetas, options->nparams);
	
	fprintf(stderr, "in range: %g..%g\n", options->range_min, options->range_max);
	fprintf(stderr, "best thetas: \t"); 
	print_vector_quiet(options->thetas, options->nthetas);
	
	gsl_matrix_free(grad_ranges);
}


/*void evaluate_region(gsl_matrix* new_x, gsl_vector* new_mean, gsl_vector* new_variance, eopts* options, gsl_rng * random){
	estimate_region(options, random);
	emulate_region(new_x, new_mean, new_variance, options);
	}*/

//! do the estimation and emulation for a region
void evaluate_region(emuResult *results, eopts* options, gsl_rng* random){
	estimate_region(options, random);
	emulate_region(results->new_x, results->new_mean, results->new_var, options);
}

//! see if a region is smooth
//! this is stupid at the moment
/*int is_smooth(double smooth_val, gsl_vector* xemu, gsl_vector* mean_emu, gsl_vector* var_emu, eopts* options){
	int i; 
	double mse = malloc(sizeof(double)*(options->nemu_points));
	double var = 0.0;
	for (i = 0; i < options-> nemu_points; i++){
		mse[i] = get_mse(gsl_vector_get(mean_emu, i), gsl_vector_get(var_emu, i));
	}

	var = gsl_stats_variance(mse, 1, options->nemu-points);
	free(mse);
	if(var > smooth_val){
		fprintf(stderr, "wiggly!\n");
		return(1);
	} else {
		fprintf(stderr, "smooth\n");
		return(0);
	}
	
	}*/


double get_mse( double mean, double variance){
	return(pow(fabs(mean-variance), 2.0));
}

/// 
/// All these functions only work in 1 parameter so far
///


//! calculate the log(1/dy^2)
/**
 * @param inverse_error is allocated to res->nemu_points long
 */
void make_goodness(emuResult *res, double* the_goodness){
	int i;
	for (i = 0; i < res->nemu_points; i++){
		the_goodness[i] = (1.0)/(pow(gsl_vector_get(res->new_var, i), 2.0));
		the_goodness[i] = log(the_goodness[i]);
	}
}

//! calculate the local diffs in the goodness
void make_diff(emuResult *res, double* the_goodness, double* diff_goodness){
	int i;
	for(i =0; i < res->nemu_points; i++){
		if(i > 0)
			diff_goodness[i] = fabs(the_goodness[i] - the_goodness[i-1]);
	}
}

//! compare successive differences, if they are small enough put a 1 in the cluster vector
void prepare_cluster(emuResult *res, double* diff_goodness, int* cluster, double diff_threshold){
	int i;
	for(i =0; i < res->nemu_points; i++){
		// this is a simple forwards derivative, 
		// for higher dims it'll make sense to go for an arbitrary star-type finite diffs methinks
		// i.e find the gradient at this point
		if(fabs(diff_goodness[i+1] - diff_goodness[i]) < diff_threshold)
			cluster[i] = 1;
	}
}


//! for a linked list of regions (see list-test.c)
// decided to just use an array and make it bigger if i need to
/* typedef struct region{ */
/* 	int region_start; */
/* 	int region_stop; */
/* 	int region_length; */
/* 	double emu_x_start; */
/* 	double emu_x_stop; */
/* } region; */

//! print out what's going on
void checkup(emuResult *res, double* goodness, double* diff_goodness, int*cluster){
	int i;
	for(i =0; i < res->nemu_points; i++){
		printf("%d:%g\t%g\t%g\t%g\t%g\t%d\n", i,gsl_matrix_get(res->new_x, i,0), gsl_vector_get(res->new_mean, i), gsl_vector_get(res->new_var, i), goodness[i], diff_goodness[i], cluster[i]);
	}
}


//! resizes an array of region structs from current_length to current_length+ grow_length
/**
 * resizes an array from current_length to current_length+grow_length
 * 
 * @param the_array the thing to be resized
 * @param current_length how long it is now
 * @param grow_length how much to grow it by
 * @return the new length (current_length+grow_length)
 */
int resize_region_array(region* the_array, int current_length, int grow_length){
	size_t current_size = sizeof(region)*current_length;
	size_t new_size = current_size + sizeof(region)*grow_length; // the grown array
	region* buffer = MallocChecked(current_size);
	copy_region_array(buffer, the_array, current_length);
	free(the_array);
	the_array = MallocChecked(new_size);
	copy_region_array(the_array, buffer, current_size);
	free(buffer);
	// and return the new length
	return(current_length+grow_length);
}


//! assign clusters
/**
 * Run through the cluster list and create a "region" struct for each cluster which is more than cluster_min
 * points long. These regions can be processed elsewhere to make the required cuts.
 * 
 * @param cluster_min is the minimum number of successive points to be a cluster
 * @param region_list is a pointer which will be assigned to an array of region structs, this should not be assigned already
 * @param nclusters is set on return to  how long this array is
`* @return region_list -> filled with the regions, nclusters -> the number of clusters we found
 */
void assign_clusters(emuResult *res, int *cluster, int cluster_min, region** region_array ,int* nclusters ){
	int i;
	int temp_cluster_count = 0;
	int cluster_begin = 0;
	int cluster_end = 0;
	int total_cluster_count = 0;
	int total_cluster_lengths = 0;

	// this was a humerously hard variable name to spell
	int rarray_length = 10;
	int rarray_size = rarray_length*sizeof(region);
	int rarray_grow_offset = 10;
	region temp_region; 
	region* temp_region_array = MallocChecked(rarray_size);



	for(i = 0; i < res->nemu_points -1; i++){
		if(cluster[i] == 1){
			if(cluster[i+1] == 1){
				temp_cluster_count++;
				if( i >0 && cluster[i-1] == 0){
					cluster_begin = i;
				}
			} else if(cluster[i+1] == 0){
				if(temp_cluster_count > cluster_min){
					// this is where you do "found a cluster, things"
					cluster_end = i;

					
					// ugly put this somewhere else
					if(total_cluster_count > rarray_length){
						fprintf(stderr, "reallocating rarray");
						// resize
						rarray_length = resize_region_array(temp_region_array, rarray_length, rarray_grow_offset);
					}

					// push temp_region onto the array
					temp_region.region_start = cluster_begin;
					temp_region.region_stop = cluster_end;
					temp_region.region_length = cluster_end - cluster_begin;
					temp_region_array[total_cluster_count] = temp_region;

					total_cluster_lengths += (cluster_end - cluster_begin);
					printf("cluster %d %d len=%d\n", cluster_begin, cluster_end, (cluster_end-cluster_begin));
					temp_cluster_count = 0;
					total_cluster_count++;
				}
				else {
					temp_cluster_count = 0;
					cluster_begin = 0;
					cluster_end =0;
				}
			}
		}
	}
	
	// fill in the other fields in the structure
	for(i = 0; i < total_cluster_count; i++){
		temp_region_array[i].emu_x_start = gsl_matrix_get(res->new_x, temp_region_array[i].region_start, 0);
		temp_region_array[i].emu_x_stop = gsl_matrix_get(res->new_x, temp_region_array[i].region_stop, 0);
	}


	//finish up
	*nclusters = total_cluster_count;
	*region_array = MallocChecked(sizeof(region)*total_cluster_count);
	copy_region_array(*region_array, temp_region_array, total_cluster_count);
	fprintf(stderr, "found %d clusters\n", total_cluster_count);
	fprintf(stderr, "total length %d\n", total_cluster_lengths);
	free(temp_region_array);
}
				
void copy_region_array(region* target, region* source, int length){
	int i;
	for (i = 0; i < length; i++){
		target[i].region_start = source[i].region_start;
		target[i].region_stop = source[i].region_stop;
		target[i].region_length = source[i].region_length;
		target[i].emu_x_start = source[i].emu_x_start;
		target[i].emu_x_stop = source[i].emu_x_stop;
	}	
}


//! take an emulator result and come up with a region list of clusters.
void create_clusters_1d(emuResult *res, region* region_list){
	int i;
	double *diff_goodness;
	double *the_goodness;
	int *cluster;
	double diff_thresh = 1.0;
	int cluster_min =5;
	int npoints = res->nemu_points;
	int number_clusters;
	region* local_region_list;
	region* rtemp;


	the_goodness = MallocChecked(sizeof(double)*npoints);
	diff_goodness = MallocChecked(sizeof(double)*npoints);
	cluster = MallocChecked(sizeof(int)*npoints);
	
	for(i= 0; i < npoints; i++){
		diff_goodness[i] = 0.0;
		cluster[i] = 0;
	}
	
	make_goodness(res, the_goodness); // create vec of "goodness"
	make_diff(res, the_goodness, diff_goodness); // take a crappy derivative of the goodness
	prepare_cluster(res, diff_goodness, cluster, diff_thresh); // lump regions which change slowly enough together

	checkup(res, the_goodness, diff_goodness, cluster); // print out what you've got so far
	assign_clusters(res, cluster, cluster_min, &local_region_list, &number_clusters);
	
	

	
	for(i = 0; i < number_clusters; i++){
		rtemp = &(local_region_list[i]);
		printf("%d\t%d\t%d\t%g\t%g\n", rtemp->region_start, rtemp->region_stop, rtemp->region_length, rtemp->emu_x_start, rtemp->emu_x_stop);
	}
	
	// blind to the very horror of this sorry life
	free(local_region_list);
}
