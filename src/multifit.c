
#include "multifit.h"
#include "estimate_threaded.h"



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
	int i;
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
	
	/* copied in the version from main.c which tries to sensibly allocate matricies up to 3d
	 * all typed in by hand! */
	/* // set the new_x values */
	/* for(i = 0; i < options->nemu_points; i++){ */
	/* 	// this doesn't make sense for many params! */
	/* 	for(j = 0; j < options->nparams; j++){	 */
	/* 		gsl_matrix_set(new_x, i, j,step_size*((double)i)+options->range_min);			 */
	/* 	} */
	/* } */
	initialise_new_x(new_x, options->nparams, options->nemu_points, options->range_min, options->range_max);
	
	
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


//! set the restrictions on the ranges allowed for the multivariate fitter
/** 
 * set the ranges allowed for the parameters in the covariance function
 * this should be the only point that these are messed with! 
 * \bold IMPORTANT
 * the matern and gaussian covariance functions use different
 * theta's so that which is the nugget will change when you switch
 * between them, if you don't spot this you'll have some seriously
 * weird results when you try and run one with the others parameters.
 */
void set_likelyhood_ranges(gsl_matrix* ranges_matrix, int nthetas){
	int i;
	// for either the gaussian or matern covfns this 
	// [0.1] range seems to work ok
	for(i=0;i<nthetas;i++){
		gsl_matrix_set(ranges_matrix, i, 0, 0.0);
		gsl_matrix_set(ranges_matrix, i, 1, 1.0);
	}

	// however we have to be VERY careful with the nugget, 
	// if we're using the matern then
	/* gsl_matrix_set(ranges_matrix, 3, 0, 0.0); */
	/* gsl_matrix_set(ranges_matrix, 3, 1, 0.1); */
	// if we're using the gaussian:
	gsl_matrix_set(ranges_matrix, 2, 0, 0.0);
	gsl_matrix_set(ranges_matrix, 2, 1, 0.2);

	/*
	 * this is ugly, i'm sorry about it 
	 */
}


//! estimate the hyperparams for a region
/**
 * @param random is a gsl_rng which has already been setup and seeded
 * @return options->thetas is set to the best hyperparams
 */
void estimate_region(eopts* options, gsl_rng *random){ /*  */
	int max_tries = 20;

	int number_steps = 40;
	gsl_matrix *grad_ranges = gsl_matrix_alloc(options->nthetas,2);

	// note you need to play with this if you change the emulator cov fn
	set_likelyhood_ranges(grad_ranges, options->nthetas);

	nelderMead(random, max_tries, number_steps, options->thetas, grad_ranges, options->xmodel, options->training, options->nmodel_points, options->nthetas, options->nparams);
	
	fprintf(stderr, "in range: %g..%g\n", options->range_min, options->range_max);
	fprintf(stderr, "best thetas: \t"); 
	print_vector_quiet(options->thetas, options->nthetas);
	
	gsl_matrix_free(grad_ranges);
}


/** call out to estimate_thetas_threaded 
 * and push the answer back into options->thetas 
 */
void estimate_region_threaded(eopts* options){
	
	optstruct newopts;	
	copy_eopts_to_optstruct(&newopts, options);
	estimate_thetas_threaded(options->xmodel, options->training, options->thetas, &newopts);
	// now the answer is back in options->thetas which is what we want
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

//! do the estimation and emulation for a region
void evaluate_region_threaded(emuResult *results, eopts* options, gsl_rng* random){
	estimate_region_threaded(options);
	emulate_region(results->new_x, results->new_mean, results->new_var, options);
}


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
/**
 * Takes a forward derivative on the vector the_goodness and returns the absolute value as diff_goodness. 
 * The emulator result is needed only to provide the number of points.
 *
 * \note for many dimensions this should be adjusted to something like 
 * (x(i) - x(i+1)) + (x(i) - x(i-1))/ 2 for each dimension, a grad basically
 */
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
	for(i =0; i < res->nemu_points-1; i++){
		// this is a simple forwards derivative, 
		// for higher dims it'll make sense to go for an arbitrary star-type finite diffs methinks
		// i.e find the gradient at this point
		if(fabs(diff_goodness[i+1] - diff_goodness[i]) < diff_threshold)
			cluster[i] = 1;
	}
}

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
int resize_region_array(region** the_array, int current_length, int grow_length){
	size_t current_size = current_length*sizeof(region);
	size_t new_size = current_size + sizeof(region)*grow_length; // the grown array
	region* buffer = MallocChecked(current_size);
	copy_region_array(buffer, *the_array, current_length);
	free(*the_array);
	*the_array = MallocChecked(new_size);
	copy_region_array(*the_array, buffer, current_length);
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
	int rarray_length = 1;
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
					if(total_cluster_count + 1 > rarray_length){
						fprintf(stderr, "reallocating rarray");
						// resize
						rarray_length = resize_region_array(&temp_region_array, rarray_length, rarray_grow_offset);
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
				

//! copy source->target
/**
 * copies the contents of the source region into the target region, useful because i can't make memcpy just take a 
 * block of these guys and do it for me.
 */
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
/**
 * Find subregions of a result which are relatively similar and smooth and then group them together,
 * if they are large enough they will be pushed into region_list.
 * 
 * Regions are found by: 
 * - Allocating a goodness to all of the points in the result via make_goodness
 * - Taking a simple forwards derivative of this goodness with make_diff
 * - Lumping together parts of the resulting list of derivatives which change slowly enough to be considered similar. 
 * - assign_clusters then takes the prepared list of regions checks if the clusters are long enoguh and if so adds them to the region array
 * 
 * The resulting array of "good" regions from assign_clusters is then copied into the region_list and returned.
 *
 * @param a result to be processed
 * @param region_list: will be filled with a list of regions which represent good clusters, That is clusters which have at least the minimum number of points within them. Given by cluster_min.
 * @param number_regions: set to the final number of regions found.
 */
void create_clusters_1d(emuResult *res, region** region_list, int* number_regions){
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

	//DEBUG 
	checkup(res, the_goodness, diff_goodness, cluster); // print out what you've got so far
	assign_clusters(res, cluster, cluster_min, &local_region_list, &number_clusters);
	
	for(i = 0; i < number_clusters; i++){
		rtemp = &(local_region_list[i]);
		printf("%d\t%d\t%d\t%g\t%g\n", rtemp->region_start, rtemp->region_stop, rtemp->region_length, rtemp->emu_x_start, rtemp->emu_x_stop);
	}
	
	*region_list = MallocChecked(sizeof(region)*number_clusters);
	*number_regions = number_clusters;
	copy_region_array(*region_list, local_region_list, number_clusters);
	free(local_region_list);
	free(the_goodness);
	free(diff_goodness);
	free(cluster);
}

//! take a region and find the model points which are contained within that region. 1d only
/** 
 * Iterate through a given region and find what subset of the model-points it represents. 
 * Quite primative and only works in 1d,need to use something like a voroni hull to make this work 
 * In higher dimensions.
 *
 * looking at the output of this, it seems like we can get a span of -1, 
 * i guess this means that the high and low index are essentially ontop of each other
 * a good way to reject a split methinks
 */
void assign_model_point(eopts* regionOpts, region* the_region){
	double region_min = the_region->emu_x_start;
	double region_max = the_region->emu_x_stop;
	double temp_val;
	int low_index = 0;
	int high_index = regionOpts->nmodel_points-1;

	temp_val = gsl_matrix_get(regionOpts->xmodel, low_index, 0);
	if(temp_val != region_min){
		while(temp_val < region_min){
			low_index++;
			temp_val = gsl_matrix_get(regionOpts->xmodel, low_index, 0);
		}
	}

	temp_val = gsl_matrix_get(regionOpts->xmodel, high_index, 0);
	if(temp_val != region_max){
		while(temp_val > region_max){
			high_index--;
			temp_val = gsl_matrix_get(regionOpts->xmodel, high_index, 0);
		}
	}
		
	// set the values 
	the_region->model_x_start = gsl_matrix_get(regionOpts->xmodel, low_index, 0);
	the_region->model_x_stop = gsl_matrix_get(regionOpts->xmodel, high_index, 0);
	the_region->model_x_span = high_index - low_index;
	fprintf(stderr, "low_index = %d\n", low_index);
	fprintf(stderr, "high_index = %d\n", high_index);
	fprintf(stderr, "span = %d\n", high_index - low_index);
}


//! take a set of split locations, make regions, evaluate them and push dump the results,
//! this is the final routine you'd want to call after doing a bunch of splits
void process_splits(gsl_matrix* the_splits, int nsplits, eopts* the_options, gsl_rng* random_number){
	int i;
	FILE *fptr;
	emuResult* results_array = MallocChecked(sizeof(emuResult)*nsplits);
	eopts* options_array = MallocChecked(sizeof(eopts)*nsplits);
	char buffer[256];
	char outfolder[128];
	
	// this should presumably be an option
	sprintf(outfolder, "%s", "output");
	

	// make the new options
	for(i = 0; i < nsplits;i++){
		split_region_options(&options_array[i], the_options, gsl_matrix_get(the_splits, i, 0), gsl_matrix_get(the_splits, i, 1));
		if(&(options_array[i])==NULL){
			fprintf(stderr, "bad split\n");
			exit(1);
		}
		alloc_emuRes(&results_array[i], &options_array[i]);
	}

	// run the options
	for(i = 0; i < nsplits; i++){
		evaluate_region(&results_array[i], &options_array[i], random_number);
	}
		
	// this is where you'd dump the results to file or whatever
	// infaaaact we can do that with persist
/* 	fptr = fopen("output.dat", "w"); */
/* 	for( i = 0; i < nsplits; i++){ */
/* 		dump_eopts(&options_array[i], fptr); */
/* 		dump_emuresult(&results_array[i], fptr); */
/* 	} */

	// loop over the splits and output them to the folder given by outfolder
	// note that we're only dumping the results, not the options
	// a bin dump of everything can  be turned on above
	for(i = 0; i < nsplits; i++){
		sprintf(buffer, "%s/split-region-%d.dat", outfolder, i);
		fptr = fopen(buffer, "w");
		if(fptr == NULL){
			fprintf(stderr, "could not open: %s\n", buffer);
			exit(1);
		}
		dump_result(&results_array[i], fptr);
		fclose(fptr);
	}


	for(i = 0; i < nsplits;i++){
		free_emuRes(&results_array[i]);
		free_eopts(&options_array[i]);
	}
	
	free(results_array);
	free(options_array);
}



//!  transforms local_split_ranges, the set of ranges with gaps into split_ranges the continuous set of ranges
/**
 * Takes a set of split ranges: local_split_ranges which have been pre-caculated and fills the gaps between them so that
 * the result: split_ranges is a set of continuous points which span the whole range of toplevel.
 * 
 * a bad diagram:
 * local_split_ranges:  start............|region 1 ..|...............|region 2..|......stop
 * split_ranges:        start|region 0...|region 1...|region 2.......|region 3..|region 4 stop
 * 
 * Ther are no gaps by the end.
 */
void fill_split_ranges(gsl_matrix* split_ranges, int ngoodregions, gsl_matrix * local_split_ranges, eopts* toplevel){
	int i;
	int split_count =0;

	// start from 0 
	for(i = 0; i < ngoodregions; i++){
		if(i == 0){
			//printf("first!\n");
			gsl_matrix_set(split_ranges, 0,0, toplevel->range_min);
			if(gsl_matrix_get(local_split_ranges, 0, 0) == toplevel->range_min){				
				// there is no gap
				//printf("no gap\n");
				gsl_matrix_set(split_ranges,0, 1, gsl_matrix_get(local_split_ranges, 0, 1));
				split_count++;
			} else {
				// there is a gap
				//printf("initial gap\n");
				gsl_matrix_set(split_ranges, 0,1, gsl_matrix_get(local_split_ranges, 0, 0));
				split_count++;
			}
		}
		//printf("middling: ");
			// now add the winner
			gsl_matrix_set(split_ranges, split_count, 0, gsl_matrix_get(local_split_ranges, i, 0));
			gsl_matrix_set(split_ranges, split_count, 1, gsl_matrix_get(local_split_ranges, i, 1));
			split_count++;
			// now we have to add a blank again
			if((i < ngoodregions -1) && (gsl_matrix_get(local_split_ranges, i+1, 0) > gsl_matrix_get(local_split_ranges, i, 1))){
				// there is another gap!
				gsl_matrix_set(split_ranges, split_count+1,  0, gsl_matrix_get(local_split_ranges, i, 1));
				gsl_matrix_set(split_ranges, split_count+1,  1, gsl_matrix_get(local_split_ranges, i+1, 1));
				split_count++;
			}
		
			if( i == (ngoodregions -1)){
				// the last one
				//printf("last!\n");
				if(gsl_matrix_get(local_split_ranges, i, 1) < toplevel->range_max){
					// there is a gap before the end
					//printf("final one\n");
					gsl_matrix_set(split_ranges, split_count, 0, gsl_matrix_get(local_split_ranges, i, 1));
					gsl_matrix_set(split_ranges, split_count, 1, toplevel->range_max);
					split_count++;
				}
				// do nothing otherwise
			}
			printf("split_count = %d\n", split_count);
	}
}

//! creates a new region from lower and upper in the first parameter of the xmodel parts of the parent
/**
 * Creates a new eopts result from the parent but only includes the data which lies within the range  
 * lower..upper. 
 * 
 * This method can then be used on the same parent multiple times with different lower and upper values
 * to carve up the parent into a series of child regions which can then be evaluated etc. 
 * 
 * if the split doesn't work, result => NULL
 */
void split_region_options(eopts *result, eopts *parent, double lower, double upper){
	assert(lower < upper);
	
	int i,j;
	int split_low = 0;// the index at which to split the model, training vecs etc
	int split_high=  parent->nmodel_points;
	int new_nemu_points = 0;
	int new_nmodel_points = 0;
	//int min_model_points = 3; // don't need this, already set by the clustering alg
	int offset = 0;
	int bad_flag = 0;
	double temp_val;
	
	// set the basic things
	new_nemu_points = 40; //(fudge)
	
	// grow split low to the last_xmodel value before the split
	for(i = 0; i < parent->nmodel_points; i++){
		// look at the first value in xmodel only
		if(gsl_matrix_get(parent->xmodel, i, 0) <= lower){
			split_low++; // don't split at this i
		} 
	}
	

	// shrink split high to the xmodel value to the right of the upper split
	for(i = parent->nmodel_points-1; i >= 0; i--){
		if(gsl_matrix_get(parent->xmodel, i, 0) > upper){
			split_high--;
		}
	}

	fprintf(stderr, "split_low = %d, split_high = %d\n", split_low, split_high);

	// now we can figure everything else out
	if(split_low == split_high){
		fprintf(stderr, "bad split\n");
		bad_flag = 1;
	}

	if(lower < parent->range_min || upper > parent->range_max){
		fprintf(stderr, "doing a split out of range\n");
		bad_flag = 1;
	}
	
	new_nmodel_points = (split_high - split_low);

	/*if(new_nmodel_points < min_model_points){
		fprintf(stderr, "bad split, not enough new points\n");
		//exit(1);
		bad_flag = 1;
		}*/
	
	if(bad_flag != 1){
		// set everything up
		result->nmodel_points = new_nmodel_points;
		result->nparams = parent->nparams;
		result->nthetas = parent->nthetas;
		result->range_min = lower;
		result->range_max = upper;
		result->nemu_points = new_nemu_points;
		result->xmodel = gsl_matrix_alloc(new_nmodel_points, parent->nparams);
		result->training = gsl_vector_alloc(new_nmodel_points);
		result->thetas = gsl_vector_alloc(parent->nthetas); 
		
		// copyin the new data
		for(i = 0; i < new_nmodel_points; i++){
			offset = i + split_low;
			for(j = 0; j < parent->nparams;  j++){
				temp_val = gsl_matrix_get(parent->xmodel, offset, j);
				gsl_matrix_set(result->xmodel, i, j, temp_val);
			}
			temp_val = gsl_vector_get(parent->training, offset);
			gsl_vector_set(result->training, i, temp_val);
		}
		
	} else if(bad_flag == 1){
		result = NULL;
		fprintf(stderr, "bad splits\n");
		exit(1);
	}

	
}


