#include "smasher.h"


//! just to qsort the regions in smasher
// want descending order, the first one is longest
int compare_regions(const region* a, const region* b){
	double temp = a->region_length - b->region_length;
	if(temp >0){
		return -1;
	} else if(temp < 0){
		return 1;
	} else 
		return 0;
}

//! runs the whole splitting thing
/**
 * @return 0 if failed, >0  (number of regions otherwise)
 * 
 * Takes toplevel and attempts to split it up into regions which are "good" i.e the emulator works well here. 
 * Then the rest of the interval is filled in with blank regions by fill_split_ranges so that on return split_ranges can be pushed
 * straight into process splits, and each subrange can be evaluated.
 * 
 * There is very little flow control right now,if things go wrong errors may spring up but the whole process may or may not terminate, 
 * this is stupid. I hope to make another driving routine ontop of this to automate exploration of various alternative splits and such.
 *
 * @split_ranges -> a set of ranges which fully span the range of toplevel
 * @nsplits -> how long split_ranges is
 * @toplevel -> specifies the model and how it is to be cut up.
 * @max_depth -> unused yet
 * @min_points -> sets the minimum number of model points for a region to be accepted
 * @random_number -> a gsl_rng used to evaluate regions
 *
 */
int smasher(gsl_matrix **split_ranges, int* nsplits, eopts* toplevel, int min_points, gsl_rng* random_number){
	int i;
	emuResult temp_result;
	region* region_list;
	int nregions;
	int ngoodregions = 0;
	int number_final_splits = 0;
	int retval = 0;
	int min_model_points = min_points;
	int init_split_ranges_length = 10;
	int split_ranges_length = init_split_ranges_length;
	gsl_matrix *local_split_ranges; // we'll do the whole store and grow and then finish the final set
	double toplevel_length = toplevel->range_max - toplevel->range_min;
	
	assert(min_model_points >0);

	local_split_ranges = gsl_matrix_alloc(split_ranges_length, 2);
	// eval the toplevel
	alloc_emuRes( &temp_result, toplevel);
	evaluate_region(&temp_result, toplevel,  random_number);
	create_clusters_1d(&temp_result, &region_list, &nregions);
	
	if(nregions == 0){
		fprintf(stderr, "didn't get any useful regions, fit again\n");
		retval = 0;
		return(retval);
	}

	// check to see if any of the regions span almost the entire toplevel, if so we should not split
	for(i = 0; i < nregions; i++){
		if ((region_list[i].emu_x_stop - region_list[i].emu_x_start)  > toplevel_length  * (5.0/6.0)){
			fprintf(stderr, "region %d is 5/6 of whole thing, no point dividing\n", i);
			retval = 1;
			return(retval); // only 1 region, stop
		}
	}

	
	assert(nregions > 0);
	// sort the regions on length, descending order
	// actually don't this messes up things afterwards
	//qsort(region_list, nregions, sizeof(region), (void*)compare_regions);

	for(i = 0; i < nregions;i++){
		assign_model_point(toplevel, &(region_list[i]));
		if(region_list[i].model_x_span > min_model_points){
			fprintf(stderr, "winner! Span: %d\t%g\t%g\n", region_list[i].model_x_span, region_list[i].model_x_start, region_list[i].model_x_stop);
			if(ngoodregions > split_ranges_length){
				fprintf(stderr, "too many split_ranges, need to realloc\n");
				exit(1);
			}
			// push to split_ranges
			gsl_matrix_set(local_split_ranges, ngoodregions, 0, region_list[i].model_x_start);
			gsl_matrix_set(local_split_ranges, ngoodregions, 1, region_list[i].model_x_stop);
			ngoodregions++;
		}

	}

	number_final_splits = ngoodregions; // it's at least as many as were "winners"
	
	// now we have set the local_split_ranges correctly, we need to allocate the final split_ranges and then
	// fill in the "gaps". The local_split_ranges contains only the winning splits, one which seem to be good
	// we need to also fill in the rest of the range so we can process the argument split_ranges with split_region_options
	// and use it to create the final set of results

	// first calculate howmany splits we'll need in the end
	// is the first local_split_ranges entry at the start of the total range?
	if(gsl_matrix_get(local_split_ranges, 0, 0) != toplevel->range_min){
		number_final_splits++;
	}
	// does the final winner run to the end of the total range?
	if(gsl_matrix_get(local_split_ranges, ngoodregions, 1) != toplevel->range_max){
		number_final_splits++;
	}

	// now figure out the rest
	for(i = 0; i < ngoodregions-1; i++){
		// check if adjacent ranges overlap
		if(gsl_matrix_get(local_split_ranges, i, 1) < gsl_matrix_get(local_split_ranges, i+1, 0)){
			number_final_splits++;
		}
	}
		 
	printf("number_final_splits %d\n", number_final_splits);
	*split_ranges = gsl_matrix_alloc(number_final_splits, 2);
	
	fill_split_ranges(*split_ranges, ngoodregions, local_split_ranges, toplevel);

	#ifdef EXTENDBYONE
	extend_regions(*split_ranges, number_final_splits, toplevel, 0);
	#endif

	*nsplits = number_final_splits;
	gsl_matrix_free(local_split_ranges);
	free(region_list);
	free_emuRes(&temp_result);
	return(nregions);
}

//! calculate the average model point spacing (in the dim'th param) and then extend all calc'd regions by this much
//! in either direction (if possible
void extend_regions(gsl_matrix* split_ranges, int number_splits, eopts* toplevel, int dim){
	int i;

	double average_spacing = 0.0;
	double extended_left = 0.0;
	double extended_right = 0.0;

	for(i = 0; i < toplevel->nmodel_points-1; i++)
		average_spacing += (gsl_matrix_get(toplevel->xmodel, i+1, dim) - gsl_matrix_get(toplevel->xmodel, i, dim));
	
	average_spacing /= toplevel->nmodel_points;

	for(i = 0; i < number_splits;i++){
		extended_left = gsl_matrix_get(split_ranges, i, 0) - average_spacing;
		extended_right = gsl_matrix_get(split_ranges, i, 1) + average_spacing;
		if( extended_left > toplevel->range_min)
			gsl_matrix_set(split_ranges, i, 0, extended_left);
		if( extended_right < toplevel->range_max)
			gsl_matrix_set(split_ranges, i, 1, extended_right);
	}
			
		 		
}

	
