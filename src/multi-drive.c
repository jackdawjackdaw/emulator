#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "emulator.h"
#include "estimator.h"
#include "maximise.h"
#include "multifit.h"
#include "sys/time.h"
#include "persist.h"
#include "ioread.h"
#include "useful.h"
#include "gsl/gsl_statistics.h"

// if this is defined, calculated regions are extended by one
// model point in each direction (if that's even possible) 
// this is to give a bit more buffering
#define EXTENDBYONE


void process_input_data(char** input_data, eopts* the_options);
void free_eopts(eopts* options);
unsigned long int get_seed(void);
void dump_result(emuResult *res, FILE *fptr);
void alloc_emuRes(emuResult *thing, eopts *options);
void free_eopts(eopts* options);
void free_emuRes(emuResult *thing);
void split_region_options(eopts *result, eopts *parent, double lower, double upper);
void read_input_from_file(char* filename, eopts* options);
double score_region(emuResult *res);
int smasher(gsl_matrix **split_ranges, int* nsplits, eopts* toplevel, int min_points, gsl_rng* random_number);
void print_splits(gsl_matrix* splits, int n);
void fill_split_ranges(gsl_matrix* split_ranges, int ngoodregions, gsl_matrix * local_split_ranges, eopts* toplevel);
void process_splits(gsl_matrix* the_splits, int nsplits, eopts* the_options, gsl_rng* random_number);
// expands regions by about one model point
void extend_regions(gsl_matrix* split_ranges, int number_splits, eopts* toplevel, int dim);
/* typedef struct emuResult{ */
/* 	int nemu_points; */
/* 	int nparams; */
/* 	gsl_matrix* new_x; */
/* 	gsl_vector* new_mean; */
/* 	gsl_vector* new_var; */
/* } emuResult; */

/** 
 * @file
 * @author Chris Coleman-Smith cec24@phy.duke.edu
 * @version 0.1
 * @section DESCRIPTION
 * 
 * This program attempts to fit a region by splitting that region into a set of subregions  
 * where the splitting is determined by selecting regions which are well fitted and then splitting 
 * around their ends. 
 * I.e 
 * start.......|smoooth...|......|smoooooth...|.......end
 * 
 * Then all of the resulting regions are re-fitted and dumped to a binary file. 
 * 
 * Note that only one level of subdivision takes place and there is almost no sensible flow control, so if something 
 * doesn't work then you might end up with crappy results and have to run again. 
 * In particular the parameter estimation is quite non-deterministic, the best parameters are taken from many runs of the neldermead
 * algorithm for random starting points. 
 */


int main (void){
	char inputfile[128];
	eopts the_options;
	char** input_data;
	int number_lines;
	int nregions;
	gsl_rng *random_number;
	const gsl_rng_type *T;
	FILE *fptr;

	// used if there is only to be one MEGA region
	emuResult wholeThing;
	

	int nsplits = 0;
	gsl_matrix *split_result;

	// hand pick the input
	//sprintf(inputfile, "../model-cut.dat");
	//sprintf(inputfile, "tests/reflected.dat");
	sprintf(inputfile, "STDIN");

	T = gsl_rng_default;
	random_number = gsl_rng_alloc(T);
	gsl_rng_set(random_number, get_seed());

	// manual options
	// have to set these, can't be bothered to unpack the input data too much
	// these are defaults
	the_options.nparams = 1;
	the_options.nthetas = 4;
	the_options.nemu_points = 100; 
	// these are now set by process_input_data
	//the_options.range_min = 0.0;
	//the_options.range_max = 4.0;

	// we kind of have to hope that the inputfile is sorted, could sort it...
	// have to free input data after this
	input_data =  unconstrained_read(inputfile, &number_lines);

	fprintf(stderr, "read in %d lines\n", number_lines);
	
	// could use some regex lib to cut out bad lines
	
	the_options.nmodel_points = number_lines;



	process_input_data(input_data, &the_options);

	print_matrix(the_options.xmodel, number_lines, 1);
	vector_print(the_options.training, number_lines); 
	
	//exit(1);


	nregions = smasher(&temp_result, &nsplits, &the_options, 1, random_number);

	if(nregions == 1){
		// only one region,don't need to split
		printf("only one region\n");
		fptr = fopen("output/wholeThing.txt", "w");
		alloc_emuRes(&wholeThing, &the_options);
		evaluate_region(&wholeThing, &the_options, random_number);
		dump_result(&wholeThing, fptr);		
		fclose(fptr);
	} else {
		// lots of regions
		printf("found %d region(s)\n", nregions);
		printf("made %d splits\n", nsplits); 
		print_splits(split_result, nsplits);
		process_splits(split_result, nsplits, &the_options, random_number);
		gsl_matrix_free(split_result);
	}


	gsl_rng_free(random_number);
	free_eopts(&the_options);
	free_emuRes(&wholeThing);	
	free_char_array(input_data, number_lines);
	return(0);
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
	
	
//! debugging
void print_splits(gsl_matrix* splits, int n){
	int i =0;
	for(i = 0; i < n; i++){
		printf("%d (%g..%g)\n", i, gsl_matrix_get(splits, i, 0), gsl_matrix_get(splits, i, 1));
	}
}

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

void dump_result(emuResult *res, FILE *fptr){
	int i;
	double goodness = 0.0;
	for(i = 0; i < res->nemu_points; i++){
		fprintf(fptr, "%g\t%g\t%g\t", gsl_matrix_get(res->new_x, i, 0), gsl_vector_get(res->new_mean, i), gsl_vector_get( res->new_var, i));
		goodness = 1 / pow(gsl_vector_get(res->new_var, i), 2.0);
		fprintf(fptr, "%g\n", goodness);
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
	

 
//! setup an emuResult struct from the given options
/**
 * copy the right parts of the options into the result and alloc the data structures inside.
 * after this the emuResult has to be freed with free_emuRes or you'll leak the arrays inside. 
 * 
 * Also note that the gsl allocs are not checked.
 */ 
void alloc_emuRes(emuResult *thing, eopts *options){
	int n = options->nemu_points;
	int np = options->nparams;
	thing->nemu_points = n;
	thing->new_x = gsl_matrix_alloc(n, np);
	thing->new_mean = gsl_vector_alloc(n);
	thing->new_var = gsl_vector_alloc(n);
	thing->nparams = np;
}

//! frees an emures struct
void free_emuRes(emuResult *thing){
	gsl_matrix_free(thing->new_x);
	gsl_vector_free(thing->new_mean);
	gsl_vector_free(thing->new_var);
}

//! frees an eopts struct
void free_eopts(eopts* options){
	gsl_matrix_free(options->xmodel);
	gsl_vector_free(options->training);
	gsl_vector_free(options->thetas);
}




//! splits up the input data 
/**
 * fixed to work with multiple params using strtok
 * input data is split on tabs and spaces
 * 
 * figures out the min and max of the input data and uses this to set the 
 * range_min and range_max fields in the_options
 * 
 * splits up the raw char input data into real things (xmodel etc)
 * this requires that the following fields in the_options be set correctly
 * nparams, nthetas, nmodel_points
 */
void process_input_data(char** input_data, eopts* the_options){
	int i,j;
	// set these off the first index
	double range_min = 0.0;
	double range_max = 0.0;
	char* split_string;
	double temp_value;
	
	assert(the_options->nmodel_points > 0); 
	// first allocate the buffers in options
	the_options->xmodel = gsl_matrix_alloc(the_options->nmodel_points, the_options->nparams);
	the_options->training = gsl_vector_alloc(the_options->nmodel_points);
	the_options->thetas = gsl_vector_alloc(the_options->nthetas);

	for(i = 0; i < the_options->nmodel_points; i++){
		split_string = strtok(input_data[i], "\t ");
		for( j = 0; j < the_options->nparams; j++){
			assert(split_string != NULL);
			sscanf(split_string, "%lg", &temp_value); 		
			gsl_matrix_set(the_options->xmodel, i, j, temp_value);
			split_string = strtok(NULL, "\t ");
			if(i == 0 && j == 0){
				range_max = temp_value;
				range_min = temp_value;
			}
			if(j ==0){ // only use the zero index param to set the max-min ranges
				if(temp_value > range_max){
					range_max = temp_value;
				}
				if(temp_value < range_min){
					range_min = temp_value;
				}
			}
	 }
		assert(split_string != NULL);		
		sscanf(split_string, "%lg", &temp_value);
		gsl_vector_set(the_options->training, i, temp_value);
	}
	
	fprintf(stderr, "read data in range: %g..%g\n", range_min, range_max);
	the_options->range_min = range_min;
	the_options->range_max = range_max;
	
}

//! read from a file of fixed length
/**
 * reads from a file of fixed length and fills the eopts struct
 * not a very portable function
 */
void read_input_from_file(char* filename, eopts* options){
	int i = 0; 
	int j = 0;
	double temp_value = 0.0;
	FILE *fptr;
	fptr = fopen(options->filename, "r");

	for(i =0; i < options->nmodel_points; i++){
		for(j = 0; j < options->nparams; j++){
			fscanf(fptr, "%lg", &temp_value);
			gsl_matrix_set(options->xmodel, i, j, temp_value);
		}
		fscanf(fptr, "%lg", &temp_value);
		gsl_vector_set(options->training, i, temp_value);
	}
	printf("read the following input matrix: %d x %d\n", options->nmodel_points, options->nparams);
	print_matrix(options->xmodel, options->nmodel_points, options->nparams);
	printf("the training data is:\n");
	print_vector_quiet(options->training, options->nmodel_points);
}



