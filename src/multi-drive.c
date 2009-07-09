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
int smasher(gsl_matrix *split_ranges, int* nsplits, eopts* toplevel, int max_depth, int min_points, gsl_rng* random_number);

/* typedef struct emuResult{ */
/* 	int nemu_points; */
/* 	int nparams; */
/* 	gsl_matrix* new_x; */
/* 	gsl_vector* new_mean; */
/* 	gsl_vector* new_var; */
/* } emuResult; */


int main (void){
	int i;
	FILE *fptr;
	char inputfile[128];
	eopts the_options;
	char** input_data;
	int number_lines;
	gsl_rng *random_number;
	const gsl_rng_type *T;

	emuResult wholeThing;
	/*emuResult res_dumpTest;
	emuResult region1;
	emuResult region2;
	emuResult region3;

	eopts opts_dumpTest;
	eopts region_1_options;
	eopts region_2_options;
	eopts region_3_options;*/


	int nsplits = 0;
	gsl_matrix *split_result;

	// hand pick the input
	sprintf(inputfile, "../model-cut.dat");

	T = gsl_rng_default;
	random_number = gsl_rng_alloc(T);
	gsl_rng_set(random_number, get_seed());

	
	
	// manual options
	// have to set these, can't be bothered to unpack the input data too much
	// these are defaults
	the_options.nparams = 1;
	the_options.nthetas = 4;
	the_options.nemu_points = 100;
	the_options.range_min = 0.0;
	the_options.range_max = 4.0;

	alloc_emuRes(&wholeThing, &the_options);


	// we kind of have to hope that the inputfile is sorted, could sort it...
	// have to free input data after this
	input_data =  unconstrained_read(inputfile, &number_lines);

	fprintf(stderr, "read in %d lines\n", number_lines);
	
	// could use some regex lib to cut out bad lines
	
	the_options.nmodel_points = number_lines;

	process_input_data(input_data, &the_options);

	print_matrix(the_options.xmodel, number_lines, 1);
	vector_print(the_options.training, number_lines); 
	
	nsplits = smasher(split_result, &nsplits, &the_options, 1, 1, random_number);
	printf("made %d splits\n", nsplits);
 

	//fclose(fptr);
	gsl_rng_free(random_number);
	free_eopts(&the_options);
	free_emuRes(&wholeThing);
	
	free_char_array(input_data, number_lines);
	return(0);
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
 */
int smasher(gsl_matrix *split_ranges, int* nsplits, eopts* toplevel, int max_depth, int min_points, gsl_rng* random_number){
	int i;
	emuResult temp_result;
	region* region_list;
	int nregions;
	int retval = 0;
	int min_model_points = 5;
	// eval the toplevel
	alloc_emuRes( &temp_result, toplevel);
	evaluate_region(&temp_result, toplevel,  random_number);

	create_clusters_1d(&temp_result, &region_list, &nregions);
	
	if(nregions == 0){
		fprintf(stderr, "didn't get any useful regions, fit again\n");
		retval = 0;
		return(retval);
	}
	
	assert(nregions >0);
	// sort the regions on length, descending order
	qsort(region_list, nregions, sizeof(region), (void*)compare_regions);
	
	for(i = 0; i < nregions;i++){
		assign_model_point(toplevel, &(region_list[i]));
	}
	

	free(region_list);
	return(nregions);
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
 * if the split doesn't work, result => NULL
 */
void split_region_options(eopts *result, eopts *parent, double lower, double upper){
	assert(lower < upper);
	
	int i,j;
	int split_low = 0;// the index at which to split the model, training vecs etc
	int split_high=  parent->nmodel_points;
	int new_nemu_points = 0;
	int new_nmodel_points = 0;
	int min_model_points = 3;
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

	if(new_nmodel_points < min_model_points){
		fprintf(stderr, "bad split, not enough new points\n");
		//exit(1);
		bad_flag = 1;
	}
	
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
	}

	
}
	

 
//! setup an emuResult struct from the given options 
void alloc_emuRes(emuResult *thing, eopts *options){
	int n = options->nemu_points;
	int np = options->nparams;
	thing->nemu_points = n;
	thing->new_x = gsl_matrix_alloc(n, np);
	thing->new_mean = gsl_vector_alloc(n);
	thing->new_var = gsl_vector_alloc(n);
	thing->nparams = np;
}


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
 * splits up the raw char input data into real things (xmodel etc)
 * this requires that the following fields in the_options be set correctly
 * nparams, nthetas, nmodel_points
 */
void process_input_data(char** input_data, eopts* the_options){
	int i,j;
	char* split_string;
	double temp_value;
	double junk;
	
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
	 }
		assert(split_string != NULL);		
		sscanf(split_string, "%lg", &temp_value);
		gsl_vector_set(the_options->training, i, temp_value);
	}
		
	
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



