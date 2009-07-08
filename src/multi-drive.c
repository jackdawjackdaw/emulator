
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
void smasher(gsl_matrix *split_ranges, int* nsplits, eopts* toplevel, int max_depth, int min_points, gsl_rng* random_number);
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
	
	smasher(split_result, &nsplits, &the_options, 1, 1, random_number);
	printf("made %d splits\n", nsplits);
 

	//fclose(fptr);
	gsl_rng_free(random_number);
	free_eopts(&the_options);
	free_emuRes(&wholeThing);
	
	free_char_array(input_data, 128, number_lines);
	return(0);
}

#define VERYGOOD 1E-10
// runs the whole splitting thing
void smasher(gsl_matrix *split_ranges, int* nsplits, eopts* toplevel, int max_depth, int min_points, gsl_rng* random_number){
	int i;
	int depth;
	double best_goodness;
	double temp_goodness;
	emuResult temp_result;
	// eval the toplevel
	alloc_emuRes( &temp_result, toplevel);
	evaluate_region(&temp_result, toplevel,  random_number);
	temp_goodness = score_region(&temp_result);
	fprintf(stderr, "in %g..%g goodness = %g\n", toplevel->range_min, toplevel->range_max, temp_goodness);
	if(temp_goodness < VERYGOOD){
		fprintf(stderr, "no need to split, this is a good fit!\n");
		*nsplits = 1;
		split_ranges = gsl_matrix_alloc(1, 2);
		gsl_matrix_set(split_ranges, 0,0, toplevel->range_min);
		gsl_matrix_set(split_ranges, 0,1, toplevel->range_max);
	}
}
	
double score_region(emuResult *res){
	int i;
	double goodness = 0.0;
	double *inverse_error;
	double *diff_error;
	// this probably has to be tuned 
	double diff_thresh = 1.0; // this means, aer they in 1 efolding?
	int total_cluster_count =0;
	int temp_cluster_count = 0;
	int total_cluster_lengths = 0;
	int cluster_min = 5; // how many successive emupoints -> cluster
	double *reduced_inverse_error;
	int *cluster;
	int cluster_begin = 0;
	int cluster_end = 0;
	// start by calculating 1/dy^2
	inverse_error = MallocChecked(sizeof(double)*res->nemu_points);
	diff_error = MallocChecked(sizeof(double)*res->nemu_points);
	cluster = MallocChecked(sizeof(int)*res->nemu_points);

	for(i = 0; i < res->nemu_points; i++){
		diff_error[i] = 0.0;
		cluster[i] = 0.0;
	}

	for(i= 0; i < res->nemu_points; i++){
		inverse_error[i] = (1.0)/(pow(gsl_vector_get(res->new_var, i), 2.0));
		inverse_error[i] = log(inverse_error[i]);
		//reduced_inverse_error[i] = inverse_error[i]*(gsl_vector_get(res->new_mean, i));
		if(i > 0)
			diff_error[i] = fabs(inverse_error[i]-inverse_error[i-1]);
	}
	
	// now look to see if successive diff_errors are less than the thresh and if so we have a cluster?
	for(i = 0; i< res->nemu_points-1; i++){
		if(fabs(diff_error[i+1]-diff_error[i]) < diff_thresh)
			cluster[i] = 1;
	}

	for(i =0; i < res->nemu_points; i++){
		printf("%d:%g\t%g\t%g\t%g\t%g\t%d\n", i,gsl_matrix_get(res->new_x, i,0), gsl_vector_get(res->new_mean, i), gsl_vector_get(res->new_var, i), inverse_error[i], diff_error[i], cluster[i]);
	}
	
	for(i = 0; i < res->nemu_points-1; i++){
		if(cluster[i] == 1){
			if(cluster[i+1] == 1){
				temp_cluster_count++;
				if(i>0 && cluster[i-1] ==0){
					cluster_begin = i;
				}
			} else if(cluster[i+1] == 0){
				// i.e we're at the end of a cluster
				if (temp_cluster_count > cluster_min){
					total_cluster_count++;					
					cluster_end = i;
					total_cluster_lengths += (cluster_end - cluster_begin);
					printf("cluster: %d %d len=%d\n", cluster_begin, cluster_end, (cluster_end-cluster_begin));
					temp_cluster_count = 0;
				}
				else {
					temp_cluster_count = 0;
					cluster_begin = 0;
					cluster_end = 0;
				}
			}
		}
	}
				

	printf("found %d clusters\n", total_cluster_count);

	if((double)(res->nemu_points - total_cluster_lengths)/(double)(res->nemu_points) < 0.2){
		printf("this whole thing is probably one cluster");
	}


	free(inverse_error);
	//free(reduced_inverse_error);
	return(goodness);
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



// RNG 
// tries to read from /dev/random, or otherwise uses the system time
unsigned long int get_seed(void){
	unsigned int seed;
	struct timeval tv;
	FILE *devrandom;

	if((devrandom = fopen("/dev/random", "r")) == NULL){
		gettimeofday(&tv, 0);
		seed = tv.tv_sec + tv.tv_usec;
		fprintf(stderr,"Got seed %u from gettimeofday()\n", seed);
	}
	else {
		fread(&seed, sizeof(seed), 1, devrandom);
		fprintf(stderr,"Got seed %u from /dev/random\n", seed);
		fclose(devrandom);
	}
	return(seed);
}
