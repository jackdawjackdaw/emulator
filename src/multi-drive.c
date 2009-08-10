#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "libEmu/emulator.h"
#include "libEmu/estimator.h"
#include "libEmu/maximise.h"
#include "multifit.h"
#include "sys/time.h"
#include "persist.h"
#include "ioread.h"
#include "useful.h"
#include "gsl/gsl_statistics.h"
#include "smasher.h"
#include "multihelper.h"

// if this is defined, calculated regions are extended by one
// model point in each direction (if that's even possible) 
// this is to give a bit more buffering
#define EXTENDBYONE
void process_input_data(char** input_data, eopts* the_options);
unsigned long int get_seed(void);
void read_input_from_file( eopts* options);
double score_region(emuResult *res);

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
	emuResult whole_thing;

	

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

	/* 
	 * in the 1d case
	 *  now there are only 3 hyperparams by default 
	 * hmm not anymore! 
	 *  -> vertical-scale theta0
	 *  -> ugly offset theta1
	 *  -> nugget theta2
	 *  -> length-scale theta3...theta(Nparams-2⎈)
	 */
	the_options.nthetas = 4;
	the_options.nemu_points = 100; 

	// we kind of have to hope that the inputfile is sorted, could sort it...
	// have to free input data after this
	input_data =  unconstrained_read(inputfile, &number_lines);

	fprintf(stderr, "read in %d lines\n", number_lines);
	
  // could use some regex lib to cut out bad lines
	
	the_options.nmodel_points = number_lines;

	process_input_data(input_data, &the_options);

	print_matrix(the_options.xmodel, number_lines, 1);
	vector_print(the_options.training, number_lines); 
	


	// split_result is a matrix of the splits, it's not a results struct
	// this is bad naming on my part. 
	nregions = smasher(&split_result, &nsplits, &the_options, 1, random_number);

	if(nregions == 1){
		// only one region,don't need to split
		printf("only one region\n");
		fptr = fopen("output/wholeThing.txt", "w");
		if(fptr == NULL){
			fprintf(stderr, "couldn't open output file\n"); 
			exit(1);
		}
		alloc_emuRes(&whole_thing, &the_options);
		evaluate_region(&whole_thing, &the_options, random_number);
		dump_result(&whole_thing, fptr);		
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
	free_emuRes(&whole_thing);	
	free_char_array(input_data, number_lines);
	return(0);
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
void read_input_from_file( eopts* options){
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



