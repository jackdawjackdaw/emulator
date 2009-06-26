
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "emulator.h"
#include "estimator.h"
#include "maximise.h"
#include "multifit.h"
#include "sys/time.h"



char** unconstrained_read(char* filename, int* line_count);
void free_char_array(char** array, int lx, int ly);
void process_input_data(char** input_data, eopts* the_options);
void free_eopts(eopts* options);
unsigned long int get_seed(void);
void dump_result(emuResult *res, FILE *fptr);
void alloc_emuRes(emuResult *thing, eopts *options);
void free_eopts(eopts* options);
void free_emuRes(emuResult *thing);
/* typedef struct emuResult{ */
/* 	int nemu_points; */
/* 	int nparams; */
/* 	gsl_matrix* new_x; */
/* 	gsl_vector* new_mean; */
/* 	gsl_vector* new_var; */
/* } emuResult; */


int main (void){
	int i;
	char inputfile[128];
	eopts the_options;
	char** input_data;
	int number_lines;
	gsl_rng *random_number;
	const gsl_rng_type *T;

	emuResult wholeThing;



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
	input_data =  unconstrained_read(inputfile, &number_lines);

	fprintf(stderr, "read in %d lines\n", number_lines);
	
	// could use some regex lib to cut out bad lines
	
	the_options.nmodel_points = number_lines;

	process_input_data(input_data, &the_options);

	print_matrix(the_options.xmodel, number_lines, 1);
	vector_print(the_options.training, number_lines);
	//estimate_region(&the_options, random_number);
	
	//evaluate_region(&wholeThing, &the_options, random_number);
	
	
	//dump_result(&wholeThing, stdout);

	gsl_rng_free(random_number);
	free_eopts(&the_options);
	free_emuRes(&wholeThing);
	return(0);
}


void dump_result(emuResult *res, FILE *fptr){
	int i;
	for(i = 0; i < res->nemu_points; i++){
		fprintf(fptr, "%g\t%g\t%g\n", gsl_matrix_get(res->new_x, i, 0), gsl_vector_get(res->new_mean, i), gsl_vector_get( res->new_var, i));
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
 * splits up the raw char input data into real things (xmodel etc)
 * this requires that the following fields in the_options be set correctly
 * nparams, nthetas, nmodel_points
 */
void process_input_data(char** input_data, eopts* the_options){
	int i,j;
	double temp_value;
	double junk;
	
	assert(the_options->nmodel_points > 0); 
	// first allocate the buffers in options
	the_options->xmodel = gsl_matrix_alloc(the_options->nmodel_points, the_options->nparams);
	the_options->training = gsl_vector_alloc(the_options->nmodel_points);
	the_options->thetas = gsl_vector_alloc(the_options->nthetas);

	for(i = 0; i < the_options->nmodel_points; i++){
		for( j = 0; j < the_options->nparams; j++){
			sscanf(input_data[i], "%lg", &temp_value); 		
			gsl_matrix_set(the_options->xmodel, i, j, temp_value);
	 }
		// HACK HACK HACK, only works for 1 param
		sscanf(input_data[i], "%lg %lg", &junk,  &temp_value);
		gsl_vector_set(the_options->training, i, temp_value);
	}
	// this is broken now, so i'm going to cheat


		
	
}


//! read from a file of fixed length
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
	vector_print(options->training, options->nmodel_points);
}


//! read a file of unknown length
/**
 * @return the read in data
 * @param line_count thenumebr of lines we read
 */
char** unconstrained_read(char* filename, int* line_count_final){
	void copy_char_arrays(char** dest, char** source, int lx, int ly);
	int i;
	FILE *fptr;
	char** input_data;
	char** temp_buffer;
	int init_number_lines = 20;
	int actual_number_lines = init_number_lines;
	int previous_number_lines;
	int line_width = 256; // assume that lines are not wider than this... (right?)
	int line_count = 0;
	char* is_end = 0;
	int buffer_size;
	fptr = fopen(filename, "r");
	if(fptr == NULL){
		fprintf(stderr, "could not open inputfile\n");
		exit(1);
	}

	input_data = malloc(sizeof(char*)*actual_number_lines);
	for(i = 0; i < actual_number_lines; i++){
		input_data[i] = malloc(sizeof(char)*line_width);
	}
	
	buffer_size = (sizeof(char*)*actual_number_lines)*sizeof(char)*line_width;
	fprintf(stderr, "buffer_size is %d\n", buffer_size);

	temp_buffer = malloc(sizeof(char*)*actual_number_lines);
	for(i = 0; i < actual_number_lines; i++){
		temp_buffer[i] = malloc(sizeof(char)*line_width);
	}
		

	do{
		// read line_width chars or up to EOF or EOL
		// if read EOF then is_end == NULL
		is_end = fgets(input_data[line_count], line_width, fptr);
		if(line_count > actual_number_lines-1){
			// i.e next read will drop us off the world
			fprintf(stderr, "allocating more space!\n");
			copy_char_arrays(temp_buffer, input_data, line_width, actual_number_lines);			
			// free the old space, this is a bit tricky since we are trying to use structured data
			free_char_array(input_data, line_width, actual_number_lines);		
			
			previous_number_lines = actual_number_lines;

			actual_number_lines = actual_number_lines + init_number_lines; // grow the size
			buffer_size = (sizeof(char*)*actual_number_lines)*sizeof(char)*line_width;
			// reallocate the buffer
			input_data = malloc(sizeof(char*)*actual_number_lines);
			for(i = 0; i < actual_number_lines; i++){
				input_data[i] = malloc(sizeof(char)*line_width);
			}
			
			// copy the data back in 
			copy_char_arrays(input_data, temp_buffer, line_width, previous_number_lines);
			// finally we have to free and realloc the temp buffer
			free_char_array(temp_buffer, line_width, previous_number_lines);
			// and allocate it again
			temp_buffer = malloc(sizeof(char*)*actual_number_lines);
			for(i = 0; i < actual_number_lines; i++){
				temp_buffer[i] = malloc(sizeof(char)*line_width);
			}
						
		}
	  line_count++;
	} while (is_end != NULL);
	line_count--; // (reading EOF overcounts by one)

	fprintf(stderr, "read %d\n", line_count);
	free_char_array(temp_buffer, line_width, actual_number_lines);

	// realloc temp to be just big enough
	temp_buffer = malloc(sizeof(char*)*line_count);
	for(i = 0; i < line_count; i++){
		temp_buffer[i]  = malloc(sizeof(char)*line_width);
	}
	// copy in the final data
	copy_char_arrays(temp_buffer, input_data, line_width, line_count); 

	free_char_array(input_data, line_width, actual_number_lines);
	fclose(fptr);
	*line_count_final = line_count;
	return(temp_buffer);
}

	
//! copy char_arrays of length lx,ly
void copy_char_arrays(char** dest, char** src, int lx, int ly){
	int i;
	for(i = 0; i < ly; i++){
		memcpy(dest[i], src[i], lx);
	}
}

//! free a 2d array
void free_char_array(char** array, int lx, int ly){
	int i;
	for(i = 0; i < ly;i++){
		free(array[i]);
	}
	free(array);
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
