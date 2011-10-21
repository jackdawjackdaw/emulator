#include "assert.h"
#include "resultstruct.h"
#include "string.h"

/**
 * free a resultstruct 
 */
void free_resultstruct(resultstruct *res){
	gsl_matrix_free(res->new_x);
	gsl_vector_free(res->emulated_mean);
	gsl_vector_free(res->emulated_var);
}

/** 
 * allocate the resultstruct from the options
 */
void alloc_resultstruct(resultstruct *res, optstruct* opts){
	/* fprintf(stderr, "allocating resultstruct\n"); */
	/* fprintf(stderr, "nemupts: %d\n", opts->nemulate_points);	 */
	/* fprintf(stderr, "nparams: %d\n", opts->nparams);	 */
	res->new_x = gsl_matrix_alloc(opts->nemulate_points, opts->nparams);
	res->emulated_mean = gsl_vector_alloc(opts->nemulate_points);
	res->emulated_var = gsl_vector_alloc(opts->nemulate_points);
	res->options = opts;
}

/**
 * src -> dst
 */
void copy_resultstruct(resultstruct *dst, resultstruct *src){
	gsl_matrix_memcpy(dst->new_x, src->new_x);
	gsl_vector_memcpy(dst->emulated_mean, src->emulated_mean);
	gsl_vector_memcpy(dst->emulated_var, src->emulated_var);
	dst->options = src->options;
}

//! fill a resultstruct from unconstrained input
/**
 * parse the results of unconstrained read and 
 * turn them into nemupoints rows of nparams positions in 
 * the parameter space of the model
 *
 * these are the points we seek the run the emulator at.
 */
void fill_resultstruct(resultstruct* res, optstruct* options, char** input_data, int number_lines){
	int i,j;
	double temp_value;
	char* split_string;

	// there's a bug, this can't handle empty lines at the end of the input!
	for(i = 0; i < options->nemulate_points; i++){
		split_string = strtok(input_data[i], "\t ");		
		for(j=0; j < options->nparams; j++){
			printf("%s\n", split_string);
			// split string into tab or space tokens
			// each time you do it split_string is pointed to the next block
			// it will come up null when you're done
			assert(split_string != NULL);
			sscanf(split_string, "%lg", &temp_value);
			//fprintf(stderr,"param: %s\n", split_string);
			gsl_matrix_set(res->new_x, i, j, temp_value);
			split_string = strtok(NULL, "\t ");
		}
	}

	fprintf(stderr, "fill_resultsruct: matrix: %d x %d\n", options->nemulate_points, options->nparams);
	print_matrix(res->new_x, options->nemulate_points, options->nparams);
}

