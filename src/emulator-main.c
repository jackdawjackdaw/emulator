#include "main.h"
#include "resultstruct.h"


void print_results(optstruct *options, resultstruct *results, FILE* fptr);

int main (int argc, char **argv){
	optstruct options;
	modelstruct the_model;
	resultstruct results;
	char buffer[128];
	int i;

	char input_file[128];
	char input_file_reference[128];
	char input_file_thetas[128];
	char** input_data;
	int number_lines = 0;
	int number_reference_points = 0;
	int number_covered_points = 0;


	parse_arguments(argc, argv, &options);	
	setup_cov_fn(&options);
	setup_optimization_ranges(&options);

	alloc_resultstruct(&results, &options);
	
	sprintf(input_file, "%s", "stdin");
	sprintf(input_file_thetas, "thetas.txt");

	sprintf(options.outputfile, "emulator-out.txt");

	assert(options.nthetas >0);
	assert(options.nparams >0);

	// hard coded, bad news here
	//options.nthetas = 3;

	sprintf(buffer, "nthetas = %d\n", options.nthetas);
	message(buffer, 2);
	sprintf(buffer, "nparams = %d\n", options.nparams);
	message(buffer, 2);

	input_data = unconstrained_read(input_file, &number_lines);
	sprintf(buffer, "emu-out:read in %d lines\n", number_lines);
	message(buffer, 2);

	if(options.nmodel_points != number_lines){
		fprintf(stderr, "options.nmodel_points = %d but read in %d\n", options.nmodel_points, number_lines);
		fprintf(stderr, "redfining options.nmodel_points to reflect read in value\n");
		// change the value to match what we actually read
		options.nmodel_points = number_lines;
		// 
	}

	alloc_modelstruct(&the_model, &options);
	fill_modelstruct(&the_model, &options, input_data, number_lines);
	// relies upon new_x which is not really well defined for arbitrary n_params, ned to figure out a way
	// to generate a set of space filling points in n-dims
	emulate_model_results(&the_model, &options, &results);
	print_results(&options, &results, stdout);


	free_resultstruct(&results);
	free_modelstruct(&the_model);
	free_optstruct(&options);
	free_char_array(input_data, number_lines);
	return(0);
}

//! dump the results to the file given, in ascii
void print_results(optstruct *options, resultstruct *results, FILE* fptr){
	int i, j;
	for(i = 0; i < options->nemulate_points; i++){
		for(j = 0; j < options->nparams; j++){
			fprintf(fptr, "%g\t", gsl_matrix_get(results->new_x, i,j));
		}
		fprintf(fptr, "%g\t", gsl_vector_get(results->emulated_mean, i));
		fprintf(fptr, "%g\n", gsl_vector_get(results->emulated_var, i));
	}
}
