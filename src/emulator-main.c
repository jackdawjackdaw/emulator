#include "main.h"
#include "resultstruct.h"

/** 
 * \file emulator-main.c
 * \brief runs the emulation process
 */



void print_results(optstruct *options, resultstruct *results, FILE* fptr);

/**
 * \brief the binary interface to the emulator fns
 * 
 * reads thetas.txt and the modeldata from stdin (same as estimator)
 *
 * the locations to run the emulator at are given by the file testPts.txt
 * or by any file pointed to by the enivronment variable INPUTFILE
 * 
 * emulate_model_results fills out a resultstruct with the estimated gp model being
 * evaluated at the various positions ion the parameter space. 
 * 
 * the results are then spewed out to emulator-out.txt or another filename.
 */


int main (int argc, char **argv){
	optstruct options;
	modelstruct the_model;
	resultstruct results;
	char buffer[128];
	int i;

	FILE *fptr = NULL;
	char input_file[128];
	char input_file_reference[128];
	char input_file_thetas[128];
	char *test_point_file = NULL;
	char** input_data;
	char** test_points;
	int number_lines_training = 0;
	int number_lines_test = 0;
	int number_reference_points = 0;
	int number_covered_points = 0;
	double *meanVector, *varVector;
	double temp;


	parse_arguments(argc, argv, &options);	
	setup_cov_fn(&options);


	sprintf(input_file, "%s", "stdin");
	sprintf(input_file_thetas, "thetas.txt");

	sprintf(options.outputfile, "emulator-out.txt");

	assert(options.nthetas >0);
	assert(options.nparams >0);

	// hard coded, bad news here
	//options.nthetas = 3;

	fprintf(stderr, "nthetas = %d\n", options.nthetas);
	fprintf(stderr, "nparams = %d\n", options.nparams);

	input_data = unconstrained_read(input_file, &number_lines_training);
	fprintf(stderr, "emu-out:read in %d lines\n", number_lines_training);

	if(options.nmodel_points != number_lines_training){
		fprintf(stderr, "options.nmodel_points = %d but read in %d\n", options.nmodel_points, number_lines_training);
		fprintf(stderr, "redfining options.nmodel_points to reflect read in value\n");
		// change the value to match what we actually read
		options.nmodel_points = number_lines_training;
		// 
	}

	alloc_modelstruct(&the_model, &options);
	fill_modelstruct(&the_model, &options, input_data, number_lines_training);
	setup_optimization_ranges(&options, &the_model);

	/**
	 * read the thetas
	 */
	fptr = fopen(input_file_thetas, "r");
	for(i = 0; i < options.nthetas; i++){
		fscanf(fptr, "%lf", &temp);
		gsl_vector_set(the_model.thetas, i, temp);
	}
	print_vector_quiet(the_model.thetas, options.nthetas);

	/**
	 * read the points that we want to emulate the model at 
	 */
	test_point_file = getenv("INPUTFILE");
	if(test_point_file == NULL){
		sprintf(test_point_file, "testPts.txt");
	} 
	fprintf(stderr, "emu-out:test_point_file: %s\n", test_point_file);
	test_points = unconstrained_read(test_point_file, &number_lines_test);
	fprintf(stderr, "emu-out:read in %d test points\n", number_lines_test);
	if(options.nemulate_points != number_lines_test){
		fprintf(stderr, "options.nemulate_points = %d but read %d\n", options.nemulate_points, number_lines_test);
		options.nemulate_points = number_lines_test;
	}
	alloc_resultstruct(&results, &options);
	fill_resultstruct(&results, &options, test_points, number_lines_test);


	// relies upon new_x which is not really well defined for arbitrary n_params, ned to figure out a way
	// to generate a set of space filling points in n-dims
	emulate_model_results(&the_model, &options, &results);
	/* meanVector = malloc(sizeof(double)*options.nemulate_points); */
	/* varVector = malloc(sizeof(double)*options.nemulate_points); */
	/* emulateAtPointList(&the_model, (results.new_x), &options, meanVector, varVector); */
	print_results(&options, &results, stdout);


	free_resultstruct(&results);
	free_modelstruct(&the_model);
	free_optstruct(&options);
	free_char_array(input_data, number_lines_training);
	free_char_array(test_points, number_lines_test);
	fclose(fptr);
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
