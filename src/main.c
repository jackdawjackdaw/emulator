#include "main.h"

	
/**
 * C.Coleman-Smith, cec24@phy.duk.edu
 * estimator
 * 
 * this is the main routine for the estimator code.
 *
 * purpose: to read input from stdin and calculate the optimal hyperparams for the data 
n * and covfn.
 *
 * key structures: modelstruct, optstruct
 * modelstruct -> contains the input data, such as training points, training values, thetas etc
 * optstruct -> contains various constants, number of model_points, number of hyperparameters in the model and so forth
 *  
 * see also: libEmu: estimate_threaded, maxwithlbfgs, a lot of the setup and teardown 
 * is in bin_support
 *
 */

int main (int argc, char ** argv){
	optstruct options;
	modelstruct the_model;
	char buffer[128];
	
	int i;

	char input_file[128];
	char theta_file[128];
	char** input_data;
	int number_lines = 0;

	/* after this the optstruct should be totally filled out */
	parse_arguments(argc, argv, &options);	
	setup_cov_fn(&options);
	setup_optimization_ranges(&options);
	
	/* now we can allocate the modelstruct */
	message("using lbfgs", 1);

	/** 
	 * \todo: the model structure should be serialised so that
	 * it can be re-read in the estimator code
	 */
	sprintf(theta_file, "thetas.txt");
	
	// we're going to read the input from the stdin
	sprintf(input_file, "%s",  "stdin");
	
	/** 
	 * we're going to read from the stdin and push it all into the 
	 * input_data buffer, unformatted and unprocessed.
	 */
	input_data = unconstrained_read(input_file, &number_lines); 
	sprintf(buffer, "read in %d lines\n", number_lines);
	message(buffer,2);

	assert(number_lines >0);

	if(options.nmodel_points != number_lines){
		sprintf(buffer, "options.nmodel_points = %d but read in %d\n", options.nmodel_points, number_lines);
		message(buffer, 2);
		sprintf(buffer, "redfining options.nmodel_points to reflect read in value\n");
		message(buffer, 2);
		// change the value to match what we actually read
		options.nmodel_points = number_lines;
	}

	alloc_modelstruct(&the_model, &options);

	// push the input_data into the model structure
	fill_modelstruct(&the_model, &options, input_data, number_lines);
	
	sprintf(buffer, "nthetas = %d\n", options.nthetas);
	message(buffer, 2);
	sprintf(buffer, "nparams = %d\n", options.nparams);
	message(buffer, 2);

	// estimate the hyperparameters for this model
	estimate_thetas_threaded(&the_model, &options);

	fprintf(stderr, "rescaled thetas:");
	for(i = 0; i < options.nthetas; i++){
		if(i != 1){
			fprintf(stderr, " %g", exp(gsl_vector_get(the_model.thetas, i)));
		} else{
			// we've not log scaled the nugget
			fprintf(stderr, " %g", gsl_vector_get(the_model.thetas, i));
		}
	}
	fprintf(stderr, "\n");

	// write the optimum thetas to a text file  ./thetas.txt 
	// perhaps seralising the_model and everything else so that it can be 
	// passed to the emulator would be more efficient
	write_thetas(theta_file, the_model.thetas, &options);

	free_modelstruct(&the_model);
	free_optstruct(&options);
	free_char_array(input_data, number_lines);
	return(0);
}

void write_thetas(char* theta_file, gsl_vector* thetas, optstruct *options){
	int i;
	FILE *fptr;
	fptr = fopen(theta_file, "w");
	for(i = 0; i < options->nthetas; i++){
		fprintf(fptr, "%g\t", gsl_vector_get(thetas, i));
	}
	fprintf(fptr, "\n");
	fclose(fptr);
}
						







