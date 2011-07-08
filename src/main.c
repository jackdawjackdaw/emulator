#include "main.h"

	
/**
 * 
 * estimator
 * @file main.c
 * @author C.Coleman-Smith, cec24@phy.duk.edu
 * \brief contains the main routine for the estimator code
 */



/**
 * \brief the main routine for the estimator code, run this to produce estimated hyperparameters
 *
 * produces estimated hyperparameters, thetas, which are written into ./thetas.txt for 
 * use by emulator
 * 
 * reads input data from stdin, input data should be laid out in columns. Denoting the 
 * training value (the value of the fn we wish to interpolate) as y_i and the model values
 * (the values of the various model  parameters x_i_1, x_i_2,... such that y_i = f(x_i_1, x_i_2, ...) 
 * the input should be setup as:
 * 
 * x_i_1 x_i_2 ... x_i_nparams y_i
 * 
 * where the whitespace can be tabs of spaces. 
 * 
 * if you have more than one parameter, i.e the model data is a result of a multidimensional
 * parameter set you need to adjust options.nparams at the command line (see parse_arguments)
 *
 * the covariance function is fundamental to the estimation and emulation of the data, 
 * setup_cov_fn sets a fn-ptr in the options optstruct which is then evaluated throughout 
 * the code. You need to be quite careful when changing the cov-fn, not only does this usually
 * require a new set of optimization ranges to work well but it may also have a different 
 * number of arguments
 *
 * note: the thetas for the gaussian covaraince fn (power exponential) are log-scaled, this
 * smooths the optimisation landscape for the estimation, making the estimation process much 
 * more robust. As such the optimisation ranges are for log_scaled thetas (they are negative) 
 * the theta values reported to the user should be regularly scaled. 
 * 
 * \todo: add a flag to the optstruct to notify if logscaled thetas are being used or not
 * 
 * the model is estimated with a regression mean and then the difference from the mean is
 * estimated with the gaussian process. The regression model is all setup and controlled by 
 * libEmu/regression.h, it's not controlled by the optstruct sadly. The default is for  
 * a simple linear regression model, this can be changed to a constant model, or a complicated 
 * polynomial model.
 * 
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

	/* malloc all the arrays needed inside the_model */
	alloc_modelstruct(&the_model, &options);

	/* push the input_data into the model structure */
	fill_modelstruct(&the_model, &options, input_data, number_lines);


	setup_optimization_ranges(&options, &the_model);
	
	sprintf(buffer, "nthetas = %d\n", options.nthetas);
	message(buffer, 2);
	sprintf(buffer, "nparams = %d\n", options.nparams);
	message(buffer, 2);

	/* estimate the hyperparameters for this model,  
	 * this is where all the computation takes place
	 * see libEmu/estimate_threaded.c libEmu/maxlbfgs.c and 
	 * refs therein for more information
	 */
	estimate_thetas_threaded(&the_model, &options);

	/* the thetas are log-scaled in the current power-exp covariance fn
	 * be careful */
	fprintf(stderr, "rescaled thetas:");
	for(i = 0; i < options.nthetas; i++){
		if(i != 1){
			fprintf(stderr, " %g", exp(gsl_vector_get(the_model.thetas, i)));
		} else{
			// we've not log scaled the nugget
			fprintf(stderr, " %g", gsl_vector_get(the_model.thetas, i));
		}
		fprintf(stderr, "\n");
	}

	// write the optimum thetas to a text file  ./thetas.txt 
	// perhaps seralising the_model and everything else so that it can be 
	// passed to the emulator would be more efficient
	write_thetas(theta_file, the_model.thetas, &options);

	free_modelstruct(&the_model);
	free_optstruct(&options);
	free_char_array(input_data, number_lines);
	return(0);
}

/**
 * \brief write the thetas vector to the given filename
 *
 * fprints the thetas vector to the given file, doesn't check to see
 * if the file can be opened :)
 */
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
						







