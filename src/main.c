#include "main.h"

// this lives in libEmu/emulator.c it's important!
extern emulator_opts the_emulator_options;


//! print the short-option switches
void print_usage(void){
	printf("emulator\n");
	printf("options are: \n");
	printf("t->number of thetas should be (2+nparams) for gaussian or 4 for matern\n");
	printf("p->number of params\n");
	printf("n->number of model_points\n");
	printf("m->number of emulator points\n");
	printf("a->min emulator value\n");
	printf("b->max emulator value\n");
}
	

//! parse the command line 
void parse_arguments(int argc, char** argv, optstruct* options){
	int theta_val = NTHETASDEFAULT;
	int param_val = NPARAMSDEFAULT;
	int nemulate_val = NEMULATEDEFAULT;
	double min_val = EMULATEMINDEFAULT;
	double max_val = EMULATEMAXDEFAULT;
	char file[128];
	int nmodel_points = 0;
	int c;

	// default
	sprintf(file, "input.txt");

	// short options only
	while (( c = getopt(argc, argv, "f:t:p:n:m:a:b:?")) != -1)
		switch(c)
			{
			case '?':
				print_usage();
				exit(1);
			case 'f':
				sprintf(file, "%s", optarg);
				break;
			case 't':
				theta_val = atoi(optarg);
				break;
			case 'p':
				param_val = atoi(optarg);
				break;
			case 'n':
				nmodel_points = atoi(optarg);
				break;
			case 'm':
				nemulate_val = atoi(optarg);
				break;
			case 'a':
				min_val = strtod(optarg, NULL);
				break;
			case 'b':
				max_val = strtod(optarg, NULL);
				break;								 
			default:
				abort();
			}

	//\todo something is wrong with the theta_val thing
	options->nthetas = theta_val;
	options->nparams = param_val;

	if(options->nthetas != options->nparams + 3){
		fprintf(stderr, "you have possbily selected a crazy value of nthetas...\n");
		// for the moment force them to work
		options->nthetas = options->nparams +2;
	}
		

	options->nmodel_points = nmodel_points;
	options->nemulate_points = nemulate_val;
	options->emulate_min = min_val;
	options->emulate_max = max_val;
	sprintf(options->filename, "%s", file);
}
 


int main (int argc, char ** argv){
	optstruct options;
	char* split_string;
	int i,j;
	double temp_value;
	gsl_matrix* xmodel_input;
	gsl_vector* training_vector;
	gsl_vector* thetas;	
	char input_file[128];
	char theta_file[128];
	char** input_data;
	int number_lines = 0;

	parse_arguments(argc, argv, &options);	
	

	#ifdef NELDER
	fprintf(stderr, "using nelder-mead\n");
	#elif BFGS
	fprintf(stderr, "using bfgs\n");
	#else 
	fprintf(stderr, "using lbfgs\n");
	#endif


	sprintf(theta_file, "thetas.txt");
	
	// testing
	//sprintf(input_file, "%s",  "../short.dat");	
	sprintf(input_file, "%s",  "stdin");

	sprintf(options.outputfile, "emulator-out.txt");

	assert(options.nthetas >0);
	assert(options.nparams >0);

	input_data = unconstrained_read(input_file, &number_lines); 
	fprintf(stderr, "read in %d lines\n", number_lines);

	assert(number_lines >0);

	if(options.nmodel_points != number_lines){
		fprintf(stderr, "options.nmodel_points = %d but read in %d\n", options.nmodel_points, number_lines);
		fprintf(stderr, "redfining options.nmodel_points to reflect read in value\n");
		// change the value to match what we actually read
		options.nmodel_points = number_lines;
	}

	//!!!! set the number of regression fns
	// this is regression model dependant
	// this is correct for a cubic  fit in each dimension plus a constant intercept
	options.nregression_fns = 1 + options.nparams;
	//!!!! 

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// this is key
	// fills in a structure in libEmu which 
	// sets gaussian or matern cov fn and 
	// the alpha option for the gaussian
	set_emulator_defaults(&the_emulator_options);
	// use the matern cov fn
	the_emulator_options.usematern = 0;
	the_emulator_options.alpha = 1.8;
	// show the default options in the lib
	print_emulator_options(&the_emulator_options);
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	// we only need 4, so maximisation is a little nicer
	if(the_emulator_options.usematern == 1 || the_emulator_options.usematern_three == 1 || the_emulator_options.usematern_five == 1){
		options.nthetas = 4;
		sprintf(theta_file, "thetas-matern.txt");
	} else {
		options.nthetas = options.nparams + 2;
		printf("nthetas = %d\n", options.nthetas);
	}

	xmodel_input = gsl_matrix_alloc(options.nmodel_points, options.nparams);
	training_vector = gsl_vector_alloc(options.nmodel_points);
	thetas = gsl_vector_alloc(options.nthetas);
	
	// proc the input_data
	// there's a bug, this can't handle empty lines at the end of the input!
	for(i = 0; i < options.nmodel_points; i++){
		split_string = strtok(input_data[i], "\t ");		
		for(j=0; j < options.nparams; j++){
			printf("%s\n", split_string);
			// split string into tab or space tokens
			// each time you do it split_string is pointed to the next block
			// it will come up null when you're done
			assert(split_string != NULL);
			sscanf(split_string, "%lg", &temp_value);
			//fprintf(stderr,"param: %s\n", split_string);
			gsl_matrix_set(xmodel_input, i, j, temp_value);
			split_string = strtok(NULL, "\t ");
		}
		assert(split_string != NULL);
		sscanf(split_string,"%lg", &temp_value);
		//fprintf(stderr,"train: %s\n", split_string);
		gsl_vector_set(training_vector, i, temp_value);
	}

	fprintf(stderr, "read the following input matrix: %d x %d\n", options.nmodel_points, options.nparams);
	//print_matrix(xmodel_input, options.nmodel_points, options.nparams);
	fprintf(stderr, "the training data is:\n");
	//print_vector_quiet(training_vector, options.nmodel_points);
	
	fprintf(stderr, "nthetas = %d\n", options.nthetas);
	fprintf(stderr, "nparams = %d\n", options.nparams);
	

	estimate_thetas_threaded(xmodel_input, training_vector, thetas, &options);


	if(the_emulator_options.usematern == 0){
		fprintf(stderr, "rescaled thetas:");
		for(i = 0; i < options.nthetas; i++){
			if(i != 1){
				fprintf(stderr, " %g", exp(gsl_vector_get(thetas, i)));
			} else{
				// we've not log scaled the nugget
				fprintf(stderr, " %g", gsl_vector_get(thetas, i));
			}
		}
		fprintf(stderr, "\n");
	}

	write_thetas(theta_file, thetas, &options);

	// calc the new means, new variance and dump to emulator-out.txt
	// we'll do this in the emulator code now
	//emulate_model(xmodel_input, training_vector, thetas, &options);

	gsl_vector_free(thetas);
	gsl_vector_free(training_vector);
	gsl_matrix_free(xmodel_input);
	free_char_array(input_data, number_lines);
	//exit(1);
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
						


/**
 * take the estimated parameters and turn them into an emulated set of model points 
 * which can then be output to stdio or whatever 
 */
void emulate_model(gsl_matrix* xmodel, gsl_vector* training, gsl_vector*thetas, optstruct* options){
	int i = 0;
	int j = 0; 
	int n_emu_points = options->nemulate_points;
	
	double temp_mean, temp_var;
	double kappa = 0; // for now
	gsl_matrix *new_x = gsl_matrix_alloc(n_emu_points, options->nparams);
	gsl_vector *new_mean = gsl_vector_alloc(n_emu_points);
	gsl_vector *new_variance = gsl_vector_alloc(n_emu_points);
	gsl_matrix *c_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_matrix *cinverse = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_vector *kplus = gsl_vector_alloc(options->nmodel_points);
	gsl_vector *h_vector = gsl_vector_alloc(options->nregression_fns);
	gsl_vector *beta_vector = gsl_vector_alloc(options->nregression_fns);
	gsl_matrix *h_matrix = gsl_matrix_alloc(options->nmodel_points, options->nregression_fns);
	gsl_vector_view new_x_row;

	gsl_matrix *temp_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_permutation *c_LU_permutation = gsl_permutation_alloc(options->nmodel_points);
	int lu_signum = 0;
	
	FILE *fptr;
	fptr = fopen(options->outputfile, "w");




	makeCovMatrix(c_matrix, xmodel, thetas,options->nmodel_points, options->nthetas, options->nparams);
	gsl_matrix_memcpy(temp_matrix, c_matrix);
	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
	gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse);

	// regression cpts
	makeHMatrix(h_matrix, xmodel, options->nmodel_points, options->nparams, options->nregression_fns);
	estimateBeta(beta_vector, h_matrix, cinverse, training, options->nmodel_points, options->nregression_fns);
	
	fprintf(stderr, "regression cpts: ");
	for(i = 0; i < options->nregression_fns; i++)
		fprintf(stderr, "%g ", gsl_vector_get(beta_vector, i));
	fprintf(stderr, "\n");
	
	// set the new_x values
	initialise_new_x(new_x, options->nparams, options->nemulate_points, options->emulate_min, options->emulate_max);

	for(i = 0; i < n_emu_points; i++){
		new_x_row = gsl_matrix_row(new_x, i);
		makeKVector(kplus, xmodel, &new_x_row.vector, thetas, options->nmodel_points, options->nthetas, options->nparams);
		makeHVector(h_vector, &new_x_row.vector, options->nparams);

		temp_mean = makeEmulatedMean(cinverse, training, kplus, h_vector, h_matrix, beta_vector, options->nmodel_points);

		kappa = covariance_fn(&new_x_row.vector, &new_x_row.vector, thetas, options->nthetas, options->nparams);
		temp_var = makeEmulatedVariance(cinverse, kplus, h_vector, h_matrix, kappa, options->nmodel_points, options->nregression_fns);
		gsl_vector_set(new_mean, i, temp_mean);
		gsl_vector_set(new_variance, i, temp_var);
	}
	 
	for(i = 0; i < n_emu_points; i++){
		for(j = 0; j < options->nparams; j++){
			fprintf(fptr, "%g\t", gsl_matrix_get(new_x, i, j));
		}
		fprintf(fptr, "%g\t", gsl_vector_get(new_mean, i));
		fprintf(fptr,"%g\n", gsl_vector_get(new_variance, i));
	}
		
	gsl_matrix_free(new_x);
	gsl_vector_free(new_mean);
	gsl_vector_free(new_variance);
	gsl_vector_free(h_vector);
	gsl_vector_free(beta_vector);
	gsl_matrix_free(c_matrix);
	gsl_matrix_free(cinverse);
	gsl_matrix_free(h_matrix);
	gsl_vector_free(kplus);
	gsl_matrix_free(temp_matrix);
	gsl_permutation_free(c_LU_permutation);
	fclose(fptr);
}

/* /\**  */
/*  * does the maximum likelyhood estimation on the model_input to try and find the best hyperparams */
/*  *\/ */
/* void estimate_thetas(gsl_matrix* xmodel_input, gsl_vector* training_vector, gsl_vector* thetas, optstruct* options){ */
/* 	const gsl_rng_type *T; */
/* 	gsl_rng *random_number; */
/* 	int max_tries = 10; */
/* 	int i;  */
/* 	int number_steps = 50; */
/* 	gsl_matrix *grad_ranges = gsl_matrix_alloc(options->nthetas, 2); */

/* 	T = gsl_rng_default; */
/* 	random_number = gsl_rng_alloc(T); */
/* 	gsl_rng_set(random_number, get_seed()); */
	
/* 	/\* set the ranges for the initial values of the NM lookup,  */
/* 	 * might want to adjust these as required etc, but whatever *\/ */
/* 	for(i = 0; i < options->nthetas; i++){ */
/* 		gsl_matrix_set(grad_ranges, i, 0, 0.1); */
/* 		gsl_matrix_set(grad_ranges, i, 1, 5.5); */
/* 	} */
	
	
/* 	// the nugget ranges for the gaussian */
/* 	/\* 	gsl_matrix_set(grad_ranges, 2, 0,  0.0000001); *\/ */
/* 	/\* 	gsl_matrix_set(grad_ranges, 2, 1, 0.01); *\/ */
	
/* 	// the nugget ranges for the matern model */
/* 	/\* if(options->nthetas == 4){ // matern *\/ */
/* 	/\* 		gsl_matrix_set(grad_ranges, 3, 0, 0.01); *\/ */
/* 	/\* 		gsl_matrix_set(grad_ranges, 3, 1, 0.1); *\/ */
/* 	/\* 	} *\/ */

/* 	nelderMead(random_number, max_tries, number_steps, thetas, grad_ranges, xmodel_input, training_vector, options->nmodel_points, options->nthetas, options->nparams); */

/* 	fprintf(stderr, "best_thetas: \t"); */
/* 	print_vector_quiet(thetas, options->nthetas); */

/* 	gsl_rng_free(random_number); */
/* 	gsl_matrix_free(grad_ranges); */
/* } */






