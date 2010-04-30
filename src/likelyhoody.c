#include "main.h" 


// this lives in libEmu/emulator.c it's important!
extern emulator_opts the_emulator_options;

void calculate_likelyhood_gauss(gsl_matrix* xmodel_input, gsl_vector* training_vector, gsl_matrix* likelyhood, optstruct* options);
void do_sanity_check( gsl_matrix* xmodel_input, gsl_vector * training_vector, optstruct *options);

//! print the short-option switches
void print_usage(void){
	printf("likelyhoody\n");
	printf("options are: \n");
	printf("t->number of thetas should be (2+nparams) for gaussian or 4 for matern\n");
	printf("p->number of params\n");
	printf("n->number of model_points\n");
	/* printf("m->number of emulator points\n"); */
	/* printf("a->min emulator value\n"); */
	/* printf("b->max emulator value\n"); */
}
	

//! parse the command line 
void parse_arguments(int argc, char** argv, optstruct* options){
	int theta_val = NTHETASDEFAULT;
	int param_val = NPARAMSDEFAULT;
	int nemulate_val = 64;
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

	if(options->nthetas != options->nparams + 2){
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
	gsl_matrix* likelyhood;
	gsl_vector* training_vector;
	gsl_vector* thetas;	
	char input_file[128];
	char** input_data;
	int number_lines = 0;

	FILE *fptr;

	parse_arguments(argc, argv, &options);	
	
	#ifdef NELDER
	fprintf(stderr, "using nelder-mead\n");
	#elif BFGS
	fprintf(stderr, "using bfgs\n");
	#else 
	fprintf(stderr, "using lbfgs\n");
	#endif

	sprintf(input_file, "%s",  "stdin");

	sprintf(options.outputfile, "likelyhoody-out.txt");

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
	// this is correct for the simple linear fit in each dimension plus a constant intercept
	options.nregression_fns = options.nparams + 1;
	//!!!! 


	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// this is key
	// fills in a structure in libEmu which 
	// sets gaussian or matern cov fn and 
	// the alpha option for the gaussian
	//set_emulator_defaults(&the_emulator_options);
	// use the matern cov fn
	the_emulator_options.usematern = 0;
	the_emulator_options.alpha = 2.0;
	// show the default options in the lib
	print_emulator_options(&the_emulator_options);
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	// we only need 4, so maximisation is a little nicer
	if(the_emulator_options.usematern == 1 || the_emulator_options.usematern_three == 1 || the_emulator_options.usematern_five == 1){
		options.nthetas = 4;
	}
	
	xmodel_input = gsl_matrix_alloc(options.nmodel_points, options.nparams);
	training_vector = gsl_vector_alloc(options.nmodel_points);
	thetas = gsl_vector_alloc(options.nthetas);
	
	// proc the input_data
	for(i = 0; i < options.nmodel_points; i++){
		split_string = strtok(input_data[i], "\t ");
		for(j=0; j < options.nparams; j++){
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

	// sanity check, this is to compare against mm output
	do_sanity_check( xmodel_input, training_vector, &options);
	

	likelyhood = gsl_matrix_alloc(options.nemulate_points, options.nemulate_points);

	calculate_likelyhood_gauss(xmodel_input, training_vector, likelyhood, &options);
	

	/* r likes to make contour plots of things which come
	 * as big matrices
	 * we can deal with this
	 */
	fptr = fopen("likeyhood-surf.dat", "w");
	fprintf(fptr, "#%d\n", options.nemulate_points);
	
	for(i = 0; i< options.nemulate_points; i++){
		for(j = 0; j < options.nemulate_points; j++){
			fprintf(fptr,"%g ", gsl_matrix_get(likelyhood, i, j));
		}
		fprintf(fptr, "\n");
	}


	fclose(fptr);
	gsl_vector_free(thetas);
	gsl_vector_free(training_vector);
	gsl_matrix_free(xmodel_input);
	gsl_matrix_free(likelyhood);
	free_char_array(input_data, number_lines);
	return(0);

}

/**
 * compare against some specific output from mm, only works if the xmodel_input and training_vector
 * are pulled from mathematica_model.dat
 *
 * there's a design, to all of this
 */
void do_sanity_check( gsl_matrix* xmodel_input, gsl_vector * training_vector, optstruct *options){
	gsl_vector *thetas_special;
	gsl_vector *h_vector;
	gsl_vector *beta_vector;
	gsl_matrix *c_matrix;
	gsl_matrix *cinverse;
	gsl_matrix *h_matrix;
	gsl_matrix *temp_matrix;
		
	double the_likelyhood;

	/* sanity check things  */
	gsl_permutation *c_LU_permutation = gsl_permutation_alloc(options->nmodel_points);
	int lu_signum = 0, i;
	c_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	cinverse = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	h_vector = gsl_vector_alloc(options->nregression_fns);
	beta_vector = gsl_vector_alloc(options->nregression_fns);
	h_matrix = gsl_matrix_alloc(options->nmodel_points, options->nregression_fns);
	temp_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);


	thetas_special = gsl_vector_alloc(3);
	gsl_vector_set(thetas_special, 0, 0.007);
	gsl_vector_set(thetas_special, 1, 0.00383);
	gsl_vector_set(thetas_special, 2, 0.929);
	

	// sanity check
	the_likelyhood = evalLikelyhood(thetas_special, xmodel_input, training_vector, options->nmodel_points, options->nthetas, options->nparams, options->nregression_fns);

	printf("the likelyhood = %g\n", the_likelyhood);

	makeCovMatrix(c_matrix, xmodel_input, thetas_special,options->nmodel_points, options->nthetas, options->nparams);
	gsl_matrix_memcpy(temp_matrix, c_matrix);
	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
	gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse);

	// regression cpts
	makeHMatrix(h_matrix, xmodel_input, options->nmodel_points, options->nparams, options->nregression_fns);
	estimateBeta(beta_vector, h_matrix, cinverse, training_vector, options->nmodel_points, options->nregression_fns);

	fprintf(stderr, "regression cpts: ");
	for(i = 0; i < options->nregression_fns; i++)
		fprintf(stderr, "%g ", gsl_vector_get(beta_vector, i));
	fprintf(stderr, "\n");
	
	gsl_matrix_free(c_matrix);
	gsl_matrix_free(cinverse);
	gsl_vector_free(h_vector);
	gsl_vector_free(beta_vector);
	gsl_vector_free(thetas_special);
	gsl_matrix_free(h_matrix);
	gsl_matrix_free(temp_matrix);
}

/* void calculate_lhood_general(gsl_matrix* xmodel_input, gsl_vector* training_vector, gsl_matrix* likelyhood, int varyparam1, int varyparam2, optstruct* options){ */
/* 	int i,j; */
/* 	double the_likelyhood, offset_one, offset_two; */
/* 	double initial_one, initial_two; */
/* 	gsl_vectors *thetas = gsl_vector_alloc(options->nthetas); */
/* 	gsl_matrix *local_like_matrix = gsl_matrix_alloc(options->nemulate_points, options->nemulate_points); */

/* 	if(varyparam1 != 1){ */
/* 		theta_initial_one = -15; */
/* 		theta_one_offset = 20 / (double)(options->nemulate_points); */
/* 	} else { */
/* 		theta_initial_one = 0.001; */
/* 		theta_one_offset = 2 / (double)(options->nemulate_points); */
/* 	} */
/* 	if(varyparam2 != 1){ */
/* 		theta_initial_two = -15; */
/* 		theta_two_offset = 20 / (double)(options->nemulate_points); */
/* 	} else { */
/* 		theta_initial_two = 0.001; */
/* 		theta_two_offset = 2 / (double)(options->nemulate_points); */
/* 	} */
	

/* } */

// calculate the likleyhood surface for varying the parameters in the gaussian covariance,
// C = sigma^2 * exp( - r^2 / (beta^2))
// do a square region which runs from options (a -> b)  in both directions
void calculate_likelyhood_gauss(gsl_matrix* xmodel_input, gsl_vector* training_vector, gsl_matrix* likelyhood, optstruct* options){
	int i, j;
	double the_likelyhood, theta_zero_offset, theta_one_offset, theta_initial;
	double theta_initial_one = -15;
	//double side_length = fabs(options->emulate_max - options->emulate_min);
	gsl_vector *thetas = gsl_vector_alloc(options->nthetas);
	gsl_matrix *local_like_matrix = gsl_matrix_alloc(options->nemulate_points, options->nemulate_points);
	

	// need to start down here for the log'd thetas, better shaped landscape
	theta_initial_one = -15;
	theta_initial = 0.001;

	theta_zero_offset = 20 / (double)(options->nemulate_points);
	theta_one_offset = 2 / (double)(options->nemulate_points);
	
	gsl_vector_set_zero(thetas);
	gsl_vector_set(thetas, 1, 0.003);

	for(i = 0; i < options->nemulate_points; i++){
		//		gsl_vector_set(thetas, 0, exp(theta_initial_one + (double)i*theta_zero_offset));		
		gsl_vector_set(thetas, 2, theta_initial_one + (double)i*theta_zero_offset);		
		for(j = 0; j < options->nemulate_points; j++){
			// varying the nugget over a huge range
			gsl_vector_set(thetas,1, theta_initial + (double)j*theta_one_offset);

			the_likelyhood = evalLikelyhood(thetas, xmodel_input, training_vector, options->nmodel_points, options->nthetas, options->nparams, options->nregression_fns);
			gsl_matrix_set(local_like_matrix, i, j, the_likelyhood);
			//printf("%g  %g  %g\n", gsl_vector_get(thetas, 0), gsl_vector_get(thetas,2), the_likelyhood);
		}
		gsl_vector_set(thetas, 1, theta_initial);
	}
	gsl_matrix_memcpy(likelyhood, local_like_matrix);
}




	
	


		
