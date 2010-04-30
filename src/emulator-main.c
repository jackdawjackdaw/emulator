#include "main.h"

int do_coverage_test(gsl_matrix* xmodel, gsl_vector* training, gsl_vector* thetas, gsl_vector* reference_values, gsl_matrix* reference_points, gsl_vector* emulated_mean, gsl_vector* emulated_var, int  number_reference_points, optstruct *options);
void emulate_model_at_point(gsl_matrix* xmodel, gsl_vector* training, gsl_vector* thetas, gsl_vector* point, optstruct* options, double* emulated_mean, double* emulated_var);
void calculate_errors(gsl_matrix* xmodel, gsl_vector* training, gsl_vector* thetas, gsl_vector* reference_values, gsl_matrix* reference_points, gsl_vector* calculated_errors, gsl_vector* emulated_mean, gsl_vector* emulated_var, int number_reference_points, optstruct *options);

void emulate_model(gsl_matrix* xmodel, gsl_vector* training, gsl_vector*thetas, optstruct* options);

void dump_errors(FILE* fptr, gsl_vector* reference_values, gsl_matrix* reference_points, gsl_vector* reference_errors, gsl_vector* emulated_mean, gsl_vector* emulated_var, int number_reference_points, optstruct* options);

// this lives in libEmu/emulator.c it's important!
extern emulator_opts the_emulator_options;


//! print the short-option switches
void print_usage(void){
	printf("coverage-test\n");
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


int main (int argc, char **argv){
	optstruct options;
	char* split_string;
	int i,j;
	double temp_value;

	gsl_matrix* reference_points;
	gsl_vector* reference_values;
	gsl_vector* emulated_mean;
	gsl_vector* emulated_var;
	gsl_matrix* xmodel_input;
	gsl_vector* training_vector;
	gsl_vector* thetas;	
	gsl_vector* reference_errors;
	char input_file[128];
	char input_file_reference[128];
	char input_file_thetas[128];
	char** input_data;
	int number_lines = 0;
	int number_reference_points = 0;
	int number_covered_points = 0;

	FILE *fptr;

	parse_arguments(argc, argv, &options);	
	
	/* hand code this for now */
	//sprintf(input_file_reference, "2d-gauss-full.txt");
	//sprintf(input_file, "../2d-gauss-emu-sample.txt");
	sprintf(input_file, "%s", "stdin");
	sprintf(input_file_thetas, "thetas.txt");

	sprintf(options.outputfile, "emulator-out.txt");

	assert(options.nthetas >0);
	assert(options.nparams >0);

	// hard coded, bad news here
	//options.nthetas = 3;

	// read the reference distribution for the coverage test
	/* input_data = unconstrained_read(input_file_reference, &number_lines);  */
	/* fprintf(stderr, "reference:read in %d lines\n", number_lines); */
	/* assert(number_lines >0); */
	/* number_reference_points = number_lines; */

	fprintf(stderr, "nthetas = %d\n", options.nthetas);
	fprintf(stderr, "nparams = %d\n", options.nparams);

	
	/* reference_points = gsl_matrix_alloc(number_reference_points, options.nparams); */
	/* reference_values = gsl_vector_alloc(number_reference_points); */

	/* // proc the reference data */
	/* for(i = 0; i < number_reference_points; i++){ */
	/* 	split_string = strtok(input_data[i],"\t "); */
	/* 	for(j=0; j < options.nparams; j++){ */
	/* 		assert(split_string != NULL); */
	/* 		sscanf(split_string, "%lg", &temp_value); */
	/* 		gsl_matrix_set(reference_points, i, j, temp_value); */
	/* 		split_string = strtok(NULL, "\t "); */
	/* 	} */
	/* 	assert(split_string != NULL); */
	/* 	sscanf(split_string, "%lg", &temp_value); */
	/* 	gsl_vector_set(reference_values, i, temp_value); */
	/* }	 */
	/* free_char_array(input_data, number_lines); */


	input_data = unconstrained_read(input_file, &number_lines);
	fprintf(stderr, "emu-out:read in %d lines\n", number_lines);

	if(options.nmodel_points != number_lines){
		fprintf(stderr, "options.nmodel_points = %d but read in %d\n", options.nmodel_points, number_lines);
		fprintf(stderr, "redfining options.nmodel_points to reflect read in value\n");
		// change the value to match what we actually read
		options.nmodel_points = number_lines;
		// 
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



	xmodel_input = gsl_matrix_alloc(options.nmodel_points, options.nparams);
	training_vector = gsl_vector_alloc(options.nmodel_points);
	thetas = gsl_vector_alloc(options.nthetas);

	// now read in the original file used to create the emulator
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
	free_char_array(input_data, number_lines);
	
	// now read in the thetas
	input_data = unconstrained_read(input_file_thetas, &number_lines);
	assert(number_lines == 1); // should just be the thetas
	split_string = strtok(input_data[0], "\t ");
	for(j = 0; j < options.nthetas; j++){
		assert(split_string != NULL);
		sscanf(split_string, "%lg", &temp_value);
		gsl_vector_set(thetas, j, temp_value);
		split_string = strtok(NULL, "\t ");
	}
	free_char_array(input_data, number_lines);

	for(j = 0; j < options.nthetas; j++){
		printf("%g\t", gsl_vector_get(thetas, j));
	}
	printf("\n");
					 
	emulate_model(xmodel_input, training_vector, thetas, &options);
	/**
	 * left the coverage test in here,but its really bad because it doesn't take into account
	 * how large the actual variance is at each point. the errors test is a bit more sensible
	 */

	/* fprintf(stderr, "starting coverage test\n"); */
	/* number_covered_points = do_coverage_test(xmodel_input, training_vector, thetas, reference_values, reference_points, number_reference_points, &options); */
	/* fprintf(stderr, "coverage test\n"); */
	/* printf("%d / %d\n", number_covered_points, number_reference_points); */
	/* printf("%g\n", (double)number_covered_points / (double)number_reference_points); */
	/* fptr = fopen("coverage-results.txt", "a"); */
	/* fprintf(fptr, "%d\t%g\n", options.nmodel_points, (double)number_covered_points / (double)number_reference_points); */
	/* fclose(fptr); */

	/* fprintf(stderr, "calculating errors\n"); */
	/* reference_errors = gsl_vector_alloc(number_reference_points); */
	/* emulated_mean = gsl_vector_alloc(number_reference_points); */
	/* emulated_var = gsl_vector_alloc(number_reference_points); */

	/* calculate_errors(xmodel_input, training_vector, thetas, reference_values, reference_points, reference_errors, emulated_mean, emulated_var, number_reference_points, &options); */
	/* fptr = fopen("cov-errors.txt", "w"); */
	/* dump_errors(fptr, reference_values, reference_points, reference_errors, emulated_mean, emulated_var, number_reference_points, &options); */
	/* fclose(fptr); */

	gsl_vector_free(thetas);
	gsl_vector_free(training_vector);
	/* gsl_vector_free(emulated_mean); */
	/* gsl_vector_free(emulated_var); */
	gsl_matrix_free(xmodel_input);
	/* gsl_vector_free(reference_values); */
	/* gsl_matrix_free(reference_points); */
	/* gsl_vector_free(reference_errors); */
	return(0);
}


void dump_errors(FILE* fptr, gsl_vector* reference_values, gsl_matrix* reference_points, gsl_vector* reference_errors, gsl_vector* emulated_mean, gsl_vector* emulated_var, int number_reference_points, optstruct* options){
	int i, j;
	for(i = 0; i < number_reference_points; i++){
		for(j = 0; j < options->nparams; j++)
			fprintf(fptr, "%g\t", gsl_matrix_get(reference_points, i, j));
		fprintf(fptr, "%g\t%g\t%g\t%g\n", gsl_vector_get(emulated_mean, i), gsl_vector_get(emulated_var,i),gsl_vector_get(reference_values,i), gsl_vector_get(reference_errors, i));
	}
}

/** 
 * do_coverage_test:
 * for each point in reference_values this emulates the xmodel data 
 * using the given thetas and a call to emulate_model_at_point and then
 * compares the reference_value of the point to the range given by 
 * emulated_mean +- emulated_var /2
 *
 * this is not a good test!
 */
/* int do_coverage_test(gsl_matrix* xmodel, gsl_vector* training, gsl_vector* thetas, gsl_vector* reference_values, gsl_matrix* reference_points, gsl_vector* emulated_mean, gsl_vector* emulated_var, int number_reference_points, optstruct *options){ */
/* 	int coverage_count = 0; */
/* 	int i,j; */
/* 	gsl_vector* the_point = gsl_vector_alloc(options->nparams); */
/* 	double temp_mean, temp_var; */
/* 	double temp_sd; */
/* 	double reference_value; */
/* 	double chi_sq = 0.0; */
/* 	double chi_sq_reduce = number_reference_points - options->nparams; */
	
/* 	for(i = 0; i < number_reference_points; i++){ */

/* 		for(j = 0; j < options->nparams; j++){ */
/* 			gsl_vector_set(the_point, j , gsl_matrix_get(reference_points, i, j)); */
/* 			//printf("%g\t", gsl_vector_get(the_point, i)); */
/* 		} */
/* 		//printf("\n"); */
		
/* 		reference_value = gsl_vector_get(reference_values, i); */
/* 		emulate_model_at_point(xmodel, training, thetas, the_point, options, &temp_mean, &temp_var); */
/* 		temp_sd = sqrt(temp_var/2); */
/* 		//printf("rv:%g\tmean:%g\tvar:%g\ttsd:%g\n", reference_value, temp_mean, temp_var, temp_sd); */
/* 		gsl_vector_set(emulated_mean, i, temp_mean); */
/* 		gsl_vector_set(emulated_var, i, temp_var); */

/* 		//chi_sq += (pow(temp_mean - reference_value, 2.0) / reference_value + temp_mean); */

/* 		if(reference_value <= temp_mean + (temp_sd) && reference_value >= temp_mean - (temp_sd)){ */
/* 			coverage_count++; */
/* 		} */
/* 	} */
	
/* 	//printf("chi_sq = %g\n", (chi_sq)); */
	
/* 	return(coverage_count); */
/* } */
	
/**
 * calculate the error at each reference point as given by eqn 14 from o'hagan and bastos
 * 
 * Di = (emulated_mean(x_i) - reference_value(x_i)) / sqrt[var(x_i)]
 */
/* void calculate_errors(gsl_matrix* xmodel, gsl_vector* training, gsl_vector* thetas, gsl_vector* reference_values, gsl_matrix* reference_points, gsl_vector* calculated_errors, gsl_vector* emulated_mean, gsl_vector* emulated_var, int number_reference_points, optstruct *options){ */
/* 	int i,j; */
/* 	gsl_vector* the_point = gsl_vector_alloc(options->nparams); */
/* 	gsl_vector* errors = gsl_vector_alloc(number_reference_points); */
/* 	double temp_mean,temp_var, temp_sd; */
/* 	double temp_error; */
/* 	double reference_value; */

/* 	for(i = 0; i < number_reference_points; i++){ */
/* 		for(j = 0; j < options->nparams; j++){ */
/* 			gsl_vector_set(the_point, j , gsl_matrix_get(reference_points, i, j)); */
/* 		} */
		
/* 		reference_value = gsl_vector_get(reference_values, i); */
/* 		emulate_model_at_point(xmodel, training, thetas, the_point, options, &temp_mean, &temp_var); */
/* 		temp_sd = sqrt(temp_var/2); */

/* 		temp_error = (temp_mean - reference_value)/(temp_sd); */
/* 		gsl_vector_set(errors, i, temp_error); */
/* 		gsl_vector_set(emulated_mean, i, temp_mean); */
/* 		gsl_vector_set(emulated_var, i, temp_var); */
/* 		//chi_sq += (pow(temp_mean - reference_value, 2.0) / reference_value + temp_mean); */
/* 	} */
/* 	gsl_vector_memcpy(calculated_errors, errors); */
/* } */
	

/**
 * for a given model and emulator (thetas) this will
 * give you the prediction at the vector point
 * @param point the location you want the emulator evaluated at
 * @return emulated_mean the mean at point
 * @return emualted_var the var at point
 *
 * It kind of sucks to set up all the crap to do the emulator and then do it for one call
 * and drop it, but lets not worry too much about that
 *
 * this doesn't work right now
 */
/* void emulate_model_at_point(gsl_matrix* xmodel, gsl_vector* training, gsl_vector* thetas, gsl_vector* point, optstruct* options, double* emulated_mean, double* emulated_var){ */
/* 	double temp_mean, temp_var; */
/* 	double kappa = 0; // for now */
/* 	gsl_matrix *c_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points); */
/* 	gsl_matrix *cinverse = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points); */
/* 	gsl_vector *kplus = gsl_vector_alloc(options->nmodel_points); */


/* 	gsl_matrix *temp_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points); */
/* 	gsl_permutation *c_LU_permutation = gsl_permutation_alloc(options->nmodel_points); */
/* 	int lu_signum = 0; */

/* 	makeCovMatrix(c_matrix, xmodel, thetas, options->nmodel_points, options->nthetas, options->nparams); */
/* 	gsl_matrix_memcpy(temp_matrix, c_matrix); */
/* 	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum); */
/* 	gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse); */

/* 	makeKVector(kplus, xmodel, point, thetas, options->nmodel_points, options->nthetas, options->nparams); */
/* 	temp_mean = makeEmulatedMean(cinverse, training, kplus, options->nmodel_points); */
/* 	kappa = covariance_fn(point, point, thetas, options->nthetas, options->nparams); */
/* 	temp_var = makeEmulatedVariance(cinverse, kplus, kappa, options->nmodel_points); */




/* 	// set the final values */
/* 	*emulated_mean = temp_mean; */
/* 	*emulated_var = temp_var; */


	

/* 	gsl_matrix_free(c_matrix); */
/* 	gsl_matrix_free(cinverse); */
/* 	gsl_vector_free(kplus); */
/* 	gsl_matrix_free(temp_matrix); */
/* 	gsl_permutation_free(c_LU_permutation); */
/* } */

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
