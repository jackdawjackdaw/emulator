#include "main.h"

int do_coverage_test(gsl_matrix* xmodel, gsl_vector* training, gsl_vector* thetas, gsl_vector* reference_values, gsl_matrix* reference_points, gsl_vector* emulated_mean, gsl_vector* emulated_var, int  number_reference_points, optstruct *options);
void emulate_model_at_point(gsl_matrix* xmodel, gsl_vector* training, gsl_vector* thetas, gsl_vector* point, optstruct* options, double* emulated_mean, double* emulated_var);
void calculate_errors(gsl_matrix* xmodel, gsl_vector* training, gsl_vector* thetas, gsl_vector* reference_values, gsl_matrix* reference_points, gsl_vector* calculated_errors, gsl_vector* emulated_mean, gsl_vector* emulated_var, int number_reference_points, optstruct *options);


void emulate_model(modelstruct* the_model, optstruct* options);

void dump_errors(FILE* fptr, gsl_vector* reference_values, gsl_matrix* reference_points, gsl_vector* reference_errors, gsl_vector* emulated_mean, gsl_vector* emulated_var, int number_reference_points, optstruct* options);



int main (int argc, char **argv){
	optstruct options;
	modelstruct the_model;
	char buffer[128];
	int i;


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
	setup_cov_fn(&options);
	setup_optimization_ranges(&options);
	
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
	emulate_model(&the_model, &options);

	free_modelstruct(&the_model);
	free_optstruct(&options);
	free_char_array(input_data, number_lines);
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
void emulate_model(modelstruct* the_model, optstruct* options){
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




	makeCovMatrix(c_matrix, the_model->xmodel, the_model->thetas,options->nmodel_points, options->nthetas, options->nparams, options->covariance_fn);
	gsl_matrix_memcpy(temp_matrix, c_matrix);
	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
	gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse);

	// regression cpts
	makeHMatrix(h_matrix, the_model->xmodel, options->nmodel_points, options->nparams, options->nregression_fns);
	estimateBeta(beta_vector, h_matrix, cinverse, the_model->training_vector, options->nmodel_points, options->nregression_fns);
	
	fprintf(stderr, "regression cpts: ");
	for(i = 0; i < options->nregression_fns; i++)
		fprintf(stderr, "%g ", gsl_vector_get(beta_vector, i));
	fprintf(stderr, "\n");
	
	// set the new_x values
	initialise_new_x(new_x, options->nparams, options->nemulate_points, options->emulate_min, options->emulate_max);

	for(i = 0; i < n_emu_points; i++){
		new_x_row = gsl_matrix_row(new_x, i);
		makeKVector(kplus, the_model->xmodel, &new_x_row.vector, the_model->thetas, options->nmodel_points, options->nthetas, options->nparams, options->covariance_fn);
		makeHVector(h_vector, &new_x_row.vector, options->nparams);

		temp_mean = makeEmulatedMean(cinverse, the_model->training_vector, kplus, h_vector, h_matrix, beta_vector, options->nmodel_points);

		kappa = options->covariance_fn(&new_x_row.vector, &new_x_row.vector, the_model->thetas, options->nthetas, options->nparams, options->cov_fn_alpha);
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
