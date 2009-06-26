
typedef struct eopts{
	int nmodel_points;
	int nemu_points;
	int nparams;
	int nthetas;
	double range_min;
	double range_max;
	gsl_matrix *xmodel;
	gsl_vector *training;
	gsl_vector *thetas;
} eopts;

// this will only work in 1d for now
//! emulate the given region, returning the new_x, emulated_mean and emulated_variances
/**
 * emulate a given region, the eopts struct has to be set correctly and also 
 * the hyperparameters have to have been set
 * @return new_x are the emulated points, emulated_mean is the mean at these points
 * and emualted_variance is the variance at these points */
void emulate_region(gsl_vector *new_x, gsl_vector* emulated_mean, gsl_vector* emulated_variance , eopts* options){
	int i, j;
	double kappa = 0.0;
	double step_size = (options->range_max - options->range_min) /((double)(options->nemu_points));
	gsl_matrix *c_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_matrix *cinverse = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_vector *kplus = gsl_vector_alloc(options->nmodel_points);
	gsl_vector_view new_x_row;

	gsl_matrix *temp_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_permutation *c_LU_permutation = gsl_permutation_alloc(options->nmodel_points);
	int lu_signum = 0;

	makeCovMatrix(c_matrix, options->xmodel, options->thetas,options->nmodel_points, options->nthetas, options->nparams);
	gsl_matrix_memcpy(temp_matrix, c_matrix);
	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
	gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse);
	
	// set the new_x values
	for(i = 0; i < options->nemu_points; i++){
		// this doesn't make sense for many params!
		for(j = 0; j < options->nparams; j++){
			gsl_matrix_set(new_x, i, j,step_size*((double)i)+options->range_min);
		}
	}
	
	for(i = 0; i < options->nemu_points; i++){
		new_x_row = gsl_matrix_row(new_x, i);
		makeKVector(kplus, options->xmodel, &new_x_row.vector, options->thetas, options->nmodel_points, options->nthetas, options->nparams);
		temp_mean = makeEmulatedMean(cinverse, options->training, kplus, options->nmodel_points);
		kappa = covariance_fn(&new_x_row.vector, &new_x_row.vector, options->thetas, options->nthetas, options->nparams);
		temp_var = makeEmulatedVariance(cinverse, kplus, kappa, options->nmodel_points);
		gsl_vector_set(emulated_mean, i, temp_mean);
		gsl_vector_set(emulated_variance, i, temp_var);
	}

	gsl_matrix_free(c_matrix);
	gsl_matrix_free(cinverse);
	gsl_vector_free(kplus);
	gsl_matrix_free(temp_matrix);
	gsl_permutation_free(c_LU_permutation);
}

//! estimate the hyperparams for a region
/**
 * @param random is a gsl_rng which has already been setup and seeded
 * @return options->thetas is set to the best hyperparams
 */
void estimate_region(eopts* options, gsl_rng *random){
	int max_tries = 20;
	int i;
	int number_steps = 40;
	gsl_matrix *grad_ranges = gsl_matrix_alloc(options->nthetas);
  
	for(i = 0; i < options->nthetas; i++){
		gsl_matrix_set(grad_ranges, i, 0, 0.0);
		gsl_matrix_set(grad_ranges, i, 1, 1.0);
	}
	
	nelderMead(random, max_tries, number_steps, options->thetas, grad_rangers, options->xmodel, options->training, null, options->nmodel_points, options->nthetas, options->nparams);
	
	fprintf(stderr, "in range: %g..%g\n", options->range_min, options->range_max);
	fprintf(stderr, "best thetas: \t"); 
	print_vector_quet(options->thetas, options->nthetas);
	
	gsl_matrix_free(grad_ranges);
}

//! do the estimation and emulation for a region
void evaluate_region(gsl_vector* new_x, gsl_vector* new_mean, gsl_vector* new_variance, eopts* options, gsl_rng * random){
	estimate_region(options, random);
	emulate_region(new_x, new_mean, new_variance, options);
}

//! see if a region is smooth
int is_smooth(double smooth_val, gsl_vector* xemu, gsl_vector* mean_emu, gsl_vector* var_emu, eopts* options){
	int i; 
	double mse = malloc(sizeof(double)*(options->nemu_points));
	double var = 0.0;
	for (i = 0; i < options-> nemu_points; i++){
		mse[i] = get_mse(gsl_vector_get(mean_emu, i), gsl_vector_get(var_emu, i));
	}

	var = gsl_stats_variance(mse, 1, options->nemu-points);
	free(mse);
	if(var > smooth_val){
		fprintf(stderr, "wiggly!\n");
		return(1);
	} else {
		fprintf(stderr, "smooth\n");
		return(0);
	}
	
}


double get_mse( double mean, double variance){
	return(pow(fabs(mean-variance), 2.0));
}
			
