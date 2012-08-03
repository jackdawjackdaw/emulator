#include "assert.h"
#include "optstruct.h"
#include "string.h"
#include "libEmu/emulator.h"
#include "libEmu/regression.h"
#include "libEmu/maxmultimin.h"

void free_optstruct(optstruct *opts){
	gsl_matrix_free(opts->grad_ranges);
}

void copy_optstruct(optstruct *dst, optstruct* src){
	dst->nthetas = src->nthetas;
	dst->nparams = src->nparams;
	dst->nmodel_points = src->nmodel_points;
	dst->nemulate_points = src->nemulate_points;
	dst->nregression_fns = src->nregression_fns;
	dst->regression_order = src->regression_order;
	dst->cov_fn_index = src->cov_fn_index;

	/* strcpy(dst->filename, src->filename); */
	/* strcpy(dst->outputfile, src->outputfile); */

	dst->fixed_nugget = src->fixed_nugget;
	dst->fixed_nugget_mode = src->fixed_nugget_mode;
	
	dst->use_data_scales = src->use_data_scales;
	// this needs to be sized by the number of hyperparams
	dst->grad_ranges = gsl_matrix_alloc(src->nthetas, 2);
	gsl_matrix_memcpy(dst->grad_ranges, src->grad_ranges);
}

/**
 * given the regression model order set the number of regression functions needed 
 * and fix the makeHVector fn called in regression.c
 * 
 */
void setup_regression(optstruct *opts)
{
	assert(opts->regression_order < 4 && opts->regression_order > -1);
	assert(opts->nparams > 0);
	int nparams = opts->nparams;
	switch(opts->regression_order){
	case 0:
		opts->nregression_fns = 1;
		makeHVector = &(makeHVector_trivial);
		break;
	case 1:
		opts->nregression_fns = 1 + nparams;
		makeHVector = &(makeHVector_linear);
		break;
	case 2:
		opts->nregression_fns = 1 + 2*nparams;
		makeHVector = &(makeHVector_quadratic);
		break;
	case 3:
		opts->nregression_fns = 1+ 3*nparams;
		makeHVector = &(makeHVector_cubic);
		break;
	default:
		opts->nregression_fns = 1;
		makeHVector = &(makeHVector_trivial);
		break;
	}
	/* printf("# set regression order to: %d\n", opts->regression_order); */
	/* printf("# set nregression_fns to: %d\n", opts->nregression_fns); */
		
}


/**
 * \brief set the covariance_fn pointer to the gaussian covariance fn
 *
 * currently we can pick between MATERN32, MATERN52 or POWEREXP
 * will freak-out / warn if nthetas are not passed in correctly
 * 
 * also sets the gradient related pointer maxlbfgs.h:(*setupdCdThetaLength)
 * from the various fns derivative_l_<cov_fn_name>
 */

void setup_cov_fn(optstruct *options)
{
	switch(options->cov_fn_index){
	case MATERN32:
		covariance_fn = &(covariance_fn_matern_three);
		makeGradMatLength = &(derivative_l_matern_three);

		if(options->nthetas != 3)
			fprintf(stderr, "# (warn) setup_cov_fn has changed nthetas, potential memory errors abound\n");
		options->nthetas = 3;
		//fprintf(stderr, "# cov_fn: MATERN32\n");
		break;
	case MATERN52:
		covariance_fn = &(covariance_fn_matern_five);
		makeGradMatLength = &(derivative_l_matern_five);

		if(options->nthetas != 3)
			fprintf(stderr, "# (warn) setup_cov_fn has changed nthetas to, potential memory errors abound\n");
		options->nthetas = 3;
		//fprintf(stderr, "# cov_fn: MATERN52\n");
		break;
	case POWEREXPCOVFN:
		// for testing
		covariance_fn = &(covariance_fn_gaussian);
		makeGradMatLength = &(derivative_l_gauss);

		if(options->nthetas != (options->nparams+2))
			fprintf(stderr, "# (warn) setup_cov_fn has changed nthetas from %d, potential memory errors\n", options->nthetas);

		options->nthetas = options->nparams+2;
		//fprintf(stderr, "# cov_fn: POWEREXP\n");
		break;
	default:
		// crap out if given a bad argument
		printf("err: cov_fn_index set to unsupported value %d\n", options->cov_fn_index);
		exit(1);
	}
	
}


/**
 * \brief set the allowed ranges for the maximisation routine lbfgs
 *
 * \todo add a flag to optstruct for log-scaled and regular thetas
 *
 * summer-2011, we should try setting lower-limits on the length-scale thetas
 * interms of the average nearest-neighbour distance in the input data.
 * i.e in a 1d model we may have samples 
 *     x --- x --- x ---  x
 * it doesn't make sense to have our length scale theta_1 < --- 
 * because we dont have information at that frequency!
 * 
 * this fills the options->grad_ranges matrix with lower and upper bounds to be used in the
 * bounded bfgs maximisation routine lbfgs (see libEmu/maxlbfgs.c) for more info on this.
 * 
 * 
 * this is not the case for the other fns, but the rest of the code may have the assumption
 * frozen into it that the ranges *are* log-scaled.
 *
 */
void setup_optimization_ranges(optstruct* options, modelstruct* the_model)
{
	int i = 0;

	/** does it make sense to have the upper limit on theta be E(10.) = 22026? probably not
	 */
	double bigRANGE = 10.0;
	double rangeMin = 0.0, rangeMax = 0.0;
	double fixedNuggetLeeWay = 0.0 ;
	double rangeMinLog = 0.0001;

	double rangeMinNugget = -5.0;//log(0.0011);
	double rangeMaxNugget = -2.0;//log(0.01); //what's a sensible upper limit here?
	/** 
	 * alloc the grad_ranges matrix in the options and 
	 * put in some sensible defaults 
	 */
	options->grad_ranges = gsl_matrix_alloc(options->nthetas, 2);

	/* 
	 * setup that we're using a log scale
	 */
	if(options->cov_fn_index == POWEREXPCOVFN){
		rangeMin = rangeMinLog;
		rangeMax = 5; // exp(5) = 148 this is big!
	} else {
		rangeMin = 0;
		rangeMax = bigRANGE;
	}

	gsl_matrix_set(options->grad_ranges, 0, 0, rangeMinLog);
	gsl_matrix_set(options->grad_ranges, 0, 1, rangeMax);


	gsl_matrix_set(options->grad_ranges, 1, 0, rangeMinNugget);
	gsl_matrix_set(options->grad_ranges, 1, 1, rangeMaxNugget);

	if(options->use_data_scales){ 
		// use length scales set by the data
		for(i = 2; i < options->nthetas; i++){

			if(options->cov_fn_index == POWEREXPCOVFN){
				rangeMin = 0.5*log(gsl_vector_get(the_model->sample_scales, i-2));
				// try stopping the max range at 25 x the lower limit...
				rangeMax = log(25*exp(rangeMin));
				
				} else {
				rangeMin = 0.5*(gsl_vector_get(the_model->sample_scales, i-2));
				// try stopping the max range at 10 x the nyquist limit...
				//rangeMax = *(gsl_vector_get(the_model->sample_scales, i-2));
			}
			
			if(rangeMin > rangeMax){
				fprintf(stderr, "#ranges failed\n");
				printf("# %d ranges: %lf %lf\n", i, rangeMin, rangeMax);
				printf("# sampleScale: %lf\n", gsl_vector_get(the_model->sample_scales, i-2));
				exit(EXIT_FAILURE);
			}
			gsl_matrix_set(options->grad_ranges, i, 0, rangeMin);
			gsl_matrix_set(options->grad_ranges, i, 1, rangeMax);
		} 

		if(isinf(rangeMin) == -1){
			rangeMin = 0.00001;
		}

	} else { // use some default scales

		for(i = 2; i < options->nthetas; i++){
			gsl_matrix_set(options->grad_ranges, i, 0, rangeMin);
			gsl_matrix_set(options->grad_ranges, i, 1, rangeMax);
		} 
	
	}


	if(options->fixed_nugget_mode == 1){
		// force the nugget to be fixed_nugget +- 20%
		fixedNuggetLeeWay = 0.20*(options->fixed_nugget);
		//gsl_matrix_set(options->grad_ranges, 1, 0, options->fixed_nugget - fixedNuggetLeeWay);
		//fix the min at the usual value
		gsl_matrix_set(options->grad_ranges, 1, 0, rangeMinNugget);
		gsl_matrix_set(options->grad_ranges, 1, 1, options->fixed_nugget + fixedNuggetLeeWay);
		printf("# (reset) %d ranges: %lf %lf (nugget)\n", 1, gsl_matrix_get(options->grad_ranges, 1,0), gsl_matrix_get(options->grad_ranges, 1,1));
	}		

	
	/** 
	 * debug info
	 * print the ranges, this is annoying
	 */
	
	double low, high;
	printf("# grad ranges (logged):\n");
	for(i = 0; i < options->nthetas; i++){
			low = gsl_matrix_get(options->grad_ranges, i, 0);
			high = gsl_matrix_get(options->grad_ranges, i, 1);
		
		if(i == 0){
			printf("# %d ranges: %lf %lf (scale)\n", i, low, high);
		} if (i == 1){
			printf("# %d ranges: %lf %lf (nugget)\n", i, low, high);
		} else if (i > 1) {
			printf("# %d ranges: %lf %lf\n", i, low, high);
		}
	}


}

/** 
 * writes the contents of opts to fptr, 
 * each field is written on a new line, a vector takes a whole line and 
 * a matrix takes nrows lines. 
 */
void dump_optstruct(FILE *fptr, optstruct* opts){
	int i;
	fprintf(fptr, "%d\n", opts->nthetas);
	fprintf(fptr, "%d\n", opts->nparams);
	fprintf(fptr, "%d\n", opts->nmodel_points);
	fprintf(fptr, "%d\n", opts->nemulate_points);
	fprintf(fptr, "%d\n", opts->regression_order);
	fprintf(fptr, "%d\n", opts->nregression_fns);
	fprintf(fptr, "%d\n", opts->fixed_nugget_mode);
	fprintf(fptr, "%lf\n", opts->fixed_nugget);
	fprintf(fptr, "%d\n", opts->cov_fn_index);
	fprintf(fptr, "%d\n", opts->use_data_scales);
	
	for(i = 0; i < opts->nthetas; i++)
		fprintf(fptr, "%lf\t%lf\n", gsl_matrix_get(opts->grad_ranges, i, 0),
						gsl_matrix_get(opts->grad_ranges, i, 1));
	
}

/**
 * reads the fptr into an optstruct, we have to allocate the grad_ranges matrix
 * when we do this
 */
void load_optstruct(FILE *fptr, optstruct* opts){
	int i; 
	double rlow, rhigh;
	fscanf(fptr, "%d", &opts->nthetas);
	fscanf(fptr, "%d", &opts->nthetas);
	fscanf(fptr, "%d", &opts->nparams);
	fscanf(fptr, "%d", &opts->nmodel_points);
	fscanf(fptr, "%d", &opts->nemulate_points);
	fscanf(fptr, "%d", &opts->regression_order);
	fscanf(fptr, "%d", &opts->nregression_fns);
	fscanf(fptr, "%d", &opts->fixed_nugget_mode);
	fscanf(fptr, "%lf", &opts->fixed_nugget);
	fscanf(fptr, "%d", &opts->cov_fn_index);
	fscanf(fptr, "%d", &opts->use_data_scales);
	
	
	// allocate the grad_ranges matrix
	opts->grad_ranges = gsl_matrix_alloc(opts->nthetas, 2);

	for(i = 0; i < opts->nthetas; i++){
		fscanf(fptr, "%lf\t%lf", &rlow, &rhigh);
		gsl_matrix_set(opts->grad_ranges, i, 0, rlow);
		gsl_matrix_set(opts->grad_ranges, i, 1, rhigh);
	}
	
}
