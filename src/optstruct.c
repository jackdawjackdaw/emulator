
#include "optstruct.h"

// this causes confusion if it's in optstruct.h
// modules, i do not understand you
#include "libEmu/maxlbfgs.h"

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

	strcpy(dst->filename, src->filename);
	strcpy(dst->outputfile, src->outputfile);

	dst->fixed_nugget = src->fixed_nugget;
	dst->fixed_nugget_mode = src->fixed_nugget_mode;
	
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
		makeHVector = *(makeHVector_trivial);
		break;
	case 1:
		opts->nregression_fns = 1 + nparams;
		makeHVector = *(makeHVector_linear);
		break;
	case 2:
		opts->nregression_fns = 1 + 2*nparams;
		makeHVector = *(makeHVector_quadratic);
		break;
	case 3:
		opts->nregression_fns = 1+ 3*nparams;
		makeHVector = *(makeHVector_cubic);
		break;
	default:
		opts->nregression_fns = 1;
		makeHVector = *(makeHVector_trivial);
		break;
	}
	printf("# set regression order to: %d\n", opts->regression_order);
	printf("# set nregression_fns to: %d\n", opts->nregression_fns);
		
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
		covariance_fn = covariance_fn_matern_three;
		makeGradMatLength = derivative_l_matern_three;
		if(options->nthetas != 3)
			fprintf(stderr, "# (warn) setup_cov_fn has changed nthetas, potential memory errors abound\n");
		options->nthetas = 3;
		fprintf(stderr, "# cov_fn: MATERN32\n");
		break;
	case MATERN52:
		covariance_fn = covariance_fn_matern_five;
		makeGradMatLength = derivative_l_matern_five;

		if(options->nthetas != 3)
			fprintf(stderr, "# (warn) setup_cov_fn has changed nthetas to, potential memory errors abound\n");
		options->nthetas = 3;
		fprintf(stderr, "# cov_fn: MATERN52\n");
		break;
	case POWEREXPCOVFN:
		covariance_fn = covariance_fn_gaussian;
		makeGradMatLength = derivative_l_gauss;

		if(options->nthetas != (options->nparams+2))
			fprintf(stderr, "# (warn) setup_cov_fn has changed nthetas from %d, potential memory errors\n", options->nthetas);

		options->nthetas = options->nparams+2;
		fprintf(stderr, "# cov_fn: POWEREXP\n");
		break;
	case POWEREXPALPHA:
		covariance_fn = covariance_fn_gaussian_alpha;
		makeGradMatLength = derivative_l_gauss_alpha;

		if(options->nthetas != (options->nparams+3))
			fprintf(stderr, "# (warn) setup_cov_fn has changed nthetas from %d, potential memory errors\n", options->nthetas);

		options->nthetas = options->nparams+3;
		fprintf(stderr, "# cov_fn: POWEREXP+ALPHA\n");
		break;
		
	default:
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
 *  
 * this is called from: emulator-main.c, main.c and rbind.c
 */
void setup_optimization_ranges(optstruct* options, modelstruct* the_model)
{
	int i = 0;
	char buffer[128];
	int thetaLengthOffset = 2;
	double low, high;
	double bigRANGE = 10.0;
	double rangeMin = 0.0, rangeMax = 0.0;
	double fixedNuggetLeeWay = 0.0 ;
	double rangeMinLog = -15;

	double rangeMinNugget = 0.0000001;
	double rangeMaxNugget = 0.01; //what's a sensible upper limit here?

	double rangeMinAlpha = 1.0;
	double rangeMaxAlpha = 2.00001;
	/** 
	 * alloc the grad_ranges matrix in the options and 
	 * put in some sensible defaults 
	 */
	options->grad_ranges = gsl_matrix_alloc(options->nthetas, 2);

	/* 
	 * setup that we're using a log scale
	 */
	if(options->cov_fn_index == POWEREXPCOVFN || options->cov_fn_index == POWEREXPALPHA ){
		rangeMin = rangeMinLog;
		rangeMax = 5;
	} else {
		rangeMin = 0;
		rangeMax = bigRANGE;
	}

	gsl_matrix_set(options->grad_ranges, 0, 0, 0.0);
	gsl_matrix_set(options->grad_ranges, 0, 1, rangeMax);


	gsl_matrix_set(options->grad_ranges, 1, 0, rangeMinNugget);
	gsl_matrix_set(options->grad_ranges, 1, 1, rangeMaxNugget);

	if(options->cov_fn_index == POWEREXPALPHA){
		gsl_matrix_set(options->grad_ranges, 2, 0, rangeMinAlpha);
		gsl_matrix_set(options->grad_ranges, 2, 1, rangeMaxAlpha);
		thetaLengthOffset = 3;
	}

	printf("covfnindex: %d\n", options->cov_fn_index);
	printf("thetaLengthOffset: %d\n", thetaLengthOffset);
	printf("nthetas: %d\n", options->nthetas);


	if(options->use_data_scales){ 
		// use length scales set by the data
		for(i = thetaLengthOffset; i < options->nthetas; i++){

			if(options->cov_fn_index == POWEREXPCOVFN || options->cov_fn_index == POWEREXPALPHA){
				rangeMin = 0.5*log(gsl_vector_get(the_model->sample_scales, i-thetaLengthOffset));
				} else {
				rangeMin = 0.5*(gsl_vector_get(the_model->sample_scales, i-thetaLengthOffset));
			}
			
			if(rangeMin > rangeMax){
				fprintf(stderr, "#ranges failed\n");
				printf("# %d ranges: %lf %lf\n", i, rangeMin, rangeMax);
				printf("# sampleScale: %lf\n", gsl_vector_get(the_model->sample_scales, i-thetaLengthOffset));
				exit(EXIT_FAILURE);
			}
			gsl_matrix_set(options->grad_ranges, i, 0, rangeMin);
			gsl_matrix_set(options->grad_ranges, i, 1, rangeMax);
		} 

		if(isinf(rangeMin) == -1){
			rangeMin = 0.00001;
		}

	} else { // use some default scales

		for(i = thetaLengthOffset; i < options->nthetas; i++){
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
	 * print the ranges, no this is annoying
	 */
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

