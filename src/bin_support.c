#include "bin_support.h"

//! print the short-option switches
void print_usage_estimator(void){
	printf("estimator\n");
	printf("options are: \n");
	printf("t->number of thetas should be (2+nparams) for gaussian or 4 for matern\n");
	printf("p->number of params\n");
	printf("n->number of model_points\n");
	printf("m->number of emulator points\n");
	printf("a->min emulator value\n");
	printf("b->max emulator value\n");
}

//! print the short-option switches
void print_usage_emulator(void){
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
				print_usage_estimator();
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
	sprintf(options->outputfile, "emulator-out.txt");

	assert(options->nthetas >0);
	assert(options->nparams >0);

	/**!!!! set the number of regression fns
	 * 
	 * this is regression model dependant
	 * 
	 * (1+ options->nparams) is correct for the simple linear fit in each dimension plus a constant intercept
	 * 
	 */
	options->nregression_fns =  1 + options->nparams;


	// use the data scales
	options->use_data_scales = 1;
}


/**
 * \brief set the covariance_fn pointer to the gaussian covariance fn
 *
 * sets the function pointer to the gaussian cov funtion. because i suck there's actually
 * a broken interface lurking here. the other cov fns may have more or fewer args than the
 * gaussian one. Be careful when changing. 
 * 
 * you *shouldn't* need to change cov fn, the gaussian one is quite suitable
 * \todo cov_fn_alpha is broken
 * \todo test changing the fn ptr
 */
void setup_cov_fn(optstruct *options){
	/*
	 * we'll use the gaussian covariance fn by default
	 */
	//message("using gaussian cov fn\n", 1); // this is annoying
	options->covariance_fn = covariance_fn_gaussian;
	options->cov_fn_alpha = 2.0;
	options->nthetas = options->nparams+2;

	/* message("using matern 3/2", 1); */
	/* options->covariance_fn = covariance_fn_matern_five; */
	/* options->nthetas = 3; */
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
 * Note that the gaussian covfn uses log scaled thetas and so the range of values is
 * log[0..upper] -> [-infty .. log(uppper)]
 * 
 * this is not the case for the other fns, but the rest of the code may have the assumption
 * frozen into it that the ranges *are* log-scaled.
 *
 *  
 * this is called from: emulator-main.c, main.c and rbind.c
 */
void setup_optimization_ranges(optstruct* options, modelstruct* the_model)
{
	int i;
	char buffer[128];
	double rangeMin, rangeMax;
	double fixedNuggetLeeWay;
	/** 
	 * alloc the grad_ranges matrix in the options and 
	 * put in some sensible defaults 
	 */
	options->grad_ranges = gsl_matrix_alloc(options->nthetas, 2);


	if(options->use_data_scales){ // use length scales set by the data
		for(i = 0; i < options->nthetas; i++){
			
			if(options->covariance_fn == covariance_fn_gaussian && i > 1){
				rangeMax = 10.0;
				// need to get the fucking xmodel too, poop
				rangeMin = 0.5*log(gsl_vector_get(the_model->sample_scales, i-2));
				if(rangeMin > rangeMax){
					fprintf(stderr, "#ranges failed\n");
					printf("# %d ranges: %lf %lf\n", i, rangeMin, rangeMax);
					printf("# sampleScale: %lf\n", gsl_vector_get(the_model->sample_scales, i-2));
					exit(EXIT_FAILURE);
				}
			} else {
				rangeMin = 0.00001;
				rangeMax = 5.0;
			}
			if(i == 0){
				printf("# %d ranges: %lf %lf (scale)\n", i, rangeMin, rangeMax); 
			} if (i == 1){
				printf("# %d ranges: %lf %lf (nugget)\n", i, rangeMin, rangeMax); 
			} else if (i > 1) {
				printf("# %d ranges: %lf %lf\n", i, rangeMin, rangeMax);
			}
			if(isinf(rangeMin) == -1){
				rangeMin = 0.00001;
			}
			gsl_matrix_set(options->grad_ranges, i, 0, rangeMin);
			gsl_matrix_set(options->grad_ranges, i, 1, rangeMax);
		}
	} else { // use some default scales
		/* these are log-scaled ranges, note the negative lower bound */
		for(i = 0; i < options->nthetas; i++){
			if(options->covariance_fn == covariance_fn_gaussian){
				gsl_matrix_set(options->grad_ranges, i, 0, -10.0);
				gsl_matrix_set(options->grad_ranges, i, 1, 5.0);	
			} else {
				/* these are the regular ranges */
				gsl_matrix_set(options->grad_ranges, i, 0, 0.00001);
				gsl_matrix_set(options->grad_ranges, i, 1, 10.0);
			}
		}
	
	}

	if(options->fixed_nugget_mode == 0){
	// and force the nugget to be small (this is still done by hand...)
		gsl_matrix_set(options->grad_ranges, 1, 0, 0.00001);
		gsl_matrix_set(options->grad_ranges, 1, 1, 0.005);
	} else {
		// force the nugget to be fixed_nugget +- 5%
		fixedNuggetLeeWay = 0.05*(options->fixed_nugget);
		gsl_matrix_set(options->grad_ranges, 1, 0, options->fixed_nugget - fixedNuggetLeeWay);
		gsl_matrix_set(options->grad_ranges, 1, 1, options->fixed_nugget + fixedNuggetLeeWay);
		printf("# %d ranges: %lf %lf (nugget)\n", i,options->fixed_nugget - fixedNuggetLeeWay
					 ,options->fixed_nugget + fixedNuggetLeeWay);
	}		

	/* for(i = 0; i < options->nthetas;i++){ */
	/* 	sprintf(buffer, "%d %g %g\n", i, gsl_matrix_get(options->grad_ranges, i, 0), gsl_matrix_get(options->grad_ranges, i, 1)); */
	/* 	message(buffer, 1); */
	/* } */

}

