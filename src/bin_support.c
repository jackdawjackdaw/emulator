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
	optoins->regression_order = 2;
	//options->nregression_fns =  1 + options->nparams;
	
	/* 
	 * use the powerexp cov fn by default 
   */
	options->cov_fn_index = POWEREXPCOVFN;
	
	// use the data scales
	options->use_data_scales = 1;
}






