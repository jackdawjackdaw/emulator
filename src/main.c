#include "stdio.h"
#include "stdlib.h"
#include "unistd.h"
#include "libEmu/estimator.h"
#include "libEmu/emulator.h"
#include "libEmu/maximise.h"
#include "ioread.h"
#include "sys/time.h"
#include "useful.h"
#include "pthread.h"

/**
 * @file
 * @author Chris Coleman-Smith cec24@phy.duke.edu
 * @version 0.2
 * @section Description
 * 
 * The main apparatus to run the emulator, reads in options from the command line and
 * reads the actual data points from the given filename
 * 
 * threaded the estimator call function
 */
 

/* 
 * 1 -> read command line parameters
 * 2 -> read and process stream data
 * 3 -> estimate thetas
 * 4 -> calculate new means and variances
 * 5 -> output
 */


/* 
 * in the 1d case
 *  now there are only 3 hyperparams by default 
 *  -> vertical-scale theta0
 *  -> nugget theta1
 *  -> length-scale theta2...theta(Nparams-2⎈)
 */
#define NTHETASDEFAULT 4
#define NPARAMSDEFAULT 1
#define NEMULATEDEFAULT 64
#define EMULATEMINDEFAULT 0.0
#define EMULATEMAXDEFAULT 4.0

//! holds command line options
/** 
 * designed to hold basic command line 
 * options
 */
typedef struct optstruct{
	int nthetas;
	int nparams;
	int nmodel_points;
	int nemulate_points;
	double emulate_min;
	double emulate_max;
	char  filename[128];
} optstruct;

struct estimate_thetas_params{
	gsl_rng* random_number;
	int max_tries;
	int number_steps;
	gsl_vector* thetas;
	gsl_matrix* grad_ranges;
	gsl_matrix* model_input;
	gsl_vector* training_vector;
	int nmodel_points;
	int nthetas;
	int nparams;
} estimate_thetas_params;
// why do you need this ^ ident?


void print_usage(void);
void parse_arguments(int argc, char** argv, optstruct* options);
void emulate_model(gsl_matrix* xmodel, gsl_vector* training, gsl_vector*thetas, optstruct* options);
void estimate_thetas(gsl_matrix* xmodel_input, gsl_vector* training_vector, gsl_vector* thetas, optstruct* options);
void read_input_bounded(gsl_matrix* model, gsl_vector* training, optstruct * options);
void read_input_fromfile(gsl_matrix *xmodel, gsl_vector *training, optstruct *options);

void estimate_thetas_threaded(gsl_matrix* xmodel_input, gsl_vector* training_vector, gsl_vector* thetas, optstruct* options);
void* estimate_thread_function(void* args);

//! print the short-option switches
void print_usage(void){
	printf("emulator\n");
	printf("options are: \n");
	printf("t->number of thetas should be (2+nparams) for gaussian or 4 for matern\n");
	printf("p->number of params\n");
	printf("n->number of model_points\n");
	printf("m->number of emulator poits\n");
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
		options->nthetas = options->nparams +3;
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
	char** input_data;
	int number_lines = 0;

	parse_arguments(argc, argv, &options);	
	

	//sprintf(input_file, "%s",  "../short.dat");	
	sprintf(input_file, "%s",  "stdin");

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
			fprintf(stderr,"param: %s\n", split_string);
			gsl_matrix_set(xmodel_input, i, j, temp_value);
			split_string = strtok(NULL, "\t ");
		}
		assert(split_string != NULL);
		sscanf(split_string,"%lg", &temp_value);
		fprintf(stderr,"train: %s\n", split_string);
		gsl_vector_set(training_vector, i, temp_value);
		}

	fprintf(stderr, "read the following input matrix: %d x %d\n", options.nmodel_points, options.nparams);
	//print_matrix(xmodel_input, options.nmodel_points, options.nparams);
	fprintf(stderr, "the training data is:\n");
	//print_vector_quiet(training_vector, options.nmodel_points);
	
	fprintf(stderr, "nthetas = %d\n", options.nthetas);
	fprintf(stderr, "nparams = %d\n", options.nparams);

	estimate_thetas_threaded(xmodel_input, training_vector, thetas, &options);

	// calc the new means, new variance and dump to stdout
	emulate_model(xmodel_input, training_vector, thetas, &options);
	gsl_vector_free(thetas);
	gsl_vector_free(training_vector);
	gsl_matrix_free(xmodel_input);
	free_char_array(input_data, number_lines);
	//exit(1);
	return(0);
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
	gsl_vector_view new_x_row;

	gsl_matrix *temp_matrix = gsl_matrix_alloc(options->nmodel_points, options->nmodel_points);
	gsl_permutation *c_LU_permutation = gsl_permutation_alloc(options->nmodel_points);
	int lu_signum = 0;

	makeCovMatrix(c_matrix, xmodel, thetas,options->nmodel_points, options->nthetas, options->nparams);
	gsl_matrix_memcpy(temp_matrix, c_matrix);
	gsl_linalg_LU_decomp(temp_matrix, c_LU_permutation, &lu_signum);
	gsl_linalg_LU_invert(temp_matrix, c_LU_permutation, cinverse);
	
	// set the new_x values
	initialise_new_x(new_x, options->nparams, options->nemulate_points, options->emulate_min, options->emulate_max);



	for(i = 0; i < n_emu_points; i++){
		new_x_row = gsl_matrix_row(new_x, i);
		makeKVector(kplus, xmodel, &new_x_row.vector, thetas, options->nmodel_points, options->nthetas, options->nparams);
		temp_mean = makeEmulatedMean(cinverse, training, kplus, options->nmodel_points);
		kappa = covariance_fn(&new_x_row.vector, &new_x_row.vector, thetas, options->nthetas, options->nparams);
		temp_var = makeEmulatedVariance(cinverse, kplus, kappa, options->nmodel_points);
		gsl_vector_set(new_mean, i, temp_mean);
		gsl_vector_set(new_variance, i, temp_var);
	}
	 
	for(i = 0; i < n_emu_points; i++){
		for(j = 0; j < options->nparams; j++){
			printf("%g\t", gsl_matrix_get(new_x, i, j));
		}
		printf("%g\t", gsl_vector_get(new_mean, i));
		printf("%g\n", gsl_vector_get(new_variance, i));
	}
		
	gsl_matrix_free(new_x);
	gsl_vector_free(new_mean);
	gsl_vector_free(new_variance);
	gsl_matrix_free(c_matrix);
	gsl_matrix_free(cinverse);
	gsl_vector_free(kplus);
	gsl_matrix_free(temp_matrix);
	gsl_permutation_free(c_LU_permutation);
}

// globals for the threads to use
pthread_mutex_t job_counter_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t results_mutex = PTHREAD_MUTEX_INITIALIZER;
/* how many lots of thread_level_tries to do */
int ntries = 4; 
/* mutex protected counter to keep track of completed jobs */
int jobnumber = 0; 
/* global spot for the best thetas to be kept in */
gsl_vector *best_thetas;
/* the best value of theta */
double best_theta_val = -1000;

#define NUMBERTHREADS 2


/* threaded estimate thetas */
void estimate_thetas_threaded(gsl_matrix* xmodel_input, gsl_vector* training_vector, gsl_vector* thetas, optstruct* options){

	/* thread data */
	int nthreads = NUMBERTHREADS;
	/* how many attempts to maximise should we make */
	/* each thread will make this number of tries and then compare its best values 
	 * to the ones in best_thetas, if it wins it will save them
	 * we only care about the *best* so it doesn't matter if we just throw 
	 * the rest out the window... 
	 */
	int thread_level_tries = 5; 
	best_thetas = gsl_vector_alloc(options->nthetas);

	pthread_t threads[NUMBERTHREADS];
	struct estimate_thetas_params params[NUMBERTHREADS];
	
	
	/* regular stuff */
	const gsl_rng_type *T;
	
	int i; 
	int number_steps = 10;
	gsl_matrix *grad_ranges = gsl_matrix_alloc(options->nthetas, 2);
	T = gsl_rng_default;

	
	/* set the ranges for the initial values of the NM lookup, 
	 * might want to adjust these as required etc, but whatever */
	for(i = 0; i < options->nthetas; i++){
		gsl_matrix_set(grad_ranges, i, 0, 0.1);
		gsl_matrix_set(grad_ranges, i, 1, 5.5);
		gsl_vector_set(best_thetas, i, 0.0);
	}


	/* setup the thread params */
	for(i = 0; i < nthreads; i++){
		// alloc a rng for each thread
		params[i].random_number = gsl_rng_alloc(T);
		// this is blocking right now (slooow)
		gsl_rng_set(params[i].random_number, get_seed_noblock());
		// not sure about this, perhaps each thread should make 10 tries
		params[i].max_tries = thread_level_tries;
		params[i].thetas = gsl_vector_alloc(options->nthetas);
		params[i].grad_ranges = gsl_matrix_alloc(options->nthetas, 2);		
		params[i].model_input = gsl_matrix_alloc(options->nmodel_points, options->nparams);
		params[i].training_vector = gsl_vector_alloc(options->nmodel_points);
		params[i].nmodel_points = options->nmodel_points;
		params[i].nthetas = options->nthetas;
		params[i].nparams = options->nparams;
		params[i].number_steps = number_steps;
		// now actually copy the stuff into the vectors / matrices
		gsl_vector_memcpy(params[i].thetas, thetas);
		gsl_matrix_memcpy(params[i].grad_ranges, grad_ranges);
		gsl_matrix_memcpy(params[i].model_input, xmodel_input);
		gsl_vector_memcpy(params[i].training_vector, training_vector);
		
	}

	// create the threads
	for(i = 0; i < nthreads; i++)
		pthread_create(&threads[i], NULL, &estimate_thread_function, &params[i]);
	
	// wait to rejoin
	for(i = 0; i < nthreads; i++)
		pthread_join(threads[i], NULL);


	fprintf(stderr, "best_thetas: \t");
	print_vector_quiet(best_thetas, options->nthetas);

	// tear down the thread params
	for(i = 0; i < nthreads; i++){
		gsl_rng_free(params[i].random_number);
		gsl_matrix_free(params[i].grad_ranges);
		gsl_matrix_free(params[i].model_input);
		gsl_vector_free(params[i].thetas);
		gsl_vector_free(params[i].training_vector);
	}

	// copy the global best_theta into the one provided 
	gsl_vector_memcpy(thetas, best_thetas);
	// now free best_thetas
	gsl_vector_free(best_thetas);
	
}	

// THIS USES GLOBAL VARIABES DEFINED ABOVE, WATCH OUT!
// what the threads actually call, put this in here so that the function
// will share the same scope as the rest of the crap here
void* estimate_thread_function(void* args){
	// cast the args back
	struct estimate_thetas_params *p = (struct estimate_thetas_params*) args;
	int next_job;
	int my_id = pthread_self();
	double my_theta_val = 0.0;
	while(1){
		/* see if we've done enough */
		pthread_mutex_lock(&job_counter_mutex);
		if(jobnumber == ntries){
			next_job = -1;
		} else {
			next_job = jobnumber;
			jobnumber++;
			printf("job: %d by %d\n", next_job, my_id); 
		}
		/* now we can unlock the job counter */
		pthread_mutex_unlock(&job_counter_mutex);
		
		/* we're done so stop */
		if(next_job == -1)
			break;
		
		/* else we do the nelder mead stuff */
		nelderMead(p->random_number, p->max_tries, p->number_steps, p->thetas, p->grad_ranges, p->model_input, p->training_vector, p->nmodel_points, p->nthetas, p->nparams);

		// kind of sneakily calling into the maximise.c api (aah well...)
		my_theta_val = evalLikelyhood(p->thetas, p->model_input, p->training_vector, p->nmodel_points, p->nthetas, p->nparams);

		pthread_mutex_lock(&results_mutex);
		printf("results locked by %d\n", my_id);
		if(my_theta_val > best_theta_val){
			// this thread has produced better thetas than previously there
			gsl_vector_memcpy(best_thetas, p->thetas); // save them
			// save the new best too
			best_theta_val = my_theta_val;
			printf("thread %d, won with %g\n", my_id, my_theta_val);
		}
		pthread_mutex_unlock(&results_mutex);
		printf("results unlocked by: %d\n", my_id);
	}
	// and relax...
	return NULL;
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



//! reads nmodel_points from the stdin
/* void read_input_bounded(gsl_matrix* model, gsl_vector* training, optstruct * options){ */
/* 	int i = 0; */
/* 	int j = 0; */
/* 	double temp_value;	  */
/* 	while(i < options->nmodel_points){ */
/* 		for(j = 0; j < options->nparams; j++){ */
/* 			scanf("%lg", &temp_value); */
/* 			//printf("%lg\n", temp_value); */
/* 			gsl_matrix_set(model, i, j, temp_value); */
/* 		} */
/* 		scanf("%lg", &temp_value); */
/* 		//printf("%lg\n", temp_value); */
/* 		gsl_vector_set(training, i, temp_value); */
/* 		i++; */
/* 	}; */

/* 	printf("read in xmodel:\n"); */
/* 	print_matrix(model, options->nmodel_points, options->nparams); */
/* 	printf("read in training_vec:\n"); */
/* 	vector_print(training, options->nmodel_points); */
/* } */



