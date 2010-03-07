#ifndef __INC_ESTIMATE_THREADED__
#define __INC_ESTIMATE_THREADED__

#include "main.h"
#include "pthread.h"


void estimate_thetas_threaded(gsl_matrix* xmodel_input, gsl_vector* training_vector, gsl_vector* thetas, optstruct* options);

// this won't work unless it has access to the nasty globals in estimate_threaded.c
void* estimate_thread_function(void* args);


//! used to pass the args into estimate_thread_function
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

// see the source for defs of the number of threads etc etc
#endif

