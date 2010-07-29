#ifndef _INC_MAXLBFGS_
#define _INC_MAXLBFGS_

#include "estimate_threaded.h"
#include "pthread.h" // for debug info
#include "lbfgs.h"
#include "../useful.h"


// for passing to the evaluate function
struct evalFnLBFGSArgs{
	int nparams;
	int nmodel_points;
	int nregression_fns;
	gsl_matrix* xmodel;
	gsl_vector* training_vector;
	gsl_matrix* h_matrix;
} evalFnLBFGSArgs;


double evalFnLBFGS(double *xinput, int nthetas, void* args);

//! wrapper for evalFnLBFGS
double evalLikelyhoodLBFGS_struct(struct estimate_thetas_params *p);

//! inits  (and fills) the contents of an evalFnLBFGSArgs struct from an estimate_thetas_params struct
void setup_evalFnLBFGSArgs(struct evalFnLBFGSArgs *arguments, struct estimate_thetas_params *params);

//! sets the vector x to some set of random initial values within the set of ranges given
void set_random_initial_value(gsl_rng* rand, gsl_vector* x, gsl_matrix* ranges,int  nthetas);

void maxWithLBFGS(gsl_rng *rand, int max_tries, int nsteps, gsl_matrix *ranges, gsl_matrix *xmodel,
									gsl_vector *trainingvector, gsl_vector *thetas, int nmodel_points, int nethas, int nparams, int nregression_fns);


#endif
