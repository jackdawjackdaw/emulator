#ifndef _INC_MAXLBFGS_
#define _INC_MAXLBFGS_

#include "pthread.h" // for debug info
#include "lbfgs.h"
#include "../useful.h"
#include "maximise.h"


// for passing to the evaluate function
struct evalFnLBFGSArgs{
	int nparams;
	int nmodel_points;
	gsl_matrix* xmodel;
	gsl_vector* training_vector;
} evalFnLBFGSArgs;


double evalFnLBFGS(double *xinput, int nthetas, void* args);



void maxWithLBFGS(gsl_rng *rand, int max_tries, int nsteps, gsl_matrix *ranges, gsl_matrix *xmodel,
									gsl_vector *trainingvector, gsl_vector *thetas, int nmodel_points, int nethas, int nparams);


#endif
