#ifndef __INC_MAXBFGS__
#define __INC_MAXBFGS__

#include "bfgs.h"
#include "../useful.h"
#include "maximise.h"


void maxWithBFGS(gsl_rng *rand, int max_tries, int nsteps, gsl_matrix *ranges, gsl_matrix* xmodel,
								 gsl_vector *trainingvector, gsl_vector* thetas, int nmodel_points, int nthetas, int nparams, int nregression_fns);


#endif
