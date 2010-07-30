#ifndef _INC_MAXLBFGS_
#define _INC_MAXLBFGS_

#include "estimate_threaded.h"
#include "pthread.h" // for debug info
#include "lbfgs.h"
#include "../useful.h"


double evalFnLBFGS(double *xinput, int nthetas, void* args);

//! wrapper for evalFnLBFGS
double evalLikelyhoodLBFGS_struct(struct estimate_thetas_params *p);

//! sets the vector x to some set of random initial values within the set of ranges given
void set_random_initial_value(gsl_rng* rand, gsl_vector* x, gsl_matrix* ranges,int  nthetas);

void maxWithLBFGS(struct estimate_thetas_params *params);


#endif
